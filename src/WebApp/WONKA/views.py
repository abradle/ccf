from django.shortcuts import render, get_object_or_404
# Import the models
from IOhandle.models import Protein,Target,Project,Molecule,ActivityPoint,Compound
from MMPMaker.models import MMPDiffMap
from WebApp.settings import DATABASES
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from WONKA.models import Residue, PharmaScore, FragScore, KeyCluster, ResShift,InternalIDLink
from OBSERVATIONS.models import Observation, Comment, UserData, ObservationChoices
from WONKA.functions import render_pharma_pdb
import bisect, urllib
from django.http import HttpResponse, HttpResponseRedirect
from Pharmacophore.models import PharmaPoint
from rdkit.ML.Cluster import Butina
from rdkit import SimDivFilters,DataStructs
from rdkit.Chem import rdMolDescriptors
import numpy
from rdkit.DataManip.Metric import rdMetricMatrixCalc
from django.core import serializers
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.contrib.auth import authenticate, login, logout
from rdkit.Chem import Draw
import StringIO
# Import smtplib for the actual sending function
import smtplib, hashlib
import struct, uuid
# Import the email modules we'll need
from email.mime.text import MIMEText
import random
from WONKA.functions import find_res_rmsd
from MMPMaker.models import Act3DFrag, ActPharmaPoint
from Viewer.models import eucl_dist
import math, ast
from django.db.models import Q
from rdkit.Chem import MCS
from djutils.decorators import async


def FindObservations(request):
    """Function to find observations given a given series of mol ids, compound ids etc"""
    # Get the host
    host = request.get_host()
    # Get the targets
    targets = Target.objects.exclude(title="DUMMY")
    # First get all the relevant IDs
    mol_ids = [int(x) for x in request.GET["mol_id"].split(",") if x]
    cmpd_ids = [Molecule.objects.get(pk=int(x)).cmpd_id.pk for x in request.GET["mol_id"].split(",") if x]
    target_id = int(request.GET["target_id"])
    water_ids = [int(x) for x in request.GET["water_id"].split(",") if x]
    ph4_ids = [int(x) for x in request.GET["ph4_id"].split(",") if x]
    prot_ids = [int(x) for x in request.GET["prot_id"].split(",") if x]
    # Now get the things that match all
    obs_c = ObservationChoices.objects.filter()
    for my_id in mol_ids:
        # Now filter this
        obs_c = obs_c.filter(mol_id=my_id)
    # Get all the observations info
    match_all_obs = Observation.objects.filter(pk__in=obs_c.values_list("obs_id", flat=True))
    mol_obs = Observation.objects.exclude(pk__in=match_all_obs).filter(pk__in=ObservationChoices.objects.filter(mol_id__in=mol_ids).values_list("obs_id", flat=True))
    cmpd_obs = Observation.objects.exclude(pk__in=mol_obs).exclude(pk__in=match_all_obs).filter(pk__in=ObservationChoices.objects.filter(cmpd_id__in=cmpd_ids).values_list("obs_id", flat=True))
    # Get the ones to show
    for ob in mol_obs:
        ob.my_src = "http://" + host + '/WONKA/request_observation/?my_pk=' + str(ob.pk) + '&me=dl.icb'
    for ob in cmpd_obs:
        ob.my_src = "http://" + host + '/WONKA/request_observation/?my_pk=' + str(ob.pk) + '&me=dl.icb'
    for ob in match_all_obs:
        ob.my_src = "http://" + host + '/WONKA/request_observation/?my_pk=' + str(ob.pk) + '&me=dl.icb'
    # Get the ones to show info for
    # Now loop through and
    context = {"mol_obs": mol_obs,
               "targets": targets,
               "match_all_obs": match_all_obs,
               "cmpd_obs": cmpd_obs}
    return render(request, 'WONKA/FindObservations.html', context)


def seen_mol(request):
    """View to like a molecule for a person"""
    # Get the 
    mol_pk = int(request.GET["mol_pk"])
    # Get the user pk
    user_pk = int(request.GET["user_pk"])
    # Now get the user 
    user = User.objects.get(pk=user_pk)
    # Now get the data
    user_d = UserData.objects.get(user_id=user)
    # Get the molecule
    mol = Molecule.objects.get(pk=mol_pk)
    # Remove it from the list
    user_d.new_mols.remove(mol)
    # Now return this
    return HttpResponse(json.dumps({"mol_pk": mol_pk}))


def get_user_info(request, user_pk):
    """Function to get the target information for the user"""
    # Get the user from this
    user = User.objects.get(pk=user_pk)
    # Now get the user data for this one
    user_d = UserData.objects.get(user_id=user)
    targets = Target.objects.exclude(title="DUMMY")
    out_d = {}
    out_d["targs"] = []
    out_d["mols"] = []
    for targ in targets:
        new_d = {"name": targ.title}
        new_mols = user_d.new_mols.filter(prot_id__target_id=targ.pk)
        new_d["tot_new"] = len(new_mols)
        new_d["pk"] = targ.pk
        new_d["new"] = bool(len(user_d.new_targs.filter(pk=targ.pk)))
        new_d["cool"] = bool(len(user_d.cool_targs.filter(pk=targ.pk)))
        mols = Molecule.objects.filter(prot_id__target_id=targ).exclude(prot_id__code__contains=targ.title)
        new_d["mols"] = len(mols)
        # Get the activity points
        acts = ActivityPoint.objects.filter(target_id=targ)
        new_d["acts"] = len(acts)
        out_d["targs"].append(new_d)
        for mol in mols:
            out_d["mols"].append({"smiles": mol.smiles,
                                  "target": mol.prot_id.target_id.title,
                                  "pk": mol.pk,
                                  "cmpd_pk": mol.cmpd_id.pk,
                                  "new": bool(len(user_d.new_mols.filter(pk=mol.pk)))
                                  })
    return HttpResponse(json.dumps(out_d))


@async
def send_email(message, subject):
    """Function to send an e-mail - when users register or whartevs """
    # Create a text/plain message
    try:
        from WebApp import personalsettings as es
    except:
        print "YOU MUST SUPPLY E-MAIL AND PASSWORD"
        return
    # me == the sender's email address
    # you == the recipient's email address
    # Send the message via our own SMTP server.
    if es.email is None or es.password is None:
        print "YOU MUST SUPPLY E-MAIL AND PASSWORD"
        return
    e_mail = es.email
    msg = "\r\n".join([
  "From: " + e_mail,
  "To: " + e_mail,
  "Subject: " + subject,
  "",
  message
    ])
    username = e_mail
    password = es.email_password
    try:
        server = smtplib.SMTP('smtp.gmail.com:587', timeout=3)
    except:
        return
    server.starttls()
    server.login(username, password)
    server.sendmail(e_mail, e_mail, msg)
    server.quit()


def index(request):
    """View  to get the main page"""
    targets = Target.objects.exclude(title="DUMMY")#.filter(title="CDK2")
    projects = Project.objects.all()
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets, "projects": None}
    return render(request, 'WONKA/index.html', context)


def cov_col(my_col, max_num=2.0):
    """Function to take a float and convert it to a HEX color on a heat map"""
    # First cutoff at -2 and +2
    my_col = min(my_col, max_num)
    my_col = max(my_col, -max_num)
    my_col += max_num
    max_num += max_num
    out_t = ((max_num - my_col) * 255 / max_num, my_col * 255 / max_num, 0)
    out_h = '#%02x%02x%02x' % out_t
    return out_h


def get_oommppaa_infos(request, cluster_pk):
    """View to take a cluster PK and return OOMMPPAA data for that target from like points nearby"""
    # Get the cluster
    dist = 1.5#float(request.GET["dist"])
    kc = KeyCluster.objects.get(pk=cluster_pk)
    # Get all like points from OOMMPPAA
    if kc.type == "MMPFrag object":
        # Filter down to just this type
        dist = 2.5
        act3dfrags = Act3DFrag.objects.filter(act_id__target_id=kc.target_id, mmp_frag_id__smistore=kc.function)#mmp_dict[kc.function]) LEAVE FOR NOW
        print len(act3dfrags)
        # Now filter these down to within 1A of the cluster centre
        out_vals = [x.actpharmapoint_set.all()[0].act_change for x in act3dfrags if eucl_dist(x.mmp_frag_id, kc) < dist]
    elif kc.type == "PharmaPoint object":
        # Filter  down to just this type
        actpps = ActPharmaPoint.objects.filter(in_diff_15=True, num_diff_15__lte=6, target_id=kc.target_id, pharma_id__smiles=kc.function)#pp_dict[kc.function]) LEAVE FOR NOW
        print len(actpps)
        # Now filter these down to distance
        out_vals = [x.act_change for x in actpps if eucl_dist(x.pharma_id, kc) < dist]
        print [eucl_dist(x.pharma_id, kc) for x in actpps]
    if len(out_vals) == 0:
        # If it's an empty list
        return HttpResponse(json.dumps({"outcome": "error", "pk": cluster_pk}))
    # Return summary information about this group
    out_vals = sorted(out_vals, reverse=True)
    out_d = {"mean": sum(out_vals) / len(out_vals),
             "min": min(out_vals),
             "max": max(out_vals),
             "stdev": numpy.std(out_vals),
             "num": len(out_vals),
             "vals": out_vals,
             # Get the colors for these buttons
             "colors": [cov_col(x) for x in out_vals],
             "pk": cluster_pk
             }
    return HttpResponse(json.dumps(out_d))


def get_cluster(request, target_id):
    """Function to get all the clusters given a lam and a type"""
    # Get the lam
    lam = float(request.GET["lam"])
    # Get the method
    method = str(request.GET["method"]) + " object"
    kcs = KeyCluster.objects.filter(type=method)
    kcs = kcs.filter(lam=lam)
    kcs = kcs.filter(target_id=target_id)
    if "smiles" in request.GET:
        my_mol = Chem.MolFromSmiles(str(request.GET["smiles"]))
        kcs = kcs.filter(function=Chem.MolToSmiles(my_mol, isomericSmiles=True))
    # Get the actual clusters
    # Now get meta data about that cluster
    out_d = {}
    if method == "ActPharmaPoint" + " object":
        for kc in kcs:
            out_d[kc.function + "__" + str(kc.pk)] = [x.act_change for x in kc.actpharma_id.all()]
    elif method == "Act3DFrag" + " object":
        for kc in kcs:
            out_d[kc.function + "__" + str(kc.pk)] = [x.actpharmapoint_set.all()[0].act_change for x in kc.act3dfrag_id.all() if x.actpharmapoint_set.all()]
    else:
        print "METHOD NOT RECOGNISED"
    return HttpResponse(json.dumps(out_d))

@async
def send_email_user(message, subject, to_email, img_data=None):
    """Function to send an e-mail - when users register or whartevs """
    import base64
    import sys
    try:
        from WebApp import personalsettings as es
    except:
        print "YOU MUST SUPPLY E-MAIL AND PASSWORD"
        return
    e_mail = es.email
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEText import MIMEText
    from email.MIMEImage import MIMEImage
    msgRoot = MIMEMultipart('related')
    msgRoot['Subject'] = subject
    msgRoot['From'] = e_mail
    msgRoot['To'] = to_email
    msgRoot.preamble = 'This is a multi-part message in MIME format.'
    msgAlternative = MIMEMultipart('alternative')
    msgRoot.attach(msgAlternative)
    msgText = MIMEText(message)
    msgAlternative.attach(msgText)
    msgText = MIMEText(message, 'html')
    msgAlternative.attach(msgText)
    if img_data:
        image = MIMEImage(base64.b64decode(img_data), _subtype="png", name="WONKA_SCREEN_GRAB.png")
        msgRoot.attach(image)
    smtp = smtplib.SMTP('smtp.gmail.com:587', timeout=3)
    smtp.starttls()
    smtp.login(e_mail, es.email_password)
    smtp.sendmail(e_mail, to_email, msgRoot.as_string())
    smtp.quit()


def add_extra_info(request, my_obs):
    """Function to add the observation information"""
    # Now add the external information
    oc = ObservationChoices()
    oc.obs_id = my_obs
    oc.save()
    my_prots = request.POST["prots"].split(",")
    [oc.prot_id.add(Protein.objects.get(molecule__pk=p)) for p in my_prots if p]
    my_mols = request.POST["mols"].split(",")
    [oc.mol_id.add(Molecule.objects.get(pk=m)) for m in my_mols if m]
    [oc.cmpd_id.add(Molecule.objects.get(pk=m).cmpd_id) for m in my_mols if m]
    my_waters = request.POST["waters"].split(",")
    my_ph4s = request.POST["ph4s"].split(",")
    [oc.kc_id.add(KeyCluster.objects.get(pk=kc_pk)) for kc_pk in my_ph4s if kc_pk]
    [oc.kc_id.add(KeyCluster.objects.get(pk=kc_pk)) for kc_pk in my_waters if kc_pk]
    oc.save()
    print "SUCCESS"
    return oc


def register_observation(request):
    """Function to register a given observation"""
    if request.method == 'POST':
        if "comments" in request.POST and 'my_html' in request.POST and 'out_blob' in request.POST and 'user' in request.POST:
            new_obs = Observation()
            user = User.objects.get(pk=int(request.POST['user']))
            # HTML text
            new_obs.html_text = request.POST['my_html']
            # BLOB
            new_obs.icm_blob = request.POST['out_blob']
            new_obs.image_page = request.POST['out_blob_im']
            # The author
            new_obs.author = user.first_name
            # The user it is linked to
            new_obs.user_id = user
            # Any comments
            new_obs.comments = request.POST['comments']
            # The target
            targ = Target.objects.get(pk=int(request.POST['target']))
            new_obs.target_id = targ
            # New ob
            new_obs.likes = 1
            # Save the uuid
            new_obs.uuid = str(uuid.uuid4())
            # Now save it
            new_obs.save()
            extra_info = add_extra_info(request, new_obs)
            # Now send an e-mail to say it's been done!
            send_email(user.first_name + " said - " + request.POST['comments'], "New observation - " +targ.title + " by " + user.first_name + "  " + user.last_name)
            # Create the URL to send out
            host = request.get_host()
            # Define the src for it
            my_url = "http://"+ host +'/WONKA/' + new_obs.uuid + '/ShowObs/'
            ## Now send an e-mail to everybody else
            users = User.objects.filter(pk__in=UserData.objects.filter(cool_targs=targ).values_list("user_id", flat=True))
            print "USERS", users
            # Now  send the e-mails to the users
            out_m = "There is a new observation for a target you are watching: " + targ.title
            out_m += "<br> The comment was made by " + user.first_name
            out_m += "<br>This can be found at the <b><a href='"+my_url+"MYUSERTEXT'>following link</a><b>"
            out_m += "<br>This observation said: " + new_obs.comments
            out_m += "<br>"
            subject = "New observation for " + targ.title + " by " + user.first_name
            # Loop  through the relevant users
            im_src = "http://" + host + '/WONKA/request_observation/?my_pk=' + str(new_obs.pk) + '&opt=png'
            out_m += '<br> The following molecules are involved:<br>'
            for mol in extra_info.cmpd_id.all():
                out_m += '<img src="http://oommppaa.sgc.ox.ac.uk/Viewer/loader/?function=2DMOL&choice='+urllib.quote(str(mol.smiles))+'"> '
            my_mess = out_m + '<br> Here is the screenshot (also attached) <br><img src="' + im_src + '"><br>'
            for us in users.filter():
                send_email_user(my_mess.replace('MYUSERTEXT', '?user_pk=' + urllib.quote(str(us.password))), subject, us.email, new_obs.image_page)
            return HttpResponse(str(new_obs.uuid))
    # nothing went well
    return HttpResponse('FAIL!')


def subscribe_targ(request, target_id, user_id):
    """Function to subscribe a user to a target"""
    user = User.objects.get(pk=user_id)
    user_d = UserData.objects.get(user_id=user)
    target = Target.objects.get(pk=target_id)
    print "Subscribing " + user.first_name + " to " + target.title
    user_d.cool_targs.add(target)
    print "Subscribed"
    return HttpResponse('Success')


def unsubscribe_targ(request, target_id, user_id):
    """Function to unsubscribe a user to a target"""
    user = User.objects.get(pk=user_id)
    user_d = UserData.objects.get(user_id=user)
    target = Target.objects.get(pk=target_id)
    print "Unsubscribing " + user.first_name + " to " + target.title
    user_d.cool_targs.remove(target)
    print "Unsubscribed"
    return HttpResponse('Success')


def subscribe_to_disqus(request):
    """View to subscribe users to DISQUS discussions automatically
    Takes a thread id and a user email"""
    from django.conf import settings
    import urllib2
    # Thread id
    thread = request.GET["thread"]
    # User email
    email = request.GET["email"]
    # Get the public key
    DISQUS_PUBLIC_KEY = getattr(settings, 'DISQUS_PUBLIC_KEY', None)
    # Creat the curl
    # Run it on the server
    payload = {"api_key": DISQUS_PUBLIC_KEY,
               "thread": thread,
               "email": email}
    url = "https://disqus.com/api/3.0/threads/subscribe.json"
    data = urllib.urlencode(payload)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req).read()
    return HttpResponse(response)


def ShowObservations(request, target_id):
    """Function to show relevant observations"""
    host = request.get_host()
    if "AUTHOR" in request.GET:
        my_auth = int(request.GET["AUTHOR"])
        my_obs = Observation.objects.filter(target_id=target_id, user_id=my_auth)
    else:
        my_obs = Observation.objects.filter(target_id=target_id)
    targets = Target.objects.exclude(title="DUMMY")
    len(my_obs)
    # Get the ones to show
    obs = my_obs
    for ob in obs:
        ob.my_src = "http://" + host + '/WONKA/request_observation/?my_pk=' + str(ob.pk) + '&me=dl.icb'
    # Get the ones to show info for
    # Now loop through and
    context = {"obs": obs, "targets": targets}
    return render(request, 'WONKA/ShowObservations.html', context)


def update_old_obs(uuid, mol_pks, prot_codes, kc_pks):
    """Function to update historic observations"""
    # First get the observation
    o = Observation.objects.get(uuid=uuid)
    oc = ObservationChoices()
    oc.obs_id = o
    oc.save()
    print oc.pk
    # Now set the molecue
    [oc.mol_id.add(x) for x in mol_pks]
    # Now update the compounds
    [oc.cmpd_id.add(x) for x in Molecule.objects.filter(pk__in=mol_pks).values_list("cmpd_id", flat=True)]
    # Now set the proteins
    [oc.prot_id.add(x) for x in Protein.objects.filter(code__in=prot_codes)]
    # Now set the waters etc
    [oc.kc_id.add(x) for x in kc_pks]
    oc.save()


def ShowObs(request, uuid):
    """View to look at a single observation"""
    from django.conf import settings
    host = request.get_host()
    obs = Observation.objects.get(uuid=uuid)
    targets = Target.objects.exclude(title="DUMMY")
    if "user_pk" in request.GET:
        user = User.objects.get(password=request.GET["user_pk"])
        user.backend='django.contrib.auth.backends.ModelBackend'
        if user is not None:
            if user.is_active:
                login(request, user)
    # Define the src for it
    obs.my_src = "http://"+host+'/WONKA/request_observation/?my_pk='+str(obs.pk)+'&me=dl.icb'
    #  Get the target 
    target = obs.target_id
    # Get the the unique id for this one
    uniq_id = hashlib.sha224(str(obs.comments)+str(obs.target_id.title)+str(obs.author)).hexdigest()
    DISQUS_PUBLIC_KEY = getattr(settings, 'DISQUS_PUBLIC_KEY', None)
    DISQUS_SHORT_NAME = getattr(settings, 'DISQUS_SHORT_NAME', None)
    # Now loop through and
    context = {"short_name": DISQUS_SHORT_NAME,
               "public_key": DISQUS_PUBLIC_KEY,
               "obs": obs, "comments": [],
               "targets": targets, "target": target, "uniq_id": uniq_id}
    return render(request, 'WONKA/ShowObs.html', context)


def add_comment(request):
    """View to add a comment"""
    if request.method == 'POST':
        if "comments" in request.POST and 'observation' in request.POST and 'user' in request.POST:
          new_c = Comment()
          new_c.user_id = User.objects.get(pk=request.POST["user"])
          new_c.comment = request.POST["comments"]
          new_c.obs_id = Observation.objects.get(pk=request.POST["observation"])
          new_c.save()
          return HttpResponse('success')
    # nothing went well
    return HttpResponse('FAIL!!!!!')


def reg_dat(ud):
    """Register all targets and molecules"""
    targs = Target.objects.exclude(title="DUMMY")
    for target in targs:
        ud.new_targs.add(target)
        new_mols = Molecule.objects.filter(prot_id__target_id=target).exclude(prot_id__code__contains=target.title)
        for mol in new_mols:
            ud.new_mols.add(mol)


def register_user(request):
    """View to register a new user"""
    new_user = User()
    host = request.get_host()
    new_user.email = request.POST['EMail']
    new_user.username = request.POST['EMail'][:30]
    new_user.first_name = request.POST['FirstName']
    new_user.last_name = request.POST['LastName']
    new_user.password = request.POST['Password']
    new_user.set_password(new_user.password)
    try:
        new_user.validate_unique()
        new_user.save()
        # Now add the new mols and stuff
        my_new_user = authenticate(username=request.POST['EMail'][:30], password=request.POST['Password'])
        # Save the user data
        ud = UserData.objects.get_or_create(user_id=new_user)[0]
        reg_dat(ud)
        send_email("ROCK ON!!!", "New user - " + new_user.first_name + " " + new_user.last_name)
        subject = "Welcome to WONKA!"
        # Now  send the e-mails to the users
        out_m = "Welcome to WONKA! " + new_user.first_name
        out_m += "<br> You can register for targets on the <a href='http://"+host+"/WONKA'> main page </a>"
        out_m += "<br> Happy WONKAING!!!"
        send_email_user(out_m, subject, new_user.email)
        login(request, my_new_user)
        return HttpResponse('success')

    except ValidationError:
        # Now return
        return HttpResponse('User already present')


def login_user(request):
    """View to login a new user"""
    user = authenticate(username=request.POST['EMail'][:30], password=request.POST['Password'])
    if user is not None:
        if user.is_active:
            login(request, user)
            send_email("ROCK ON!!!", "User login - " + user.first_name + " " + user.last_name)
            # Redirect to a success page.
            return HttpResponse('success')
        else:
            # Return a 'disabled account' error message
            return HttpResponse('Account disabled')
    else:
        # Return an 'invalid login' error message.
        return HttpResponse('Invalid username or password')


def logout_user(request):
    """View to logout a new user"""
    logout(request)
    return HttpResponse('You have now logged out!')


def fp_view(request, target_id):
    """Function to get the fingerprints for a target"""
    target = Target.objects.get(pk=target_id)
    my_prot = Protein.objects.filter(target_id=target_id, pdb_info__isnull=False)[0]
    context = {"my_temp": my_prot.code, "target": target}
    return render(request, 'WONKA/fpview.html', context)


def mol_view(request):
    """Function to view a 2D depiction of a molecule -> as PNG"""
    my_choice = request.GET['choice'].split("_")[0]
    try:
        mol = Chem.MolFromSmiles(str(InternalIDLink.objects.filter(internal_id=my_choice)[0].mol_id.smiles))
    except IndexError:
        mol = Chem.MolFromSmiles(str(Molecule.objects.get(pk=my_choice).smiles))
    image = Draw.MolToImage(mol)
    output = StringIO.StringIO()
    image.save(output, format="PNG")
    contents = output.getvalue()
    return HttpResponse(contents)


def request_observation(request):
    """View to request an observation"""
    import binascii
    observation = Observation.objects.get(pk=int(request.GET['my_pk']))
    if "opt" in request.GET:
        return HttpResponse(binascii.a2b_base64(observation.image_page))
    else:
        return HttpResponse(binascii.a2b_base64(observation.icm_blob))


def dmat_sim(fps, ntopick):
    """Function to pick a series of mols based on their fingerprints"""
    ds = []
    for i in range(1,len(fps)):
        ds.extend(DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i],returnDistance=True))
    mmp = SimDivFilters.MaxMinPicker()
    ids = mmp.Pick(numpy.array(ds),len(fps),ntopick)
    return ids


def annotate_mols(new_mols, i, host, button_list, color_list, target_id):
    """Function to annotate molecules wih the relevant colour and URL information"""
    for mol in new_mols:
        try:
          mol.acts = ActivityPoint.objects.filter(cmpd_id=mol.cmpd_id, target_id=target_id)[0].activity 
        except IndexError:
          mol.acts = []
        mol.button = button_list[i]
        mol.bg = color_list[i]
        mol.code = mol.prot_id.code
        mol.src = '/Viewer/loader/?function=2DMOL&choice='+urllib.quote(mol.smiles, '')
        mol.url_centre = "http://"+host+'/Viewer/loader/?function=VIEWMOLPK&choice='+str(mol.pk)
        mol.url_mol = "http://"+host+'/Viewer/loader/?function=VIEWMOLPK&choice='+str(mol.pk)
    return new_mols


def get_map(request, mol_id):
    """Function to return a map -> specific to the SGC """
    # Get the molecule
    mol = Molecule.objects.get(pk=mol_id)
    prot = mol.prot_id
    # Get the target
    target = prot.target_id.title
    # Get the modelid
    model = prot.code.split("_")[0]
    # Get the CHAIN
    chain = prot.code.split("_")[1].upper()
    # file = m001.chainA.aligned.map
    file_path = "/opt/edensity/"+target+"/"+model+".chain"+chain+".aligned.map"
    # Now serve this file
    return HttpResponse(open(file_path, "rb").read())


def rgb(minimum, maximum, value):
    """Function to return a colour based on a heatscale"""
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    g = 255 - b - r
    return (r, g, b)


def cluster_fps(ms, cutoff=0.2):
    """RDKit function to cluster a set of molecules"""
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina
    # first generate the distance matrix:
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 3, 1024) for x in ms]
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])
    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs


def get_frags(request, target_id):
    """View to return all cmps less than a certain mol weight"""
    # Get the target
    target = Target.objects.get(pk=target_id)
    # Get all of them
    tot_mols = Compound.objects.filter(heavy_atom_count__lte=18, pk__in=list(set(Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title).values_list("cmpd_id", flat=True))))
    return HttpResponse(json.dumps([x.pk for x in tot_mols]))


def getMCSSmiles(mol, mcs):
    """Function get the smiles of an MCS from an MCS smarts and a molecule"""
    mcsp = Chem.MolFromSmarts(mcs.smarts)
    match = mol.GetSubstructMatch(mcsp)
    return Chem.MolFragmentToSmiles(mol, atomsToUse=match,
                                    isomericSmiles=True,
                                    canonical=False)


def view_clusts(request, target_id):
    """Function to return a JSON of the appropriate molecule clusters"""
    target = Target.objects.get(pk=target_id)
    tot_mols = Compound.objects.filter(pk__in=list(set(Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title).values_list("cmpd_id", flat=True))))
    len(tot_mols)
    if "CUTOFF" in request.GET:
        cutoff = float(request.GET["CUTOFF"])
    else:
        cutoff = 0.7
    clusters = cluster_fps([Chem.MolFromSmiles(str(x.smiles)) for x in tot_mols], cutoff=cutoff)
    out_d = {}
    # Loop through the clusters - assign an MCS and a mols
    for i, c in enumerate(clusters):
        out_d[i] = {}
        out_d[i]["mols"] = [tot_mols[x].pk for x in c]
        if len(c) > 1:
            # Now get the MCS
            res = MCS.FindMCS([Chem.MolFromSmiles(str(tot_mols[x].smiles)) for x in c])
            out_d[i]["mcs"] = res.smarts
        else:
            out_d[i]["mcs"] = "NONE"
    # Now return these clusters
    return HttpResponse(json.dumps(out_d))


def SplitSet(request, target_id):
    """View to split a dataset based on compound clustering etc."""
    # Return
    return render(request, 'WONKA/SplitSet.html', {"target_id": target_id})


def find_min_max(coord_list):
    """Function to return the min and the max of a set of coords"""
    x_coords = [x.x_com for x in coord_list]
    y_coords = [x.y_com for x in coord_list]
    z_coords = [x.z_com for x in coord_list]
    return max(x_coords), min(x_coords), max(y_coords), min(y_coords), max(z_coords), min(z_coords)

 
def SummariseBack(request, obs_id):
    """View to show a summarise view for an observation (returns to the previous environment)"""
    host = request.get_host()
    obs = Observation.objects.get(pk=obs_id)
    my_text = obs.html_text.split('<h4>Hi', 1)[0]
    my_text = my_text.replace("var bodycontent = body[0];", "var bodycontent = body[0];on_load_fun();")
    src = "http://"+host+'/WONKA/request_observation/?my_pk='+str(obs.pk)+'&me=dl.icb'
    return render(request, 'WONKA/Summarise.html', {"html_text": my_text, "my_src": src})


def get_ph4(ph4_lam, target_id, pharma_num_map, ph4_list, pharma_lookup, my_sites, color_list, my_mols, host):
    """Function to get the PH4 clusters"""
    print "GETTING ALL PH4S"
    all_ph4s = KeyCluster.objects.filter(lam=ph4_lam, target_id=target_id).filter(type="PharmaPoint object").exclude(function__in=["MagicSix", "WeakDonor", "Carbonyl"]).exclude(function__contains="Attach").order_by().order_by("size").reverse()
    # These are x,y and z limits either side - for each site
    # Now 
    my_d = {}
    for i, p in enumerate(all_ph4s):
        # Now assign the colour based on feature
        if p.function in pharma_num_map:
            my_ind = pharma_num_map[p.function]
        else:
            my_ind = 0
        my_d[color_list[my_ind]+"_"+str(i)] = [m.mol_id.pk for m in  p.pharma_id.all()]
        p.bg = ph4_list[my_ind]#color_list[my_ind]
        p.button = ph4_list[my_ind]#button_list[my_ind]
        if p.function in pharma_lookup:
            p.smiles = pharma_lookup[p.function]
        else:
            p.smiles = p.function
        p.my_list = Molecule.objects.filter(pk__in=[x.pk for x in my_mols]).filter(pk__in=p.pharma_id.filter().values_list("mol_id__pk",flat=True)).values_list("pk", flat=True)
        p.rel_size = float(len(p.my_list)) / float(len(my_mols))
        p.per_size = int(float(len(p.my_list)) / float(len(my_mols)) * 100)
        p.url_centre = "http://" + host + '/WONKA/show_point/?x_com='+str(p.x_com)+'&y_com='+str(p.y_com)+'&z_com='+str(p.z_com)
        p.url_ph4 = "http://" + host + '/WONKA/show_point/?x_com='+str(p.x_com)+'&y_com='+str(p.y_com)+'&z_com='+str(p.z_com)
        p.sites = []
        for my_site in my_sites:
            site_mols = my_site.mol_id.filter().values_list("pk", flat=True)
            if len([x for x in site_mols if x in p.my_list]) > 0:
                p.sites.append(my_site.pk) 
        p.sites = " ".join(["SITE"+ str(x) for x in  p.sites])
    return all_ph4s, my_d


def get_waters(x_max, waters, color_list, my_mols, button_list, host, my_sites, delta):
    """Function to get the water clusters"""
    print "GETTING TOP WATERS"
    # Filter
    if x_max:
#### Put in code here that finds the top waters -> e.g. from pharmacophore points  =NEEEDS FIXING
        big_waters = waters.filter(size__gte=len(my_mols) * 3 / 10)
        big_waters = big_waters.order_by("size").reverse()
        top_waters = []
        for i, w in enumerate(big_waters):
            w.bg = color_list[i % len(color_list)]
            w.button = button_list[i % len(button_list)]
            # Get all the mols that don't have a water in this cluster
            w.my_list = Molecule.objects.filter(pk__in=[x.pk for x in my_mols]).filter(pk__in=[x.pk for x in my_mols]).exclude(prot_id__in=w.water_id.all().values_list("prot_id__pk",flat=True)).values_list("pk", flat=True)
            # Now show the size the way it is
            w.rel_size = float(len(my_mols) - len(w.my_list)) / float(len(my_mols))
            w.per_size = int(w.rel_size * 100)
            w.url_centre = "http://" + host + '/WONKA/show_point/?x_com=' + str(w.x_com) + '&y_com=' + str(w.y_com) + '&z_com=' + str(w.z_com)
            w.sites = []
            for my_site in my_sites:
                near_waters = KeyCluster.objects.filter(pk=w.pk).filter(x_com__lte=my_site.x_max + delta, x_com__gte=my_site.x_min - delta)
                near_waters = near_waters.filter(y_com__lte=my_site.y_max + delta, y_com__gte=my_site.y_min - delta)
                near_waters = near_waters.filter(z_com__lte=my_site.z_max + delta, z_com__gte=my_site.z_min - delta)
                if len(near_waters) > 0:
                    w.sites.append(my_site.pk)
                if w.sites:
                    top_waters.append(w)
    else:
#### Put in code here that finds the top waters -> e.g. from pharmacophore points
        w_pk_list = []
        big_waters = waters.filter(size__gte=len(my_mols) * 3 / 10)
        big_waters = big_waters.order_by("size").reverse()
        top_waters = []
        for i, w in enumerate(big_waters):
            w.bg = color_list[i % len(color_list)]
            w.button = button_list[i % len(button_list)]
            # Get all the mols that don't have a water in this cluster
            w.my_list = Molecule.objects.filter(pk__in=[x.pk for x in my_mols]).exclude(prot_id__in=w.water_id.all().values_list("prot_id__pk",flat=True)).values_list("pk", flat=True)
            # Now show the size the way it is
            w.rel_size = float(len(my_mols) - len(w.my_list)) / float(len(my_mols))
            w.per_size = int(w.rel_size * 100)
            w.url_centre = "http://" + host + '/WONKA/show_point/?x_com=' + str(w.x_com) + '&y_com=' + str(w.y_com) + '&z_com=' + str(w.z_com)
            w.sites = []
            for my_site in my_sites:
                near_waters = KeyCluster.objects.filter(pk=w.pk).filter(x_com__lte=my_site.x_max + delta, x_com__gte=my_site.x_min - delta)
                near_waters = near_waters.filter(y_com__lte=my_site.y_max + delta, y_com__gte=my_site.y_min - delta)
                near_waters = near_waters.filter(z_com__lte=my_site.z_max + delta, z_com__gte=my_site.z_min - delta)
                if len(near_waters) > 0:
                    w.sites.append(my_site.pk)
                if w.sites:
                    if w.pk not in w_pk_list: 
                        top_waters.append(w)
                        w_pk_list.append(w.pk)
            w.sites = " ".join(["SITE"+ str(x) for x in w.sites])
    return top_waters


def find_res_shift(x_min, x_max, y_min, y_max, z_min, z_max, target_id, my_sites, res_two_three_dict, my_mols, color_list, button_list):
    """Function to find the relavant residue shifts"""
    print "FINDING MAX SHIFTS"
    max_shift = []
    # Get the delta value
    delta = 5.0
    # Filter residues to the ones within 1.0 A of any molecule AND then sort by size
    tot_res = Residue.objects.filter(target_id=target_id)
    if x_max:
        criterion1 = Q(x_max__gte=x_max + delta)
        criterion2 = Q(x_max__gte=x_min + delta)
        near_res = tot_res.exclude(criterion1 & criterion2)
        criterion1 = Q(x_min__lte=x_max - delta)
        criterion2 = Q(x_min__lte=x_min - delta)
        near_res = near_res.exclude(criterion1 & criterion2)
        criterion1 = Q(y_max__gte=y_max + delta)
        criterion2 = Q(y_max__gte=y_min + delta)
        near_res = near_res.exclude(criterion1 & criterion2)
        # Now do y_min
        criterion1 = Q(y_min__lte=y_max - delta)
        criterion2 = Q(y_min__lte=y_min - delta)
        near_res = near_res.exclude(criterion1 & criterion2)
        # Now do Z
        # First Z_max
        criterion1 = Q(z_max__gte=z_max + delta)
        criterion2 = Q(z_max__gte=z_min + delta)
        near_res = near_res.exclude(criterion1 & criterion2)
        # Now Z min
        criterion1 = Q(z_min__lte=z_max - delta)
        criterion2 = Q(z_min__lte=z_min - delta)
        near_res = near_res.exclude(criterion1 & criterion2)
        near_res = set(near_res.filter().values_list("res_name", "res_num"))
    else:
        tot_near_res = []
        tot_res_d = {}
        for my_site in my_sites:
            criterion1 = Q(x_max__gte=my_site.x_max + delta)
            criterion2 = Q(x_max__gte=my_site.x_min + delta)
            near_res = tot_res.exclude(criterion1 & criterion2)
            criterion1 = Q(x_min__lte=my_site.x_max - delta)
            criterion2 = Q(x_min__lte=my_site.x_min - delta)
            near_res = near_res.exclude(criterion1 & criterion2)
            criterion1 = Q(y_max__gte=my_site.y_max + delta)
            criterion2 = Q(y_max__gte=my_site.y_min + delta)
            near_res = near_res.exclude(criterion1 & criterion2)
            # Now do y_min
            criterion1 = Q(y_min__lte=my_site.y_max - delta)
            criterion2 = Q(y_min__lte=my_site.y_min - delta)
            near_res = near_res.exclude(criterion1 & criterion2)
            # Now do Z
            # First Z_max
            criterion1 = Q(z_max__gte=my_site.z_max + delta)
            criterion2 = Q(z_max__gte=my_site.z_min + delta)
            near_res = near_res.exclude(criterion1 & criterion2)
            # Now Z min
            criterion1 = Q(z_min__lte=my_site.z_max - delta)
            criterion2 = Q(z_min__lte=my_site.z_min - delta)
            near_res = near_res.exclude(criterion1 & criterion2)
            # Now we get the near res for this site
            near_res = set(near_res.filter().values_list("res_name", "res_num"))
            for res in near_res:
                if res in tot_res_d:
                    tot_res_d[res].append(my_site.pk)
                else:
                    tot_res_d[res] = [my_site.pk]
            tot_near_res.extend(list(near_res))
        near_res = tot_near_res
    print "Getting clusters"
    my_res = ResShift.objects.filter(target_id=target_id, res_name__in=[x[0] for x in near_res], res_num__in=[x[1] for x in near_res])
    # Only find those close to the BOX / main
    out_res_d = {}
    for i, val in enumerate(sorted(my_res.values_list("max_shift", "res_name", "pk", "res_num"),reverse=True)):
        my_mol = Molecule()
        # Define the site the residues are in
        res_hash = (val[1], val[3])
        if res_hash in tot_res_d:
            my_mol.sites = " ".join(["SITE"+ str(x) for x in tot_res_d[res_hash]])
        #my_mol.my_list = [(x[0]) for x in sorted(ResShift.objects.filter(target_id=target).values_list("max_shift"),reverse=True)[:5]]
        if val[1] in res_two_three_dict:
            this_res_name = res_two_three_dict[val[1]]
        else:
            this_res_name = "UNI"
        my_mol.res = "^" + this_res_name + str(val[3])
        out_res_d[my_mol.res] = {}
        my_mol.my_name = val[1] + ": " + str(val[3])
        my_mol.shift = val[0]
        my_mol.button = button_list[i % len(button_list)]
        my_mol.bg = color_list[i % len(color_list)]
        my_mol.res_cl = {}
        # Now get how the molecules rank on this residue move
        # instead we want to go trhrough molecules
        my_mol.my_list = []
        # Now colour the clusters
        for item in my_mols:
            this_res = tot_res.filter(res_name=val[1], res_num=val[3],
                                  prot_id__molecule=item)
            if len(this_res) ==0:
                new_mol = Molecule()
                # Get the PK from here
                new_mol.pk = item.pk
                new_mol.shift = 0.0
                new_mol.colour = ""
                out_res_d[my_mol.res][item.prot_id.code] = ""
                my_mol.my_list.append(new_mol)
            elif len(this_res) == 1:
                this_res = this_res[0]
                new_mol = Molecule()
                # Get the PK from here
                new_mol.pk = item.pk
                new_mol.shift = this_res.max_shift
                new_mol.clus_id = "RESCL" + str(this_res.clust_id) + "_" +  val[1] + "_" + str(val[3])
                my_mol.res_cl["RESCL" + str(this_res.clust_id) + "_" +  val[1] + "_" + str(val[3])] = [color_list[this_res.clust_id % len(color_list)],  button_list[this_res.clust_id % len(button_list)]]
                new_mol.colour = color_list[this_res.clust_id % len(color_list)]
                out_res_d[my_mol.res][this_res.prot_id.code] = button_list[this_res.clust_id % len(button_list)]
                my_mol.my_list.append(new_mol)
            else:
                print "ERROR MORE THAN ONE MOLS"
        # Now append this guy to the list
        max_shift.append(my_mol)
    return json.dumps(out_res_d), max_shift


from threading import Thread

class ThreadWithReturnValue(Thread):
    """Function takenfrom stack overflow - > generate a thread with return value"""
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs, Verbose)
        self._return = None
    def run(self):
        if self._Thread__target is not None:
            self._return = self._Thread__target(*self._Thread__args,
                                                **self._Thread__kwargs)
    def join(self):
        Thread.join(self)
        return self._return


def return_wonka_summarise_context(request, target_id):
    """Function to find the WONKA summarise context -> to be multithreaded"""
    # Now get the host - for making URLS
    host = request.get_host()
    x_max = False
    targets = Target.objects.exclude(title="DUMMY")
    ### These configure the different options for naming and colouring
    # Dict to convert 3 letter codes to one letter codes
    res_two_three_dict = {'ALA': "A",
                      'ARG': "R",
                      'ASN': "N",
                      'ASP': "D",
                      'CYS': "C", 
                      'GLU': "E",
                      'GLN': "Q",
                      'GLY': "G",
                      'HIS': "H",
                      'HIP': "H",
                      "HIE": "H",
                      'ILE': "I",
                      'LEU': "L",
                      'LYS': "K",
                      'MET': "M",
                      'PHE': "F",
                      'PRO': "P",
                      'SER': "S",
                      'THR': "T",
                      'TRP': "W",
                      'TYR': "Y",
                      'VAL': "V"}
    button_list = {1: "green", 2: "purple", 3: "orange", 4: "pink", 0: "grey", 5: "yellow", 6: "lightyellow", 7: "cyan", 8: "black", 9: "brown"}
    ph4_list = {1: "red", 2: "blue", 3: "brown", 4: "cyan", 0: "grey", 5: "green"}
    color_list = {1: "#36A603", 2: "#672997", 3: "#9B3900", 4: "pink", 0: "#656565", 5: "#969608", 6: "#d9d97f", 7: "#00bdbd", 8: "#000000", 9: "#564f00"}
    smiles_lookup = {"[*:1]C": "methyl", "[*:1]O": "hydroxyl","[*:1]N": "amine", "[*:1]S(N)(=O)=O": "sulphonamide",'[*:1][N+](=O)[O-]': "Nitro",'[*:1]C(N)=O':'terminal amide'}
    pharma_lookup = {"SingleAtomAcceptor": "Acceptor","SingleAtomDonor": "Donor", "RH6_6": "Six ring", "ThreeWayAttach": "Hydrophobic","Methyl": "Methyl"}
    pharma_num_map = {u'PosN': 0,
     u'ZnBinder5': 0,
     u'RH6_6': 3,
     u'Methyl': 3,
     u'SingleAtomDonor': 2,
     u'RH5_5': 3,
     u'iPropyl': 3,
     u'Carbonyl': 5,
     u'BasicGroup': 0,
     u'SingleAtomAcceptor': 1,
     u'TriFluro': 5,
     u'SingleF': 5,
     u'SmallHalogen': 5,
     u'SingleCl': 5,
     u'BigHalogen': 5,
     u'AcidicGroup': 0,
     u'Arom6': 4,
     u'Hydrophobe': 3,
     u'Arom5': 4,
     u'Cyano': 5,
     u'Imidazole': 5}
    ######
    # When the clusters have a target_id
    target = Target.objects.get(pk=target_id)
    # Now pull the relevant option
    tot_mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)
    if "MY_CMPS" in request.GET:
        my_cmpds = [int(x) for x in request.GET["MY_CMPS"].split(",")]
        tot_mols = tot_mols.filter(cmpd_id__in=my_cmpds)
    mol_lam = 5.0
    mol_clusts = KeyCluster.objects.filter(lam=mol_lam, target_id=target_id, type="Molecule object").order_by("size").reverse()
    if "my_site" in request.GET:
        my_m_c = mol_clusts.get(pk=int(request.GET["my_site"]))
        mol_clusts = mol_clusts.filter(pk=int(request.GET["my_site"]))
        tot_mols = tot_mols.filter(pk__in=my_m_c.mol_id.filter().values_list("pk", flat=True))
    # First make sure two identical compounds are not shown for the same molecule cluster
    my_mols = []
    my_id = 0
    for clu_num, mol_clu in enumerate(mol_clusts):
        clu_cmps = []
        clu_mols = mol_clu.mol_id.filter(pk__in=tot_mols)
        if clu_num == my_id:
            if len(clu_mols) > 0:
                target.mol = clu_mols[0].prot_id.code
            else:
                my_id += 1
        for mol in clu_mols:
            if mol.cmpd_id.pk not in clu_cmps:
                mol.rmsd = 0.0
                mol.category = mol.prot_id.code
                my_mols.append(mol)
                clu_cmps.append(mol.cmpd_id.pk)
            else:
                # Calc the RMSD
                for mol_2 in my_mols:
                    if mol.cmpd_id.pk == mol_2.cmpd_id.pk:
                        mol_2.rmsd = find_res_rmsd(Chem.MolFromMolBlock(str(mol.sdf_info)), Chem.MolFromMolBlock(str(mol_2.sdf_info)), str(mol.cmpd_id.pk))
    ### NOW ADD THE COMPOUND SIMILARITY DATA
    clusters = cluster_fps([Chem.MolFromSmiles(str(x.cmpd_id.smiles)) for x in my_mols], cutoff=0.7)
    for i, clus in enumerate(clusters):
        for mol_in_clu in clus:
            mol_out = my_mols[mol_in_clu]
            mol_out.cluster = i
    # Get the delta value - from where a point can project
    delta = 1.5
    # Define the Lamda
    water_lam = 1.5
    my_lam = 2.0
    kcs = KeyCluster.objects.filter(lam=my_lam, target_id=target_id)
    # Filter based on distances
    try:
        maps = request.GET['map']
    except KeyError:
        maps = None
    # If you get a dict of map coords - given by the box
    if maps:
        my_d = ast.literal_eval(maps)
        points = [float(x) for x in my_d["my_var"]]
        # These are x,y and z limits either side
        x_max = max(points[3], points[0])
        x_min = min(points[3], points[0])
        y_max = max(points[4], points[1])
        y_min = min(points[4], points[1])
        z_max = max(points[5], points[2])
        z_min = min(points[5], points[2])
        delta = 0.0
        kcs = kcs.filter(x_com__lte=x_max, x_com__gte=x_min)
        kcs = kcs.filter(y_com__lte=y_max, y_com__gte=y_min)
        kcs = kcs.filter(z_com__lte=z_max, z_com__gte=z_min)
    else:
        # In the event their are no molecule
        x_min = None
        x_max = None
        y_min = None
        y_max = None
        z_min = None
        z_max = None
        pass
    my_sites = KeyCluster.objects.filter(lam=mol_lam, target_id=target_id, type="Molecule object")
    my_frags = kcs.filter(type="MMPFrag object")
    frags = my_frags.filter(function__in=smiles_lookup)
#    rings = my_frags.filter(ringname__isnull=False)
    # Get the pharmacophores, waters and residues
    ph4_lam = 2.0
    waters = KeyCluster.objects.filter(lam=water_lam, target_id=target_id, type="Water object")
    tot_res = Residue.objects.filter(target_id=target_id).filter(prot_id__molecule__pk__in=[x.pk for x in my_mols])
    # A dictionary to hold all the lists against features
#    feat_list = {"key_frag": {}, "key_ph4": {}, "key_waters": {},
#               "wie_feat": {}, "max_shift": {}, "key_rings": {},
#               "all_ph4s": {}}
    feat_list = {# So get the waters, max shifts and pharmacophores. We need to do this site by site and then to 
                 }
    #water_thread = ThreadWithReturnValue(target=get_waters, args=(x_max, waters, color_list, my_mols, button_list, host, my_sites, delta))
    #water_thread.start()
    top_waters = get_waters(x_max, waters, color_list, my_mols, button_list, host, my_sites, delta)
    print "UPDATING MOLS"
    for mol in my_mols:
        i_ids = InternalIDLink.objects.filter(mol_id=mol).values_list("internal_id", flat=True)
        if i_ids:
            mol.my_id = i_ids[0]
            try:
                mol.acts = str(ActivityPoint.objects.filter(internal_id=mol.my_id, target_id=target_id)[0].activity)
            except IndexError:
                mol.acts = None
        else:
            mol_acts = ActivityPoint.objects.filter(cmpd_id=mol.cmpd_id, target_id=target_id)
            if mol_acts:
                mol.acts = str(mol_acts[0].activity)
            else:
                mol.acts = None
        mol.code = mol.prot_id.code
        mol.src = '/Viewer/loader/?function=2DMOL&choice='+urllib.quote(mol.prot_id.code, '')
    # Now get the clusters
    top_clusts = mol_clusts.order_by("size")
    for i, clus in enumerate(top_clusts):
        clus.bg = color_list[i % len(color_list)]
        clus.button = button_list[i % len(button_list)]
        clus.url_centre = "http://"+host+'/WONKA/show_point/?x_com='+str(clus.x_com)+'&y_com='+str(clus.y_com)+'&z_com='+str(clus.z_com)
        clus.url_ph4 = "http://"+host+'/WONKA/show_point/?x_com='+str(clus.x_com)+'&y_com='+str(clus.y_com)+'&z_com='+str(clus.z_com)
        clus.my_list = clus.mol_id.filter().values_list("pk", flat=True)
    # Get the thread return values
    all_ph4s, feat_list["all_ph4s"] = get_ph4(ph4_lam, target_id, pharma_num_map, ph4_list, pharma_lookup, my_sites, color_list, my_mols, host)
    out_res_d, max_shift = find_res_shift(x_min, x_max, y_min, y_max, z_min, z_max, target_id, my_sites, res_two_three_dict, my_mols, color_list, button_list)
    print "COMPACTING JSONs"
    feat_list["all_mols"] = [x.pk for x in my_mols]
    tot_list = json.dumps(feat_list)
    #ph4_pks = json.dumps([x.pk for x in all_ph4s])
    #frag_pks = json.dumps([x.pk for x in frags])
    context = {"target": target, "targets": targets,
               #"key_frag": top_frags, "key_ph4": top_ph4s,
               "key_waters": top_waters,
               "out_res_d": out_res_d,
               #"frag_pks": frag_pks,
               #"ph4_pks": ph4_pks,
               "max_shift": max_shift,
               "tot_mols": my_mols,
               "my_sites": my_sites,
               #"key_rings": top_rings,
               "tot_list": tot_list,
               # OOMMPPAA specific
               #"all_frags": frags.order_by("size").reverse(),
               #"coll_waters": coll_top_waters 
               "mol_clusts": top_clusts,
               "all_ph4s": all_ph4s
               }
    return context


def Summarise(request, target_id):
    """View to summarise the data for a given target"""
    context = return_wonka_summarise_context(request, target_id)
    # Now check to see if stuff is in
    host = request.get_host()
    ll_mols = []
    if "ALT_MOLS" in request.GET:
        # Now append these to the contex
        for mol_pk in request.GET["ALT_MOLS"].split(","):
            ll_mols.append(Molecule.objects.get(pk=int(mol_pk)))
    if "GET_LL" in request.GET:
        ll_mol_pks = Molecule.objects.filter(prot_id__target_id__title=request.GET["GET_LL"])
        for mol_pk in ll_mol_pks:
            ll_mols.append(mol_pk)
    if ll_mols:
        for mol in ll_mols:
            i_ids = InternalIDLink.objects.filter(mol_id=mol).values_list("internal_id", flat=True)
            if i_ids:
                mol.my_id = i_ids[0]
                try:
                    mol.acts = str(ActivityPoint.objects.filter(internal_id=mol.my_id, target_id=target_id)[0].activity)
                except IndexError:
                    mol.acts = None
            else:
                mol_acts = ActivityPoint.objects.filter(cmpd_id=mol.cmpd_id, target_id=target_id)
                if mol_acts:
                    mol.acts = str(mol_acts[0].activity)
                else:
                    mol.acts = None
            mol.code = mol.prot_id.code
            mol.src = "http://"+host+'/Viewer/loader/?function=VIEWCMPDPK&choice='+str(mol.cmpd_id.pk)
        context["method"] = request.GET["method"]
    context["ll_mols"] = ll_mols
    return render(request, 'WONKA/Summarise.html', context)


def SummariseORIG(request, target_id):
    """View to summarise the data for a given target (original -> activeICM version)"""
    context = return_wonka_summarise_context(request, target_id)
    # Now check to see if stuff is in
    host = request.get_host()
    ll_mols = []
    if "ALT_MOLS" in request.GET:
        # Now append these to the contex
        for mol_pk in request.GET["ALT_MOLS"].split(","):
            ll_mols.append(Molecule.objects.get(pk=int(mol_pk)))
        for mol in ll_mols:
            i_ids = InternalIDLink.objects.filter(mol_id=mol).values_list("internal_id", flat=True)
            if i_ids:
                mol.my_id = i_ids[0]
                try:
                    mol.acts = str(ActivityPoint.objects.filter(internal_id=mol.my_id, target_id=target_id)[0].activity)
                except IndexError:
                    mol.acts = None
            else:
                mol_acts = ActivityPoint.objects.filter(cmpd_id=mol.cmpd_id, target_id=target_id)
                if mol_acts:
                    mol.acts = str(mol_acts[0].activity)
                else:
                    mol.acts = None
            mol.code = mol.prot_id.code
            mol.src = "http://"+host+'/Viewer/loader/?function=VIEWCMPDPK&choice='+str(mol.cmpd_id.pk)
        context["method"] = request.GET["method"]
    context["ll_mols"] = ll_mols
    return render(request, 'WONKA/SummariseSTASH.html', context)  


def NewComplex(request, target_id):
    """Trial view - first draft of WONKA"""
    host = request.get_host()
    target = get_object_or_404(Target, pk=target_id)
    target.mol = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)[0].prot_id.code
    targets = Target.objects.exclude(title="DUMMY")
    my_mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)
    # 
    othertargs = Target.objects.filter(title="PHIPA")
    for targ in othertargs:
      other_mols = Molecule.objects.filter(prot_id__target_id=targ).exclude(prot_id__code__startswith=targ.title)
      for mol in other_mols:
        mol.code = mol.prot_id.code
        mol.src = "http://" + host + '/Viewer/loader/?function=2DMOL&choice=' + urllib.quote(mol.cmpd_id.smiles, "")
      targ.mols = other_mols
    for mol in my_mols:
        mol.code = mol.prot_id.code
        mol.src = "http://" + host + '/Viewer/loader/?function=2DMOL&choice=' + mol.cmpd_id.smiles
    context = {"target": target, "targets": targets, "mols": my_mols, "othertargs": othertargs}
    return render(request, 'WONKA/NewComplex.html', context)


def get_interactions(request):
    """Function to get the interactions for a molecule"""
    dist_dict = {"SingleAtomAcceptor_SingleAtomDonor": {"dist": 4.0}, # H-bonding
                    "SingleAtomAcceptor_WeakDonor": {"dist": 3.0}, # Weak H-bond
                    "Halogen_SingleAtomAcceptor": {"dist": 4.0}, # Halogen bonding
                    "AcidicGroup_BasicGroup": {"dist": 4.0}, # Acid-base
                    "Arom5_Arom6": {"dist": 5.5},"Arom6_Arom6": {"dist": 5.5},"Arom5_Arom5": {"dist": 5.5},# Aromatic-aromatic interactions
                    "Arom6_Carbonyl": {"dist": 4.5}, "Arom5_Carbonyl": {"dist": 4.5},# Carbonyl-aromatic interactions - CARBONLY from PROTEIN ONLY!!!!
                    "Hydrophobe_Hydrophobe": {"dist": 4.5}}#Hydrophobic interactions
    mol_pk = request.GET['obs_id']
    my_dist = request.GET['dist']
    host = request.get_host()
    mol = Molecule.objects.get(pk=mol_pk)
    out_l = []
    prot = mol.prot_id
    # Get the interactions
    interactions = ProbeBit.objects.filter(prot_id=prot, mol_id=mol, dist__lte=my_dist)
    i = -1
    for my_int in interactions:
        if my_int.type not in dist_dict:
            continue
        if my_int.dist > dist_dict[my_int.type]["dist"]:
            continue
        print "HERE"
        i += 1
        out_l.append({})
        f = my_int.probe_source_id
        out_l[i]["url_1"] = "http://"+host+'/WONKA/show_point/?x_com='+str(f.x_com)+'&y_com='+str(f.y_com)+'&z_com='+str(f.z_com)
        f = my_int.probe_dest_id
        out_l[i]["url_2"] = "http://"+host+'/WONKA/show_point/?x_com='+str(f.x_com)+'&y_com='+str(f.y_com)+'&z_com='+str(f.z_com)
        out_l[i]["dist"] = my_int.dist
        out_l[i]["type"] = my_int.type
        out_l[i]["angle_1"] = my_int.angle_1
        out_l[i]["angle_2"] = my_int.angle_2
    return HttpResponse(json.dumps(out_l))


def show_point(request):
    """Function to show a points position on click"""
    x_com = request.GET['x_com']
    y_com = request.GET['y_com']
    z_com = request.GET['z_com']
    out = render_pharma_pdb(1, 1.0, 1.0, "SingleAtomAcceptor", float(x_com), float(y_com), float(z_com))
    print out
    return HttpResponse(out)


def get_grob(request, prot_id):
    """Function to get X-ray density into the viewer"""
    p = Protein.objects.get(pk=prot_id)
    try:
        my_val = open(str(p.cif_info.file))
    except IOError:
        my_val = str(p.cif_info)
    return HttpResponse(my_val)


def get_pot_ph4s(request):
    """Function to find the five closest potential ph4 differences for this compound"""
    # First get the molecule
    mol_pk = request.GET['mol_pk']
    mol = Molecule.objects.get(pk=mol_pk)
    # Get the lam
    lam = request.GET['lam']
    # Next get the clusters
    kcs = KeyCluster.objects.filter(lam=lam, target_id=mol.prot_id.target_id)
    mol_cs = KeyCluster.objects.filter(pharma_id__in=PharmaPoint.objects.filter(mol_id=mol)).values_list("pk", flat=True)
    print mol_cs
    kcs = kcs.exclude(pk__in=mol_cs)
    ph4s = kcs.filter(type="PharmaPoint object").exclude(function__in=["MagicSix", "WeakDonor"]).exclude(function__contains="Attach")
    # Finally find the five most populous ones not in this
    key_ph4s = ph4s.order_by("size").reverse()[:5]
    return HttpResponse(json.dumps({"my_pks": [x.pk for x in key_ph4s]}))


def ScoreMol(request):
    """Function to score a given mols points and place them in %ile lists
    Takes a molecule and dict of points and scores"""
    target_id = request.GET['target']
    pdb_code = request.GET['code']
    mol = Molecule.objects.filter(prot_id__code=pdb_code)[0]
    opt = request.GET['opt']
    # Get the target
    target = Target.objects.get(pk=target_id)
    # Get the mols
    tot_mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__contains=target.title).exclude(pk=mol.pk)
    if opt in ["ph4", "out"]:
        # Do this for PH4 points
        # Get the scores for all the other molecules
        scores = sorted(PharmaScore.objects.filter(pharma_id__mol_id__in=tot_mols).values_list("score", flat=True))
        # Get the scores for my molecule
        ph4_scores = PharmaScore.objects.filter(pharma_id__mol_id=mol)
        len_scores = len(scores)
    pharma_list = []
    if opt == "ph4":
        for i, p in enumerate(ph4_scores):
            # Now find the index of where it sits
            ind = bisect.bisect_left(scores, p.score)
            # Now calc the %ile
            percentile = float(ind + 1) / float(len_scores)
            # Put it in the Bfactor column
            pharma_list.append({"type": p.pharma_id.smiles, "x_com": p.pharma_id.x_com, "y_com": p.pharma_id.y_com, "z_com": p.pharma_id.z_com, "score": percentile})
            #out_f.write(render_pharma_pdb(i+1, percentile, percentile, p.pharma_id.smiles, p.pharma_id.x_com, p.pharma_id.y_com, p.pharma_id.z_com))
        return HttpResponse(json.dumps(pharma_list, sort_keys=False))
    elif opt == "out":
        out_f = StringIO.StringIO()
        for i, p in enumerate(ph4_scores):
            # Now find the index of where it sits
            ind = bisect.bisect_left(scores, p.score)
            # Now calc the %ile
            percentile = float(ind + 1) / float(len_scores)
            # Put it in the Bfactor column
            out_f.write(render_pharma_pdb(i + 1, percentile, percentile,
                                          p.pharma_id.smiles, p.pharma_id.x_com,
                                          p.pharma_id.y_com, p.pharma_id.z_com))
        contents = out_f.getvalue()
        out_f.close()
        return HttpResponse(contents)
    elif opt == "frag":
        # Now return a photo-montage of the fragments - with scores
        # Get the scores for all the other molecules
        scores = sorted(FragScore.objects.filter(frag_id__mol_id__in=tot_mols).values_list("score", flat=True))
        # Get the scores for my molecule
        frag_scores = FragScore.objects.filter(frag_id__mol_id=mol)
        len_scores = len(scores)
        out_f = StringIO.StringIO()
        for i, p in enumerate(frag_scores):
            # Now find the index of where it sits
            ind = bisect.bisect_left(scores, p.score)
            # Now calc the %ile
            percentile = float(ind + 1) / float(len_scores)
            # Put it in the Bfactor column
            out_f.write(render_pharma_pdb(i + 1, percentile, percentile,
                                          "SingleAtomAcceptor", p.frag_id.x_com,
                                          p.frag_id.y_com, p.frag_id.z_com))
        contents = out_f.getvalue()
        out_f.close()
        return HttpResponse(contents)      
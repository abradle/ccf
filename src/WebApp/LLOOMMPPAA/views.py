from django.shortcuts import render, get_object_or_404
# Import the models
from IOhandle.models import Protein,Target,Project,Molecule,Compound
from MMPMaker.models import MMPDiffMap
from LLOOMMPPAA.models import GPVal, PossReact, ReactionProcess
from WebApp.settings import DATABASES
import json
from LLOOMMPPAA.analysis import make_dendrogram
from LLOOMMPPAA.functions import find_diff_div
from django.http import HttpResponse
from WONKA.models import Residue
from LLOOMMPPAA.models import ReactantLib, Reaction, ProductLib, RunAnalysis, DivList, HistBinVal
import collections, numpy
from rdkit import Chem
from django.core.exceptions import ValidationError
from LLOOMMPPAA.reactions import run_react_proc, init_react_proc
from django.shortcuts import redirect
from StringIO import StringIO
from LLOOMMPPAA.functions import div_pick, get_fps, best_from_wonka
from LLOOMMPPAA.models import LLConf
from PLIFS.models import Interaction
from rdkit.Chem import AllChem
from LLOOMMPPAA.functions import div_pick, get_fps
from LLOOMMPPAA.models import LLConf
from PLIFS.models import Interaction


def index(request):
    """View for the index page"""
    targets = Target.objects.exclude(title="DUMMY")
    projects = Project.objects.all()
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets, "projects": None}
    return render(request, 'LLOOMMPPAA/index.html', context)


def div_viewer(request, target_id):
    """View to find diverse molecules for a target from a set"""
    if 'num' in request.GET:
        tot_num = request.GET["num"]
    else:
        tot_num = 30
    context = {"target_id": target_id, "tot_num": tot_num}
    return render(request, 'LLOOMMPPAA/div_viewer.html', context)
    out_d = find_diff_div(target_id, tot_num)[0]
    host = request.get_host()
    # Now go through and make these options
    for item in out_d:
        context[item] = []
        for my_cmp in list(out_d[item]):
            my_mol = Molecule()
            my_mol.src = "http://" + host + '/Viewer/loader/?choice=' + str(my_cmp) + '&function=VIEWCMPDPK'
            context[item].append(my_mol)
    return render(request, 'LLOOMMPPAA/div_viewer.html', context)


def res_finder(request, target_id):
    """View to find the residues required for my list"""
    res_two_three_dict = {'ALA': "A",
                      'ARG': "R",
                      'ASN': "N",
                      'ASP': "D",
                      'CYS': "C", 
                      'GLU': "E",
                      'GLN': "Q",
                      'GLY': "G",
                      'HIS': "H",
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
    my_res = list(Residue.objects.filter(pk__in=list(set(Residue.objects.filter(probe__wonka_probebit_dest__target_id=target_id).order_by("pk").values_list("pk", flat=True)))).values("pk", "res_name", "res_num"))
    for res in my_res:
        res["icm_id"] = "^" + res_two_three_dict[res["res_name"]] + str(res["res_num"])
    # Now make this information JSON
    out_d = {"my_res": my_res}
    # Now return this object
    return HttpResponse(json.dumps(out_d))


def int_picker(request, target_id):
    """View to find compounds that all have a given interaction"""
    # Get a list of the appropriate interactions - these are all those actually having a intereaction
    target = Target.objects.get(pk=target_id)
    target.prot = Protein.objects.filter(target_id=target_id).exclude(pdb_info="")[0].code
    context = {"target_id": target_id, "target": target}
    return render(request, 'LLOOMMPPAA/int_picker.html', context)


def get_sdf_from_list(request):
    """ """
    my_mols = [int(x) for x in request.GET["MOLS"].split(",")]
    dj_mols = Molecule.objects.filter(pk__in=my_mols)
    out_sd = "$$$$\n".join([mol.sdf_info for mol in dj_mols])
    return HttpResponse(out_sd)


def DownloadLibs(request, job_id):
    """View to give the chemist an easy download of their stuff."""
    # Get thies process
    r = ReactionProcess.objects.get(pk=job_id)
    # Get the corresponding analysis run number
    my_cont = {
     "ra_id": [x.pk for x in RunAnalysis.objects.filter(react_proc_id=r)],
     }
    return render(request, 'LLOOMMPPAA/DownloadLibs.html', my_cont)


def FollowUpMaker(request, target_id):
    """View to make follow ups - finds all the potential addition points and
    allows for upload of SDF files to follow up on particular series"""
    # Get the potential reactions
    from LLOOMMPPAA.reactions import find_follow_ups, define_reacts
    define_reacts()
    find_follow_ups(target_id)
    poss_react = PossReact.objects.filter(mol_id__prot_id__target_id=target_id)
    prot_code = Protein.objects.filter(target_id=target_id)[0].code
    return render(request, 'LLOOMMPPAA/FollowUpMaker.html', {"poss_react": poss_react, "prot_code": prot_code})


def FollowUpFrag(request, frag_id, choice):
    """View to provide the view for a specific follow up"""
    # Get the reaction we're talking about
    p = PossReact.objects.filter(pk=frag_id)
    if choice == "2":
        # Get the frag - no stereo - because we don't need it because we're doing 3D superposition
        my_context = p[0].retained_frag_context
        r_frag = p[0].retained_frag
    elif choice == "1":
        my_context = p[0].replaced_frag_context
        r_frag = p[0].replaced_frag
    # Proteins to chose
    prot_ids = [p[0].mol_id.prot_id]
    # Libraries to chose
    r_lib = ReactantLib.objects.filter()
    # Libraries to chose
    p_lib = ProductLib.objects.filter()
    # Set running
    print "MY CONTEXT", my_context
    context = {"poss_react": p,
               "prot_ids": prot_ids,
               "r_lib": r_lib,
               "choice": choice,
               "my_context": my_context,
               "react_frag":  r_frag,
               "p_lib": p_lib,
               "mol_id": p[0].mol_id.pk,
               "react_id": p[0].react_id.pk}
    # Produces a URL to then check up on the progress
    return render(request, 'LLOOMMPPAA/FollowUpFrag.html', context)


def get_sdf_info(request, frag_id):
    """Function to get the fragment from a reaction PK"""
    # Get them
    p = PossReact.objects.get(pk=frag_id)
    # Get the mols
    mol_one = Chem.MolFromMolBlock(str(p.ret_frag_sdf_info))
    mol_one.SetProp("_Name", "FRAG_ONE")
    mol_two = Chem.MolFromMolBlock(str(p.rep_frag_sdf_info))
    mol_two.SetProp("_Name", "FRAG_TWO")
    out_s = Chem.MolToMolBlock(mol_one) +"\n$$$$\n"+ Chem.MolToMolBlock(mol_two)
    return HttpResponse(out_s)


def FollowUpReview(request, job_id):
    """Function to review the conformations and PLIFS generated by a react_proc"""
    # Get the ReactProc
    r = ReactionProcess.objects.get(pk=job_id)
    # Get the LLConfs for this one
    if "shape_dist" in request.GET:
        s_d = request.GET["shape_dist"]
    else:
        s_d = 1.0
    # Return the dict of compounds
    my_cmps = list(set(r.ll_conf.filter(clash__gte=-1.0, rmsd__gte=0.5, shape_dist__lte=s_d).values_list("mol_id__cmpd_id__smiles", flat=True)))
    out_d = {}
    for my_cmp in my_cmps:
        out_d[my_cmp] = r.ll_conf.filter(clash__gte=-1.0, rmsd__gte=0.5, shape_dist__lte=s_d, mol_id__cmpd_id__smiles=my_cmp)
    # Return a list of proteins
    prots = r.prot_id.all()
    # Return the different diverse lists -> based on analyses SMILES string. Use this to annotate the mols
    my_div_list = []
    # Return the difference maps -> to review the metrics
    my_diff_maps = [u'MACCS', u'RANDOM', u'USRCAT',
                   u'MORGAN', u'PLIFS_default', 'FOCUSSED']
    my_types = [u'PosN', u'RH6_6', u'iPropyl', u'TriFluro', u'SmallHalogen',
                u'SingleF', u'Imidazole', u'Methyl', u'SingleAtomDonor',
                u'Cyano', u'Arom6', u'ChainTwoWayAttach', u'Arom5',
                u'BasicGroup', u'WeakDonor', u'ThreeWayAttach', u'BigHalogen',
                u'AcidicGroup', 'ALL', u'Carbonyl', u'SingleAtomAcceptor',
                u'SingleCl', u'MagicSix']
    numbers = ["5", "10", "15", "20", "25", "30", "50", "100"]
    my_hist_types = [u'num_hba',
 u'num_hbd',
 u'num_rot_bonds',
 u'num_rings',
 u'num_morg',
 u'num_maccs',
 u'num_chir']
    my_cont = {"cmpd_dict": out_d,
               "numbers": numbers,
               "prots": prots,
               "p_code": r.mol_id.prot_id.code,
               "m_pk": r.mol_id.pk,
               "my_hist_types": my_hist_types,
               "types": my_types,
               "methods": my_diff_maps,
               "ra_id": [x.pk for x in RunAnalysis.objects.filter(react_proc_id=r)],
               "target_id": r.mol_id.prot_id.target_id.pk,
               }
    return render(request, 'LLOOMMPPAA/FollowUpReview.html', my_cont)


def MakeDendrogram(request, ra_id):
    """View to return a dendrogram for an analysis run"""
    react_anal_id = RunAnalysis.objects.get(pk=ra_id)
    return HttpResponse(make_dendrogram(react_anal_id), content_type="image/png")


def GetMols(request, ra_id):
    """View to return a list of mols for a given analysis and method"""
    # Get the RA
    anal_id = RunAnalysis.objects.get(pk=ra_id)
    # Get the method
    method_id = request.GET["METHOD"]
    div_id = request.GET["DIV_ID"]
    #num_samps = request.GET["NUM_SAMPS"]
    # Get the mols that relate to this method
    # Use the maxmin picker to pick the mols
    #d_l = DivList.objects.filter(run_id=ra_id, num_samps=num_samps, method_id=method_id)
    ntopick = 10
    target_id = anal_id.react_proc_id.mol_id.prot_id.target_id.pk
    mols = LLConf.objects.filter(clash__gte=anal_id.clash, rmsd__gte=anal_id.rmsd, shape_dist__lte=anal_id.shape_dist, reactionprocess=anal_id.react_proc_id).values_list("mol_id__cmpd_id__pk", "mol_id__pk") 
    len(mols)
    # Now get the FPS
    tot_bits = None
    tot_bits = sorted(list(set(Interaction.objects.filter(scheme_id__scheme_name=scheme, target_id=target_id, mol_id__cmpd_id__in=[x[0] for x in mols]).values_list("interaction_id__probe_dest_id__pk", flat=True))))
    # Fill the vector
    vect_l, opt = get_fps(method_id, [x[0] for x in mols], target_id, tot_bits)
    print vect_l
    div_picks = div_pick(vect_l, ntopick, opt)
    # Now define the list to make
    out_list = []
    for i, x in enumerate(div_picks):
        out_list.append((mols[x][0], i, mols[x][1], method_id))
    return HttpResponse(json.dumps({"MOLS": out_list, "DIV_ID": div_id}))


def ring_test(smiles):
    """Function to test if a mol has a ring"""
    # Get the mol
    rdmol = Chem.MolFromSmiles(str(smiles))
    if rdmol:
        return rdmol.GetRingInfo().NumRings() > 0
    else:
        return False


def find_mols(request, ra_id):
    """Function to get the molecules for a given analysis"""
    # Get the method
    method_id = request.GET["METHOD"]
    #num_samps = request.GET["NUM_SAMPS"]
    if "FILTER" in request.GET:
        my_filter = request.GET["FILTER"]
    else:
        my_filter = None
    # Get the mols that relate to this method
    # Use the maxmin picker to pick the mols
    #d_l = DivList.objects.filter(run_id=ra_id, num_samps=num_samps, method_id=method_id)
    if "NO_TO_PICK" in request.GET:
        ntopick = int(request.GET["NO_TO_PICK"])
    else:
        ntopick = 30
    div_picks, mols = pick_this_div(method_id, ntopick, my_filter)
    return [mols[x][1] for x in div_picks]


def get_sdf_file(ra_id, method_id="PLIFS_DEFAILT", ntopick=30, my_filter=None):
    """Function to get an SDF file for an RA ID"""
    div_picks, mols = pick_this_div(ra_id, method_id, ntopick, my_filter)
    my_mols = Molecule.objects.filter(pk__in=[mols[x][1] for x in div_picks])
    out_sd = Chem.SDWriter("OUTPUT."+str(ra_id)+".sdf")
    for m in my_mols:
        out_sd.write(Chem.MolFromMolBlock(str(m.sdf_info)))


def pick_this_div(ra_id, method_id, ntopick=30, my_filter=None):
    """Function to get the indices to select molecules - returns the Compound PKs of the selected molecules""" 
    # Get the RA
    anal_id = RunAnalysis.objects.get(pk=ra_id)
    target_id = anal_id.react_proc_id.mol_id.prot_id.target_id.pk
    if method_id == "FOCUSSED":
        my_picks = best_from_wonka(target_id, anal_id)
        # Now sort on heavy atoms per int
        my_picks.sort(key=lambda x: x[4], reverse=False)
        # Don't pick the same cmpd twice
        counter = 0
        done_list = []
        div_picks = []
        for i, x in enumerate(my_picks):
            mol_id = x[5]
            mol = Molecule.objects.get(pk=mol_id)
            if mol.cmpd_id_id in done_list:
                continue
            div_picks.append(mol_id)
            done_list.append(mol.cmpd_id_id)
            counter += 1
            if counter == ntopick:
                break
    else:
        mols = LLConf.objects.filter(clash__gte=anal_id.clash, rmsd__gte=anal_id.rmsd, shape_dist__lte=anal_id.shape_dist, reactionprocess=anal_id.react_proc_id).order_by("mol_id__cmpd_id__pk").distinct("mol_id__cmpd_id__pk").values_list("mol_id__cmpd_id__pk", "mol_id__pk")
        if my_filter == "RING":
            # Filter out all the rings 
            mols = [x for x in mols if ring_test(str(Compound.objects.get(pk=x[0]).smiles))]
        len(mols)
        # Now get the FPS
        tot_bits = None
        if method_id[:5] == "PLIFS":
            scheme = method_id.split("__")[0].split("_")[1]
            tot_bits = sorted(list(set(Interaction.objects.filter(scheme_id__scheme_name=scheme, target_id=target_id, mol_id__cmpd_id__in=[x[0] for x in mols]).values_list("interaction_id__probe_dest_id__pk", flat=True))))
        # Fill the vector
        vect_l, opt = get_fps(method_id, list(set([x[0] for x in mols])), target_id, tot_bits)
        if method_id[:5] == "PLIFS":
            if len(method_id.split("__")) > 1:
                opt = method_id.split("__")[1]
    return div_pick(vect_l, ntopick, opt), mols


def get_mol_sdf(request, ra_id):
    """View to return a list of mols for a given analysis and method"""
    #
    div_picks = find_mols(request, ra_id)
    # Now define the list to make
    sio = StringIO()
    w = Chem.SDWriter(sio)
    for mol_id in div_picks:
        m = Chem.MolFromMolBlock(str(Molecule.objects.get(pk=mol_id).sdf_info))
        w.write(m)
    w.flush()
    return HttpResponse(sio.getvalue())


def get_mol_ids(request, ra_id):
    """Function to get molids for the diverse picked molecules"""
    div_picks = find_mols(request, ra_id)
    return HttpResponse(json.dumps(div_picks))


def compare_maps(ra_id, method_id, type_id, method_comp=None, type_comp=None):
    """Function to compare maps / or just print off a given map"""
    # Get the map
    map_one = GPVal.objects.filter(my_anal_id=ra_id,
                           type_id=type_id,
                           method_id=method_id)
    if method_comp and type_comp:
        map_two = GPVal.objects.filter(my_anal_id=ra_id,
                               type_id=type_comp,
                               method_id=method_comp)
        # Now do the comparison
        for gpval in map_one:
            comp_gp = map_two.filter(gp_id=gpval.gp_id)
            if comp_gp:
                if gpval.value != 0.0 and comp_gp[0].value != 0.0:
                    gpval.out_val = gpval.value / comp_gp[0].value
                else:
                    gpval.out_val = 0.0
            else:
                gpval.out_val = gpval.value
    # Now render this data
    out_m = ""
    for my_p in map_one:
        if method_comp and type_comp:
            my_mol = Chem.MolFromPDBBlock(str(my_p.pdb_info))
            if my_p.out_val != 0.0:
                atm = my_mol.GetAtomWithIdx(0)
                atm.GetPDBResidueInfo().SetTempFactor(my_p.out_val)
                out_m += Chem.MolToPDBBlock(my_mol)
        else:
            if my_p.value != 0.0:
                out_m += my_p.pdb_info
    return out_m


def compare_hists(ra_id, method_id, type_id, method_comp=None, type_comp=None):
    """Function to compare maps / or just print off a given map
    Type_id is the value (like number of HBA)
    Method is the method we are using - e.g. PLIFS"""
    # Get the bins
    map_one = HistBinVal.objects.filter(my_anal_id=ra_id,
                           type_id=type_id,
                           method_id=method_id)
    if method_comp and type_comp:
        map_two = HistBinVal.objects.filter(my_anal_id=ra_id,
                               type_id=type_comp,
                               method_id=method_comp)
        # Now do the comparison
        for gpval in map_one:
            comp_gp = map_two.filter(bin_id__min_char=gpval.bin_id.min_char,bin_id__max_char=gpval.bin_id.max_char)
            if comp_gp:
                if gpval.value != 0.0 and comp_gp[0].value != 0.0:
                    gpval.out_val = gpval.value / comp_gp[0].value
                else:
                    gpval.out_val = 0.0
            else:
                print "NO COMPARISON"
                gpval.out_val = gpval.value
    # Now render this data
    out_l = []
    for my_p in map_one:
        out_l.append([my_p.bin_id.min_char, my_p.bin_id.max_char, my_p.out_val])
    return json.dumps(out_l)


def CompareHists(request, ra_id):
    """View to compare maps """
    if "METHOD_ID" in request.GET:
        method_id = request.GET["METHOD_ID"]
    else:
        return HttpResponse("ERROR - MUST SPECIFIY METHOD")
    if "TYPE_ID" in request.GET:
        type_id = request.GET["TYPE_ID"]
    else:
        return HttpResponse("ERROR - MUST SPECIFY TYPE")
    if "METHOD_COMP" in request.GET:
        method_comp = request.GET["METHOD_COMP"]
    else:
        method_comp = None
    if "TYPE_COMP" in request.GET:
        type_comp = request.GET["TYPE_COMP"]
    else:
        type_comp = None
    return HttpResponse(compare_hists(ra_id, method_id, type_id, method_comp, type_comp))


def CompareMaps(request, ra_id):
    """View to compare maps """
    if "METHOD_ID" in request.GET:
        method_id = request.GET["METHOD_ID"]
    else:
        return HttpResponse("ERROR - MUST SPECIFIY METHOD")
    if "TYPE_ID" in request.GET:
        type_id = request.GET["TYPE_ID"]
    else:
        return HttpResponse("ERROR - MUST SPECIFY TYPE")
    if "METHOD_COMP" in request.GET:
        method_comp = request.GET["METHOD_COMP"]
    else:
        method_comp = None
    if "TYPE_COMP" in request.GET:
        type_comp = request.GET["TYPE_COMP"]
    else:
        type_comp = None
    return HttpResponse(compare_maps(ra_id, method_id, type_id, method_comp, type_comp))


def CheckProgress(request, job_id):
    """View to check the progress of a given LLOOMMPPAA run"""
    # Get the job
    react_proc = ReactionProcess.objects.filter(pk=job_id)
    # LIST OF AVAILABLE STAGES
    tot_stages = ["ADD REACTANTS", "ADD PRODUCTS", "RUN REACTION", "GENERATE CONFS", "MAKING PROTEIN PROBES",
                  "MAKING MOLECULE PROBES", "MAKING LLOOMMPPAA PROBES", "FINDING INTERACTIONS"]
    # Define the length of this
    num_stages = len(tot_stages)
    switch = 0
    to_do_stages = []
    done_stages = []
    for stage in tot_stages:
        if stage == react_proc[0].proc_stage:
            switch = 1
            continue
        if switch == 1:
            to_do_stages.append(stage)
        if switch == 0:
            done_stages.append(stage)
    context = {"react_proc": react_proc, "tot_stages": tot_stages, "num_stages": num_stages, "done_stages": done_stages, "to_do_stages": to_do_stages}
    return render(request, 'LLOOMMPPAA/CheckProgress.html', context)


def SubmitJob(request):
    """View to submit a job of performing LLOOMMPPAA analysis on a library"""
    # Get the fragment to chose for the fragment to chose
    request_d = dict(request.POST)
    # Define a new ReactionProcess
    new_rp = ReactionProcess()
    # These need to be passed into the request object to
    # The context
    new_rp.context = request_d["context"][0]
    # The fragment to react
    new_rp.react_frag = request_d["react_frag"][0]
    # The molecule being used
    new_rp.mol_id = Molecule.objects.get(pk=request_d["mol_id"][0])
    # The name of the reaction
    new_rp.react_id = Reaction.objects.get(pk=request_d["react_id"][0])
    # Register the stage
    new_rp.proc_stage = "INIT"
    new_rp.stage_completion = 0
    if "project_name" in request.POST:
        new_rp.title = request_d["project_name"][0]
    else:
        # Now fail
        return HttpResponse("FAIL - NO PROJECT NAME")
    # Now save this reaction process
    try:
        new_rp.validate_unique()
        new_rp.save()
    except ValidationError:
        return HttpResponse("FAIL - NOT A NEW REACTION")
    # Get these two
    react = new_rp.reactant_queue.through
    prods = new_rp.product_queue.through
    # Get the library
    if "library_options" in request.POST or "product_options" in request.POST:
        # Now add the reactants or the prodcuts
        if "library_options" in request.POST:
            [new_rp.reactants_id.add(ReactantLib.objects.get(pk=int(x))) for x in request_d["library_options"]]
            new_rp.save()
            # Now add this to the process
            # Add these - do we want to do it now?????
            # The list of reactants left
            new_rp.is_reaction = True
            # The list of products left
            new_rp.save()
        else:
            new_rp.is_reaction = False
        if "product_options" in request.POST:
            # Now add this to  the process
            # The total library of products
            [new_rp.products_id.add(ProductLib.objects.get(pk=int(x))) for x in request_d["product_options"]]
            new_rp.save()
            # The list of products left
            new_rp.save()
    else:
        new_rp.delete()
        return HttpResponse("FAIL - NO PROTEIN SPECIFIED FOR CLASHES")
    # Get the protein
    if "protein_options" in request.POST:
        protein_options = request_d["protein_options"]
        prots = Protein.objects.filter(pk__in=[int(x) for x in protein_options])
        # The proteins used as a conflict
        [new_rp.prot_id.add(p) for p in prots]
        new_rp.save()
    else:
        new_rp.delete()
        return HttpResponse("FAIL - NO PROTEIN SPECIFIED FOR CLASHES")
    new_rp.stage_completion = 100
    new_rp.save()
    init_react_proc(new_rp)
    return redirect('/LLOOMMPPAA/' + str(new_rp.pk) + '/CheckProgress')


def TargetSummary(request, target_id):
    """View to assess a given targets current processing jobs and jump to analysi"""
    # Get all the jobs
    my_diff_maps = [u'MACCS', u'RANDOM', u'USRCAT',
                   u'MORGAN', u'PLIFS_default', "PLIFS_default__MAX", 'FOCUSSED']
    # Yes 
    my_types = [u'PosN', u'RH6_6', u'iPropyl', u'TriFluro', u'SmallHalogen',
                u'SingleF', u'Imidazole', u'Methyl', u'SingleAtomDonor',
                u'Cyano', u'Arom6', u'ChainTwoWayAttach', u'Arom5',
                u'BasicGroup', u'WeakDonor', u'ThreeWayAttach', u'BigHalogen',
                u'AcidicGroup', 'ALL', u'Carbonyl', u'SingleAtomAcceptor',
                u'SingleCl', u'MagicSix']
    numbers = ["5", "10", "15", "20", "25", "30", "50", "100"]
    context = {"JOBS": ReactionProcess.objects.filter(mol_id__prot_id__target_id=target_id),
               "target_id": target_id,
               "types": my_types,
               "methods": my_diff_maps,
                "numbers": numbers,
               }
    # Deliver them to the front end
    return render(request, 'LLOOMMPPAA/TargetSummary.html', context)

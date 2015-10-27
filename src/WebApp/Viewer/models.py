import matplotlib
## Required for the Django app
matplotlib.use('Agg')
from django.db import models
from IOhandle.models import Protein,Molecule,Target,ActivityPoint,Compound
from Pharmacophore.models import PharmaPoint
from MMPMaker.models import MMPDiffMap,MMPComparison,ActMapPoint,MMPFrag,MMP,ActPharmaPoint
from MMPMaker.functions import find_points
#from MakeMap import main as makemap
from math import exp,sqrt,pow
import ast,string,random,StringIO,gzip,sys,numpy,math,os
from base64 import encodestring
from rdkit.Chem import RDConfig,ChemicalFeatures,AllChem,Draw,MCS
from rdkit.Chem.Draw import SimilarityMaps
from rdkit import DataStructs,Chem
from pylab import cm, colorbar, subplot, get_cmap, ones, imshow, arange
try:
    import Image, ImageDraw, ImageFont
except ImportError:
    from PIL import Image, ImageDraw, ImageFont
import json


def get_hit_dict(option=None, maps=None, out_put=None, target_id=None, extra=None):
    """Function to get a hit dictionary from a target id"""
    print "GETTING THE HIT DICT"
    target = Target.objects.get(pk=target_id)
      # Get the molecules -> these are the hits
    mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)
    m_dict = {}
    # We want all the modifications off these hits -> i.e. MMPS from them
    for m in mols:
        mol_list = {}
        sps = SynthPoint.objects.filter(target_id=target_id, mmp_frag_id__mol_id=m)
        # Now loop through the synth points
        for s in sps:
        # Get the molecules and append this to a list -> we want the biggest matching part
            for x in s.mol_id.all():
                if x.cmpd_id.smiles in mol_list:
                    if s.mmp_frag_id.core_size < mol_list[x.cmpd_id.smiles][0]:
                        mol_list[x.cmpd_id.smiles] = [s.mmp_frag_id.core_size, s.pk]
                    else:
                        pass
                else:
                    mol_list[x.cmpd_id.smiles] = [s.mmp_frag_id.core_size, s.pk]
        if mol_list:
            m_dict[m.smiles] = [(x, mol_list[x][1]) for x in mol_list]
    # So now we have a dict with key -> smiles and value a list of compound smiles and their pks
    return json.dumps(m_dict, sort_keys=False)


def select_ph4(feature, group=None):
    """Function to """
    my_d = {"Donor": ['SingleAtomDonor'],
     "Acceptor": ['SingleAtomAcceptor'],
     "Arom": ['Arom6', 'Arom5'],
     "Arom Six": ['Arom6'],
     "Arom Five": ['Arom5'],
     "Hydrophobic": ['ThreeWayAttach', 'iPropyl', 'tButyl', 'ChainTwoWayAttach', 'RH6_6', 'RH3_3', 'RH4_4', 'RH5_5'],
     "Three ring": ['RH3_3'],
     "Four ring": ['RH4_4'],
     "Five ring": ['RH5_5'],
     "Six ring": ['RH6_6'],
     "Funky": ["Imidazole","Guanidine","Nitro2","ZnBinder5","ZnBinder1"],
     "Acid/Base": ["AcidicGroup","BasicGroup"],
     "MagicSix": ["SingleF","SingleCl","Methyl","Methoxy","TriFluro","Cyano"],
     "Halogen": ["BigHalogen","SmallHalogen"],
     "Fluoride": ["SingleF"],"Chloride":["SingleCl"],"Methyl":["Methyl"],"Methoxy":["Methoxy"],
     "TriFluoro": ["TriFluro"],"Cyano":["Cyano"]}
    if not group:
        return my_d[feature]
    else:
        if hasattr(group[0], "smiles"):
            return group.filter(smiles__in=my_d[feature])
        else:
            return group.filter(type__in=my_d[feature])


PHARMA_OPT_LIST = {"Standard": ["Donor","Acceptor","Arom","Hydrophobic","MagicSix","Halogen"],
                   "Magic Six": ["Fluoride","Chloride","Methyl","Methoxy","TriFluoro","Cyano"],
                   "Rings": ["Three ring","Four ring","Five ring","Six ring","Arom Five","Arom Six"]}


def gen_box_plot(option=None, maps=None, out_put=None, target_id=None, extra=None):
    """Function to make a square plot to show exploration for a given feature"""
    # Get the ranges of the box 
    target = Target.objects.get(pk=target_id)
    bins = None
    # Now take this in as an argument
    opt_list = PHARMA_OPT_LIST[extra]
    spheres = ["one radius", "two radius", "three radius"]
    # Define  where the pharmapoints can come from
    if out_put == "xtal":
        prot_list = list(Protein.objects.filter(target_id=target_id).exclude(code__contains=target.title))
        #orig_sele = PharmaPoint.objects.filter(mol_id__prot_id__in=prot_list)
        orig_sele = Probe.objects.filter(mol_id__prot_id__in=prot_list)

        desired_val = 20.0
    elif out_put == "oommppaa":
        prot_list = list(Protein.objects.filter(code=target.title+"ChEMBL"))
#        orig_sele = PharmaPoint.objects.filter(mol_id__prot_id__in=prot_list)
        orig_sele = PharmaPoint.objects.filter(mol_id__prot_id__in=prot_list)

        desired_val = 20.0
    elif out_put == "lloommppaa":
        #FUNCTION TO GET THE INFORMATION FOR LLOOMMPPAA
        pk_list = Molecule.objects.filter(prot_id__code=target.title+"SYNTH").values_list("pk", flat=True)
#        orig_sele = PharmaPoint.objects.filter(mol_id__pk__in=pk_list)
        orig_sele = Probe.objects.filter(mol_id__pk__in=pk_list)

        desired_val = 20.0
    elif out_put == "subset":
        pk_list = SubSet.objects.get(pk=option).mol_id.all().values_list("pk", flat=True)
 #       orig_sele = PharmaPoint.objects.filter(mol_id__pk__in=pk_list)
        orig_sele = Probe.objects.filter(mol_id__pk__in=pk_list)

        desired_val = 20.0
    elif out_put[:4] == "comp":
        if target.title == "BRD4":
            ls_bins = ["0-2","2-4","4-10"]
            bins = {"0-2": (0,2.0), "2-4":(2.0,4.0), "4-10":(4,10)}            
        else:
            ls_bins = ["0-0.5","0.5-2","2-10"]
            bins = {"0-0.5": (0,0.5), "0.5-2":(0.5,2.0), "2-10":(2,10)}
        # First get the activity reducing points
        aps = ActPharmaPoint.objects.filter(in_diff_15=True, num_diff_15__lte=3, target_id=target_id)
        pks = aps.values_list("pharma_id__pk", flat=True)
        # Now get the bins
        bin_dict = {}
        for bin in bins:
            bins[bin] = {"act": aps.filter(act_change__gt=bins[bin][0],act_change__lte=bins[bin][1]).values_list("pharma_id__pk", flat=True),
                         "inact":aps.filter(act_change__lt=-bins[bin][0],act_change__gte=-bins[bin][1]).values_list("pharma_id__pk", flat=True)}
            bin_dict[bin] = {"act": {}, "inact": {}}
        #orig_sele = PharmaPoint.objects.filter(pk__in=pks, smiles__in=select_ph4(opt_list[int(out_put[4:])]))
        orig_sele = Probe.objects.filter(pk__in=pks, type__in=select_ph4(opt_list[int(out_put[4:])]))

        print orig_sele
        desired_val = 3.0
    # Now get the coords
    my_d = ast.literal_eval(maps)
    points = [float(x) for x in my_d["my_var"]]
    # These are x,y and z limits either side
    x_max = max(points[3], points[0])
    x_min = min(points[3], points[0])
    y_max = max(points[4], points[1])
    y_min = min(points[4], points[1])
    z_max = max(points[5], points[2])
    z_min = min(points[5], points[2])
    previous_sele = []
    min_x = max_x = (x_max + x_min) / 2.0
    min_y = max_y = (y_max + y_min) / 2.0
    min_z = max_z = (z_max + z_min) / 2.0
    # Fix this to be a set distance
    x_step = 1.0#((x_max - x_min) * 0.5 / float(option))
    y_step = 1.0#((y_max - y_min) * 0.5 / float(option))
    z_step = 1.0#((z_max - z_min) * 0.5 / float(option))

    count_ph4 = {}
    sphere_dict = {}
    # These could be the averages for this position
    for sphere in spheres:
        sphere_dict[sphere] = {"desired": desired_val}
    out_list = {}
    last_val = {}
    for my_type in opt_list:
        last_val[my_type] =0
        count_ph4[my_type] = len(select_ph4(my_type, orig_sele))
        out_list[my_type] = {"desired": desired_val }#float("{0:.2f}".format(float(desired_val)*100)) }/ count_ph4[my_type]}
    # Define the "three" spheres
    res_d = []

    for sphere in spheres:
        # Work out the area
        min_x = min_x - x_step
        max_x = max_x + x_step
        min_y = min_y - y_step
        max_y = max_y + y_step
        min_z = min_z - z_step
        max_z = max_z + z_step
        # Find the number of features in this space
        current_sele = orig_sele
        current_sele = current_sele.filter(x_com__lte=max_x, x_com__gte=min_x)
        current_sele = current_sele.filter(y_com__lte=max_y, y_com__gte=min_y)
        current_sele = current_sele.filter(z_com__lte=max_z, z_com__gte=min_z)
        if bins:
          # Now place all the bin information in for this shell
            for bin in bins:
                bin_dict[bin]["act"][sphere] = float("{0:.2f}".format(float(len([x for x in current_sele if x.pk in bins[bin]["act"] ]))))
                bin_dict[bin]["inact"][sphere] = float("{0:.2f}".format(float(len([x for x in current_sele if x.pk in bins[bin]["inact"] ]))))
        else:
          # Now filter by type:
            for my_type in opt_list:
                #curr_val = float("{0:.2f}".format(float(len([x for x in current_sele if x.smiles in select_ph4(my_type)]))))
                curr_val = float("{0:.2f}".format(float(len([x for x in current_sele if x.type in select_ph4(my_type)]))))

                # Now fill the sphere dict
                sphere_dict[sphere][my_type] = curr_val-last_val[my_type]
                # Convert numbers to % of total
                out_list[my_type][sphere] = curr_val-last_val[my_type]#*100) / count_ph4[my_type]))
                last_val[my_type] = curr_val
    if out_put in ["oommppaa","lloommppaa","subset","xtal"]:
        for opt in opt_list:
            res_d.append({"desired": out_list[opt]["desired"], "feature": opt.lower(), "one radius": out_list[opt]["one radius"], "two radius": out_list[opt]["two radius"], "three radius": out_list[opt]["three radius"]})
    elif out_put == "NOT POSSIBLE":
        for opt in spheres:
            res_d.append({"desired": sphere_dict[opt]["desired"], "feature": opt, "acceptor": sphere_dict[opt]["Acceptor"], "donor": sphere_dict[opt]["Donor"], "arom": sphere_dict[opt]["Arom"]})         
    elif out_put[:4] == "comp":
        for bin in ls_bins:
            my_d = {"bin":str(bin)}
            for sphere in spheres:
                my_d["inact"+" "+sphere] =bin_dict[bin]["inact"][sphere]
                my_d["act"+" "+sphere] =bin_dict[bin]["act"][sphere]
            res_d.append(my_d)
    my_json = res_d
    return json.dumps(my_json, sort_keys=False)


def get_graph_base(option=None, maps=None, out_put=None, target_id=None, extra=None):
    """Function to return the base JSON for a graph - given some options"""
    if int(option) in [1,3]:
        my_list = [{
        "balloonText": "<b>[[title]]</b><br><span style='font-size:14px'>[[category]]: <b>[[value]]</b></span>",
        "fillAlphas": 0.8,
        "labelText": "[[value]]",
        "lineAlpha": 0.3,
        "title": "ONE",
        "type": "column",
    "color": "#000000",
        "valueField": "one radius"
    }, {
        "balloonText": "<b>[[title]]</b><br><span style='font-size:14px'>[[category]]: <b>[[value]]</b></span>",
        "fillAlphas": 0.8,
        "labelText": "[[value]]",
        "lineAlpha": 0.3,
        "title": "TWO",
        "type": "column",
    "color": "#000000",
        "valueField": "two radius"
    },{
        "balloonText": "<b>[[title]]</b><br><span style='font-size:14px'>[[category]]: <b>[[value]]</b></span>",
        "fillAlphas": 0.8,
        "labelText": "[[value]]",
        "lineAlpha": 0.3,
        "title": "THREE",
        "type": "column",
    "color": "#000000",
        "valueField": "three radius"
    },
    {
        "balloonText": "<span style='font-size:13px;'>[[title]] in [[category]]:<b>[[value]]</b></span>",
        "bullet": "round",
        "bulletBorderAlpha": 1,
        "bulletColor": "#FFFFFF",
        "useLineColorForBulletBorder": True,
        "fillAlphas": 0,
        "lineThickness": 2,
        "lineAlpha": 1,
        "bulletSize": 7,
        "title": "Desired",
        "valueField": "desired"
    }]
    elif int(option) ==2:
        my_list = [{
        "balloonText": "<b>[[title]]</b><br><span style='font-size:14px'>[[category]]: <b>[[value]]</b></span>",
        "fillAlphas": 0.8,
        "labelText": "[[value]]",
        "lineAlpha": 0.3,
        "title": "Arom",
        "type": "column",
    "color": "#000000",
        "valueField": "arom"
    },{
        "balloonText": "<b>[[title]]</b><br><span style='font-size:14px'>[[category]]: <b>[[value]]</b></span>",
        "fillAlphas": 0.8,
        "labelText": "[[value]]",
        "lineAlpha": 0.3,
        "title": "Acceptor",
        "type": "column",
    "color": "#000000",
        "valueField": "acceptor"
    }, {
        "balloonText": "<b>[[title]]</b><br><span style='font-size:14px'>[[category]]: <b>[[value]]</b></span>",
        "fillAlphas": 0.8,
        "labelText": "[[value]]",
        "lineAlpha": 0.3,
        "title": "Donor",
        "type": "column",
    "color": "#000000",
        "valueField": "donor"
    },
    {
        "balloonText": "<span style='font-size:13px;'>[[title]] in [[category]]:<b>[[value]]</b></span>",
        "bullet": "round",
        "bulletBorderAlpha": 1,
        "bulletColor": "#FFFFFF",
        "useLineColorForBulletBorder": True,
        "fillAlphas": 0,
        "lineThickness": 2,
        "lineAlpha": 1,
        "bulletSize": 7,
        "title": "Desired",
        "valueField": "desired"
    }]
    print int(option)
    if int(option) ==1:
        m_theme = "light"
    elif int(option) == 2:
        m_theme = "dark"
    elif int(option) == 3:
        m_theme = "light"
    chart_basis = {
    "type": "serial",
  "theme": m_theme,
  "handDrawn": False,
    "handDrawScatter":3,
    "legend": {
        "horizontalGap": 10,
        "maxColumns": 1,
        "position": "right",
    "useGraphSettings": True,
    "markerSize": 10
    },
    "dataProvider": {},
    "valueAxes": [{
        # DECIDE WHETHER WE WANT STACKS OR NOT
        #"stackType": "regular",
        "axisAlpha": 0.3,
        "gridAlpha": 0
    }],
    "graphs": my_list,
    "categoryField": "feature",
    "categoryAxis": {
        "gridPosition": "start",
        "axisAlpha": 0,
        "gridAlpha": 0,
        "position": "left"
    },
  "exportConfig":{
    "menuTop":"20px",
        "menuRight":"20px",
        "menuItems": [{
        "icon": '/lib/3/images/export.png',
        "format": 'png'
        }]
    }
    }
    return json.dumps(chart_basis)


def eucl_dist(mol_1, mol_2):
    """Function to find the euclidean distance between two points.
    Takes two objects with attribute x_com, y_com and z_com
    Returns a float"""
    x_squared = pow((mol_1.x_com - mol_2.x_com), 2)
    y_squared = pow((mol_1.y_com - mol_2.y_com), 2)
    z_squared = pow((mol_1.z_com - mol_2.z_com), 2)
    dist = sqrt(x_squared + y_squared + z_squared)
    return dist


class ViewHandler():
    """Class to deal with functions for viewing molecules/proteins in 
    different ways.
    Consists of a dict relating calls from the web to functions here.
    Arguments to functions are handled in views.py"""
    # These are the functions to be used here for IOhandling
    def __init__(self):
        # Routines for loading files in as PDBs
        # All these routines take an option of PDB_code(either 1, a list or all) and return a PDB file
        self.pdbroutine_list = {"CHECKPOINTS": check_points, "VIEWWATER": view_waters,
                                "VIEWPROTEIN": view_protein,
                                "VIEWALLMOLS": view_all_mols,
                                "VIEWMOL": view_mol, "CHECKMOL": check_mol,
                                "GETMMP": get_mmp, "GETMAP": get_map,
                                "MOLS": get_mol, "2DMOL": view_2dmol,
                                "MAKESIM": make_similarity_map, "VIEWMOLPK": view_mol_pk,
                                "VIEWCMPDPK": view_cmpd_pk, "VIEWALLRES": view_all_res,
                                "BOXPLOT": gen_box_plot,
                                "GETBASIS": get_graph_base, "GETHITS": get_hit_dict}


def make_sdf(my_pts, my_name="EPMAP", my_atm_num=11):
    """Function to make a line (C-C bond) for two points"""
    # Make an empty molecule
    mol = Chem.MolFromSmiles("C")
    em = AllChem.EditableMol(mol)
    em.RemoveAtom(0)
    for pt in range(len(my_pts)):
        # Firstly make all the atoms (Na or could be set by Ph4 type
        em.AddAtom(Chem.Atom(my_atm_num))
    gm = em.GetMol()
    Chem.SanitizeMol(gm)
    AllChem.EmbedMolecule(gm)
    cnf = gm.GetConformer()
    gp = AllChem.rdGeometry.Point3D()
    # Now add the points
    for i, pt in enumerate(my_pts):
        gp.x = pt.x_com
        gp.y = pt.y_com
        gp.z = pt.z_com
        cnf.SetAtomPosition(i, gp)
    rdm = cnf.GetOwningMol()
    rdm.SetProp("_Name", my_name)
    return Chem.MolToMolBlock(rdm)


def view_mol_pk(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to get the 3D info for a mol from a pk"""
    return Molecule.objects.get(pk=option).sdf_info


def view_cmpd_pk(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to get a 2D depiction based on a compound pk"""
    return view_2dmol(Compound.objects.get(pk=int(option)).smiles, extra="SIMPLE")


def make_similarity_maps(mol, weights, colorMap=cm.PiYG, scale=-1, size=(250, 250), sigma=None,coordScale=1.5, step=0.01, colors='k', contourLines=10, alpha=0.5, **kwargs):
    """Function to calculate similarity maps
    Heavily based on the similarity map function in the RDKit. A few changes
    to deal with exceptions and change rendering,
    Takes an RDKit molecule and a list of atom-based weights.
    Returns an image."""
    if mol.GetNumAtoms() < 2:
        raise ValueError("too few atoms")
    fig = Draw.MolToMPL(mol, coordScale=coordScale, size=size, **kwargs)
    if sigma is None:
        if mol.GetNumBonds() > 0:
            bond = mol.GetBondWithIdx(0)
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            sigma = 0.3 * math.sqrt(sum([(mol._atomPs[idx1][i] - mol._atomPs[idx2][i]) ** 2 for i in range(2)]))
        else:
            sigma = 0.3 * math.sqrt(sum([(mol._atomPs[0][i] - mol._atomPs[1][i]) ** 2 for i in range(2)]))
        sigma = round(sigma, 2)
    x, y, z = Draw.calcAtomGaussians(mol, sigma, weights=weights, step=step)
    # scaling
    if scale <= 0.0:
        maxScale = max(math.fabs(numpy.min(z)), math.fabs(numpy.max(z)))
    else:
        maxScale = scale
    # coloring
    cax = fig.axes[0].imshow(z, cmap=colorMap, interpolation='bilinear', origin='lower', extent=(0,1,0,1), vmin=-maxScale, vmax=maxScale)
    cbar = fig.colorbar(cax, shrink=.75, pad=.02,ticks=[-maxScale, 0, maxScale], orientation='vertical')
    cbar.ax.set_yticklabels(['', '', ''])  # contour lines
    fig.axes[0].contour(x, y, z, contourLines, colors=colors, alpha=alpha, **kwargs)
    return fig


def view_2dmol(option, maps=None, out_put=None, target_id=None, legends=None, extra=None):
    """Function to render a mol image from a smiles.
    The input (option) could be 1) a list of smiles 2) a smiles 3) pdb_code
    Returns a molecule image as data"""
    option = str(option)
    print option
    try:
        option = ast.literal_eval(option)
    except:
        pass
    if type(option) is list:
        mols = [Chem.MolFromSmiles(str(x)) for x in option]
        p = Chem.MolFromSmarts(MCS.FindMCS(mols).smarts)
        AllChem.Compute2DCoords(p)
        [AllChem.GenerateDepictionMatching2DStructure(x, p) for x in mols]
        image = Draw.MolsToGridImage(mols, 2, legends=legends)
    # If it's a PDB code
    elif Chem.MolFromSmiles(str(option)) is None and type(option) is str:
        mol = Chem.MolFromSmiles(str(Molecule.objects.filter(prot_id__code=option)[0].cmpd_id.smiles))
        AllChem.GenerateDepictionMatching3DStructure(mol, Chem.MolFromMolBlock(str(Molecule.objects.filter(prot_id__code=option)[0].sdf_info)))
        image = Draw.MolToImage(mol)
    # If it's a 
    elif type(option) is str:
        mol = Chem.MolFromSmiles(str(option))
        if extra == "SIMPLE":
            image = Draw.MolToImage(mol)
        elif extra is None:
            h_map = None
            image = Draw.MolToImage(mol, size=(100,100), fitImage=True, highlightMap=h_map)
        else:
            sub = Chem.MolFromSmiles(str(extra))
            h_map = get_h_map(mol, sub)
            image = Draw.MolToImage(mol, highlightMap=h_map)
        if maps is None:
            pass
        else:
            maps = float(maps)
            draw = ImageDraw.Draw(image)
            dim = (20, 0) + (20, image.size[1])
            draw.line(dim, fill=(255-int(255*maps),int(255*maps),0), width=10)
    else:
        print "NOT VALID TYPE"
        return "NOT VALID TYPE"
    output = StringIO.StringIO()
    image.save(output, format="PNG")
    contents = output.getvalue()
    return contents


def view_all_res(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to view all of one type of residue for a target"""
    from WONKA.models import Residue
    res_two_three_dict = {'A': "ALA",
                      'R': "ARG",
                      'N': "ASN",
                      'D': "ASP",
                      'C': "CYS", 
                      'E': "GLU",
                      'Q': "GLN",
                      'G': "GLY",
                      'H': "HIS",
                      'I': "ILE",
                      'L': "LEU",
                      'K': "LYS",
                      'M': "MET",
                      'F': "PHE",
                      'P': "PRO",
                      'S': "SER",
                      'T': "THR",
                      'W': "TRP",
                      'Y': "TYR",
                      'V': "VAL"}
    r = Residue.objects.filter(prot_id__code=option, res_name=res_two_three_dict[extra[0]], res_num=extra[1:])
    if r:
        return r[0].pdb_info
    else:
        return ""


def view_protein(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to return the PDB file for a protein given it's code
    Takes a code for a Protein object
    Returns an open file handle"""
    my_protein = Protein.objects.get(code=option)
    try:
        my_val = open(str(my_protein.pdb_info.file))
    except IOError:
        my_val = str(my_protein.pdb_info)
    return my_val


def view_waters(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to get the Waters for a protein"""
    from WONKA.models import Water
    return make_sdf(Water.objects.filter(prot_id__code=option),my_name="EPMAP",my_atm_num=11)


def view_all_mols(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to show the SD block of all the molecules for a target
    Takes a Target as input
    Returns an SD block of all the molecules for that target.
    """
    my_mols = Molecule.objects.filter(prot_id__target_id=int(option))
    new_mol = ""
    for mol in my_mols:
        new_mol += (str(mol.sdf_info))+"\n\n$$$$\n"
    return new_mol


def get_mmp(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to find the MMP comparison
    Takes a pk of the MMP comparison
    Returns the SD block of that data"""
    mmp = MMPComparison.objects.get(pk=option)
    return mmp.sdf_info


def get_map(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to show the information for a given map
    Takes the pk for the map
    Returns the PDB information for that map as a GZIP"""
    mymap = MMPDiffMap.objects.get(pk=option)
    out = StringIO.StringIO()
    f = gzip.GzipFile(fileobj=out, mode='w')
    f.write(str(mymap.pdb_info))
    f.close()
    contents = out.getvalue()
    return contents


def check_mol(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to check whether an input is either a valid smiles or a valid 
    protein code
    Takes a string and a Target
    Returns an answer to be used by jquery"""
    my_mols = Molecule.objects.filter(prot_id__code__icontains=option)
    target = Target.objects.get(pk=target_id)
    if len(my_mols) == 0:
        tmpmol = Chem.MolFromSmiles(str(option))
        if tmpmol is None:
            return "None molecule"
        # Now do a similarity search on this against all the molecules
        cmps = [Chem.MolFromSmiles(str(x)) for x in Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title).exclude(cmpd_id__smiles="DUMMY").values_list("cmpd_id__smiles", flat=True)]
        fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024) for x in cmps]
        sims = DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprintAsBitVect(tmpmol, 2, nBits=1024), fps)
        ind = max(enumerate(sims), key=lambda x: x[1])[0]
        mycmp = cmps[ind]
        my_mols = Molecule.objects.filter(cmpd_id__smiles=Chem.MolToSmiles(mycmp, isomericSmiles=True))
    # Now return the appropriate PDBcode
    return my_mols[0].prot_id.code


def check_points(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to check how many points correspond to a selection
    Takes the min and max activity change and the min and max number of pharmacophore point differences in option
    Returns the number of points"""
    spl_opt = option.split(",")
    points = ActPharmaPoint.objects.filter(target_id=target_id, in_diff_15=True, num_diff_15__gte=spl_opt[0], num_diff_15__lte=spl_opt[1], act_change__gte=spl_opt[2], act_change__lte=spl_opt[3])
    return len(points)


def view_mol(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to render the 3D coordinates of a molecule
    Takes a PDB code as input
    Returns an SD block"""
    my_mols = Molecule.objects.filter(prot_id__code=option)
    new_mol = ""
    for mol in my_mols:
        new_mol += (str(mol.sdf_info)) + "\n\n$$$$\n"
    return new_mol


def make2dcanv(mmpcomps):
    """Function to make a 2D canv from MMPComps"""
    rdmols = []
    acts = []
    subs = []
    smsubs = []
    donemols = []
    for m in mmpcomps:
        # Ensure that this is not just the same comparison
        mol1 = Chem.MolFromMolBlock((str(m.xtal_mol.sdf_info)))
        mol2 = Chem.MolFromMolBlock((str(m.chembl_mol.sdf_info)))
        if [m.xtal_mol.cmpd_id.pk,m.chembl_mol.pk] in donemols or [m.chembl_mol.cmpd_id.pk,m.xtal_mol.pk] in donemols:
            # Don't do the same comparison twice
            continue
        else:
            donemols.append([m.xtal_mol.cmpd_id.pk,m.chembl_mol.cmpd_id.pk])
        # Set the molecule name for the 3D display
        acts.append(render_act(m.xtal_act))
        acts.append(render_act(m.chembl_act))
        # Generate the two-d depictions after canonicalising the smiles
        mol1 = Chem.MolFromSmiles(Chem.MolToSmiles(mol1, isomericSmiles=True))
        mol2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol2, isomericSmiles=True))
        smp = MCS.FindMCS([mol1,mol2],completeRingsOnly=True, ringMatchesRingOnly=True,timeout=0.5).smarts
        p = Chem.MolFromSmarts(smp)
        subs.append(p)
        smsubs.append(smp)
        AllChem.Compute2DCoords(p)
        AllChem.GenerateDepictionMatching2DStructure(mol1,p, acceptFailure=True)
        AllChem.GenerateDepictionMatching2DStructure(mol2,p, acceptFailure=True)
        rdmols.extend([mol1, mol2])
    # So now we have the mols in a list with actvity information in a list
    # Order this list of molecules based on scaffold (p)
    # Get a list of the indices of rdmols to rearrange
    myinds = sorted(range(len(smsubs)),key=lambda x:smsubs[x])
    nmols = []
    nacts = []
    nsubs = []
    # Now rearrange everthing to suit
    for ind_m in myinds:
        nmols.extend([rdmols[ind_m*2],rdmols[ind_m*2+1]])
        nacts.extend([acts[ind_m*2],acts[ind_m*2+1]])
        nsubs.append(subs[ind_m])
    image = draw_acts(nmols,nacts,nsubs)
    output = StringIO.StringIO()
    image.save(output, format="PNG")
    contents = output.getvalue()
    return contents


def makesdinf(mmpcomps):
    """Function to make 3D molecule from MMPComps"""
    new_mol = ""
    i =0
    for mol in mmpcomps:
        i +=1
        rd_mol = Chem.MolFromMolBlock((str(mol.xtal_mol.sdf_info)))
        rd_mol.SetProp("_Name","MOL"+str(i))
        new_mol += Chem.MolToMolBlock(rd_mol)+"\n\n$$$$\n"
        i +=1
        rd_mol = Chem.MolFromMolBlock((str(mol.chembl_mol.sdf_info)))
        rd_mol.SetProp("_Name","MOL"+str(i))
        new_mol += Chem.MolToMolBlock(rd_mol)+"\n\n$$$$\n"
    return new_mol


def find_actmapoints(option, maps):
    """Function to return a series of activity map points from a map and a JSON dict"""
    my_d = ast.literal_eval(option)
    if "my_res" in my_d:
        maps = MMPDiffMap.objects.filter(pk__in=[int(x) for x in maps.split(",")])
        act = maps[0]
        inact = maps[1]
        shape = maps[2]
        mols = my_d["my_res"].split("|")
        mol_d = {}
        # Now make a dictionary of the inds for the maps
        for mol in mols:
            if mol.split("_")[1].split(".//")[0] not in ["act","inact","shape"]:
                continue
            mol_d[mol.split("_")[1].split(".//")[0]] = get_ls_inds(mol.split("_")[1].split(".//")[1])
        my_points = []
        if "act" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=act, point_id__in=mol_d["act"]))
        if "inact" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=inact, point_id__in=mol_d["inact"]))
        if "shape" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=shape, point_id__in=mol_d["shape"]))
        # Find things with a CofM within these ranges
        return my_points
    elif "my_var" in my_d:
        points = [float(x) for x in my_d["my_var"]]
        return points[3], points[0], points[4], points[1], points[5], points[2]


def get_ls_inds(my_ls):
    """Function to find the indexes for points based on 1:4,5,7:9 type notation. Returns a list"""
    ls = my_ls.split(",")
    me = [range(int(x.split(":")[0]), int(x.split(":")[1]) + 1) for x in ls if ":" in x]
    import itertools
    me_out = list(itertools.chain(*me))
    me_out.extend([int(x) for x in ls if ":" not in x])
    return me_out


def get_mol(option, maps, out_put, target_id=None, extra=None):
    """Function to get the molecules in an MMP comparison
    Option is a dict => values to find
    Maps is a list of map pks to cross reference back"""
    # First get my maps
    my_points = find_actmapoints(option, maps)
    # Now get all the mols that are associated to this
    mmpcomps = MMPComparison.objects.filter(actmappoint__in=my_points).distinct()
    #Now make the molecule SD information
    # Now make the image to look at
    if out_put == "images":
        return make2dcanv(mmpcomps)
    elif out_put == "sds":
        return makesdinf(mmpcomps)


def find_act_map_points(option, maps=None, out_put=None, target_id=None, extra=None):
    """Function to return a series of activity map points from a map and a JSON dict
    Option is a json dict of points and maps is the indices of maps to use
    Returns the points that apply."""
    my_d = ast.literal_eval(option)
    print my_d
    if "my_res" in my_d:
        maps = MMPDiffMap.objects.filter(pk__in=[int(x) for x in maps.split(",")])
        act = maps[0]
        inact = maps[1]
        shape = maps[2]
        mols = my_d["my_res"].split("|")
        mol_d = {}
        # Now make a dictionary of the inds for the maps
        for mol in mols:
            if mol.split("_")[1].split(".//")[0] not in ["act", "inact", "shape"]:
                continue
            mol_d[mol.split("_")[1].split(".//")[0]] = get_list_inds(mol.split("_")[1].split(".//")[1])
        my_points = []
        if "act" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=act, point_id__in=mol_d["act"]))
        if "inact" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=inact, point_id__in=mol_d["inact"]))
        if "shape" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=shape, point_id__in=mol_d["shape"]))
        # Find things with a CofM within these ranges
        return my_points
    elif "my_var" in my_d:
        points = [float(x) for x in my_d["my_var"]]
        return points[3], points[0], points[4], points[1], points[5], points[2]


def render_act(act):
    """Function to render activity data
    Takes a float of activity data and the operator for that data 
    Returns the data rendered to 2DP with activity associated"""
    # Deal with inactive data
    if act.operator == "<":
        return act.source + " Inactive: under {0:.2f}".format(act.activity)
    else:
        return act.source + " Activity: {0:.2f}".format(act.activity)


def get_h_map(mol, sub):
    """Function to create and return a highlight map
    Takes two rdmols - one for the molecule and the second for the substruct"""
    highlightAtoms = mol.GetSubstructMatch(sub)
    highlightMap = {}
    for idx in highlightAtoms:
        highlightMap[idx] = (0.7, 0.7, 0.7)
    return highlightMap


def draw_acts(mols, legends, subs, molsPerRow=2, subImgSize=(200, 200), **kwargs):
    """Function to draw a grid of activity data. 
    Based very strongly on the RDKit Draw.MolsToGridImage()
    Takes a list of RDKit molecules, assoicated legends, the associated 
    substituions, the numnber of mols per row and the size
    Returns a PIL image"""
    try:
        import Image, ImageDraw
    except ImportError:
        from PIL import Image, ImageDraw
    nRows = len(mols) // molsPerRow
    # Add the extra one if they are divisble
    if len(mols) % molsPerRow:
        nRows += 1
    res = Image.new("RGB", (molsPerRow * subImgSize[1], nRows * subImgSize[1]), (255, 255, 255))
    for i, mol in enumerate(mols):
        h_map = get_h_map(mol, subs[i // 2])
        # Define the row and the column
        row = i // molsPerRow
        col = i % molsPerRow
        # Now make hte molecular image
        molimage = Draw.MolToImage(mol, subImgSize, legend=legends[i], highlightMap=h_map)
        # Now paste this image
        res.paste(molimage, (col * subImgSize[0], row * subImgSize[1]))
    draw = ImageDraw.Draw(res)
# draw.line((res.size[0]/2, 0,res.size[0]/2,res.size[1]), width=2,fill=(0,0,0))
    for i in range(len(mols) / 2):
        if i == 0:
            draw.line((0, 0, res.size[0], 0), width=2, fill=(0, 0, 0))
        else:
            draw.line((0, res.size[1] * i / nRows, res.size[0], res.size[1] * i / nRows ), width=2, fill=(0, 0, 0))
    draw.line((0, res.size[1], res.size[0], res.size[1]), width=2, fill=(0, 0, 0))
    return res


def get_empty_image(error_message="Invalid Molecule"):
    """Function to return an empyty image if molecule rendering fails
    Takes an error message
    Returns an empty image with that message."""
    try:
        import Image, ImageDraw, ImageFont
    except ImportError:
        from PIL import Image, ImageDraw, ImageFont
    out_map = Image.new("RGB", (240, 240), (255, 255, 255))
    draw = ImageDraw.Draw(out_map)
    draw.text((out_map.size[0] / 2, out_map.size[1] / 2), error_message, (0, 0, 0))
    output = StringIO.StringIO()
    out_map.save(output, format="PNG")
    contents = output.getvalue()
    return contents


def make_similarity_map(option, maps=1, out_put=None, target_id=None, extra=None):
    """Function to make a sim map for a given PDB code relating to a molecule
    Takes option as the code of the molecules protein, choice indicates the 
    option for the map, target_id is the Target
    Returns a png image  as data"""
    # Find the molecule from the option i
    dj_mols = Molecule.objects.filter(prot_id__code=option)
    if len(dj_mols) != 0:
        dj_mol = dj_mols[0]
    else:
        # Assume it is a smiles string
        tmpmol = Chem.MolFromSmiles(str(option))
        if tmpmol is None:
            # If it is not then return an empty image
            #Except if the values are all nought!
            return get_empty_image()
        # Now do a similarity search on this against all the molecules
        cmps = [Chem.MolFromSmiles(str(x)) for x in Molecule.objects.exclude(cmpd_id__smiles="DUMMY").filter(prot_id__target_id=target_id).values_list("cmpd_id__smiles", flat=True)]
        if len(cmps) == 0:
            return get_empty_image()
        fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024) for x in cmps]
        sims = DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprintAsBitVect(tmpmol, 2, nBits=1024), fps)
        ind = max(enumerate(sims), key=lambda x: x[1])[0]
        mycmp = cmps[ind]
        dj_mols = Molecule.objects.filter(prot_id__target_id=target_id).filter(cmpd_id__smiles=Chem.MolToSmiles(mycmp, isomericSmiles=True))
        dj_mol = dj_mols[0]
    # Locate the target from this
    target_id = dj_mol.prot_id.target_id
    # Now iterate over the atoms in this molecule. Everyone starts of with a score of 100. Everyone is in a pharmacophore group. If it contributes to a positive pharmacophore group you add 10. If it contributes to a negative pharmacophore group you minus 10.
    my_mol = Chem.MolFromMolBlock(str(dj_mol.sdf_info))
    atms = my_mol.GetAtoms()
    weights = {a.GetIdx(): 0.1 for a in atms}
    # Find all the atoms matching the pharmacophore substructures for this type
    fdefName = os.path.join(os.path.join(os.path.split(sys.argv[0])[0], 'data/media'), 'BaseFeatures.fdef')
    if not os.path.isfile(fdefName):
        fdefName = os.path.join('data', 'media', 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    features = factory.GetFeaturesForMol(my_mol)
    # These maps indicate the properties to do with act, inact maps
    if int(maps) == 1:
## These are hard coded -> which is find e
        my_maps = MMPDiffMap.objects.filter(type__in=["act", "inact"], target_id=target_id).filter(activity_type=str("['IC50', 'Ki', 'TM']"))
        act_points = my_maps[0].actmappoint_set.all()
        inact_points = my_maps[1].actmappoint_set.all()
    # This map indicates pharmacophoric interests
    elif int(maps) == 3:
        my_points = PharmaPoint.objects.filter(mol_id__prot_id__target_id=target_id)
    if int(maps) == 1:
        for x in weights:
            weights[x] = 0.1
        # Now loop through the features
        for feat in features:
            fp = PharmaPoint()
            fp.x_com = feat.GetPos().x
            fp.y_com = feat.GetPos().y
            fp.z_com = feat.GetPos().z
            fp.type = feat.GetType()
            fp.ids = feat.GetAtomIds()
            # Now do the same again for the inactivity points
            for inact in inact_points:
                #
                if inact.mmp_comparsion_id.activity_change < 0.5 or inact.mmp_comparsion_id.num_change > 3:
                    continue
                pp = inact.pharma_id
                # Calculate the euclidean distance between this featur position and the other pharmapoistion
                if eucl_dist(pp, fp) > 4:
                    continue
                # If they are the same type
                if str(pp.uniq_id.smiles) == str(fp.type):
                    for x in fp.ids:
                        weights[x] -= exp(-eucl_dist(pp, fp)) * 300
    elif int(maps) == 3:
        # A map showing the pharmacophoric fit to the data -> i.e. conserved pharmacophoric features
        for x in weights:
            weights[x] = 0.1
        # Now loop through the features
        for feat in features:
            fp = PharmaPoint()
            fp.x_com = feat.GetPos().x
            fp.y_com = feat.GetPos().y
            fp.z_com = feat.GetPos().z
            fp.type = feat.GetType()
            fp.ids = feat.GetAtomIds()
            # Go through the activity points
            # Now loop through the features
            for pp in my_points:
                # Calculate the euclidean distance between this featur position and the other pharmapoistion
                if eucl_dist(pp, fp) > 4.0:
                    continue
                if str(pp.uniq_id.smiles) == str(fp.type):
                    for x in fp.ids:
                        weights[x] += exp(-eucl_dist(pp, fp)) * 100
    # Now do the explicity H's for choice 1
    if int(maps) == 1:
        # Add explicit H's
        hmy_mol = AllChem.AddHs(my_mol)
        AllChem.ConstrainedEmbed(hmy_mol, my_mol)
        my_mol = hmy_mol
        emy_mol = AllChem.EditableMol(my_mol)
        # And get this Conf
        conf = my_mol.GetConformer()
        # Loop through the explicit H's in the model to produce the extension list
        atoms = my_mol.GetAtoms()
        for atm in atoms:
            # Find the H's
            if atm.GetAtomicNum() == 1:
                # Now add this to the weight list
                ap = PharmaPoint()
                ap.id = atm.GetIdx()
                emy_mol.ReplaceAtom(ap.id, Chem.Atom(53))
                weights[ap.id] = 0.1
                ap.x_com = conf.GetAtomPosition(ap.id).x
                ap.y_com = conf.GetAtomPosition(ap.id).y
                ap.z_com = conf.GetAtomPosition(ap.id).z
                for act in act_points:
                    pp = act.pharma_id
                    weights[ap.id] += exp(-eucl_dist(pp, ap)) * 100
        my_mol = emy_mol.GetMol()
        Chem.SanitizeMol(my_mol)

    # Now compute 2D coords
    AllChem.Compute2DCoords(my_mol)
    output = StringIO.StringIO()
    try:
        if int(maps) == 1:
            out_map = make_similarity_maps(my_mol, weights, colorMap=cm.RdBu, alpha=0.00)
            out_map.savefig(output, format="PNG", bbox_inches='tight', dpi=35)
        else:
            out_map = make_similarity_maps(my_mol, weights, colorMap=cm.jet, alpha=0.00)
            out_map.savefig(output, format="PNG", bbox_inches='tight', dpi=35)
        contents = output.getvalue()
        return contents
    except ValueError:
        #Except if the values are all nought!
        return get_empty_image()


def get_list_inds(my_ls):
    """Function to find the indexes for points based on 1:4,5,7:9 type notation.
    Takes a string.
    Returns a list"""
    ls = my_ls.split(",")
    me = [range(int(x.split(":")[0]), int(x.split(":")[1]) + 1) for x in ls if ":" in x]
    import itertools
    me_out = list(itertools.chain(*me))
    me_out.extend([int(x) for x in ls if ":" not in x])
    return me_out

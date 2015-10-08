# Create your views here.
from IOhandle.models import Target, Protein, Molecule
from SUPERSTAR.models import PlifVisGrid, PlifProbeScore, PlifProbeGridScoreNew
from SUPERSTAR.functions import make_json
from PLIFS.models import PlifProbe
from django.http import HttpResponse
from django.shortcuts import render
from Viewer.models import make_similarity_maps
from rdkit import Chem
from rdkit.Chem import AllChem
import ast, math, StringIO
from pylab import cm
from WONKA.models import InternalIDLink


def ShowPlifs(request, target_id):
    """View to show PLIFS given a target id"""
    target = Target.objects.get(pk=target_id)
    if 'grid_space' in request.GET:
        grid_space = float(request.GET['grid_space'])
    else:
        grid_space = 0.5
    my_prot = Protein.objects.filter(target_id=target_id, pdb_info__isnull=False)[0]
    context = {"my_temp": my_prot.code, "target": target, "grid_space": grid_space}
    return render(request, 'SUPERSTAR/ShowPlifs.html', context)


def GetPlifs(request):
    """View to return the appropriate JSON for the target fingerprints"""
    targ = Target.objects.get(pk=request.GET['target_id'])
    grid_space = float(request.GET['grid_space'])
    my_fps = PlifVisGrid.objects.filter(target_id=targ, grid_space=grid_space)
    if len(my_fps) == 0:
        tot_list = make_json(targ)
    else:
        tot_list = str(my_fps[0].json_text)
    return HttpResponse(tot_list)


def loadmol(request):
    """Function to load a molecule"""
    # First get the mol
    my_choice = request.GET['choice'].split("_")[0]
    dj_mol = InternalIDLink.objects.filter(internal_id=my_choice)[0].mol_id
    # Get the molecule
    return HttpResponse(dj_mol.sdf_info)


def get_plif_json(request):
    """Function to be callable by David's stuff to return a JSON of ints
    Can return three error codes
    1) No TARGET asked for
    2) Target not available
    3) No data available for target
    If successful returns a Json"""
    # Get the target
    if "TARGET" in request.GET:
        my_targ = request.GET["TARGET"]
    else:
        return "NO KEY ERROR"
    # Now get the target
    dj_targ = Target.objects.filter(title=my_targ)
    if not dj_targ:
        return "NO TARGET ERROR"
    # Now get the data
    plif_vis = PlifVisGrid.objects.filter(target_id=dj_targ[0], grid_space=0.5)
    if not plif_vis:
        return "NO DATA FOR TARGET ERRROR"
    # Now get the JSON data
    tot_list = str(plif_vis.json_text)
    return HttpResponse(tot_list)


def get_sim_map(request):
    """Return the similarity map for this molecule based on scores"""
    # First get the mol
    my_choice = request.GET['choice'].split("_")[0]
    # Get the grid space
    grid_space = float(request.GET['grid_space'])
    # Get the django mol
    dj_mol = InternalIDLink.objects.filter(internal_id=my_choice)[0].mol_id
    # Get the PLIFS probes
    plif_ps = PlifProbe.objects.filter(mol_id=dj_mol)
    # Get the molecule
    rdmol = Chem.MolFromMolBlock(str(dj_mol.sdf_info))
    # Score dictionary - give everything 1
    score_d = {}
    for x in rdmol.GetAtoms():
        score_d[x.GetIdx()] = math.log10(1.0)
    for plif in plif_ps:
        atom_ids = ast.literal_eval(plif.atom_ids)
        my_score = PlifProbeGridScoreNew.objects.filter(plif_probe=plif, grid_space=grid_space)
        if my_score:
            score = my_score[0].score + 0.01
        else:
            continue
        for my_id in atom_ids:
            score_d[my_id] += math.log10(score)
    for x in rdmol.GetAtoms():
        score_d[x.GetIdx()] = min(score_d[x.GetIdx()], 1.0)
        score_d[x.GetIdx()] = max(score_d[x.GetIdx()], -1.0)
    # Now compute
    print score_d
    AllChem.Compute2DCoords(rdmol)
    output = StringIO.StringIO()
    out_map = make_similarity_maps(rdmol, score_d, colorMap=cm.RdBu, alpha=0.00)
    out_map.savefig(output, format="PNG", bbox_inches='tight', dpi=35)
    contents = output.getvalue()
    return HttpResponse(contents)
import os,sys
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "WebApp.settings")
from WebApp import settings
from LLOOMMPPAA.models import *
from IOhandle.models import Molecule, Target
from rdkit import Chem
from LLOOMMPPAA.functions import div_pick, get_fps
from PLIFS.models import *

def get_sdf_file(target_id=2, method_id="PLIFS_DEFAILT", ntopick=30, my_filter=None):
    """Function to get an SDF file for an RA ID"""
    # First get the RA_ID
    ra_id = RunAnalysis.objects.all()[0].pk
    # Now get the picks
    div_picks, mols = pick_this_div_new(ra_id, method_id, ntopick, my_filter)
    target = Target.objects.get(pk=target_id)
    my_mols = Molecule.objects.filter(pk__in=[mols[x][1] for x in div_picks])
    # Now write this out
    out_sd = Chem.SDWriter("/data/OUTPUT.sdf")
    for m in my_mols:
        out_sd.write(Chem.MolFromMolBlock(str(m.sdf_info)))

def pick_this_div_new(ra_id, method_id, ntopick=30, my_filter=None):
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
        if len(mols) < 30:
            print "USING CLASH OF -1.5 A"
            mols = LLConf.objects.filter(clash__gte=-1.5, rmsd__gte=anal_id.rmsd, shape_dist__lte=anal_id.shape_dist, reactionprocess=anal_id.react_proc_id).order_by("mol_id__cmpd_id__pk").distinct("mol_id__cmpd_id__pk").values_list("mol_id__cmpd_id__pk", "mol_id__pk")
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
get_sdf_file()
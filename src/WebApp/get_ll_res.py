from IOhandle.models import  *
from LLOOMMPPAA.models import *
from LLOOMMPPAA.functions import *
from PLIFS.models import *


def pick_div(ntopick=30, react_proc=3):
    """Return the inidces of a list of mols (Django Compounds) to pick """
    # Now get the picks -- use default values
    clash = -1.0
    rmsd = 0.5
    shape_dist = 1.0
    react_proc = ReactionProcess.objects.filter(pk=react_proc)
    mols = LLConf.objects.filter(clash__gte=clash, rmsd__gte=rmsd, shape_dist__lte=shape_dist,reactionprocess=react_proc).order_by("mol_id__cmpd_id__pk").distinct("mol_id__cmpd_id__pk").values_list("mol_id__cmpd_id__pk", "mol_id__pk")
    tot_bits = sorted(list(set(Interaction.objects.filter(target_id=react_proc.mol_id.target_id, mol_id__cmpd_id__in=[x[0] for x in mols]).values_list("interaction_id__probe_dest_id__pk", flat=True))))
    if not tot_bits:
        print "NO INTERACTIONS FOUND"
        return
    # Fill the vector
    vect_l, opt = get_fps("PLIFS", [x[0] for x in mols], react_proc.mol_id.target_id.pk, tot_bits)
    print vect_l
    # This is the PKey of the Compounds found
    picks =  div_pick(vect_l, ntopick, opt)
    print picks
    cmpds = Compound.objects.filter(pk__in=picks)
    # Now you can make an SDFile of these Compounds
    out_sd = Chem.SDWriter("OUTPUT.results.sdf")
    for cmpd in cmpds:
        out_sd.write(Chem.MolFromSmiles(str(cmpd.smiles)))

### IN YOUR CASE YOU CAN USE
#pick_div(ntopick=num_to_pick, react_proc=pk_to_use)
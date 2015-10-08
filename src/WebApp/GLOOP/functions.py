from GLOOP.models import Transformation
from MMPMaker.models import MMPFrag, MMP, Act3DFrag
from MMPMaker.functions import cansmirk
from IOhandle.models import ActivityPoint
from IOhandle.functions import centre_of_mass_from_SD_block
from django.core.exceptions import ValidationError
from rdkit import Chem
from rdkit.Chem import AllChem
import sys


def return_act_change(act_one, act_two):
    """Function to return the activity difference between"""
    # If they are not the same type -> ignore
    if  act_one.operator == "<" and act_two.operator == "<":
        # Both compoounds are inactive so ignore this one
        act_change = "INACTIVE__INACTIVE"
        diff_act = 0.0
    elif  act_one.operator == "<":
        # Check if the inactive compound has a higher measured activity
        if act_one.activity > act_two.activity:
            # This is now invalid is the inactvie comopound has a greater activity
            diff_act = 0.0
            act_change = "INACTIVE__INACTIVE"
        else:
            diff_act = act_two.activity - act_one.activity
            act_change = "INACTIVE__ACTIVE"
    elif act_two.operator == "<":
        # Check if the inactivity compound has a higher measured activity
        if act_two.activity > act_one.activity:
            # This is also invalid
            diff_act = 0.0
            act_change = "INACTIVE__INACTIVE"
        else:
            diff_act = act_two.activity - act_one.activity
            act_change = "ACTIVE__INACTIVE"
    else:
        # If its just a normal one
        diff_act = act_two.activity - act_one.activity
        act_change = "ACTIVE__ACTIVE"
    # Act change
    act_type = act_one.source + "__" + act_two.source
    return diff_act, act_type, act_change


def get_frag_change(in_sdf):
    """Function to get a molecules HeavyAtomCount and Molecular Weight"""
    m = Chem.MolFromMolBlock(in_sdf)
    return Chem.Lipinski.HeavyAtomCount(m), float(Chem.rdMolDescriptors.CalcExactMolWt(m))


def get_frag_diff(sdf_one, sdf_two):
    """Function to find the difference between two fragments"""
    fr_change_ha_one, fr_change_molwt_one = get_frag_change(sdf_one)
    fr_change_ha_two, fr_change_molwt_two = get_frag_change(sdf_two)
    return fr_change_ha_one - fr_change_ha_two, fr_change_molwt_one - fr_change_molwt_two


def get_ha_ratio(context, ha_change):
    """Functon to get the heavy atom ratio between the trasnformation and the mol"""
    mol = Chem.MolFromSmiles(str(context))
    return float(ha_change) / float(Chem.Lipinski.HeavyAtomCount(mol))


def store_transformation(before_frag_sdf, after_frag_sdf,
                         before_frag_smiles, after_frag_smiles,
                         mol_one, act_one, mol_two, act_two,
                         method_used, context):
    """Function to store transformations"""
    #### FILL IN THESE FIELDS
    new_trans = Transformation()
    # The change in activity
    new_trans.activity_change, new_trans.activity_type, new_trans.change_type = return_act_change(act_one, act_two)
    # Two foreign key relations to the molecules it relates to
    new_trans.mol_one = mol_one
    new_trans.mol_two = mol_two
    # Two foreign key relations to the activity points it relates to it relates to
    new_trans.act_one = act_one
    new_trans.act_two = act_two
    # Now assign the target
    new_trans.target_id = mol_one.prot_id.target_id
    # The SMARTS of the transformation
    new_trans.trans_smirks, new_trans.context = cansmirk(before_frag_smiles, after_frag_smiles, context)
    # The before frag smiles
    new_trans.before_smiles = before_frag_smiles
    # The before frag SDF
    new_trans.before_sdf = before_frag_sdf
    # The after frag smiles
    new_trans.after_smiles = after_frag_smiles
    # The after SDF
    new_trans.after_sdf = after_frag_sdf
    # Now get the centre of mass
    mol = Chem.MolFromMolBlock(before_frag_sdf)
    ot_mol = Chem.MolFromMolBlock(after_frag_sdf)
    if not mol:
        print "THIS IS NOT CORRECT", before_frag_sdf
        return
    if not ot_mol:
        print "THIS IS NOT CORRECT", after_frag_sdf
        return
    new_trans.x_com, new_trans.y_com, new_trans.z_com = centre_of_mass_from_SD_block(before_frag_sdf)
    new_trans.x_com_dest, new_trans.y_com_dest, new_trans.z_com_dest = centre_of_mass_from_SD_block(after_frag_sdf)
    if not new_trans.x_com:
        print "NOT A FRAG", before_frag_sdf
        return
    if not new_trans.x_com_dest:
        print "NOT A FRAG", after_frag_sdf
        return
    # The size of the fragment that replaces
    new_trans.ha_change, new_trans.mol_wt_change = get_frag_change(after_frag_sdf)
    # The difference in size between the two fragments
    new_trans.ha_ratio = get_ha_ratio(new_trans.context, new_trans.ha_change)
    # The difference in size between the two fragments
    new_trans.ha_diff, new_trans.mol_wt_diff = get_frag_diff(before_frag_sdf, after_frag_sdf)
    # How the molecules were derived - XTAL_XTAL, ACT_ACT, XTAL_ACT, ACT_XTAL all valid
    new_trans.mols_derived = method_used
    # Now store the shape difference between them
    if mol_one.sdf_info != "DUMMY" and mol_two.sdf_info != "DUMMY":
        rdmol1 = Chem.MolFromMolBlock(str(mol_one.sdf_info))
        rdmol2 = Chem.MolFromMolBlock(str(mol_two.sdf_info))
    new_trans.shape_dist = AllChem.ShapeTanimotoDist(rdmol1, rdmol2, ignoreHs=True)
    try:
        new_trans.validate_unique()
        new_trans.save()
    except ValidationError:
        pass


def get_act_points_from_mols(mol_one, mol_two, cmpd_id=None):
    """Function to get the relevant activity data point from a pair of molecule"""
    # Get Activity points from mol_one
    act_points_one = ActivityPoint.objects.filter(cmpd_id=mol_one.cmpd_id,
                                                  target_id=mol_one.prot_id.target_id,
                                                  source__in=["IC50", "Ki"]).order_by("activity").reverse()
    # Get Activity points from mol_two
    if not cmpd_id:
        act_points_two = ActivityPoint.objects.filter(cmpd_id=mol_two.cmpd_id,
                                                  target_id=mol_two.prot_id.target_id,
                                                  source__in=["IC50", "Ki"]).order_by("activity").reverse()
    else:
        act_points_two = ActivityPoint.objects.filter(cmpd_id=cmpd_id,
                                                  target_id=mol_two.prot_id.target_id,
                                                  source__in=["IC50", "Ki"]).order_by("activity").reverse()        
    # If one or the other doesn't have one
    if not act_points_one or not act_points_two:
        return None, None
    # Now check they are the same type - I want to prioritise IC50 -> IC50 

        # If the user has added multiple points. The the biggest one is taken
    act_point_one = act_points_one[0]
        # If the user has added multiple points. The the biggest one is taken
    act_point_two = act_points_two[0]
    # Now return these
    return act_point_one, act_point_two


def find_act_xtal(target_id):
    # Activity 3D frags
    act_3d_frags = Act3DFrag.objects.filter(mol_id__prot_id__target_id=target_id).iterator()
    tot = Act3DFrag.objects.filter(mol_id__prot_id__target_id=target_id).count()
    old = -1
    # Now ACT on XTAL and XTAL ACT
    for i, act_frag in enumerate(act_3d_frags):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rRegistering ACT to XTAL comps %d%% complete..." % old)
            sys.stdout.flush()
        # After frag SDF
        after_frag_sdf = act_frag.sdf_info
        # Context 
        context = act_frag.mmp_frag_id.mmp_link.context
        after_frag_smiles = act_frag.mmp_frag_id.smistore
        mol_two = act_frag.mol_id
        mol_one = act_frag.over_mol_id
        # Get the activity points
        act_one, act_two = get_act_points_from_mols(mol_one, mol_two)
        if not act_one:
            continue
        # Before frag SDF
        my_mmps = act_frag.over_mol_id.mmpfrag_set.filter(mmp_link=act_frag.mmp_frag_id.mmp_link)
        # After frag SDF
        if len(my_mmps) == 0:
            print "NO MATCH ERRRORRR!!!!"
        elif len(my_mmps) > 1:
            print "MULTIPLE OPTIONS -> LETS PICK BOTH"
        for my_mmp in my_mmps:
            before_frag_sdf = str(my_mmp.sdf_info)
            before_frag_smiles = str(my_mmp.smistore)
            # Store the frontwards o
            store_transformation(before_frag_sdf, after_frag_sdf,
                     before_frag_smiles, after_frag_smiles,
                     mol_one, act_one, mol_two, act_two,
                     "ACT_XTAL", context)
            # Now make the comparison both ways
            store_transformation(after_frag_sdf, before_frag_sdf,
                     after_frag_smiles, before_frag_smiles, 
                     mol_two, act_two, mol_one, act_one,
                     "XTAL_ACT", context )


def define_trans_groups(act_points):
    """Function to define a given group of transformations, e.g. IC50 transformations
    Takes a list of activity points"""
    # Loop through the transformations
    # 
    tg = TransGroups()
    

def find_transformations(target_id):
    """Find all the transformations for a given target"""
    # Get all the context
    mmps = MMP.objects.all().iterator()
    # Now circle through them
    tot = MMP.objects.all().count()
    old = -1
    # Now ACT on XTAL and XTAL ACT
    for i, my_mmp in enumerate(mmps):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rRegistering XTAL to XTAL comps %d%% complete..." % old)
            sys.stdout.flush()
        print my_mmp.context
        if len(my_mmp.context) < 9:
            continue
        # Now find ALL the fragments with 3D information
        mol_frags = MMPFrag.objects.filter(mmp_link=my_mmp).exclude(mol_id__sdf_info="DUMMY").exclude(mol_id__isnull=True).iterator()
        # Now loop through them to find the points
        # First crystal on crystal
        for mol_frag in mol_frags:
            for other_mol_frag in mol_frags:
                if mol_frag.pk == other_mol_frag.pk:
                    continue
                # Get the activity points
                act_one, act_two = get_act_points_from_mols(mol_frag.mol_id, other_mol_frag.mol_id)
                if not act_one:
                    continue
                # Now make the comparison
                store_transformation(mol_frag.sdf_info, other_mol_frag.sdf_info,
                         mol_frag.smistore, other_mol_frag.smistore,
                         mol_frag.mol_id, act_one, other_mol_frag.mol_id, act_two,
                         "XTAL_XTAL", my_mmp.context)
    # Now do act on act
    find_act_xtal(target_id)
    # Now do act on xtal for double and triple cuts######
    for i, my_mmp in enumerate(mmps):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rRegistering ACT to XTAL comps %d%% complete..." % old)
            sys.stdout.flush()
        act_frags = MMPFrag.objects.filter(mmp_link=my_mmp).filter(mol_id__sdf_info="DUMMY").filter(mol_id__isnull=True).iterator()
        mol_frags = MMPFrag.objects.filter(mmp_link=my_mmp).exclude(mol_id__sdf_info="DUMMY").exclude(mol_id__isnull=True).iterator()
        for mol_frag in mol_frags:
            for act_frag in act_frags:
                act_one, act_two = get_act_points_from_mols(mol_frag.mol_id, act_frag.activitypoint.cmpd_id)
                # Now make the comparison
                store_transformation(mol_frag.sdf_info, act_frag.sdf_info,
                         mol_frag.smistore, act_frag.smistore,
                         mol_frag.mol_id, act_one, act_frag.mol_id, act_two,
                         "ACT_XTAL", my_mmp.context)
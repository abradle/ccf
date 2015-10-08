from MMPMaker.models import MMP,MMPFrag,MMPComparison
from IOhandle.models import Molecule,Target,Protein,ActivityPoint
import WebApp.settings as settings
import sys
import random
import string
import os
import ast
import tempfile, pickle
import subprocess
import re
from collections import Counter
from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem, Draw, MCS, ChemicalFeatures
#from LLOOMMPPAA.models import ContextTable, CoreTable, CmpdSmisp
from LLOOMMPPAA.models import DivList, LLConf,  CmpdScore, RunAnalysis
from IOhandle.functions import add_new_comp, make_smarts_from_frag, Compound
from MMPMaker.functions import make_h_frag, find_core
from MMPMaker.helpers import adjust_arom_Ns, gen_confs_lloommppaa
from math import *
from Pharmacophore.models import PharmaPoint, find_pharmacophore_points
from MMPMaker.models import ActPharmaPoint
from OOMMPPAA.cluster import dpmeans
from Pharmacophore.models import PharmaPoint
from MMPMaker.models import ActPharmaPoint
from OOMMPPAA.cluster import dpmeans
from rdkit.Chem import Lipinski, ChemicalFeatures
from django.db import connections, connection
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from WONKA.models import InternalIDLink
from rdkit import SimDivFilters, DataStructs
import numpy, math, random
### Not on the server yet - just for LLOOMMPPAA
try:
    from PLIFS.models import Interaction, InteractionScheme
except:
    pass
try:
    from LLOOMMPPAA.sascore import calculateScore
except:
    pass
try:
    from usrcat.toolkits.rd import generate_moments
except:
    pass
from rdkit.Chem import Descriptors


def check_same(my_point, my_ints):
    """Function to check if two points are the same"""
    for my_int in my_ints:
        if my_int.function == my_point.type:
            if eucl_dist([my_point.x_com, my_point.y_com, my_point.z_com], [my_int.x_com, my_int.y_com, my_int.z_com]) < 2.0:
                return True
    return False


def best_from_wonka(target_id, react_anal_id):
    """Function to find a focussed set of compounds from a target"""
    from WONKA.models import KeyCluster
    target = Target.objects.get(pk=target_id)
    vect_type = "default"
    # First get the interactions from WONKA
    my_ints = KeyCluster.objects.filter(lam=2.0, target_id=target_id).filter(type="PharmaPoint object").exclude(function__in=["MagicSix", "WeakDonor", "Carbonyl"]).exclude(function__contains="Attach").order_by().order_by("size").reverse()
    print "FOUND ALL INTERACTIONS"
    # Now find ligands that only have these interactions
    # Get the ligands
    my_mols = react_anal_id.react_proc_id.ll_conf.filter(clash__gte=react_anal_id.clash, rmsd__gte=react_anal_id.rmsd, shape_dist__lte=react_anal_id.shape_dist).values_list("mol_id", flat=True)
    # Now iterate through them
    out_ans = []
    print "ASSESSING MOLS"
    tot = len(my_mols)
    old = -1
    for i, my_pk in enumerate(my_mols):
        # Get the complexity object
        mol = Molecule.objects.get(pk=my_pk)
        if i * 100 / tot > old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
        # Get the molecular interactions
        mol_ints = mol.plifprobe_set.all()
        matches = [x for x in mol_ints if check_same(x, my_ints)]
        # Now append a dict of 
        if matches and mol_ints:
            out_ans.append([len(mol_ints), len(matches), float(len(matches)) / float(len(mol_ints)), float(mol.cmpd_id.mol_wt) / float(len(matches)), float(mol.cmpd_id.heavy_atom_count) / float(len(matches)), mol.pk])
    return out_ans


def write_ints_pdb(mol_id, int_scheme="defaults"):
    """Function to write all interactions for a molecule to a file"""
    # First get the interactions
    my_ints = Interaction.objects.filter(mol_id=mol_id, scheme_id__scheme_name=int_scheme)
    # Next write them out
    mol = Chem.MolFromSmiles("C")
    mol_txt = ""
    for item in my_ints:
        pp = item.interaction_id.probe_source_id
        pp_2 = item.interaction_id.probe_dest_id
        # Make an editable molecule
        em = AllChem.EditableMol(mol)
        em.RemoveAtom(0)
        sv = em.AddAtom(Chem.Atom(6))
        sv_2 = em.AddAtom(Chem.Atom(6))
        em.AddBond(sv, sv_2)
        gm = em.GetMol()
        Chem.SanitizeMol(gm)
        AllChem.EmbedMolecule(gm)
        cnf = gm.GetConformer()
        sp = AllChem.rdGeometry.Point3D()
        # Now add the coords
        sp.x = pp.x_com
        sp.y = pp.y_com
        sp.z = pp.z_com
        cnf.SetAtomPosition(sv, sp)
        sp = AllChem.rdGeometry.Point3D()
        # Now add the coords
        sp.x = pp_2.x_com
        sp.y = pp_2.y_com
        sp.z = pp_2.z_com
        cnf.SetAtomPosition(sv_2, sp)
        # Now get the molecule
        out_mol = cnf.GetOwningMol()
        out_mol = Chem.MolFromPDBBlock(Chem.MolToPDBBlock(out_mol))
        # Now add this value
        mol_txt += Chem.MolToPDBBlock(out_mol) + "\n"
    out_f = open(str(mol_id.pk) + "__" + int_scheme +".pdb", "w")
    out_f.write(mol_txt)
    out_f.close()
    # Now write out all the extensions
    # Next write them out
    mol = Chem.MolFromSmiles("C")
    mol_txt = ""
    for item in my_ints:
        pp = item.interaction_id.probe_source_id
        pp_2 = item.interaction_id.probe_dest_id
        if pp.type == "SingleAtomDonor":
            pass
        elif  pp_2.type == "SingleAtomDonor":
            pp = pp_2
        else:
            continue
        # Make an editable molecule
        em = AllChem.EditableMol(mol)
        em.RemoveAtom(0)
        sv = em.AddAtom(Chem.Atom(6))
        sv_2 = em.AddAtom(Chem.Atom(6))
        em.AddBond(sv, sv_2)
        gm = em.GetMol()
        Chem.SanitizeMol(gm)
        AllChem.EmbedMolecule(gm)
        cnf = gm.GetConformer()
        sp = AllChem.rdGeometry.Point3D()
        # Now add the coords
        sp.x = pp.x_extent
        sp.y = pp.y_extent
        sp.z = pp.z_extent
        cnf.SetAtomPosition(sv, sp)
        # Now add the coords
        if pp.x_extent_2:
            sp = AllChem.rdGeometry.Point3D()
            sp.x = pp.x_extent_2
            sp.y = pp.y_extent_2
            sp.z = pp.z_extent_2
            cnf.SetAtomPosition(sv_2, sp)
        else:
            cnf.SetAtomPosition(sv, sp)
        # Now get the molecule
        out_mol = cnf.GetOwningMol()
        out_mol = Chem.MolFromPDBBlock(Chem.MolToPDBBlock(out_mol))
        # Now add this value
        mol_txt += Chem.MolToPDBBlock(out_mol) + "\n"
    out_f = open(str(mol_id.pk) + "__" + int_scheme +"EXTS.pdb", "w")
    out_f.write(mol_txt)
    out_f.close()


def find_clashes(target_id=2):
    """Helper function to take a Ligands native protein-ligand ints and find the closest contact"""
    # Get the molecules
    from MMPMaker.helpers import check_for_clashes
    targ = Target.objects.get(pk=target_id)
    my_mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__contains=targ.title)
    # Go through the mols
    clashes = []
    for mol in my_mols:
        # Get the protein coords
        prot_pos = get_coords(targ, [mol.prot_id.code])
        rdmol = Chem.MolFromMolBlock(str(mol.sdf_info))
        my_clash = check_for_clashes(rdmol, prot_pos)
        print mol.prot_id.code, my_clash
        clashes.append(my_clash)
    # So add -0.5 Angstrom on to that or take -1.0 as the worst because otherwise it's getting too much
    return max(min(clashes) - 0.5, - 1.0)


def calc_cmpd_complex():
    """Function to calculate the complexity of all compounds"""
    # Take all the compounds
    cmpds = Compound.objects.exclude(smiles="DUMMY")
    tot = len(cmpds)
    old = -1
    for i, cmpd in enumerate(cmpds):
        # Get the complexity object
        if i * 100 / tot > old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
        rdmol = Chem.MolFromSmiles(str(cmpd.smiles))
        if not rdmol:
            continue
        complexity = Cmpd_Complexity.objects.get_or_create(cmpd_id=cmpd)[0]
        # Synthetic accesibilty of a compound
        complexity.sa_score = calculateScore(rdmol)
        # The number of heavy atoms
        complexity.size = Descriptors.HeavyAtomCount(rdmol)
        # The number of on bits for MF
        complexity.num_morg = len(AllChem.GetMorganFingerprint(rdmol, 2).GetNonzeroElements())
        # The number of on bits for MACCS keys
        complexity.num_maccs = MACCSkeys.GenMACCSKeys(rdmol).GetNumOnBits()
        # The number of chiral centers
        complexity.num_chir = len(Chem.FindMolChiralCenters(rdmol))
        # The number of rings
        complexity.num_rings = AllChem.CalcNumRings(rdmol)
        # The number of HBA
        complexity.num_hba = AllChem.CalcNumLipinskiHBA(rdmol)
        # The number of HBD
        complexity.num_hbd = AllChem.CalcNumLipinskiHBD(rdmol)
        # The number of rot bonds
        complexity.num_rot_bonds = AllChem.CalcNumRotatableBonds(rdmol)
        # The number of het atoms
        complexity.num_het_atms = AllChem.CalcNumHeteroatoms(rdmol)
        # Now save this
        complexity.save()
    old = 100
    sys.stdout.write("\r%d%% complete..." % old)
    sys.stdout.flush()


def inv_man_dist(fp_1, fp_2):
    """Function to find the inverse manhattan distance between two float vectors"""
    tot_dist = 0.0
    for i, fp_bit_1 in enumerate(fp_1):
        tot_dist += math.fabs((fp_1[i] - fp_2[i]))
    return tot_dist#(1.0 / (1.0 + (1.0 / float(len(fp_1))) * tot_dist))


def bulk_inv_man_dist(fp, fps):
    """Function to find the inverse manhattan distance between one molecule
    and many other molecules"""
    out_vect = []
    for fp_vect in fps:
        out_vect.append(inv_man_dist(fp, fp_vect))
    return out_vect


def pick_random(my_list, ntopick):
    """Function to pick a series of random elements from a list"""
    ids_out = []
    fps_out = []
    i = 0
    while i < ntopick:
        my_choice = random.randrange(0, len(my_list))
        if my_choice in ids_out:
            continue
        ids_out.append(my_choice)
        fps_out.append(my_list[my_choice])
        i += 1
    return fps_out, ids_out


def bulk_shape_dist(mol, other_mols):
    """Function to calculate the shape difference between this molecule and a list of other molecules
    Takes an RDMol
    and a list of RMols
    Returns a list of values"""
    return [AllChem.ShapeProtrudeDist(mol, other_mol) for other_mol in other_mols]


def div_pick(fps, ntopick, opt=None):
    """Function to pick a diverse set of compounds based on fingerprints
    Taken from GLs IPython notebook"""
    ds = []
    mmp = SimDivFilters.MaxMinPicker()
    rand_seed = random.randint(1, len(fps))
    if opt == "USRCAT":
        def inv_man_dist(i, j, fps=fps):
            """Function to find the manhattan distance between two float vectors"""
            tot_dist = 0
            fp_1 = fps[i]
            fp_2 = fps[j]
            for i, fp_bit_1 in enumerate(fp_1):
                tot_dist += math.fabs((fp_1[i] - fp_2[i]))
            return tot_dist
        ids = mmp.LazyPick(inv_man_dist, len(fps) - 1, ntopick, seed=rand_seed)
        return ids
    elif opt == "SHAPE_DIST":
        def shape_dist(i, j, fps=fps):
            """Function to calculate the shape difference between this molecule and a list of other molecules
            Takes an RDMol
            and a list of RMols
            Returns a list of values"""
            mol = fps[i]
            other_mol = fps[j]
            return AllChem.ShapeProtrudeDist(mol, other_mol)
        ids = mmp.LazyPick(shape_dist, len(fps)-1, ntopick, seed=rand_seed)
        return ids
    elif opt == "RANDOM":
        fps_out, ids_out = pick_random(fps, ntopick)
        return ids_out
    # If this is trying to optimise the range of fingerprints
    # PICK 10000 subsets
    # Chose the one with most onbits
    elif opt == "MAX":
        print "HERE"
        best_count = 0
        best_list = []
        for i in range(10000):
            tot_list = []
            # Pick a random set of fps
            rand_fps, rand_ids = pick_random(fps, ntopick)
            for bit_s in rand_fps:
                tot_list.extend(list(bit_s.GetOnBits()))
            my_count = len(set(tot_list))
            # If this is better than the other list
            if my_count > best_count:
                best_count = my_count
                best_list = rand_ids
                print best_count
        if not best_list:
            return rand_ids
        else:
            return best_list
    elif not opt:
        def tani_dist(i, j, fps=fps):
            """Function to calculate the tanimoto dis-similarity between molecules"""
            fp_one = fps[i]
            fp_two = fps[j]
            my_t_dist = DataStructs.TanimotoSimilarity(fp_one, fp_two)
            return 1.0 - my_t_dist
        ids = mmp.LazyPick(tani_dist, len(fps) - 1, ntopick, seed=rand_seed)
        return ids


def get_counts_for_each(my_d):
    """Function to take in the dictionary and return counts for each compound set"""
    ele_set = set()
    samp_set = set()
    out_d = {}
    for test in my_d:
        out_d[test] = {}
        for num_samps in my_d[test]:
            # Now count the number of each ele
            samp_set.add(num_samps)
            tot_l = []
            # Now go through the iterations
            for my_iter in my_d[test][num_samps]:
                tot_l.extend(my_d[test][num_samps][my_iter])
            # Tot_ans is a dict key = comp, val = count
            tot_ans = dict(Counter(tot_l))
            out_d[test][num_samps] = tot_ans
            # Go through and add these to the total set
            for my_ans in tot_ans:
                ele_set.add(my_ans)
    #  Now return this dict
    return out_d, list(samp_set), list(ele_set)


def do_diff_analysis(target_id, num_iters=50, num_trial=30, clash=-1.0, rmsd=0.5, shape_dist=1.0, react_id=None, redo=False):
    """Function to do the different analyses - and bootstrap the findings"""
    from LLOOMMPPAA.analysis import find_populated_interactions
    # Get the target
    target = Target.objects.get(pk=target_id)
    print "PROCESSING ANALYSES FOR TARGET - ", target.title
    # Now do the analysis
    num_to_do = []
    ra = RunAnalysis.objects.get_or_create(clash=clash, rmsd=rmsd, shape_dist=shape_dist, num_iters=num_iters, num_trials=num_trial, react_proc_id=react_id)
    print ra
    if ra[1] == False and redo != True:
        print ra[1]
        print redo
        print "DONE ALREADY"
        return
    out_ans = find_diff_div(target_id, tot_num_list=[num_trial], num_trials=num_iters, clash=clash, rmsd=rmsd, shape_dist=shape_dist, react_id=react_id)
    # So now store these answer
    out_d, samp_set, ele_set = get_counts_for_each(out_ans)
    # Now loop through the methods
    for method in out_d:
        # The number of samples
        for num_samp in out_d[method]:
            ra = RunAnalysis.objects.get(clash=clash, rmsd=rmsd, shape_dist=shape_dist, num_iters=num_iters, num_trials=num_samp, react_proc_id=react_id)
            div_list = DivList.objects.get_or_create(method_id=method, num_samps=num_samp, react_proc_id=react_id, run_id=ra)[0]
            # The compound itself
            for cmpd in out_d[method][num_samp]:
                # Either get or create the score
                dj_cmpd = Compound.objects.get(pk=cmpd)
                exist_score = div_list.score_id.filter(cmpd_id=dj_cmpd)
                if exist_score:
                    exist_score = div_list.score_id.filter(cmpd_id=dj_cmpd)[0]
                    exist_score.score = out_d[method][num_samp][cmpd]
                    exist_score.save()
                else:
                    dj_score = CmpdScore()
                    dj_score.cmpd_id = dj_cmpd
                    dj_score.score = out_d[method][num_samp][cmpd]
                    dj_score.save()
                    # Add this score to the list
                    div_list.score_id.add(dj_score)
                    div_list.save()


def get_fps(vect_type, dj_pks, target_id, tot_bits):
    vect_l = []
    for dj_pk in dj_pks:
      # Append a Morgan fingerprint to this list
        dj_mol = Compound.objects.get(pk=dj_pk)
        rdmol = Chem.MolFromSmiles(str(dj_mol.smiles))
        conf = AllChem.EmbedMolecule(rdmol)
        if not rdmol or conf == -1:
            print "NONE MOL"
            continue
        # Now assign the various fingerprints
        if vect_type == "MORGAN":
            vect_l.append(AllChem.GetMorganFingerprint(rdmol, 2))
        elif vect_type == "MACCS":
            vect_l.append(AllChem.GetMACCSKeysFingerprint(rdmol))
        elif vect_type == "RANDOM":
            vect_l.append(AllChem.GetMACCSKeysFingerprint(rdmol))
        elif vect_type == "USRCAT":
            my_moments = generate_moments(rdmol)
            vect_l.append(my_moments[0])
        elif vect_type == "SHAPE_DIST":
            vect_l.append(rdmol)
        elif vect_type[:5] == "PLIFS":
            pbs = Interaction.objects.filter(target_id=target_id, mol_id__cmpd_id=dj_mol, scheme_id__scheme_name=vect_type.split("_")[1]).values_list("interaction_id__probe_dest_id__pk", flat=True)
            plif_v = DataStructs.ExplicitBitVect(len(tot_bits), False)
            for i, bit in enumerate(tot_bits):
                if bit in pbs:
                    plif_v.SetBit(i)
                else:
                    continue
            vect_l.append(plif_v)
    # Now use the diversity picker to produce the diverse molecules
    if vect_type == "USRCAT":
        opt = "USRCAT"
    elif vect_type == "SHAPE_DIST":
        opt = "SHAPE_DIST"
    elif vect_type == "RANDOM":
        opt = "RANDOM"
    elif len(vect_type.split("__")) > 1:
        if vect_type.split("__") == "MAX":
            opt = "MAX"
        else:
            print "ERROR UNKNOWN DIV TYPE"
            raise ValueError
    else:
        opt = None
    return vect_l, opt


def find_diff_div(target_id, opt=None, mols=None, tot_num_list=[10], num_trials=1, clash=-1.0, rmsd=0.5, shape_dist=1.0, react_id=None):
    """Function to take a set of MOLs and pick a diverse range of them"""
    # Get the target
    target = Target.objects.get(pk=target_id)
    # First get the LLOOMMPPAA mols for this target and react_id
    if not mols:
        if react_id:
            mols = list(set(LLConf.objects.filter(clash__gte=clash, rmsd__gte=rmsd, shape_dist__lte=shape_dist, reactionprocess=react_id).values_list("mol_id__cmpd_id__pk", flat=True)))
        else:
            mols = list(set(LLConf.objects.filter(clash__gte=clash, rmsd__gte=rmsd, shape_dist__lte=shape_dist).values_list("mol_id__cmpd_id__pk", flat=True)))
    # Check that it's ok for all of the selections
    print "YOU HAVE THIS MANY MOLS", len(mols)
    tot_num_list_old = tot_num_list
    tot_num_list = [x for x in tot_num_list if x <= len(mols)]
    banned_list = [x for x in tot_num_list_old if x not in tot_num_list]
    # Now print this out if it exists
    if banned_list:
        print "The following group sizes are two big: ", banned_list
    # Now loop through them
    tot_vects = []
    # Now add the SCHEMES to this
    [(tot_vects.append("PLIFS_" + int_sch.scheme_name), tot_vects.append("PLIFS_" + int_sch.scheme_name + "__MAX")) for int_sch in InteractionScheme.objects.all()]
    out_d = {}
    for vect_type in tot_vects:
        print "Finding div mols for: ", vect_type
        # Get the tot_bits
        tot_bits = None
        if vect_type[:5] == "PLIFS":
            scheme = vect_type.split("_")[1]
            tot_bits = sorted(list(set(Interaction.objects.filter(scheme_id__scheme_name=scheme, target_id=target_id, mol_id__cmpd_id__in=mols).values_list("interaction_id__probe_dest_id__pk", flat=True))))
            print tot_bits
        # Fill the vector
        vect_l, opt = get_fps(vect_type, mols, target_id, tot_bits)
        for tot_num in tot_num_list:
            if vect_type in out_d:
                out_d[vect_type][tot_num] = {}
            else:
                out_d[vect_type] = {tot_num: {}}
            for trial_num in range(num_trials):
                sys.stdout.write("\rDoing trial %d of %d..." % (trial_num, num_trials))
                sys.stdout.flush()
                div_picks = div_pick(vect_l, int(tot_num), opt)
                out_d[vect_type][tot_num][trial_num] = []
                for pick in div_picks:
                    # Now get the correct compound id - using the pick as the index
                    out_d[vect_type][tot_num][trial_num].append(mols[pick])
    # Now return again
    print "Found diverse molecules"
    return out_d


def make_map(my_pts):
    """Function to make a map from a series of points"""
    mol = Chem.MolFromSmiles("C")
    em = AllChem.EditableMol(mol)
    em.RemoveAtom(0)
    for pt in range(len(my_pts)):
        # Firstly make all the atoms (Na or could be set by Ph4 type
        em.AddAtom(Chem.Atom(11))
    gm = em.GetMol()
    Chem.SanitizeMol(gm)
    AllChem.EmbedMolecule(gm)
    cnf = gm.GetConformer()
    gp = AllChem.rdGeometry.Point3D()
    for i, pt in enumerate(my_pts):
        gp.x = pt.x_com
        gp.y = pt.y_com
        gp.z = pt.z_com
        cnf.SetAtomPosition(i, gp)
    return Chem.MolToMolBlock(cnf.GetOwningMol())


def view_2dmol(option, file_path, legends=None):
    """Function to render a mol image from a smiles
    Input is a smiles string, a path to a file to save it and a possible legends
    Returns a file name"""
    option = str(option)
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
    elif type(option) is str:
        mol = Chem.MolFromSmiles(str(option))
        image = Draw.MolToImage(mol)
    else:
        print "NOT VALID TYPE"
        return "NOT VALID TYPE"
    image.save(file_path, format="PNG")
    return os.path.split(file_path)[1]


def find_cmpds_xtalise(target_id):
    """Function to, for a given target, suggest compounds to crystalise and to
    get activity data for.
    Takes a target id
    Returns None"""
    # Find the single split contexts found in active compounds BUT not crystal compounds
    print "FINDING MMPS..."
    mmps = MMP.objects.exclude(context__contains="[*:2]").exclude(mmpfrag__mol_id__prot_id__target_id=target_id).filter(mmpfrag__act_id__target_id=target_id).distinct()
    tot_mmps = len(mmps)
    old_val = -1
    target = Target.objects.get(pk=target_id)
    # Now find the crystal structure MMP with the nearest similarity to this one
    my_mols = Molecule.objects.exclude(prot_id__code__startswith=target.title).filter(prot_id__target_id=target_id)
    len(my_mols)
    # Get all the fingerprints for xtal structures of this mol
    Mol_fps = [MACCSkeys.GenMACCSKeys(Chem.MolFromMolBlock(str(x.sdf_info))) for x in my_mols]
    # Now find the number of active compounds these relate to
    for i, m in enumerate(mmps):
        new_point = ExpPoint()
        new_point.MMP_link = m
        new_point.target_id = target
        new_point.save()
        # For those fine
        if i * 100 / tot_mmps != old_val:
            old_val = i * 100 / tot_mmps
            sys.stdout.write(str(old_val) + "% COMPLETED....\n")
            sys.stdout.flush()
        # Add the links in
        [new_point.MMPFragLink.add(mmpf) for mmpf in MMPFrag.objects.filter(mmp_link=m, act_id__target_id=target_id)]
        my_mps = MMPFrag.objects.filter(mmp_link=m, act_id__target_id=target_id, core_size__lte=6)
        new_point.NumHeavyAtomGTSix = len(my_mps)
        ### Define the f_name
        f_name = ''.join(random.sample(string.letters * 10, 10))
        new_point.ContImage = view_2dmol(str(m.context), os.path.join(settings.STATIC_ROOT,"MMPMaker","Images",f_name)+".png",legends=None)
        # Now finde the molecules linked to this molecule
        mol_list = [str(x.act_id.cmpd_id.smiles) for x in my_mps]

        # Now define the second filename
        f_name = ''.join(random.sample(string.letters * 10, 10))
        if len(mol_list) == 1:
            new_point.MolImage = view_2dmol(mol_list[0], os.path.join(settings.STATIC_ROOT, "MMPMaker", "Images", f_name) + ".png", legends=None)
        elif len(mol_list) > 1:
            new_point.MolImage = view_2dmol(mol_list, os.path.join(settings.STATIC_ROOT, "MMPMaker", "Images", f_name) + ".png", legends=None)
        else:
            new_point.save()
            continue
        # Work out all the similarities
        my_mol = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(str(m.context)))
        sim = AllChem.DataStructs.BulkTanimotoSimilarity(my_mol, Mol_fps)
        # Find the max
        max_value = max(sim)
        max_index = sim.index(max_value)
        # Find this nearest compound
        new_point.NearestCmpd = my_mols[max_index].prot_id.code
        new_point.NearestCmpdDist = max_value
        # Now overlay my molecule and save these coordinates
#### code to add the 3D coordinates of this newly found compound
#    overlay_mol = Chem.MolFromMolBlock(str(my_mols[max_index].sdf_info))
##    # Now make the molecule with this MCS
#    out_m = Chem.MolFromSmiles(str(m.context))
#    my_mcs = Chem.MolFromSmarts(MCS.FindMCS([overlay_mol,out_m],matchValences=True, ringMatchesRingOnly=True, completeRingsOnly=True).smarts)
#    core_mol =AllChem.ReplaceSidechains(overlay_mol,my_mcs)
##    # Now generate the  3D co-ordinates of this molecule
#    try:
#      AllChem.ConstrainedEmbed(out_m, core_mol )
#    except:
#      print m.context
#      print Chem.MolToSmiles(core_mol)
#      sys.exit()
#    new_point.SDFOverlay = Chem.MolToMolBlock(out_m)
#    new_point.save()
    return


def get_atom(atom_in):
    """Function to get an atom VDW radius from a PDB atomtype"""
    atom_dict = {"N": 7, "C": 6, "O": 8, "S": 16, "H": 1}
    if atom_in[0] in atom_dict:
        return atom_dict[atom_in[0]]
    else:
        print "UNRECOGNISED ATOM TYPE GIVING DEFAULT VALUE - CARBON"
        print "ATOM TYPE ", atom_in
        return atom_dict["C"]


def get_coords(target, prot_list=None, comp_coords=None):
    """Function to return the coords for a protein as a list of points
    Takes a target and an optional list of proteins. If none then it just calculates it for the target
    Returns a list of points"""
    # Only find coor
    if prot_list is None:
        prot_pos = []
        prot = str(Protein.objects.get(code=target.title + "TEMP", target_id=target).pdb_info)
        for line in prot.split("\n"):
            if line[:4] != "ATOM":
                continue
            new_point = Protein()
            new_point.atom_num = int(get_atom(line[12:16].strip()))
            new_point.res_name = line[17:20].strip()
            new_point.res_num = line[22:26].strip()
            new_point.x = float(line[30:38])
            new_point.y = float(line[38:46])
            new_point.z = float(line[46:54])
            new_point.res = ""
            if comp_coords:
                if min([eucl_dist([new_point.x, new_point.y, new_point.z], coord) for coord in comp_coords]) < 10.0:
                    prot_pos.append(new_point)
            else:
                prot_pos.append(new_point)
        return prot_pos
    elif type(prot_list) is list:
        prots = Protein.objects.filter(code__in=prot_list, target_id=target).values_list("pdb_info",flat=True)
        tot_pos = []
        for prot in prots:
            prot_pos = []
            for line in prot.split("\n"):
                if line[:4] != "ATOM":
                    continue
                new_point = Protein()
                new_point.atom_num = int(get_atom(line[12:16].strip()))
                new_point.x = float(line[30:38])
                new_point.y = float(line[38:46])
                new_point.z = float(line[46:54])
                new_point.res_name = line[17:20].strip()
                new_point.res_num = line[22:26].strip()
                if comp_coords:
                    if min([eucl_dist([new_point.x, new_point.y, new_point.z], coord) for coord in comp_coords]) < 10.0:
                        prot_pos.append(new_point)
                else:
                    prot_pos.append(new_point)
            # Add in this list
            tot_pos.append(prot_pos)
        return tot_pos
    else:
        print "Prot list is not list..."
        return None


def check_target_rmsd_spesh(target_id, mol_id):
    """Function to find the RMSD values for a target using different numbers
    of LLOOMMPPAA confs"""
    target = Target.objects.get(pk=target_id)
    mols = Molecule.objects.filter(prot_id__target_id=target_id, pk__gte=mol_id)
    out_f = open(target.title + "_rmsds_TOT_SPESH.csv", "w")
    out_f.write("MOL,FRAG,NUM_ROT,HEAV_AT_FRAG,HEAV_AT_NON_MATCH,NUM_ROT_NON_MATCH,300_MAX,300_MIN,300_CONFS,200_MAX,200_MIN,200_CONFS,100_MAX,100_MIN,100_CONFS,50_MAX,50_MIN,50_CONFS,20_MAX,20_MIN,20_CONFS,10_MAX,10_MIN,10_CONFS,2_MAX,2_MIN,2_CONFS\n")
    for mol in mols:
        out_f = check_cmpd_rmsd(mol.pk, out_f)


def check_target_rmsd(target_id):
    """Function to find the RMSD values for a target using different numbers
    of LLOOMMPPAA confs"""
    target = Target.objects.get(pk=target_id)
    mols = Molecule.objects.filter(prot_id__target_id=target_id)
    out_f = open(target.title + "_rmsds_TOT.csv", "w")
    out_f.write("MOL,FRAG,NUM_ROT,HEAV_AT_FRAG,HEAV_AT_NON_MATCH,NUM_ROT_NON_MATCH,300_MAX,300_MIN,300_CONFS,200_MAX,200_MIN,200_CONFS,100_MAX,100_MIN,100_CONFS,50_MAX,50_MIN,50_CONFS,20_MAX,20_MIN,20_CONFS,10_MAX,10_MIN,10_CONFS,2_MAX,2_MIN,2_CONFS\n")
    for mol in mols:
        out_f = check_cmpd_rmsd(mol.pk, out_f)


def check_cmpd_rmsd(mol_id, out_f):
    """Function to test the compound RMSD for a given conformation"""
    from MMPMaker.helpers import find_rmsd
    # Get the molecule
    my_mol = Molecule.objects.get(pk=mol_id)
    rdmol = Chem.MolFromMolBlock((my_mol.sdf_info))
    my_frags = MMPFrag.objects.filter(mol_id=my_mol).exclude(smistore__contains="[*:2]").exclude(smistore__contains="[*:3]")
    # Get the confotmers
    num_rot = Chem.Lipinski.NumRotatableBonds(rdmol)
    for my_f in my_frags:
        heav_at = Chem.Lipinski.HeavyAtomCount(Chem.MolFromSmiles(str(my_f.smistore)))
        heav_at_non = Chem.Lipinski.HeavyAtomCount(Chem.MolFromSmiles(str(my_f.mmp_link.context)))
        num_rot_non = Chem.Lipinski.NumRotatableBonds(Chem.MolFromSmiles(str(my_f.mmp_link.context)))
        # Skip less than three
        if heav_at < 3:
            continue
        out_f.write(str(mol_id) + "," + str(my_f.pk) + "," + str(num_rot) + "," + str(heav_at) + "," + str(heav_at_non) + "," + str(num_rot_non))
        context = Chem.MolFromSmiles(str(my_f.smistore))
        context = Chem.MolFromSmarts(make_smarts_from_frag(Chem.MolToSmiles(context, isomericSmiles=True)))
        print my_f.smistore
        coreMol = find_core(rdmol, rdmol, context, option="SIMPLE")
        print Chem.MolToSmiles(coreMol)
        print "300"
        out_confs = gen_confs_lloommppaa(rdmol, coreMol, num_confs=300, ff="MMFF", max_iters=300)
        my_rmsds = [find_rmsd(Chem.MolFromMolBlock(x[0]), rdmol) for x in out_confs]
        if not my_rmsds:
            out_f.write("\n")
            continue
        num_confs = len([x for x in out_confs if x[3] > 0.35])
        print "MIN", min(my_rmsds)
        print "MAX", max(my_rmsds)
        out_f.write(","+str(max(my_rmsds)) + "," + str(min(my_rmsds)) + "," + str(num_confs))
        print "200"
        out_confs = gen_confs_lloommppaa(rdmol, coreMol, num_confs=200, ff="MMFF", max_iters=200)
        my_rmsds = [find_rmsd(Chem.MolFromMolBlock(x[0]), rdmol) for x in out_confs]
        num_confs = len([x for x in out_confs if x[3] > 0.35])
        print "MIN", min(my_rmsds)
        print "MAX", max(my_rmsds)
        out_f.write(","+str(max(my_rmsds)) + "," + str(min(my_rmsds)) + "," + str(num_confs))
        print "100"
        out_confs = gen_confs_lloommppaa(rdmol, coreMol, num_confs=100, ff="MMFF", max_iters=100)
        my_rmsds = [find_rmsd(Chem.MolFromMolBlock(x[0]), rdmol) for x in out_confs]
        num_confs = len([x for x in out_confs if x[3] > 0.35])
        print "MIN", min(my_rmsds)
        print "MAX", max(my_rmsds)
        out_f.write(","+str(max(my_rmsds)) + "," + str(min(my_rmsds)) + "," + str(num_confs))
        print "50"
        out_confs = gen_confs_lloommppaa(rdmol, coreMol, num_confs=50, ff="MMFF", max_iters=50)
        my_rmsds = [find_rmsd(Chem.MolFromMolBlock(x[0]), rdmol) for x in out_confs]
        num_confs = len([x for x in out_confs if x[3] > 0.35])
        print "MIN", min(my_rmsds)
        print "MAX", max(my_rmsds)
        out_f.write(","+str(max(my_rmsds)) + "," + str(min(my_rmsds)) + "," + str(num_confs))
        print "20"
        out_confs = gen_confs_lloommppaa(rdmol, coreMol, num_confs=20, ff="MMFF", max_iters=20)
        my_rmsds = [find_rmsd(Chem.MolFromMolBlock(x[0]), rdmol) for x in out_confs]
        num_confs = len([x for x in out_confs if x[3] > 0.35])
        print "MIN", min(my_rmsds)
        print "MAX", max(my_rmsds)
        out_f.write("," + str(max(my_rmsds)) + "," + str(min(my_rmsds)) + "," + str(num_confs))
        print "10"
        out_confs = gen_confs_lloommppaa(rdmol, coreMol, num_confs=10, ff="MMFF", max_iters=10)
        my_rmsds = [find_rmsd(Chem.MolFromMolBlock(x[0]), rdmol) for x in out_confs]
        num_confs = len([x for x in out_confs if x[3] > 0.35])
        print "MIN", min(my_rmsds)
        print "MAX", max(my_rmsds)
        out_f.write(","+str(max(my_rmsds))+"," + str(min(my_rmsds)) + "," + str(num_confs))
        print "2"
        out_confs = gen_confs_lloommppaa(rdmol, coreMol, num_confs=2, ff="MMFF", max_iters=2)
        my_rmsds = [find_rmsd(Chem.MolFromMolBlock(x[0]), rdmol) for x in out_confs]
        num_confs = len([x for x in out_confs if x[3] > 0.35])
        print "MIN", min(my_rmsds)
        print "MAX", max(my_rmsds)
        out_f.write(","+ str(max(my_rmsds)) + ","+str(min(my_rmsds))+","+str(num_confs))
        out_f.write("\n")
    return out_f


def make_cmpd_conds(m, mol1, mol2, coreMol, target, sp, new_prot, comp_ref, ll_smiles, mmp=None, context=None, sdf_flag=False, title="NONE", prot_pos=[], react_proc=None):
    """Function to make all the LLOOMMPPAA conformations
    m - the django molecule
    mol1 - rdmol of the xtal structure
    mol2 - rdmol of the lloommppaa structure
    coreMol - the shared core with 3D coords
    target - the target
    sp - the synthpoint for this connection
    new_prot - the lloommppaa prot
    comp_ref -  the django cmpd
    ll_smiles - the smiles for the lloommppaa mol
    mmp - the mmp list of attributes
    context - the context
    sdf_flag - whether this is from an OAKLEY style or normal style
    """
    from LLOOMMPPAA.models import LLConf
    # If it is a H-change then the coreMol is this one
    if mmp:
        if mmp[0].sdf_info == "":
            coreMol = mol1
    # Check there is only one option. If there are two then we have a problem. With the current implementation we can only
    if context:
        if len(mol2.GetSubstructMatches(context)) > 1:
            print "ERROR TWO POSSIBLE MATCHES"
            print "ERROR " + Chem.MolToSmiles(context)
            print "ERROR " + Chem.MolToSmiles(mol2)
            return
    #Now the function to do the constraining
    try:
        mol2 = AllChem.ConstrainedEmbed(mol2, coreMol)
    except ValueError:
        # If it is a nitrile it needs extra sanitisation before embedding.
        Chem.SanitizeMol(coreMol)
        try:
            mol2 = AllChem.ConstrainedEmbed(mol2, coreMol)
        except ValueError:
            try:
              # Setting the mol as a molblock and reading back in helps
                mol2 = AllChem.ConstrainedEmbed(Chem.MolFromMolBlock(Chem.MolToMolBlock(mol2)), Chem.MolFromMolBlock(Chem.MolToMolBlock(coreMol)))
            except ValueError:
              # Excpet when it doesn't so print out this
                print "ERROR MOL ONE:\n", Chem.MolToSmiles(mol1)
                print "ERROR MOL TWO:\n", Chem.MolToSmiles(mol2)
                print "ERROR CORE MOL:\n", Chem.MolToSmiles(coreMol)
                if context:
                    print "ERROR CONTEXT:\n", Chem.MolToSmiles(context)
                print "ERROR CANNOT CONSTRAINED EMBED MOLECULE"
                return
        # If the molecule is the same as the core then why bother
        if Chem.MolToSmiles(mol2, isomericSmiles=True) == Chem.MolToSmiles(coreMol, isomericSmiles=True):
            return
    # Now add in here the conformation generation and filtering
    out_confs = gen_confs_lloommppaa(mol2, coreMol, num_confs=200, ff="MMFF", prot_pos=prot_pos, max_iters=200)
    from rdkit.Chem.rdShapeHelpers import ShapeProtrudeDist
    if out_confs:
        # Problem here is that a new one will be made every time we run the analysis because sdf_info will be different...
        # However you would want the molecule to be made if it sat in very different areas
        print "Scoring conformers..."
#        my_f = "ll_" + str(title) + "_MOD_CL_15.sdf"
#        print my_f
#        out_sd = Chem.SDWriter(my_f)
        i = 0
        for conf, clash, strain, rmsd in out_confs:
            i += 1
            # Now store these conformations protrude dist
            shape_dist = ShapeProtrudeDist(mol1, Chem.MolFromMolBlock(str(conf)))
            newmol = Molecule.objects.get_or_create(prot_id=new_prot, smiles=comp_ref.smiles, sdf_info=conf, cmpd_id=comp_ref, lig_id="UNL", occupancy=0.0, chain_id="X")[0]
            ll_conf = LLConf()
            ll_conf.shape_dist = shape_dist
            ll_conf.target_id = target
            ll_conf.mol_id = newmol
            ll_conf.mol_ref = m
            ll_conf.prot_id = m.prot_id
            ll_conf.clash = clash
            ll_conf.strain = strain
            ll_conf.conf_num = i
            # The min RMSD of this conf
            ll_conf.rmsd = rmsd
            try:
                ll_conf.validate_unique()
                ll_conf.save()
                # Now add it to the reaction process object - if it exists
                if react_proc:
                    react_proc.ll_conf.add(ll_conf)
            except ValidationError:
                pass
            # Add this to the list - for the synthpoint
            #sp.mol_id.add(newmol)
            # Now score this molecule - clash informations
            # Now write these to a file of conformations
            #out_sd.write(Chem.MolFromMolBlock(conf))
            # score_mol(newmol, target.pk, sp, mmp[0].mmp_link,  mmp[0].mol_id, div, clash)
        # Now remove this from that quueue
        react_proc.product_queue.remove(comp_ref)
        react_proc.save()


def make_synth_points(target_id=3, mmpdbs=["MMPDB"], sdf_flag=False, react_proc=None):
    """Function to find synthetic attachment points located in a
    large database.
    Takes a target id, an optional list of databases (defined in the dbs.py file)
    and a flag to determine wheter this is MMP search or simple SDF superstructure.
    If sdf_flag is true - mmpdbs is a dict - mapping internal id to SDF library.
    Returns None"""
    # Do the same on all molecules by removing H's in turn
    # Need to delete all of them at the start
    #try:
    #    SynthPoint.objects.filter(target_id=target_id).delete()
    #except:
    #    for s in SynthPoint.objects.filter(target_id=target_id):
    #        s.mol_id.all().delete()
    #        ep = s.extrapoint_set.all()
    #        for e in ep:
    #            e.delete()
    #Make the pharmacophore factory
    print "\nSetting up environment"
    if react_proc:
        # Loop through the crystal structures
        react_proc.proc_stage = "GENERATE CONFS"
        react_proc.stage_completion = 0
        react_proc.save()
    fdefName = os.path.join(os.path.join(os.path.split(sys.argv[0])[0], 'data/media'), 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    # Get the target object
    target = Target.objects.get(pk=target_id)
    # Get the protein coordinates
    # Make a new one
    new_prot = Protein.objects.get_or_create(code=target.title + "SYNTH", target_id=target)[0]
    if not react_proc:
        mols = Molecule.objects.filter(prot_id=new_prot)
        #   for m in mols:
        #       m.delete()
        # Get all the molecules related to this target - with xtal strcutre
        mols = Molecule.objects.filter(prot_id__target_id=target).exclude(prot_id__code__contains=target.title)
        if not sdf_flag:
            Mol_fps = [MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(str(x.cmpd_id.smiles))) for x in ActivityPoint.objects.filter(target_id=target_id)]
            Mol_fps.extend([MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(str(x.cmpd_id.smiles))) for x in Molecule.objects.filter(prot_id__target_id=target_id)])
        cont_list = []
        tot_list = []
        # First loop through these and find all the relevant H-changes and contexts
        print "Finding contexts and H-changes...."
        tot = len(mols)
    old = -1
    # If we are just taking substructures from the SDF file
    if sdf_flag:
        # Loop through the crystal structures
        for m in mols:
            # Get the internal ID for this
            iid = InternalIDLink.objects.filter(mol_id=m)[0].internal_id
            # Look to see if there is a database - IF NOT CONTINUE
            db_path = 'C:\\Users\\abradley\\Downloads\\' + iid + ".sdf"
            if os.path.isfile(db_path):
                db = Chem.SDMolSupplier(str(db_path))
            else:
                print "NO SUCH FILE", db_path
                continue
            mol1 = Chem.MolFromMolBlock(str(m.sdf_info))
            sp_dict = {}
            # Loop through the database
            # Filter on the bound conformation - to stop clashes
            conf = mol1.GetConformer()
            comp_coords = []
            for atom in mol1.GetAtoms():
                atom_id = atom.GetIdx()
                position = conf.GetAtomPosition(atom_id)
                comp_coords.append((position.x, position.y, position.z))
            prot_pos = get_coords(target, [m.prot_id.code], comp_coords)
            tot = len(db)
            for counter, mol2 in enumerate(db):
                print counter, "of ", tot
#                if counter < 1715:
#                    continue
                # Get a maximum substruct
#                context = Chem.MolFromSmarts(MCS.FindMCS([mol2, mol1], matchValences=True, completeRingsOnly=True).smarts)
#                print Chem.MolToSmiles(context)
                context = Chem.MolFromSmarts("CNC(=O)C(C)C")
                # Use this to get the coreMol
                coreMol = find_core(mol1, mol2, context, option="SIMPLE")

                # Make a comp_ref
                comp_ref = add_new_comp(mol2)
                # Get the smiles of the molecule
                ll_smiles = Chem.MolToSmiles(mol2, isomericSmiles=True)
                # Now do this
                print "MAKING CONFORMATIONS"
                make_cmpd_conds(m, mol1, mol2, coreMol, target, None, new_prot,
                                comp_ref, ll_smiles, mmp=None, context=None,
                                sdf_flag=True, title=iid+"_TO_MODEL_"+str(comp_ref.pk)+".sdf",
                                prot_pos=prot_pos)
    elif react_proc:
        # Loop through the crystal structures
        print "GENERATING CONFS"
        react_proc.proc_stage = "GENERATE CONFS"
        react_proc.stage_completion = 0
        react_proc.save()
        m = react_proc.mol_id
        iid = InternalIDLink.objects.filter(mol_id=m)[0].internal_id
        # Get the internal ID for this
        mol1 = Chem.MolFromMolBlock(str(m.sdf_info))
        # Filter on the bound conformation - to stop clashes
        conf = mol1.GetConformer()
        comp_coords = []
        for atom in mol1.GetAtoms():
            atom_id = atom.GetIdx()
            position = conf.GetAtomPosition(atom_id)
            comp_coords.append((position.x, position.y, position.z))
        prot_pos = get_coords(target, [prot_id.code for prot_id in react_proc.prot_id.all()], comp_coords)
        tot_cmps = react_proc.product_queue.all()
        tot = len(tot_cmps)
        print "GOT ", tot, "COMPOUNDS"
        for counter, comp_ref in enumerate(tot_cmps):
            sys.stdout.write("\rGenerating conformations for compund %d of %d..." % (counter, tot))
            sys.stdout.flush()
            react_proc.stage_completion = int((float(counter) / float(tot)) * 100)
            react_proc.save()
            # Get the molecule
            mol2 = Chem.MolFromSmiles(str(comp_ref.smiles))
            print counter, "of ", tot
            context = Chem.MolFromSmarts(str(react_proc.context))
            # Use this to get the coreMol
            coreMol = find_core(mol1, mol2, context, option="SIMPLE")
            # Get the smiles of the molecule
            ll_smiles = Chem.MolToSmiles(mol2, isomericSmiles=True)
            # Now do this
            make_cmpd_conds(m, mol1, mol2, coreMol, target, None, new_prot, comp_ref, ll_smiles, mmp=None, context=None, sdf_flag=True, title=iid+"_TO_MODEL_"+str(comp_ref.pk)+".sdf", prot_pos=prot_pos, react_proc=react_proc)
    else:
        # Loop through the molecules - getting the H-addition positions
        for mol_counter, m in enumerate(mols):
            # Print the progress
            if mol_counter * 100 / tot != old:
                old = mol_counter * 100 / tot
                sys.stdout.write("\r%d%% complete..." % old)
                sys.stdout.flush()
            # First get my pre-existing MMPFrags for all molecules
            smi_list = [[x, None] for x in MMPFrag.objects.filter(mol_id=m).exclude(smistore__contains="[*:2]").exclude(smistore="*[H]")]
            # Now loop over the H's and replace with [*:1] to get my H-index frag
            rdm = Chem.MolFromMolBlock(str(m.sdf_info))
            rdm = AllChem.AddHs(rdm)
            atms = rdm.GetAtoms()
            for atm in atms:
                if atm.GetSmarts() == "[H]":
                    em = AllChem.EditableMol(rdm)
                    em.ReplaceAtom(atm.GetIdx(), Chem.Atom("*"))
                    # Now make the smiles
                    nm = em.GetMol()
                    smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(nm, isomericSmiles=True).replace("*", "*:1")),isomericSmiles=True)
                    cont = MMP.objects.get_or_create(context=Chem.MolToSmiles(Chem.MolFromSmiles(smi), isomericSmiles=True))[0]
                    cont.compound_id.add(m.cmpd_id)
                    cont.save()
                    mmpid = make_h_frag(m, cont, str(m.cmpd_id.smiles))
                    smi_list.append((mmpid, atm.GetIdx()))
            # Now loop through and only take the context once
            for mysmi in smi_list:
                # Check the context contains the conserved substruct for this target
                if m.smiles in cons_subs:
                    if not Chem.MolFromSmiles(str(mysmi[0].mmp_link.context)).HasSubstructMatch(Chem.MolFromSmiles(cons_subs[m.smiles])):
                        continue
                # Make sure the mmpid is not None
                if mysmi[0] is None:
                    continue
                # Make sure it is => 5 heavy atoms
                if Lipinski.HeavyAtomCount(Chem.MolFromSmiles(mysmi[0].mmp_link.context)) <= 5:
                    continue
                # Make sure the context is not already in the conect list
                if mysmi[0].mmp_link.context in cont_list:
                    continue
                else:
                    # Append to the context list and to the total list
                    cont_list.append(mysmi[0].mmp_link.context)
                    tot_list.append([mysmi[0], m, mysmi[1]])
        # Final print out
        sys.stdout.write("\r%d%% complete..." % old)
        sys.stdout.flush()
        # Now search this list in the database
        # Set up the progress counter
        mmp_tot = len(tot_list)
        mmp_old = -1
        print "Searching specified databases for MMPS from ", str(mmp_tot), ""
        print " ".join(mmpdbs)
        # For each of the context's produced above - generate synthpoints and MMPs
        for mmp_counter, mmp in enumerate(tot_list):
            # Print the progress
            print str(mmp[0].mmp_link.context)
            if mmp_counter * 100 / mmp_tot != mmp_old:
                mmp_old = mmp_counter * 100 / mmp_tot
                sys.stdout.write("\r Synth point creation %d%% complete..." % mmp_old)
                sys.stdout.flush()
            # Iterate through the dbs
            cmps = []
            for db in mmpdbs:
                out_cmps = search_mmp_db(context=str(mmp[0].mmp_link.context), dbname=db)
                if out_cmps is None:
                    continue
                cmps.extend(out_cmps)
            # Now deal with these compounds
            if len(cmps) == 0:
                continue
            print str(mmp[0].mmp_link.context), ": found " + str(len(cmps)) + " molecules"
            print "WORKING OUT DIVERSITY"
            old = 1.0
            tot = 0.0
            # Now work out the most diverse mol
            if len(Mol_fps) > 0:
                for c in cmps:
                    my_mol = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(str(c)))
                    sim = max(AllChem.DataStructs.BulkTanimotoSimilarity(my_mol, Mol_fps))
                    if sim < old:
                        old = sim
                    tot += sim
                # Now work out the average diversity
                tot = tot / len(cmps)
            # Get an RDKit molecule of the thing we're
            oldrdm = Chem.MolFromMolBlock(str(mmp[0].mol_id.sdf_info))
            rdm = AllChem.AddHs(oldrdm, addCoords=True)
            try:
                AllChem.ConstrainedEmbed(rdm, oldrdm)
                H_allow = True
            except ValueError:
                # Wont work for symmetric molecules where fragment is the symmetric
                # part. (Triangle bound smooth, use isotopes to fix??)
                rdm = Chem.MolFromMolBlock(str(mmp[0].mol_id.sdf_info))
                rdm = AllChem.AddHs(rdm, addCoords=True)
                H_allow = False
                pass
            # Find the atom in the substructure match but a neighbour to one that isn't
            context = Chem.MolFromSmiles(str(mmp[0].mmp_link.context))
            context = Chem.MolFromSmarts(make_smarts_from_frag(Chem.MolToSmiles(context, isomericSmiles=True)))
            # Get the conformer for the molecule
            cnf = rdm.GetConformer()
            # Code to pick up the start and end point
            if mmp[2] is None:
                # If it is an MMPFrag and not just an H pick up the coords
                eid = None
                sid = None
                smatch = rdm.GetSubstructMatch(context)
                not_match = [x.GetIdx() for x in rdm.GetAtoms() if x.GetIdx() not in smatch]
                for x in not_match:
                    if rdm.GetAtomWithIdx(x).GetSmarts() == "[H]":
                        continue
                    for y in rdm.GetAtomWithIdx(x).GetNeighbors():
                        if y.GetIdx() in smatch:
                            # These are the synthpoints
                            if eid not in [None, x] or sid not in [None, y.GetIdx()]:
                                print "ERROR MULTIPLE OPTIONS"
                                print eid
                                print sid
                                print x
                                print y.GetIdx()
                            else:
                                eid = x
                                sid = y.GetIdx()
                # Now assign these
                epv = cnf.GetAtomPosition(eid)
                spv = cnf.GetAtomPosition(sid)
            elif H_allow == True:
                epv = cnf.GetAtomPosition(mmp[2])
                spv = cnf.GetAtomPosition(rdm.GetAtomWithIdx(mmp[2]).GetNeighbors()[0].GetIdx())
            else:
                continue
            # Add the synth point - i.e. the position this is found at.
            sp = add_synth_point(spv, epv, target, mmp[0], len(cmps), old, tot)
            # Now make the molecule
            mol1 = Chem.MolFromMolBlock(str(mmp[1].sdf_info))
            # Now make conformations for the compounds returned
            print "MAKING CONFORMATIONS"
            for cmp_counter, c in enumerate(cmps):
                print "COMPOUND ", cmp_counter + 1, "of ", len(cmps)
                # Work out it's diversity from the other exisiting tested molecules
                my_mol = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(str(c)))
                div = max(AllChem.DataStructs.BulkTanimotoSimilarity(my_mol, Mol_fps))
                # Generate a new molecule with a conformation
                cm = Chem.MolFromSmiles(str(c))
                comp_ref = add_new_comp(cm)
                mol2 = cm
                coreMol = find_core(mol1, mol2, context)
                #make_cmpd_conds(mmp, mol2, coreMol, context=None, sdf_flag=False)
                make_cmpd_conds(m, mol1, mol2, coreMol, target, None, new_prot, comp_ref, c, mmp, context=None, sdf_flag=False)


def get_smarts(mol, feat):
    """Function to get the smiles for molecule for reduced graphs
    Takes an RDKit molecule and the feature
    Returns the smiles"""
    atms = feat.GetAtomIds()
    if len(atms) > 1:
        smi = ""
        for i, a in enumerate(atms):
            atm = mol.GetAtomWithIdx(a)
            smi += atm.GetSmarts()
            if i == 0:
                smi += "1"
        smi += "1"
    else:
        smi = mol.GetAtomWithIdx(atms[0]).GetSmarts()
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        mol = adjust_arom_Ns(mol)
        if mol is None:
            if sorted(smi) == sorted("c1coccc1"):
                mol = Chem.MolFromSmiles("C1CC=COC=1")
            elif sorted(smi) == sorted("c1[nH]cc[nH]1"):
                mol = Chem.MolFromSmiles("c1ncc[nH]1")
            elif sorted(smi) == sorted("c1ccnc[nH]1"):
                mol = Chem.MolFromSmiles("c1ccncn1")
            elif len("".join(re.findall("[a-z]+", smi))) == 7 and len(smi) == 9 and len([x for x in smi if x == "1"]) == 2:
                #And then just return as there is no easy way of forming a canonical form of all the 7 rings
                return smi
            else:
                print "STILL BROKEN", smi, feat.GetFamily()
                return smi
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def find_reduced_graphs(mol):
    """Function to find reduced graphs for a molecule.
    Takes an RDKit molecule
    Returns None"""
    rdmol = Chem.MolFromMolBlock(str(mol.sdf_info))
    factory = ChemicalFeatures.BuildFeatureFactoryFromString(feats)
    myfs = factory.GetFeaturesForMol(rdmol)
    nodedict = {}
    atomlist = {}
    counter = 0
    ringlist = []
    for f in myfs:
        # Get the rings and acceptors and donors and save them
        if f.GetFamily() == "Aliphatic":
            counter += 1
            nodedict[counter] = f
            ringlist.append(counter)
            for i in f.GetAtomIds():
                if i in atomlist:
                    atomlist[i].append(counter)
                else:
                    atomlist[i] = [counter]
        elif f.GetFamily() == "Aromatic":
            counter += 1
            ringlist.append(counter)
            nodedict[counter] = f
            for i in f.GetAtomIds():
                if i in atomlist:
                    atomlist[i].append(counter)
                else:
                    atomlist[i] = [counter]
        elif f.GetFamily() == "Acceptor":
            counter += 1
            nodedict[counter] = f
            for i in f.GetAtomIds():
                if i in atomlist:
                    atomlist[i].append(counter)
                else:
                    atomlist[i] = [counter]
        elif f.GetFamily() == "Donor":
            counter += 1
            nodedict[counter] = f
            for i in f.GetAtomIds():
                if i in atomlist:
                    atomlist[i].append(counter)
                else:
                    atomlist[i] = [counter]
    old_list = []
    mylist = [atomlist[x] for x in atomlist]
    for atom in mylist:
        if len([x for x in atom if x in ringlist]) > 1:
            continue
        union = set(atom)
        for atom2 in mylist:
            if [x for x in atom if x in atom2]:
                # Don't do the union if they both contain ring nodes from different rings
                tmpunion = union.union(union, atom2)
                if len([x for x in tmpunion if x in ringlist]) > 1:
                    continue
                else:
                    union = tmpunion
        old_list.append(tuple(union))
    # Now take only unique elements
    out_list = set(old_list)
    # Now go through these combinations
    for item in out_list:
        name = ""
        main = None
        for x in item:
            if nodedict[x].GetFamily() in ["Aromatic", "Aliphatic"]:
                name = nodedict[x].GetFamily() + name
                main = nodedict[x]
            elif nodedict[x].GetFamily() == "Acceptor":
                if "Acceptor" in name:
                    pass
                elif "Donor" in name and "AcceptorDonor" not in name:
                    name = name.replace("Donor", "AcceptorDonor")
                else:
                    name = name + nodedict[x].GetFamily()
            elif nodedict[x].GetFamily() == "Donor":
                if "Donor" in name:
                    pass
                elif "Acceptor" in name and "AcceptorDonor" not in name:
                    name = name.replace("Acceptor", "AcceptorDonor")
                else:
                    name = name + nodedict[x].GetFamily()
        if main is None:
            main = nodedict[x]
        new_ph4 = RedGraph()
        # Record the coords
        pos = main.GetPos()
        new_ph4.x_com, new_ph4.y_com, new_ph4.z_com = pos.x, pos.y, pos.z
        # Relate back to the molecule
        new_ph4.mol_id = mol
        # Relate back to the smiles of the substructure
        new_ph4.func_group = get_smarts(rdmol, main)
        new_ph4.type = "BLAH"
        # Relate back to ring size (do this through the name)
        new_ph4.ring_size = len(main.GetAtomIds())
        new_ph4.feat = name
        new_ph4.save()
    return


def make_reduced_graphs(target_id):
    """Function to make reduced graphs for all the molecules in a target
    Takes a target id
    Returns None"""
    target = Target.objects.get(pk=target_id)
    mols = Molecule.objects.filter(prot_id__target_id=target).exclude(prot_id__code__contains=target.title)
    rgs = RedGraph.objects.filter(mol_id__in=mols)
    for r in rgs:
        r.delete()
    for mol in mols:
        find_reduced_graphs(mol)
    return None


def eucl_dist(mol_1, mol_2):
    """Function to find the Euclidean distance between two molecules in 3D
    Takes two len=3 tuples
    Returns a float"""
    return sqrt(pow((mol_1[0]-mol_2[0]),2)+pow((mol_1[1]-mol_2[1]),2)+pow((mol_1[2]-mol_2[2]),2))


def make_actpharma_clusters(target_id):
    """Function to`` make the clusters for the ActPharmaPoint for a target
    Takes a target id
    Returns None"""
    type_list = ["Arom", "Acceptor", "Donor"]
    lam = 4.0
    opt_list = ["ACT", "INACT"]
    clus_dict = {"ACT": {}, "INACT": {}}
    for o in opt_list:
        for t in type_list:
            clus_dict[o][t] = []
    # Refresh them
    MMPClus.objects.filter(target_id=target_id).delete()
    # Take all the objects
    for t in type_list:
        for opt in opt_list:
            # Get the points
            if opt == "INACT":
                myps = ActPharmaPoint.objects.filter(target_id=target_id, pharma_id__smiles__contains=t, in_diff_15=True, num_diff_15__lte=4,act_change__lt=0.0)
            else:
                myps = ActPharmaPoint.objects.filter(target_id=target_id, pharma_id__smiles__contains=t, in_diff_15=True, num_diff_15__lte=4,act_change__gt=0.0)
            dp = dpmeans([(x.pharma_id.x_com, x.pharma_id.y_com, x.pharma_id.z_com) for x in myps],lam,0,False)
            dp.run()
            # Now assign the points in the lists to their clusters
            out_d = {}
            for i, val in enumerate(dp.dataClusterId):
                if val in out_d:
                    out_d[val].append(myps[i])
                else:
                    out_d[val] = [myps[i]]
            for val in out_d:
                # Now make the cluster
                nclus = MMPClus()
                nclus.sdf_info = make_map([x.pharma_id for x in out_d[val]])
                nclus.type = t
                acts = [x.act_change for x in out_d[val]]
                nclus.mean_change = sum(acts) / len(acts)
                if nclus.mean_change < 0.0:
                    nclus.max_change = min(acts)
                else:
                    nclus.max_change = max(acts)
                nclus.target_id = Target.objects.get(pk=target_id)
                nclus.num_points = len(acts)
                # Now find the number of different scaffolds this is from
                nclus.num_scaffs = len(set([x.act_3d_frag_id.mol_id.pk for x in out_d[val]]))
                nclus.save()
                # Add the centre of mass of the cluster, referring to the pk
                clus_dict[opt][t].append([(dp.clusters[val]), nclus.pk])
                # Now find the MMPCompTD objects it relates to
                for x in out_d[val]:
                    nclus.mmp_comp_id.add(MMPComparison.objects.filter(xtal_act__target_id=target_id, xtal_act__cmpd_id__pk=x.act_3d_frag_id.over_mol_id.cmpd_id.pk, chembl_act__cmpd_id__pk=x.act_3d_frag_id.mol_id.cmpd_id.pk)[0])
                    nclus.save()
    # Now find the conflicts and add them as links. Conflicts are clusters of the same type opposite sign for which the cluster centre is within 1.5 A
    for t in type_list:
        # for the types take good and bad and find ones that are close by
        print t
        for act in clus_dict["ACT"][t]:
            for inact in clus_dict["INACT"][t]:
                if eucl_dist(act[0], inact[0]) < 2.5:
                    # join them
                    mmpclus = MMPClus.objects.get(pk=act[1])
                    mmpclus.mmp_clus_id = MMPClus.objects.get(pk=inact[1])
                    mmpclus.save()


def search_mmp_db(context='[*:1]C.[*:2]C(=O)NCC(=O)N1CCCc2ccccc21', dbname="MMPDB"):
    """code to connect to an external MMPDB and search for MMPs
    1) Should take a timeout limit
    2) Number of heavy atoms max
    3) Max number of results"""
    # Take the context and find its ID in the DB
    try:
        # Only get those with more than 3 heavy atoms
        c = ContextTable.objects.using(dbname).get(context_smi=context.replace("[1*]", "[*:1]"),context_size__gt=3)
    except ObjectDoesNotExist:
        return None
    # Look now for the compound ID's related to this context # get a list of smisp then do an in search
    cores = CoreTable.objects.using(dbname).filter(context_id=c.context_id).values_list("cmpd_id", flat=True)
    # Search in the appropriate table column, only things that are less than 10 heavy atoms greater
    comps = CmpdSmisp.objects.using(dbname).filter(cmpd_size__lt=c.context_size + 10,cmpd_id__in=cores).values_list("SMILES", flat=True)
    # Return all appropriate compounds
    return comps


def search_rdkit_db(substruct,lib="ZINC",limit=None,opt=None,extra=2):
    """code to connect to an external Postgres RDKit cartridge database and do a substructure search. Returns all related molecules as a list"""
    # Search for smiles
    print "Looking for:",substruct
    if opt is None:
        my_type = ""
        tmp_mol = Chem.MolFromSmiles(substruct)
        if tmp_mol is None:
            print "INVALID SMILES"
            return None
    # Search for smarts
    elif opt == "SMARTS":
        my_type = "::qmol"
        tmp_mol = Chem.MolFromSmarts(substruct)
        if tmp_mol is None:
            print "INVALID SMARTS"
            return None
    cursor = connections[lib].cursor()
    numhevayatoms_extra = Lipinski.HeavyAtomCount(tmp_mol)+extra
    if limit == None:
        my_s = "select * from mols where m@>'"+substruct+"'"+my_type+" and mol_numheavyatoms(m) < "+str(numhevayatoms_extra)+";"
        print my_s
        cursor.execute(my_s)
    else:
        my_s = "select * from mols where m@>'" + substruct+"'"+my_type+" and mol_numheavyatoms(m) < "+str(numhevayatoms_extra)+" limit "+str(limit)+";"
        print my_s
        cursor.execute(my_s)
    return cursor.fetchall()


def search_db(smi, libs=[("ZINC", "RDKit")], limit=10, opt="SMARTS", extra=2):
    """Function to take a substruct and search an or several external database(s).
    Lib is a list of DB's to query, as a tuple with the type that the are as the second entry
    Limit is the maximum number to return
    Opt is whether to search for smarts or smiles
    Extra is how much bigger the query molecule can be than the entry molecule
    Returns a dictionary of molecules. DB name: [compounds]"""
    smarts = make_smarts_from_frag(smi)
    tot_cmps = {}
    for lib in libs:
        print "NOW MAKING SERVER REQUEST..."
        if lib[1] == "RDKit":
            cmpds = search_rdkit_db(smarts, lib[0], limit, opt, extra)
            # Now add this to the dictionary
            tot_cmps[lib[1]] = cmpds
        elif lib[1] == "MMP":
            out = search_mmp_db(smi.replace("[1*]", "[*:1]"), lib[0])
            if out is None:
                continue
            else:
                tot_cmps[lib[1]] = out[:10]
        else:
            print "Type not recognised", lib[1]
    print "COMPOUNDS RETRIEVED..."
    return tot_cmps


def find_actpharma_clusts(target_id):
    """Function to make the clusters for the ActPharmaPoint for a target
    and write them to a file
    Takes a target_id
    Returns None"""
    # List of types to consider
    type_list = ["Arom", "Acceptor", "Donor"]
    # List of lambda values to use
    lam_lis = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    opt_list = ["ACT", "INACT"]
    clus_dict = {"ACT": {}, "INACT": {}}
    for o in opt_list:
        for t in type_list:
            clus_dict[o][t] = []
    # Refresh them
    MMPClus.objects.filter(target_id=target_id).delete()
    # Take all the objects
    out_sd = Chem.SDWriter("clusts.sdf")
    for lam in lam_lis:
        for t in type_list:
            for opt in opt_list:
                # Get the points
                if opt == "INACT":
                    myps = ActPharmaPoint.objects.filter(target_id=target_id, pharma_id__smiles__contains=t, in_diff_15=True, num_diff_15__lte=4, act_change__lt=0.0)
                else:
                    myps = ActPharmaPoint.objects.filter(target_id=target_id, pharma_id__smiles__contains=t, in_diff_15=True, num_diff_15__lte=4, act_change__gt=0.0)
                dp = dpmeans([(x.pharma_id.x_com, x.pharma_id.y_com, x.pharma_id.z_com) for x in myps], lam, 0, False)
                dp.run()
                # Now get the clusters
                my_pts = []
                for cluster in dp.clusters:
                    m = MMPClus()
                    m.x_com = dp.clusters[cluster][0]
                    m.y_com = dp.clusters[cluster][1]
                    m.z_com = dp.clusters[cluster][2]
                    my_pts.append(m)
                mol = Chem.MolFromMolBlock(make_map(my_pts))
                mol.SetProp("_Name", str(lam) + "_" + t + "_" + opt)
                out_sd.write(mol)


feats = """# $Id$
#
# RDKit base fdef file.
# Created by Greg Landrum
#

AtomType NDonor [N&!H0&v3,N&!H0&+1&v4,n&H1&+0]
AtomType AmideN [$(N-C(=O))]
AtomType SulfonamideN [$([N;H0]S(=O)(=O))]
AtomType NDonor [$([Nv3](-C)(-C)-C)]

AtomType NDonor [$(n[n;H1]),$(nc[n;H1])]

AtomType ChalcDonor [O,S;H1;+0]
DefineFeature SingleAtomDonor [{NDonor},{ChalcDonor}]
  Family Donor
  Weights 1
EndFeature

# aromatic N, but not indole or pyrole or fusing two rings
AtomType NAcceptor [n;+0;!X3;!$([n;H1](cc)cc)]
AtomType NAcceptor [$([N;H0]#[C&v4])]
# tertiary nitrogen adjacent to aromatic carbon
AtomType NAcceptor [N&v3;H0;$(Nc)]

# removes thioether and nitro oxygen
AtomType ChalcAcceptor [O;H0;v2;!$(O=N-*)]
Atomtype ChalcAcceptor [O;-;!$(*-N=O)]

# Removed aromatic sulfur from ChalcAcceptor definition
Atomtype ChalcAcceptor [o;+0]

# Hydroxyls and acids
AtomType Hydroxyl [O;H1;v2]

# F is an acceptor so long as the C has no other halogen neighbors. This is maybe
# a bit too general, but the idea is to eliminate things like CF3
AtomType HalogenAcceptor [F;$(F-[#6]);!$(FC[F,Cl,Br,I])]

DefineFeature SingleAtomAcceptor [{Hydroxyl},{ChalcAcceptor},{NAcceptor},{HalogenAcceptor}]
  Family Acceptor
  Weights 1
EndFeature

# this one is delightfully easy:
DefineFeature AcidicGroup [C,S](=[O,S,P])-[O;H1,H0&-1]
  Family NegIonizable
  Weights 1.0,1.0,1.0
EndFeature

AtomType Carbon_NotDouble [C;!$(C=*)]
AtomType BasicNH2 [$([N;H2&+0][{Carbon_NotDouble}])]
AtomType BasicNH1 [$([N;H1&+0]([{Carbon_NotDouble}])[{Carbon_NotDouble}])]
AtomType PosNH3 [$([N;H3&+1][{Carbon_NotDouble}])]
AtomType PosNH2 [$([N;H2&+1]([{Carbon_NotDouble}])[{Carbon_NotDouble}])]
AtomType PosNH1 [$([N;H1&+1]([{Carbon_NotDouble}])([{Carbon_NotDouble}])[{Carbon_NotDouble}])]
AtomType BasicNH0 [$([N;H0&+0]([{Carbon_NotDouble}])([{Carbon_NotDouble}])[{Carbon_NotDouble}])]
AtomType QuatN [$([N;H0&+1]([{Carbon_NotDouble}])([{Carbon_NotDouble}])([{Carbon_NotDouble}])[{Carbon_NotDouble}])]


DefineFeature BasicGroup [{BasicNH2},{BasicNH1},{BasicNH0};!$(N[a])]
  Family PosIonizable
  Weights 1.0
EndFeature

# 14.11.2007 (GL): add !$([N+]-[O-]) constraint so we don't match
# nitro (or similar) groups
DefineFeature PosN [#7;+;!$([N+]-[O-])]
 Family PosIonizable
 Weights 1.0
EndFeature

# imidazole group can be positively charged (too promiscuous?)
DefineFeature Imidazole c1ncnc1
  Family PosIonizable
  Weights 1.0,1.0,1.0,1.0,1.0
EndFeature
# guanidine group is positively charged (too promiscuous?)
DefineFeature Guanidine NC(=N)N
  Family PosIonizable
  Weights 1.0,1.0,1.0,1.0
EndFeature

# the LigZn binder features were adapted from combichem.fdl
DefineFeature ZnBinder1 [S;D1]-[#6]
  Family ZnBinder
  Weights 1,0
EndFeature
DefineFeature ZnBinder2 [#6]-C(=O)-C-[S;D1]
  Family ZnBinder
  Weights 0,0,1,0,1
EndFeature
DefineFeature ZnBinder3 [#6]-C(=O)-C-C-[S;D1]
  Family ZnBinder
  Weights 0,0,1,0,0,1
EndFeature

DefineFeature ZnBinder4 [#6]-C(=O)-N-[O;D1]
  Family ZnBinder
  Weights 0,0,1,0,1
EndFeature
DefineFeature ZnBinder5 [#6]-C(=O)-[O;D1]
  Family ZnBinder
  Weights 0,0,1,1
EndFeature
DefineFeature ZnBinder6 [#6]-P(=O)(-O)-[C,O,N]-[C,H]
  Family ZnBinder
  Weights 0,0,1,1,0,0
EndFeature



# aromatic rings of various sizes:
#
# Note that with the aromatics, it's important to include the ring-size queries along with
# the aromaticity query for two reasons:
#   1) Much of the current feature-location code assumes that the feature point is
#      equidistant from the atoms defining it. Larger definitions like: a1aaaaaaaa1 will actually
#      match things like 'o1c2cccc2ccc1', which have an aromatic unit spread across multiple simple
#      rings and so don't fit that requirement.
#   2) It's *way* faster.
#

#
# 21.1.2008 (GL): update ring membership tests to reflect corrected meaning of
# "r" in SMARTS parser
#
AtomType AromR4 [a;r4,!R1&r3]
DefineFeature Arom4 [{AromR4}]1:[{AromR4}]:[{AromR4}]:[{AromR4}]:1
 Family Aromatic
 Weights 1.0,1.0,1.0,1.0
EndFeature
AtomType AromR5 [a;r5,!R1&r4,!R1&r3]
DefineFeature Arom5 [{AromR5}]1~[{AromR5}]~[{AromR5}]~[{AromR5}]~[{AromR5}]~1
 Family Aromatic
 Weights 1.0,1.0,1.0,1.0,1.0
EndFeature
AtomType AromR6 [a;r6,!R1&r5,!R1&r4,!R1&r3]
DefineFeature Arom6 [{AromR6}]1~[{AromR6}]~[{AromR6}]~[{AromR6}]~[{AromR6}]~[{AromR6}]:1
 Family Aromatic
 Weights 1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
AtomType AromR7 [a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]
DefineFeature Arom7 [{AromR7}]1~[{AromR7}]~[{AromR7}]~[{AromR7}]~[{AromR7}]~[{AromR7}]~[{AromR7}]:1
 Family Aromatic
 Weights 1.0,1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
AtomType AromR8 [a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]
DefineFeature Arom8 [{AromR8}]1:[{AromR8}]:[{AromR8}]:[{AromR8}]:[{AromR8}]:[{AromR8}]:[{AromR8}]:[{AromR8}]:1
 Family Aromatic
 Weights 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
EndFeature

# hydrophobic features
# any carbon that is not bonded to a polar atom is considered a hydrophobe
#
# 23.11.2007 (GL): match any bond (not just single bonds); add #6 at
#  beginning to make it more efficient
AtomType Carbon_Polar [#6;$([#6]~[#7,#8,#9])]
# 23.11.2007 (GL): don't match charged carbon
AtomType Carbon_NonPolar [#6;+0;!{Carbon_Polar}]

DefineFeature ThreeWayAttach [D3,D4;{Carbon_NonPolar}]
  Family Hydrophobe
  Weights 1.0
EndFeature

DefineFeature ChainTwoWayAttach [R0;D2;{Carbon_NonPolar}]
  Family Hydrophobe
  Weights 1.0
EndFeature

# hydrophobic atom
AtomType Hphobe [c,s,S&H0&v2,Br,I,{Carbon_NonPolar}]
AtomType RingHphobe [R;{Hphobe}]

# nitro groups in the RD code are always: *-[N+](=O)[O-]
DefineFeature Nitro2 [N;D3;+](=O)[O-]
  Family LumpedHydrophobe
  Weights 1.0,1.0,1.0
EndFeature

#
# 21.1.2008 (GL): update ring membership tests to reflect corrected meaning of
# "r" in SMARTS parser
#



AtomType Ring6 [r6,!R1&r5,!R1&r4,!R1&r3]
DefineFeature RH6_6 [{Ring6};{RingHphobe}]1[{Ring6};{RingHphobe}][{Ring6};{RingHphobe}][{Ring6};{RingHphobe}][{Ring6};{RingHphobe}][{Ring6};{RingHphobe}]1
  Family LumpedHydrophobe
  Weights 1.0,1.0,1.0,1.0,1.0,1.0
EndFeature

AtomType Ring5 [r5,!R1&r4,!R1&r3]
DefineFeature RH5_5 [{Ring5};{RingHphobe}]1[{Ring5};{RingHphobe}][{Ring5};{RingHphobe}][{Ring5};{RingHphobe}][{Ring5};{RingHphobe}]1
  Family LumpedHydrophobe
  Weights 1.0,1.0,1.0,1.0,1.0
EndFeature

AtomType Ring4 [r4,!R1&r3]
DefineFeature RH4_4 [{Ring4};{RingHphobe}]1[{Ring4};{RingHphobe}][{Ring4};{RingHphobe}][{Ring4};{RingHphobe}]1
  Family LumpedHydrophobe
  Weights 1.0,1.0,1.0,1.0
EndFeature

AtomType Ring3 [r3]
DefineFeature RH3_3 [{Ring3};{RingHphobe}]1[{Ring3};{RingHphobe}][{Ring3};{RingHphobe}]1
  Family LumpedHydrophobe
  Weights 1.0,1.0,1.0
EndFeature

DefineFeature tButyl [C;!R](-[CH3])(-[CH3])-[CH3]
  Family LumpedHydrophobe
  Weights 1.0,0.0,0.0,0.0
EndFeature

DefineFeature iPropyl [CH;!R](-[CH3])-[CH3]
  Family LumpedHydrophobe
  Weights 1.0,1.0,1.0
EndFeature

## Aliphatic rings
AtomType AliphR4 [A;r4]
DefineFeature Aliph4 [{AliphR4}]1~[{AliphR4}]~[{AliphR4}]~[{AliphR4}]~1
 Family Aliphatic
 Weights 1.0,1.0,1.0,1.0
EndFeature
AtomType AliphR5 [A;r5]
DefineFeature Aliph [{AliphR5}]1~[{AliphR5}]~[{AliphR5}]~[{AliphR5}]~[{AliphR5}]~1
 Family Aliphatic
 Weights 1.0,1.0,1.0,1.0,1.0
EndFeature
AtomType AliphR6 [A;r6,!R1&r5,!R1&r4,!R1&r3]
DefineFeature Aliph6 [{AliphR6}]1~[{AliphR6}]~[{AliphR6}]~[{AliphR6}]~[{AliphR6}]~[{AliphR6}]~1
 Family Aliphatic
 Weights 1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
AtomType AliphR7 [A;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]
DefineFeature Aliph7 [{AliphR7}]1~[{AliphR7}]~[{AliphR7}]~[{AliphR7}]~[{AliphR7}]~[{AliphR7}]~[{AliphR7}]~1
 Family Aliphatic
 Weights 1.0,1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
AtomType AliphR8 [A;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]
DefineFeature Aliph8 [{AliphR8}]1~[{AliphR8}]~[{AliphR8}]~[{AliphR8}]~[{AliphR8}]~[{AliphR8}]~[{AliphR8}]~[{AliphR8}]~1
 Family Aliphatic
 Weights 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
"""

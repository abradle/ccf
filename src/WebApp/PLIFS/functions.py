from PLIFS.models import PlifProbe, PlifProbeBit
from IOhandle.models import Target, Molecule, Protein
from Pharmacophore.models import PharmaPoint
from WONKA.models import PharmaScore, Residue
from WONKA.functions import get_all_coords
from Viewer.models import eucl_dist
from rdkit import Chem
import math, os, ast, sys
from LLOOMMPPAA.models import LLConf
from rdkit.Chem import rdGeometry
from rdkit.Chem import AllChem
from math import exp
from django.core.exceptions import ValidationError
# Make the pharmacophore features for a molecule
from rdkit.Chem import ChemicalFeatures
from PLIFS.models import PlifProtein, PlifProbe, PlifProbeBit, InteractionScheme, Interaction
from django.db import reset_queries


def find_interactions(target_id, t_count, tot_prots, scheme_name="default", mol_ids=None, react_id=None):
    """Function to take all the potential interactions and score them based on what is allowed"""
    # First get the interactions
    print "FINDING INTERACTIONS"
    if mol_ids:
        tot = PlifProbeBit.objects.filter(target_id=target_id, mol_id__in=mol_ids).count
        plif_ints = PlifProbeBit.objects.filter(target_id=target_id, mol_id__in=mol_ids).iterator()
    if react_id:
        mol_ids = LLConf.objects.filter(rmsd__gte=0.45, reactionprocess=react_id).values_list("mol_id__pk", flat=True)
        print "FOUND MOL IDS"
        tot = PlifProbeBit.objects.filter(target_id=target_id, mol_id__pk__in=mol_ids).count()
        plif_ints = PlifProbeBit.objects.filter(target_id=target_id, mol_id__pk__in=mol_ids).iterator()
    else:
        tot = PlifProbeBit.objects.filter(target_id=target_id).count()
        plif_ints = PlifProbeBit.objects.filter(target_id=target_id).iterator()
    print "FOUND PLIF INTERACTIONS"
    # Second get the interaction dict
    my_scheme = InteractionScheme.objects.get(scheme_name=scheme_name)
    scheme = ast.literal_eval(my_scheme.scheme_dict)
    # Third go through interactions of this type
    old = -1
    if tot == 0:
        target = Target.objects.get(pk=target_id)
        print "\nNo compounds for target: ", target.title
        return
    tot_ints = []
    tot_counter = 0
    counter = -1
    for interaction in plif_ints:
        counter += 1
        if counter * 100 / tot != old:
            old = counter * 100 / tot
            sys.stdout.write("\rChecking %s interactions for target %d of %d: %d%% complete..." % (scheme_name, t_count + 1, tot_prots, old))
            sys.stdout.flush()
            if react_id:
                # Now update this
                react_id.stage_completion = float(old) / float(react_id.tot_sch - react_id.sch + 1)
                react_id.save()
        if interaction.type in scheme:
            pass
        else:
            continue
        # Fourth - if allowed - determine the interactions present
        my_dist = scheme[interaction.type]["dist"]
        if interaction.dist > my_dist:
            continue
        # For the donor acceptor interaction
        if interaction.type == "SingleAtomAcceptor_SingleAtomDonor" and scheme[interaction.type]["angle_1"]:
            angle_1_min = int(scheme[interaction.type]["angle_1"].split(":")[0])
            angle_1_max = int(scheme[interaction.type]["angle_1"].split(":")[1])
            if interaction.angle_1 < angle_1_max and interaction.angle_1 > angle_1_min:
                pass
            # Check the second angle
            elif interaction.angle_2 < angle_1_max and interaction.angle_2 > angle_1_min:
                pass
            # Check the third angle (e.g. protonated lysine)
            elif interaction.angle_3 < angle_1_max and interaction.angle_3 > angle_1_min:
                pass
            else:
                continue
        else:
            # If  there is an angle constraint
            if scheme[interaction.type]["angle_1"] != None:
                angle_1_min = int(scheme[interaction.type]["angle_1"].split(":")[0])
                angle_1_max = int(scheme[interaction.type]["angle_1"].split(":")[1])
                # Now do the test
                if interaction.angle_1:
                    if interaction.angle_1 < angle_1_max and interaction.angle_1 > angle_1_min:
                        pass
                    else:
                        continue
            if scheme[interaction.type]["angle_2"] != None:
                angle_2_min = int(scheme[interaction.type]["angle_2"].split(":")[0])
                angle_2_max = int(scheme[interaction.type]["angle_2"].split(":")[1])
                # Now do the test
                if interaction.angle_2:
                    if interaction.angle_2 < angle_2_max and interaction.angle_2 > angle_2_min:
                        pass
                    else:
                        continue
        # Save it
        my_int = Interaction()
        my_int.scheme_id = my_scheme
        my_int.interaction_id = interaction
        my_int.target_id = interaction.target_id
        my_int.mol_id = interaction.mol_id
        my_int.prot_id_id = interaction.prot_id_id
        try:
            my_int.validate_unique()
            my_int.save()
        except ValidationError:
            continue
        #tot_ints.append(my_int)
        #tot_counter += 1
        #if tot_counter == 5000:
        #    Interaction.objects.bulk_create(tot_ints)
        #    tot_counter = 0
        #    tot_ints = []
        #    reset_queries()
    # Now print the final one
    Interaction.objects.bulk_create(tot_ints)
    reset_queries
    old = 100
    sys.stdout.write("\rChecking %s interactions for target %d of %d: %d%% complete..." % (scheme_name, t_count + 1, tot_prots, old))
    sys.stdout.flush()


def give_p_and_theta(ph4, probe):
    """Function to find p and theta values
    Input is two three by three tuples. First element of each is the ring centre"""
    plane_one_point_1 = rdGeometry.Point3D()
    plane_one_point_2 = rdGeometry.Point3D()
    plane_one_point_3 = rdGeometry.Point3D()
    plane_two_point_1 = rdGeometry.Point3D()
    plane_two_point_2 = rdGeometry.Point3D()
    plane_two_point_3 = rdGeometry.Point3D()
    # Point one
    plane_one_point_1.x = ph4.x_com
    plane_one_point_1.y = ph4.y_com
    plane_one_point_1.z = ph4.z_com
    # Point one
    plane_two_point_1.x = probe.x_com
    plane_two_point_1.y = probe.y_com
    plane_two_point_1.z = probe.z_com
    ### NEED TO DEFINE A PLANE
    # Point two
    plane_one_point_2.x = ph4.x_extent
    plane_one_point_2.y = ph4.y_extent
    plane_one_point_2.z = ph4.z_extent
    # Point two
    plane_two_point_2.x = probe.x_extent
    plane_two_point_2.y = probe.y_extent
    plane_two_point_2.z = probe.z_extent
    # Point three
    plane_one_point_3.x = probe.x_extent_2
    plane_one_point_3.y = probe.y_extent_2
    plane_one_point_3.z = probe.z_extent_2
    # Point three
    plane_two_point_3.x = probe.x_extent_2
    plane_two_point_3.y = probe.y_extent_2
    plane_two_point_3.z = probe.z_extent_2
    # Make the first plane normal
    plane_one_vector_1 = plane_one_point_1.DirectionVector(plane_one_point_2)
    plane_one_vector_2 = plane_one_point_1.DirectionVector(plane_one_point_3)
    norm_plane_one = plane_one_vector_1.CrossProduct(plane_one_vector_2)
    norm_plane_one.Normalize()
    # Make the second plane normal
    plane_two_vector_1 = plane_two_point_1.DirectionVector(plane_two_point_2)
    plane_two_vector_2 = plane_two_point_1.DirectionVector(plane_two_point_3)
    norm_plane_two = plane_two_vector_1.CrossProduct(plane_two_vector_2)
    norm_plane_two.Normalize()
    # Define P - the angle between the normals
    angle_1 = math.degrees(norm_plane_one.AngleTo(norm_plane_two))
    # Define THETA - the angle from the normal to the second vector
    # Define the vector between the ring centres
    centre_of_mass_vect = plane_one_point_1.DirectionVector(plane_two_point_1)
    centre_of_mass_vect.Normalize()
    # Define angle between this vector and the plane
    angle_2 = math.degrees(norm_plane_two.AngleTo(centre_of_mass_vect))
    if angle_1 > 90:
        # Bring it down
        angle_1 = 180 - angle_1
    if angle_2 > 90:
        angle_2 = 180 - angle_2
    # Now define it
    if angle_1 < 30:
        if angle_2 < 30:
            definition = "ff"
        elif angle_2 < 60:
            definition = "ft"
        elif angle_2 <= 90:
            definition = "fe"
        else:
            print "ERROR", angle_2
            definition = "ERROR"
    elif angle_1 < 60:
        if angle_2 < 30:
            definition = "of"
        elif angle_2 < 60:
            definition = "ot"
        elif angle_2 <= 90:
            definition = "oe"
        else:
            print "ERROR", angle_2
            definition = "ERROR"
    elif angle_1 <= 90:
        if angle_2 < 30:
            definition = "ef"
        elif angle_2 < 60:
            definition = "et"
        elif angle_2 <= 90:
            definition = "ee"
        else:
            print "ERROR", angle_2
            definition = "ERROR"
    else:
        print "ERROR", angle_1
        definition = "ERROR"
    return angle_1, angle_2, definition


def calc_int(ph4, probe, max_dist, my_ph4s, target):
    """Function to caluclate and store and in interaction"""
    # This is the centroid distance for aromatics and the Donor -> Acceptor distance for H-bonds
    dist = eucl_dist(ph4, probe)
    if dist > max_dist:
        return
    pb = PlifProbeBit()
    pb.probe_source_id = ph4
    pb.probe_dest_id = probe
    pb.prot_id = probe.prot_id
    pb.mol_id = ph4.mol_id
    pb.definition = pb.type
    try:
        pb.validate_unique()
    except ValidationError:
        return
    if probe.type in my_ph4s["Hydrophobic"]:
        pb.type = "Hydrophobe_Hydrophobe"
    elif ph4.type in my_ph4s["MagicSix"]:
        pb.type = "Magic"
    else:
        my_l = [ph4.type, probe.type]
        my_l.sort()
        pb.type = "_".join(my_l)
    # Check for validation
    pb.target_id = target
    pb.mol_id = ph4.mol_id
    pb.probe_id = probe
    pb.dist = dist
    # Now work out the angles if needed
    point_1 = rdGeometry.Point3D()
    point_1.x = ph4.x_com
    point_1.y = ph4.y_com
    point_1.z = ph4.z_com
    point_2 = rdGeometry.Point3D()
    point_2.x = probe.x_com
    point_2.y = probe.y_com
    point_2.z = probe.z_com
    # Donor -> acceptor
    if probe.x_extent and not ph4.x_extent:
        point_3 = rdGeometry.Point3D()
        point_3.x = probe.x_extent
        point_3.y = probe.y_extent
        point_3.z = probe.z_extent
        probe_vect = point_3.DirectionVector(point_2)
        other_vect = point_3.DirectionVector(point_1)
        pb.angle_1 = math.degrees(probe_vect.AngleTo(other_vect))
    # Acceptor -> donor
    elif not probe.x_extent and ph4.x_extent:
        point_3 = rdGeometry.Point3D()
        point_3.x = ph4.x_extent
        point_3.y = ph4.y_extent
        point_3.z = ph4.z_extent
        ph4_vect = point_3.DirectionVector(point_1)
        other_vect = point_3.DirectionVector(point_2)
        pb.angle_1 = math.degrees(ph4_vect.AngleTo(other_vect))
    # Donor -> acceptor if two attachment points
    if probe.x_extent_2 and not ph4.x_extent:
        point_3 = rdGeometry.Point3D()
        point_3.x = probe.x_extent_2
        point_3.y = probe.y_extent_2
        point_3.z = probe.z_extent_2
        probe_vect = point_3.DirectionVector(point_2)
        other_vect = point_3.DirectionVector(point_1)
        pb.angle_2 = math.degrees(probe_vect.AngleTo(other_vect))
    # Acceptor -> donor
    elif not probe.x_extent and ph4.x_extent_2:
        point_3 = rdGeometry.Point3D()
        point_3.x = ph4.x_extent_2
        point_3.y = ph4.y_extent_2
        point_3.z = ph4.z_extent_2
        ph4_vect = point_3.DirectionVector(point_1)
        other_vect = point_3.DirectionVector(point_2)
        pb.angle_2 = math.degrees(ph4_vect.AngleTo(other_vect))
    # Donor -> acceptor if three attachment points
    if probe.x_extent_3 and not ph4.x_extent:
        point_3 = rdGeometry.Point3D()
        point_3.x = probe.x_extent_3
        point_3.y = probe.y_extent_3
        point_3.z = probe.z_extent_3
        probe_vect = point_3.DirectionVector(point_2)
        other_vect = point_3.DirectionVector(point_1)
        pb.angle_3 = math.degrees(probe_vect.AngleTo(other_vect))
    # Acceptor -> donor
    elif not probe.x_extent and ph4.x_extent_3:
        point_3 = rdGeometry.Point3D()
        point_3.x = ph4.x_extent_3
        point_3.y = ph4.y_extent_3
        point_3.z = ph4.z_extent_3
        ph4_vect = point_3.DirectionVector(point_1)
        other_vect = point_3.DirectionVector(point_2)
        pb.angle_3 = math.degrees(ph4_vect.AngleTo(other_vect))
    # Finally if they both have an extent it's a PI-PI interaction or a carbonyl pi interaction
    if probe.x_extent and ph4.x_extent:
        if "Carbonyl" in [probe.type, ph4.type]:
          ### Need to work out the carbonyl pi interaction
            pass
        else:
            angles = give_p_and_theta(ph4, probe)
            pb.angle_1 = angles[0]
            pb.angle_2 = angles[1]
    return pb


def score_point(target_point, score_points):
    """Function to score a given ph4 point.
Takes a target point as input and a list of score points"""
    out_score = 0
    for score_p in score_points:
        out_score += exp(-eucl_dist(target_point, score_p))
    return out_score


def create_ph4_scores(target_id, opt="XTAL"):
    """Creates a dict of scores for each pharma type for a target. If a point doesn't exist
    it then gives None. If it does it gives a ranked list of scores for each.
    Takes a target_id and an optional option (just XTAL or OOMMPPAA too)"""
    target = Target.objects.get(pk=target_id)
    tot_mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__contains=target.title)
    tot_set = PharmaPoint.objects.filter(mol_id__in=tot_mols)
    # Evaluate the set
    len(tot_set)
    # Now score the mols
    len_mols = len(tot_mols)
    for i, m in enumerate(tot_mols):
        print (i + 1) * 100 / len_mols
        # Set up the score set (set of ph4s to score against)
        score_ph4s = tot_set.exclude(mol_id=m)
        # And the test set (set of ph4s to test)
        test_ph4s = tot_set.filter(mol_id=m)
        # Now test these
        for ph4 in test_ph4s:
            ph4_score = PharmaScore()
            ph4_score.pharma_id = ph4
            ph4_score.type = opt
            try:
                ph4_score.validate_unique()
            except ValidationError:
                ph4_score = PharmaScore.objects.get(pharma_id=ph4, type="XTAL")
            # Now get the score for this one
            ph4_score.score = score_point(ph4, score_ph4s.filter(smiles=ph4.smiles))
            ph4_score.save()
            # Now return all the scores based on type and overall


def score_target(target_id, opt=None, clash=-1.0, rmsd=0.5, shape_dist=0.2, react_id=None):
    """Function to score actual crystal structures against all of the protein structures
    for a target."""
    target = Target.objects.get(pk=target_id)
    # List of pharamcohpores
    my_ph4s = {"Hydrophobic": ["RingTwoWayAttach", "ThreeWayAttach", "Hydrophobe",'iPropyl','tButyl','ChainTwoWayAttach','RH6_6','RH3_3','RH4_4','RH5_5'],
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
    # Contrast Ph4 this is a dict from the proteins perspective. Which groups can each feature interact with
    # Append the magic six to each
    contrast_ph4 = {"SingleAtomAcceptor": ["SingleAtomDonor", "BigHalogen", "SmallHalogen"], # H-bonding and Halogen bonding
                    "SingleAtomDonor": ["SingleAtomAcceptor", "SmallHalogen", "BigHalogen"],# H-bonding and Halogen bonding
                    "AcidicGroup": ["BasicGroup"], "BasicGroup": ["AcidicGroup"], # Acid-base
                    "Arom6": ["Arom6", "Arom5", "Carbonyl"], "Arom5": ["Arom6", "Arom5", "Carbonyl"],# Aromatic-aromatic interactions
                    #Hydrophobic interactions
                    "ChainTwoWayAttach": ["Hydrophobe", "ThreeWayAttach", "ChainTwoWayAttach", "RingTwoWayAttach", "BigHalogen", "SmallHalogen"],
                    "RingTwoWayAttach":  ["Hydrophobe", "ThreeWayAttach", "ChainTwoWayAttach", "RingTwoWayAttach", "BigHalogen", "SmallHalogen"],
                    "ThreeWayAttach":  ["Hydrophobe", "ThreeWayAttach", "ChainTwoWayAttach", "RingTwoWayAttach", "BigHalogen", "SmallHalogen"],
                    "Hydrophobe":  ["Hydrophobe", "ThreeWayAttach", "ChainTwoWayAttach", "RingTwoWayAttach", "BigHalogen", "SmallHalogen"]}# Hydrophobes
    # Ignore if above this dist
    max_dist = 6.0
    # First of all score the magic six against the residues
    # Get the pharmapoints - for the molecules
    if opt == None:
        tot_mols = Molecule.objects.filter(prot_id__target_id=target).exclude(prot_id__code__contains=target.title)
        tot_ph4s = PlifProbe.objects.filter(mol_id__in=tot_mols)
    elif opt == "LLOOMMPPAA":
        tot_ph4s = PlifProbe.objects.filter(mol_id__prot_id__code=target.title + "SYNTH")
    # Now find their proximity to the residues and store in the DB
    # Get the probes
    # Then score the molecule pharmacophores against the protein pharmacophores
    # These should be max 6 A away from a mol on 
    my_coords = [(plif.x_com, plif.y_com, plif.z_com) for plif in tot_ph4s]
    x_max = max([x[0] for x in my_coords])
    x_min = min([x[0] for x in my_coords])
    y_max = max([x[1] for x in my_coords])
    y_min = min([x[1] for x in my_coords])
    z_max = max([x[2] for x in my_coords])
    z_min = min([x[2] for x in my_coords])
    res_ph4 = PlifProbe.objects.filter(type__in=contrast_ph4, prot_id__isnull=False)
    extra = 6.0
    res_ph4 = res_ph4.filter(x_com__lte=x_max+extra, y_com__lte=y_max+extra, z_com__lte=z_max+extra)
    res_ph4 = res_ph4.filter(x_com__gte=x_min-extra, y_com__gte=y_min-extra, z_com__gte=z_min-extra)
    pp_list = []
    pp_counter = 0
    # Compare all against all
    if opt == None:
        old = -1
        tot = len(tot_mols)
        # Loop through the mols
        for counter, mol in enumerate(tot_mols):
            if counter * 100 / tot != old:
                old = counter * 100 / tot
                sys.stdout.write("\rMaking plif interactions %d%% complete..." % old)
                sys.stdout.flush()
            # Filter down to this mols PH4
            mol_ph4 = tot_ph4s.filter(mol_id=mol)
            # And to its own proteins PH4s
            my_res_ph4 = res_ph4.filter(prot_id=mol.prot_id)
            # Now loop through the residue probes
            for probe in my_res_ph4:
                # Get the potential counterparts probes
                my_mol_ph4s = mol_ph4.filter(type__in=contrast_ph4[probe.type])
                # Loop through these
                for ph4 in my_mol_ph4s:
                    # Calculate the interactions
                    my_pp = calc_int(ph4, probe, max_dist, my_ph4s, target)
                    if my_pp:
                        pp_list.append(my_pp)
                        pp_counter += 1
                        if pp_counter == 5000:
                            PlifProbeBit.objects.bulk_create(pp_list)
                            pp_list = []
                            pp_counter = 0
        PlifProbeBit.objects.bulk_create(pp_list)
        old = counter * 100 / tot
        sys.stdout.write("\rMaking plif interactions %d%% complete..." % old)
        sys.stdout.flush()
        print "\nPlif creation complete"
    # Just compare to the "NATIVE" structure
    elif opt == "LLOOMMPPAA":
        # Get the proteins that have LLOOMMPPAAs = doing nothing at the moment
        if Protein.objects.filter(code=target.title + "SYNTH"):
            # Get all the llconfs - within a certain RMSD, Strain and Clash
            if react_id:
                ll_confs = LLConf.objects.filter(target_id=target, clash__gte=clash, rmsd__gte=rmsd, shape_dist__lte=shape_dist, reactionprocess=react_id).order_by("pk")
                react_id.proc_stage = "FINDING INTERACTIONS"
                react_id.save()
            else:
                ll_confs = LLConf.objects.filter(target_id=target, clash__gte=clash, rmsd__gte=rmsd, shape_dist__lte=shape_dist)
            tot_confs = ll_confs.count()
            ll_confs = ll_confs.iterator()
            #  Loop through them
            done_list = []
            for conf_c, ll_conf in enumerate(ll_confs):
                if ll_conf.mol_id.pk in done_list:
                    continue
                if conf_c < react_id.stage_completion * tot_confs / 100:
                    continue
                if react_id:
                    react_id.stage_completion = int((float(conf_c) / float(tot_confs)) * 100)
                    react_id.save()
                done_list.append(ll_conf.mol_id.pk)
#### IF WE HAVE A MOL2 PROTEIN FOR THE NATIVE PROTEIN USE THAT 
                if res_ph4.filter(prot_id=ll_conf.mol_ref.prot_id, res_id__res_name="MOL2"):
                    print "GETTING MOL2 PLIFS"
                    my_ph4 = res_ph4.filter(prot_id=ll_conf.mol_ref.prot_id, res_id__res_name="MOL2").iterator()
                    tot = res_ph4.filter(prot_id=ll_conf.mol_ref.prot_id, res_id__res_name="MOL2").count()
#### IF NOT JUST USE THE STANDARD ONES
                else:
                    my_ph4 = res_ph4.filter(prot_id=ll_conf.mol_ref.prot_id).iterator()
                    tot = res_ph4.filter(prot_id=ll_conf.mol_ref.prot_id).count()
                if react_id:
                    react_id.proc_stage = "FINDING INTERACTIONS"
                    react_id.stage_completion = conf_c * 100 / tot_confs
                    react_id.save()
                old = -1
                for i, probe in enumerate(my_ph4):
                    # Now filter on type - and the protein needs to be the same
                    # Do the percent clock
                    perc_comp = int((float(i) / float(tot)) * 100)
                    if perc_comp != old:
                        sys.stdout.write("\rMaking plif interactions for conf %d of %d: %d%% complete..." % (conf_c, tot_confs, old))
                        sys.stdout.flush()
                        old = perc_comp
                    mol_ph4 = tot_ph4s.filter(type__in=contrast_ph4[probe.type], mol_id=ll_conf.mol_id)
                    for ph4 in mol_ph4:
                        my_pp = calc_int(ph4, probe, max_dist, my_ph4s, target)
                        if my_pp:
                            pp_list.append(my_pp)
                            pp_counter += 1
                            if pp_counter == 1000:
                                PlifProbeBit.objects.bulk_create(pp_list)
                                pp_list = []
                                pp_counter = 0
                                reset_queries()
            PlifProbeBit.objects.bulk_create(pp_list)
            sys.stdout.write("\rMaking plif interactions complete...")
            sys.stdout.flush()
            react_id.proc_stage = "MAKE INT SCHEMES"
            react_id.stage_completion = -1
            react_id.save()
            print "\nPlif creation complete"


def find_attached_h(atom_id, rdmol, my_mol=None):
    """Function to take in an Atom ID and a molecule
    and return the attached H"""
    # Get all the bonds
    bonds = rdmol.GetAtomWithIdx(atom_id).GetBonds()
    # List for the output to go
    attached_hs = []
    met_c = 0
    his_c = 0
    for bond in bonds:
        end_id = bond.GetEndAtomIdx()
        start_id = bond.GetBeginAtomIdx()
        # Check it's not the start atom
        if end_id != atom_id:
            if bond.GetEndAtom().GetSmarts() == "[H]":
                attached_hs.append(bond.GetEndAtom())
            elif bond.GetEndAtom().GetSmarts() == "C":
                met_c += 1
            elif bond.GetEndAtom().GetSmarts() == "c":
                his_c += 1
        elif start_id != atom_id:
            if bond.GetBeginAtom().GetSmarts() == "[H]":
                attached_hs.append(bond.GetBeginAtom())
            elif bond.GetBeginAtom().GetSmarts() == "C":
                met_c += 1
            elif bond.GetBeginAtom().GetSmarts() == "c":
                his_c += 1
    if attached_hs: 
        return attached_hs
    elif met_c == 3:
        #print "METHYLATED AMINE"
        return "METHYLATED"
    elif his_c == 2:
        #print "WRONG HISTIDINE"
        return "HISTIDINE"
    elif my_mol:
        return my_mol
    else:
        print atom_id


def assign_extents(new_fragment, conf, feat, rdmol, opt=None):
    """Function to assign x, y and z extents to a pharmacophore feature"""
    # Now do the x_extent for others
        # If it is a donor or an aromatic ring - get another atom(s) to give vector(s)
    if new_fragment.type in ["Carbonyl"]:
        attached = feat.GetAtomIds()[0]
        attached_2 = feat.GetAtomIds()[1]
        # Store both sides of the carbonyl
        if rdmol.GetAtomWithIdx(attached).GetAtomicNum() == 6:
            c_id = attached
            o_id = attached_2
        # Assign the co-ordinates of the C
        new_fragment.x_com = conf.GetAtomPosition(c_id).x
        new_fragment.y_com = conf.GetAtomPosition(c_id).y
        new_fragment.z_com = conf.GetAtomPosition(c_id).z
        # Adssign the coordinates of the o
        new_fragment.x_extent = conf.GetAtomPosition(o_id).x
        new_fragment.y_extent = conf.GetAtomPosition(o_id).y
        new_fragment.z_extent = conf.GetAtomPosition(o_id).z
    elif new_fragment.type in ["SingleAtomDonor", "WeakDonor"]:
        attached_2 = find_attached_h(feat.GetAtomIds()[0], rdmol, opt)
        if attached_2 in ["METHYLATED", "HISTIDINE", "MOL"]:
            # Skip these as they have no free H
            return None
        if not attached_2:
            print feat.GetAtomIds()
            return rdmol
        attached = attached_2[0]
        new_fragment.x_extent = conf.GetAtomPosition(attached.GetIdx()).x
        new_fragment.y_extent = conf.GetAtomPosition(attached.GetIdx()).y
        new_fragment.z_extent = conf.GetAtomPosition(attached.GetIdx()).z
        if len(attached_2) > 1:
            attached = attached_2[1]
            new_fragment.x_extent_2 = conf.GetAtomPosition(attached.GetIdx()).x
            new_fragment.y_extent_2 = conf.GetAtomPosition(attached.GetIdx()).y
            new_fragment.z_extent_2 = conf.GetAtomPosition(attached.GetIdx()).z
        if len(attached_2) > 2:
            attached = attached_2[2]
            new_fragment.x_extent_3 = conf.GetAtomPosition(attached.GetIdx()).x
            new_fragment.y_extent_3 = conf.GetAtomPosition(attached.GetIdx()).y
            new_fragment.z_extent_3 = conf.GetAtomPosition(attached.GetIdx()).z
    elif new_fragment.type in ["Arom5", "Arom6"]:
        # Store these ring positions
        attached = feat.GetAtomIds()[0]
        attached_2 = feat.GetAtomIds()[1]
        new_fragment.x_extent = conf.GetAtomPosition(attached).x
        new_fragment.y_extent = conf.GetAtomPosition(attached).y
        new_fragment.z_extent = conf.GetAtomPosition(attached).z
        new_fragment.x_extent_2 = conf.GetAtomPosition(attached_2).x
        new_fragment.y_extent_2 = conf.GetAtomPosition(attached_2).y
        new_fragment.z_extent_2 = conf.GetAtomPosition(attached_2).z
    return new_fragment


def init_probe(feat, target, mol_id=None, res_id=None, prot_id=None):
    """Function to initiate a PlifProbe"""
    new_fragment = PlifProbe()
    new_fragment.type = str(feat.GetType())
    if mol_id:
        new_fragment.mol_id = mol_id
        new_fragment.unique_id = "MOL" + str(mol_id.pk)
    elif res_id:
        new_fragment.prot_id = prot_id
        new_fragment.unique_id = "RES" + str(res_id.pk)
        new_fragment.res_id
    # Now set the target
    new_fragment.target_id = target
    # Set the coords
    coord = feat.GetPos()
    # Assign the co-ordinates
    new_fragment.x_com = coord.x
    new_fragment.y_com = coord.y
    new_fragment.z_com = coord.z
    return new_fragment


def make_pharma_probes(target_id, opt=None, clash=-1.0, rmsd=0.5, shape_dist=0.2, react_id=None):
    """Function to make a protein pharamcophore"""
    target = Target.objects.get(pk=target_id)
    tot_res = Residue.objects.filter(target_id=target_id)
    if not opt:
        # Do the proteins - but only a limited distance
        # Make the limited description for these proteins
        if react_id:
            react_id.proc_stage = "MAKING PROTEIN PROBES"
            react_id.stage_completion = 0
            react_id.save()
        factory = ChemicalFeatures.BuildFeatureFactoryFromString(prot_feats)
        prots = Protein.objects.filter(target_id=target, pdb_info__isnull=False).exclude(code__contains=target.title)
        tot_prots = len(prots)
        # Loop through the mols
        prot_old = -1
#        for counter, p in enumerate(prots):
#            # Trim the coords to X A from the ligand
#            delta = 10.0
#            prot_res = []
#            if react_id:
#                my_prg = int((float(counter) / float(tot_prots)) * 100)
#                if my_prg != prot_old:
#                    react_id.stage_completion = my_prg
#                    prot_old = my_prg
#                    react_id.save()
#            for mol in Molecule.objects.filter(prot_id=p):
#                rdmol = Chem.MolFromMolBlock(str(mol.sdf_info))
##                x_coords, y_coords, z_coords = get_all_coords(rdmol)
##                x_max = max(x_coords)
##                y_max = max(y_coords)
##                z_max = max(z_coords)
##                x_min = min(x_coords)
##                y_min = min(y_coords)
##                z_min = min(z_coords)
##                near_res = tot_res.filter(prot_id=p).exclude(x_max__gte=x_max + delta, x_min__gte=x_min + delta).exclude(x_max__lte=x_max - delta, x_min__lte=x_min - delta)
##                near_res = near_res.exclude(y_max__gte=y_max + delta, y_min__gte=y_min + delta).exclude(y_max__lte=y_max - delta, y_min__lte=y_min - delta)
##                near_res = near_res.exclude(z_max__gte=z_max + delta, z_min__gte=z_min + delta).exclude(z_max__lte=z_max - delta, z_min__lte=z_min - delta)
##                prot_res.extend(near_res)
#            old = -1
#            tot = len(prot_res)
#            #for res_count, res in enumerate(prot_res):
        for p in prots:
            print p.code
            ## Mol2 required for hydrogens currently
            if PlifProtein.objects.filter(prot_id=p):
                rdmol = Chem.MolFromMol2Block(str(PlifProtein.objects.filter(prot_id=p)[0].mol2_data), removeHs=False)
            else:
                rdmol = Chem.MolFromPDBBlock(str(p.pdb_info))
                # Add the Hs => but this is pretty ruddy dodgy
                rdmol = AllChem.AddHs(rdmol, addCoords=True)
            if rdmol is None:
                print "ERROR WITH PROTEIN ", p.pk
                continue

            conf = rdmol.GetConformer()
            feats = factory.GetFeaturesForMol(rdmol)
            for i, feat in enumerate(feats):
                # Assign the residue
                info = rdmol.GetAtomWithIdx(feat.GetAtomIds()[0]).GetPDBResidueInfo()
                if info:
                    res_num = info.GetResidueNumber()
                    res_name = info.GetResidueName()
                else:
                    res_num = i
                    res_name = "MOL2"
                my_res = Residue.objects.get_or_create(target_id=target, prot_id=p, res_name=res_name, res_num=res_num)[0]
                new_fragment = init_probe(feat, target, mol_id=None, res_id=my_res, prot_id=p)
                new_fragment.res_id = my_res
                new_fragment.unique_id = "RES" + str(my_res.pk)
                ### Now we've set x_com
                # Now validate this is in 1) a new place 2) a new feature 3) on a unique residue
                new_fragment = assign_extents(new_fragment, conf, feat, rdmol)
                if not new_fragment:
                    continue
                # Now assign the atom_ids
                new_fragment.atom_ids = str(feat.GetAtomIds())
                try:
                    new_fragment.validate_unique()
                except ValidationError:
                    new_fragment = PlifProbe.objects.get(x_com=new_fragment.x_com,y_com=new_fragment.y_com,z_com=new_fragment.z_com,unique_id=new_fragment.unique_id, type=new_fragment.type)
                    new_fragment.atom_ids = str(feat.GetAtomIds())
                    new_fragment.save(update_fields=["atom_ids"])
                    continue
                new_fragment.save()
        sys.stdout.write("\rMaking Protein probes for protein complete")
        sys.stdout.flush()
        print "Made protein PLIFS"
    # Different ph4s for these molecules
    fdefName = os.path.join(os.path.join(os.path.split(sys.argv[0])[0], 'data/media'), 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    # Get the molecules that relate to this
    if not opt:
        if react_id:
            react_id.proc_stage = "MAKING MOLECULE PROBES"
            react_id.stage_completion = 0
            react_id.save()
        mols = Molecule.objects.filter(prot_id__in=prots)
    elif opt == "LLOOMMPPAA":
        # Take  the molecules from LLOOMMPPAA for testing
        if react_id:
            ll_confs = LLConf.objects.filter(target_id=target, clash__gte=clash, rmsd__gte=rmsd, shape_dist__lte=shape_dist, reactionprocess=react_id)
            react_id.proc_stage = "MAKING LLOOMMPPAA PROBES"
            react_id.stage_completion = 0
            react_id.save()
        else:
            ll_confs = LLConf.objects.filter(target_id=target, clash__gte=clash, rmsd__gte=rmsd, shape_dist__lte=shape_dist)
        mols = [ll.mol_id for ll in ll_confs]
    old = -1
    tot = len(mols)
    # List and counter to trigger the bulk create
    new_counter = 0
    new_fragments = []
    done_list = []
    mol_old = -1
    print "Making probes for " + str(tot) + " molecules"
    for counter, m in enumerate(mols):
        if m.pk in done_list:
            continue
        done_list.append(m.pk)
        sys.stdout.write("\rMaking molecule probes %d of %d complete..." % (counter, tot))
        sys.stdout.flush()
        if react_id:
            my_prg = int((float(counter) / float(tot)) * 100)
            if my_prg != mol_old:
                react_id.stage_completion = my_prg
                mol_old = my_prg
                react_id.save()
        # Get the molecule
        orig_mol = Chem.MolFromMolBlock(str(m.sdf_info))
        # Add the Hydrogens
        rdmol = AllChem.AddHs(orig_mol, addCoords=True)
        # Get the conformer
        conf = rdmol.GetConformer()
        for feat in factory.GetFeaturesForMol(rdmol):
            # Initialise this new object
            new_fragment = init_probe(feat, target, mol_id=m)
            new_fragment = assign_extents(new_fragment, conf, feat, rdmol, "MOL")
            if not new_fragment:
                continue
            # Now assign the atom_ids
            new_fragment.atom_ids = str(feat.GetAtomIds())
            try:
                new_fragment.validate_unique()
                new_fragment.save()
            except ValidationError:
                continue
                
            new_fragments.append(new_fragment)
            new_counter += 1
            if new_counter == 5000:
                print "COMMITING"
                #PlifProbe.objects.bulk_create(new_fragments)
                new_fragments = []
                done_list = []
                new_counter = 0 
    #PlifProbe.objects.bulk_create(new_fragments)
    sys.stdout.write("\rMaking molecule probes %d%% complete..." % 100)
    sys.stdout.flush()
    print "Made molecule probes"
    if not opt and react_id:
        react_id.proc_stage = "MAKING LLOOMMPPAA PROBES"
        react_id.stage_completion = -1
        react_id.save()
    elif react_id:
        react_id.proc_stage = "FINDING INTERACTIONS"
        react_id.stage_completion = -1
        react_id.save()


def restore_probs(target_id, opt=None, clash=-1.0, rmsd=0.5, shape_dist=1.0, react_id=None):
    """Function to restore the pharmacophore probes"""
    # Destroy the probes
    #[p.delete() for p in PlifProbe.objects.filter(target_id=target_id)]
    if react_id:
        # Make the ligand and protein probes
        if react_id.proc_stage == "MAKING PROTEIN PROBES":
            make_pharma_probes(target_id, react_id=react_id)
        # Now make the ligand probes on the LLOOMMPPAA conformations
        if react_id.proc_stage == "MAKING LLOOMMPPAA PROBES":
            make_pharma_probes(target_id, "LLOOMMPPAA", clash=clash, rmsd=rmsd, shape_dist=shape_dist, react_id=react_id)
        # Now find the interactions
        if react_id.proc_stage == "FINDING INTERACTIONS":
            score_target(target_id, "LLOOMMPPAA", clash=clash, rmsd=rmsd, shape_dist=shape_dist, react_id=react_id)
    elif opt:
        # Make the ligand and protein probes
        make_pharma_probes(target_id)
        # Now make the ligand probes on the LLOOMMPPAA conformations
        make_pharma_probes(target_id, "LLOOMMPPAA", clash=clash, rmsd=rmsd, shape_dist=shape_dist)
        # Now find the interactions
        score_target(target_id, "LLOOMMPPAA", clash=clash, rmsd=rmsd, shape_dist=shape_dist)
        #make_json(Target.objects.get(pk=target_id))
    else:
        # Makes the PLIF probes
        make_pharma_probes(target_id)
        # Make the PLIFProbBits
        score_target(target_id)


def make_ints(target_id, opt=None):
    """Function to make the interactions for a given target
    opt - 'LLOOMMPPAA' or None"""
    from PLIFS.interaction_defs import init_schemes
    target = Target.objects.get(pk=target_id)
    # Recognise the OPT
    if opt is None:
        prob_opt = None
        mol_ids = Molecule.objects.filter(prot_id__target_id=target).exclude(prot_id__code__contains=target.title)
    elif opt == "LLOOMMPPAA":
        prob_opt = "LLOOMMPPAA"
        mol_ids = Molecule.objects.filter(prot_id__code=target.title + "SYNTH")
    elif opt == "ALL":
        prob_opt = "LLOOMMPPAA"
        mol_ids = Molecule.objects.filter(prot_id__code=target.title + "SYNTH") | Molecule.objects.filter(prot_id__target_id=target).exclude(prot_id__code__contains=target.title)
    else:
        print "Unrecognised opt ALL"
        return
    # Make the plifs and plifprobebits
    restore_probs(target_id, prob_opt)
    # Make the ints
    init_schemes(mol_ids=mol_ids, refresh=True)


prot_feats = """# $Id$
#
# RDKit base fdef file.
# Created by Greg Landrum
#

AtomType NDonor [N&!H0&v3,N&!H0&+1&v4,n&H1&+0]
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

# Hydroxyls and acids
AtomType CHalcacceptor [O,S;H1;+0]

DefineFeature SingleAtomAcceptor [{CHalcacceptor},{NAcceptor}]
  Family Acceptor
  Weights 1
EndFeature

# this one is delightfully easy:
DefineFeature AcidicGroup [C](=[O])-[O;H1,H0&-1]
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


# hydrophobic features
# any carbon that is not bonded to a polar atom is considered a hydrophobe
#
# 23.11.2007 (GL): match any bond (not just single bonds); add #6 at
#  beginning to make it more efficient
AtomType Carbon_Polar [#6;$([#6]~[#7,#8,#9])]
# 23.11.2007 (GL): don't match charged carbon
AtomType Carbon_NonPolar [#6;+0;!{Carbon_Polar}]
# hydrophobic atom


DefineFeature Hydrophobe [D1;{Carbon_NonPolar}]
  Family Hydrophobe
  Weights 1.0
EndFeature

DefineFeature RingTwoWayAttach [R;D2;{Carbon_NonPolar}]
  Family Hydrophobe
  Weights 1.0
EndFeature

DefineFeature ThreeWayAttach [D3,D4;{Carbon_NonPolar}]
  Family Hydrophobe
  Weights 1.0
EndFeature

DefineFeature ChainTwoWayAttach [R0;D2;{Carbon_NonPolar}]
  Family Hydrophobe
  Weights 1.0
EndFeature


"""

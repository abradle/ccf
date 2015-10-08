from Pharmacophore.models import PharmaPoint
from math import exp
from MMPMaker.models import MMPFrag
from Viewer.models import eucl_dist
from IOhandle.models import Target, Molecule, Protein, Compound
from WONKA.models import InternalIDLink, PharmaScore, FragScore, KeyCluster, ResShift, Water, Residue
from django.core.exceptions import ValidationError, FieldError
import bisect
from rdkit import Chem
from rdkit.Chem import AllChem, rdForceFieldHelpers 
import math,ast,os,sys
from OOMMPPAA.cluster import dpmeans
from IOhandle.functions import centre_of_mass_from_SD_block, add_new_comp
from rdkit.Chem import rdGeometry
import json, numpy, csv
from rdkit.Chem import ChemicalFeatures
from Pharmacophore.models import find_pharmacophore_points
from MMPMaker.functions import make_mmp_database, find_frag_com
from loading import load_activity_data
from IOhandle.functions import desalt_compound
from OBSERVATIONS.models import UserData


def make_oommppaa_clusters(target_id):
    """Function to make Clusters of Act3DFrags and ActPharmaPoints at a range of lamda values"""
    # Get the points
    target = Target.objects.get(pk=target_id)
    from MMPMaker.models import ActPharmaPoint, Act3DFrag
    apps = ActPharmaPoint.objects.filter(in_diff_15=True, target_id=target_id,)
    # Get the types
    app_types = list(set([x.pharma_id.smiles for x in apps]))
    # Get the types - 
    a3d = Act3DFrag.objects.filter(act_id__target_id=target_id).exclude(mmp_frag_id__smistore="*[H]")
    # Go through these and give them the correct cofm (if possible)
    new_l = []
    for a in a3d:
        ms = MMPFrag.objects.filter(mmp_link=a.mmp_frag_id.mmp_link,mol_id=a.over_mol_id)
        if len(ms) == 0:
            print "NONE"
            a3d = a3d.exclude(a)
        elif len(ms) >1:
            print "DUPLICATES"
            a3d = a3d.exclude(a)
        else:
            new_m = a.mmp_frag_id
            new_m.x_com = ms[0].x_com
            new_m.y_com = ms[0].y_com
            new_m.y_com = ms[0].z_com
            new_m.save()
    a3d_types = list(set([x.mmp_frag_id.smistore for x in a3d]))
    # Loop through the lamda values
    lams = [float("%.2f" % x) for x in numpy.linspace(0.0, 3.0, 31)]
    # Make the clusters
    # Do the lam
    for lam in lams:
        print lam
        make_clusters(apps, app_types, lam, target)
        make_clusters(a3d, a3d_types, lam, target)


def return_points(rdmol, feat):
    #Returns the three positions required
    from rdkit.Chem import ChemicalFeatures as CF
    import sys
    import os
    # Make the limited description for these proteins
    conf = rdmol.GetConformer()
    pos_1 = feat.GetPos()
    pos_2 = conf.GetAtomPosition(feat.GetAtomIds()[0])
    pos_3 = conf.GetAtomPosition(feat.GetAtomIds()[1])
    return [pos_1, pos_2, pos_3]


def get_min_and_max(mols):
    """Function to get the min and max coords for a set of mols"""
    # Loop through the mols
    x_coords = []
    y_coords = []
    z_coords = []
    for rdmol in mols:
        # Gives the atoms
        atoms = rdmol.GetAtoms()
        conf = rdmol.GetConformer()
        numatoms = 0.0
        # Assume all heavy atoms have the same mass
        for atom in atoms:
            numatoms += 1.0
            coords = conf.GetAtomPosition(atom.GetIdx())
            if coords.x == 0.0 and coords.y == 0.0 and coords.z == 0.0:
                print atom.GetSmarts()
                sys.exit()
            x_coords.append(float(coords.x))
            y_coords.append(float(coords.y))
            z_coords.append(float(coords.z))
    return min(x_coords), max(x_coords), min(y_coords), max(y_coords), min(z_coords), max(z_coords)


def make_clusters(object_list, type_list, lam, target_id):
    """Function to make clusters from object lists"""
    my_method = "DPMEANS"
    my_types = set([str(x) for x in object_list])
    if len(my_types) > 1:
        print "Error - multiple types in cluster list"
        print my_types
        raise TypeError
    elif len(my_types) == 0:
        print "No types..."
        return
    else:
        my_type = list(my_types)[0]
    # Set up the counters
    old = -1
    tot = len(type_list)
    for counter, t in enumerate(type_list):
        if counter * 100 / tot != old:
            old = counter * 100 / tot
            sys.stdout.write("\rProcessing clusters %d%% complete..." % old)
        # Get the points and set them all to size =0
        existing_clusts = KeyCluster.objects.filter(lam=lam, function=t, method=my_method, type=my_type, target_id=target_id)
        # Set them all to size zero
        for kc in existing_clusts:
            kc.size = 0
            kc.save()
        if my_type == "Water object":
            myps = object_list
        elif my_type == "PharmaPoint object":
            myps = object_list.filter(smiles=t)
        elif my_type == "Act3DFrag object":
            myps = object_list.filter(mmp_frag_id__smistore=t)
        elif my_type == "ActPharmaPoint object":
            myps = object_list.filter(pharma_id__smiles=t)
        elif my_type == "MMPFrag object":
            myps = object_list.filter(smistore=t)
        elif my_type == "Molecule object":
            myps = object_list
        # Get the clusters
        # Fix the order
        len(myps)
        if my_type == "ActPharmaPoint object":
            com_list = [(x.pharma_id.x_com, x.pharma_id.y_com, x.pharma_id.z_com) for x in myps]
        elif my_type == "Act3DFrag object":
            com_list = [(x.mmp_frag_id.x_com, x.mmp_frag_id.y_com, x.mmp_frag_id.z_com) for x in myps]
        else:
            com_list = [(x.x_com, x.y_com, x.z_com) for x in myps]
        dp = dpmeans(com_list, lam, 0, False)
        dp.run()
        # Now assign the points in the lists to their clusters
        clusters = []
        # Now assign the clusters
        for i, cluster in enumerate(dp.clusters):
            # Try and find one within 0.2 A of this one
            m = KeyCluster()
            m.x_com = dp.clusters[cluster][0]
            m.y_com = dp.clusters[cluster][1]
            m.z_com = dp.clusters[cluster][2]
            # Try and find other objects that have this same cluster (half the lambda)
            my_matches = [kc for kc in existing_clusts.filter(lam=lam, method=my_method, type=my_type, function=t, target_id=target_id) if eucl_dist(m, kc) < float(lam) / 2.0]
            if my_matches:
                if len(my_matches) > 1:
                    print "MULTIPLE OPTIONS ERROR.....KeyCluster making"
                m = my_matches[0]
            # If not make a new one
            else:
                m.lam = lam
                m.method = my_method
                m.type = my_type
                m.function = t
                m.target_id = target_id
            # Now update the information appropriately
            m.cluster = i
            m.size = 0
            # Now save this point
            try:
                m.validate_unique()
                m.save()
            except ValidationError:
                # Update the centre of mass for each
                m = KeyCluster.objects.get(lam=lam, function=t, method=my_method, type=my_type, cluster=i, target_id=target_id)
                m.x_com = dp.clusters[cluster][0]
                m.y_com = dp.clusters[cluster][1]
                m.z_com = dp.clusters[cluster][2]
                m.cluster = i
                m.size = 0
                m.save()
            # Now rezero this guy
            if  my_type == "MMPFrag object":
                [m.mmp_frag_id.remove(x) for x in m.mmp_frag_id.iterator()]
            elif  my_type == "PharmaPoint object":
                [m.pharma_id.remove(x) for x in m.pharma_id.iterator()]
            elif my_type == "Water object":
                [m.water_id.remove(x) for x in m.water_id.iterator()]
            elif my_type == "Molecule object":
                [m.mol_id.remove(x) for x in m.mol_id.iterator()]
            elif my_type == "Act3DFrag object":
                [m.act3dfrag_id.remove(x) for x in m.act3dfrag_id.iterator()]
            elif my_type == "ActPharmaPoint object":
                [m.actpharma_id.remove(x) for x in m.actpharma_id.iterator()]
# Now zero the foreign keys
            clusters.append(m)
        for i, val in enumerate(dp.dataClusterId):
            my_clust = [x for x in clusters if x.cluster == val][0]
            if  my_clust.type == "MMPFrag object":
                my_clust.mmp_frag_id.add(myps[i])
                my_clust.frag_size = myps[i].core_size
            elif  my_clust.type == "PharmaPoint object":
                my_clust.pharma_id.add(myps[i])
            elif my_clust.type == "Water object":
                my_clust.water_id.add(myps[i])
            elif my_clust.type == "Molecule object":
                my_clust.mol_id.add(myps[i])
            elif my_type == "Act3DFrag object":
                my_clust.act3dfrag_id.add(myps[i])
            elif my_type == "ActPharmaPoint object":
                my_clust.actpharma_id.add(myps[i])
        # Now assign the sizes to these clusters
        for my_clust in clusters:
            if  my_clust.type == "MMPFrag object":
                # Find the the number of distinct molecules relating to this point
                my_clust.size = len(set(my_clust.mmp_frag_id.filter().values_list("mol_id__pk")))
            elif  my_clust.type == "PharmaPoint object":
                my_clust.size = len(set(my_clust.pharma_id.filter().values_list("mol_id__pk")))
            elif my_clust.type == "Water object":
                my_clust.size = len(set(my_clust.water_id.filter().values_list("prot_id__pk")))
            elif my_clust.type == "Molecule object":
                my_clust.size = len(set(my_clust.mol_id.filter().values_list("prot_id__pk")))
                # Define min and max
                my_clust.x_min, my_clust.x_max, my_clust.y_min, my_clust.y_max, my_clust.z_min, my_clust.z_max = get_min_and_max([Chem.MolFromMolBlock(x) for x in my_clust.mol_id.filter().values_list("sdf_info",flat=True)])
            elif my_type == "Act3DFrag object":
                my_clust.size = len(set(my_clust.act3dfrag_id.filter().values_list("pk")))
            elif my_type == "ActPharmaPoint object":
                my_clust.size = len(set(my_clust.actpharma_id.filter().values_list("pk")))
            my_clust.save()


def find_clusts(target_id):
    """Function to make the clusters for the PharmaPoints and Fragments for a target
    and store them to a file
    Takes a target_id
    Returns None"""
    # Get the target
    target = Target.objects.get(pk=target_id)
    # List of pharmacophores to look at
    my_pharmas = PharmaPoint.objects.filter(mol_id__prot_id__target_id=target_id).exclude(mol_id__prot_id__code__contains=target.title)
    # List if mmpfrags to consider
    my_mmpfrags = MMPFrag.objects.filter(mol_id__prot_id__target_id=target_id).exclude(mol_id__prot_id__code__contains=target.title).exclude(smistore="*[H]").exclude(x_com__isnull=True)
    # List of types to consider
    # and then we want to filter all the hydrophobes into one grouping
    pharma_type_list = list(set(my_pharmas.values_list("smiles", flat=True)))
    # List of mmpfrags to consider
    mmp_type_list = list(set(my_mmpfrags.filter().values_list("smistore", flat=True)))
    # CLuster the waters together
    # List of lambda values to use
    lam_lis = [2.0]#[1.0, 1.5, 2.0, 2.5, 3.0] #[0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 4.0, 4.5, 5.0] #3.0
    # Loop through the different lamda
    my_waters = Water.objects.filter(prot_id__target_id=target_id)
    for lam in lam_lis:
        print lam
        if len(my_mmpfrags) == 0:
            print "\nNO MMPS..."
        else:
            print "\nClustering MMPs.."
            make_clusters(my_mmpfrags, mmp_type_list, lam, target)
        if len(my_pharmas) == 0:
            print "\nNO PHARMAS..."
        else:
            print "\nClustering Pharmacophores..."
            make_clusters(my_pharmas, pharma_type_list, lam, target)
#    # Now make the area around which we cluster the points
    lam_lis = [1.5]
    if len(my_waters) == 0:
        pass
    else:
        for lam in lam_lis:
            print "\nClustering Waters..."
            make_clusters(my_waters, [1], lam, target)
    # Get the mols
    mols = Molecule.objects.filter(prot_id__target_id=target).exclude(prot_id__code__contains=target.title)
    # Get the centre of mass of these
    for mol in mols:
        mol.x_com, mol.y_com, mol.z_com = centre_of_mass_from_SD_block(str(mol.sdf_info))
    # Now cluster these dp_means - lam = 8.0
    for lam in [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]:
        print "\nClustering mols... lam=", lam
        make_clusters(mols, [1], lam, target)
    # Now clean up the clusters -> delete all empty ones
    KeyCluster.objects.filter(target_id=target, size=0).delete()


def get_ring_names(target_id):
    """Function to assign ring names to rings"""
    ring_dict = {'c1ccoc1': 'furan',
'c1ccoc1': 'pyrrole',
'c1ccsc1': 'thiophene',
'c1c[nH]cn1': 'imidazole',
'c1cn[nH]c1': 'pyrazole',
'c1cnoc1': 'isoxazole',
'c1cnsc1': 'thiazole',
'c1ccccc1': 'benzene',
'c1ccncc1': 'pyridine',
'c1cnccn1': 'pyrazine',
'c1cncnc1': 'pyrimidine',
'c1ccnnc1': 'pyridazine',
'c1cnnnc1': '1,2,3-triazine',
'c1cnncn1': '1,2,3-triazine',
'c1ncncn1': '1,3,5-triazine',
'c1cc2ccccc2o1': 'benzofuran',
'c1occ2ccccc12': 'isobenzofuran',
'C1CCCCC1': 'cyclohexane',
'C1CCOCC1': 'cyclohexane',
'C1CCNCC1': 'cyclohexane',
'c1scc2c1OCCO2': 'funky s'}
    kcs = KeyCluster.objects.filter(type='MMPFrag object',target_id=target_id)
    for k in kcs:
        my_smi = Chem.MolToSmiles(Chem.MolFromSmiles(str(k.function.replace("[*:1]","[H]").replace("[*:2]","[H]").replace("[*:3]","[H]"))),isomericSmiles=True)
        if my_smi in ring_dict:
            k.ringname = ring_dict[my_smi]
            k.save()


def get_res(pdb_block, chain=None):
    """"Script to get the residues out of the pdb information"""
    my_d = {}
    for line in pdb_block:
        if line[:4] != "ATOM":
            continue
        line = line.strip()
        if chain:
            if line[21].strip() != chain:
                continue
        res_name = line[17:20].strip()
        res_num = line[22:26].strip()
        my_res = res_name + str(res_num)
        if my_res not in my_d:
            my_d[my_res] = [line]
        else:
            my_d[my_res].append(line)
    out_d = {}
    for res in my_d:
        out_d[res] = Chem.MolFromPDBBlock("\n".join(my_d[res]))
    return out_d


def find_res_rmsd(mol1, mol2, res):
    """Calculate the differences between the atom coordinates of two identical structures"""
    if (not mol1) or (not mol2):
        print "NONE MOL", res
        return None
    # Gets atoms in mol1 (e.g. 14,5,3...) that match mol2 (1,2,3...)
    matchpatterns = mol1.GetSubstructMatches(mol2, uniquify=False)
    # Check to see if the molecules actually DO contain common substructures
    if not matchpatterns:
        # In this instance it may be only partial occupancy
        matchpatterns = mol2.GetSubstructMatches(mol1, uniquify=False)
        if not  matchpatterns:
            print "NO MATCH ", res
            return 0.0
    differences = []
    # Get the conformers to access the coords
    conf1 = mol1.GetConformer(0)
    conf2 = mol2.GetConformer(0)
    # May be more than one matching pattern. Calculate all of them.
    for matchlist in matchpatterns:
        # The total distance
        tot_dist = 0
        for idx2, idx1 in enumerate(matchlist):
            # Get the atom coords
            try:
                atm1 = conf1.GetAtomPosition(idx1)
                atm2 = conf2.GetAtomPosition(idx2)
            except:
                print res
                return None
            # Find the distance
            dist = find_dist(atm1, atm2)
            tot_dist += dist
        # Find the mean square distance
        mean_dist = float(tot_dist) / float(len(matchlist))
        # Append the root mean square distance
        differences.append(math.sqrt(mean_dist))
    # Return the differences corresponding to all of the ways of matching the molecules
    return min(differences)


def find_all_rmsd(prot_1, prot_2):
    """Function to get the residue to residue RMSD - for all residues"""
    # Define the chain like this - if autoloaded
    try:
        chain_1 = prot_1.code.split("_")[1]
    except IndexError:
        print "NO CHAIN ID", prot_1.code
        chain_1 = None
    # Get the residues for this protein
    res_prot_1 = get_res(str(prot_1.pdb_info).split("\n"), chain_1.upper())
    # Define the chain like this - if autoloaded
    try:
        chain_2 = prot_2.code.split("_")[1]
    except IndexError:
        chain_2 = None
        print "NO CHAIN ID", prot_2.code

    res_prot_2 = get_res(str(prot_2.pdb_info).split("\n"), chain_2.upper())
    # Check these dicts match
    rmsds = {}
    for res in res_prot_1:
        if res in res_prot_2:
            rmsds[res] = find_res_rmsd(res_prot_1[res], res_prot_2[res], res)
    return rmsds


def get_all_coords(rdmol):
    """Function to get all the coordinates for an RDMol
    Returns three lists"""
    conf = rdmol.GetConformer()
    x_coords = []
    y_coords = []
    z_coords = []
    for atm in rdmol.GetAtoms():
        x_coords.append(conf.GetAtomPosition(atm.GetIdx()).x)
        y_coords.append(conf.GetAtomPosition(atm.GetIdx()).y)
        z_coords.append(conf.GetAtomPosition(atm.GetIdx()).z)
    return x_coords, y_coords, z_coords


def register_res(prot, targ, prot_c, tot_prots):
    """Function to register all the resdiues for a protein"""
    try:
        chain = prot.code.split("_")[1]
    except IndexError:
        print "NO CHAIN ID", prot.code
        chain = None
    # Now get a dict of res
    out_d = get_res(str(prot.pdb_info).split("\n"), chain.upper())
    # Now go through it 
    old = -1
    tot = len(out_d)
    for i, item in enumerate(out_d):
        if not out_d[item]:
            continue
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rProtein %d of %d Registering residues %d%% complete..." % (prot_c + 1, tot_prots, old))
            sys.stdout.flush()
        res = Residue.objects.get_or_create(target_id=targ, prot_id=prot, res_num=item[3:], res_name=item[:3])[0]
        # Now get the min and max
        x_coords, y_coords, z_coords = get_all_coords(out_d[item])
        res.x_max = max(x_coords)
        res.x_min = min(x_coords)
        res.y_max = max(y_coords)
        res.y_min = min(y_coords)
        res.z_max = max(z_coords)
        res.z_min = min(z_coords)
        # Now add the pdb information - so we can use this later
        res.pdb_info = Chem.MolToPDBBlock(out_d[item])
        res.save()


def find_targ_rmsd(targ):
    """Function to find all RMSDs against all other RMSDs
    Takes a target
    Returns three dictionaries"""
    # Get all the proteins
    prots = [x for x in Protein.objects.filter(pdb_info__isnull=False, target_id=targ) if str(x.pdb_info) != ""]
    # Register all the residues
    # Print progress here for the loading of the proteins
    tot_prots = len(prots)
    for prot_counter, prot in enumerate(prots):
        register_res(prot, targ, prot_counter, tot_prots)
    print "Registered residues"
    # A dict to hold all the RMSD calcs in the end
    tot_d = {}
    # A dict to keep track of an individual proteins residues
    prot_dict = {}
    for p in prots:
        prot_dict[p.pk] = {}
    clust_d = {}
    res_set = set()
    # Loop through and calculate the RMSDs
    # Print a sensible protein counter
    for prot_c, prot_1 in enumerate(prots):
        sys.stdout.write("\rProtein %d of %d" % (prot_c + 1, tot_prots))
        sys.stdout.flush()
        for prot_2 in prots:
            # Function to find the rmsds between every residue in the protein
            out_d = find_all_rmsd(prot_1, prot_2)
            # Now tidy this up
            for res in out_d:
                res_set.add(res)
                if res in tot_d:
                    tot_d[res].append(out_d[res])
                else:
                    tot_d[res] = [out_d[res]]
                # Log this into the dictionary for all proteins
                if res in prot_dict[prot_1.pk]:
                    prot_dict[prot_1.pk][res].append(out_d[res])
                else:
                    prot_dict[prot_1.pk][res] = [out_d[res]]
                if res in prot_dict[prot_2.pk]:
                    prot_dict[prot_2.pk][res].append(out_d[res])
                else:
                    prot_dict[prot_2.pk][res] = [out_d[res]]
                # Log this into the dictionary for clustering
                if res in clust_d:
                    if prot_1.pk in clust_d[res]:
                        clust_d[res][prot_1.pk][prot_2.pk] = out_d[res]
                    else:
                        clust_d[res][prot_1.pk] = {prot_2.pk: out_d[res]}
                else:
                    clust_d[res] = {prot_1.pk: {prot_2.pk: out_d[res]}}
    # A dict to find the RMSDs beteeen 
    out_d = {}
    # Loop through the residues
    for res in list(res_set):
        out_d[res] = {}
        for prot in prots:
            if prot.pk not in clust_d[res]:
                # If it doesn't exist give it zero for each case
                out_d[res][prot.pk] = [0.0 for x in prots]
                continue
            else:
                out_d[res][prot.pk] = []
            for prot_comp in prots:
                if prot_comp.pk not in clust_d[res][prot.pk]:
                    out_d[res][prot.pk].append(0.0)
                else:
                    if clust_d[res][prot.pk][prot_comp.pk]:
                        out_d[res][prot.pk].append(clust_d[res][prot.pk][prot_comp.pk])
                    else:
                        out_d[res][prot.pk].append(0.0)
    print "Got residue RMSDs"
    return tot_d, prot_dict, out_d


def calc_com(atoms):
    """Function to calculate the unweighted centre of mass"""
    numatoms = len(atoms)
    x_coord = 0.0
    y_coord = 0.0
    z_coord = 0.0
    # Assume all heavy atoms have the same mass
    for coords in atoms:
        x_coord += float(coords[0])
        y_coord += float(coords[1])
        z_coord += float(coords[2])
    return x_coord / numatoms, y_coord / numatoms, z_coord / numatoms


def make_probe(x_com, y_com, z_com, type, target_id):
    """Function to make a new probe"""
    new_probe = Probe()
    new_probe.x_com = x_com
    new_probe.y_com = y_com
    new_probe.z_com = z_com
    new_probe.type = type
    new_probe.target_id = Target.objects.get(pk=target_id)
    new_probe.save()


def cluster_mols(target_id):
    """Function to cluster molecules for a target"""
    target = Target.objects.get(pk=target_id)
    # Get the mols
    mols = Molecule.objects.filter(prot_id__target_id=target).exclude(prot_id__code__contains=target.title)
    # Get the centre of mass of these
    for mol in mols:
        mol.x_com, mol.y_com, mol.z_com = centre_of_mass_from_SD_block(str(mol.sdf_info))
    # Now cluster these dp_means - lam = 8.0
    for lam in [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]:
        print "Clustering mols --- lam=", lam
        make_clusters(mols, [1], lam, target)


def render_pharma_pdb(index, occupancy, bfactor, atom, x_coord, y_coord, z_coord):
    type_dict = {u'SHAPE': 'Na', u'SingleAtomAcceptor': 'Li',
                 u'SingleAtomDonor': 'Na', u'RH5_5': 'K',
                 u'RH6_6': 'K', u'Arom6': 'Zn', u'Arom5': 'Zn',
                 u'ThreeWayAttach': 'K',  u'MagicSix': 'Na', u'Methyl': 'Au',
                 u'Imidazole' : 'Zn', u'SmallHalogen': 'Na', u'SingleF': 'Na', u'SingleCl': 'Na'}
    try:
        my_s = "HETATM"+" "*(5-len(str(index+1)))+str(index+1)+" "+type_dict[atom]+" "*(5-len(type_dict[atom]))+type_dict[atom]+" "*(5-len(type_dict[atom]))+" "*(4-len(str(index+1)))+str(index+1)+" "*(12-len("{0:.3f}".format(x_coord)))+"{0:.3f}".format(x_coord)+" "*(8-len("{0:.3f}".format(y_coord)))+"{0:.3f}".format(y_coord)+" "*(8-len("{0:.3f}".format(z_coord)))+"{0:.3f}".format(z_coord)+" "*(6-len("{0:.2f}".format(min(1.0,math.fabs(occupancy)))))+"{0:.2f}".format(min(1.0,math.fabs(occupancy)))+" "*(6-len("{0:.2f}".format(bfactor)))+"{0:.2f}".format(bfactor)+" "*(12-len(str(type_dict[atom])))+type_dict[atom]+"1+\n"
        #my_s = "COMPND    UNNAMED\nAUTHOR    GENERATED BY OPEN BABEL 2.3.2\n" + my_s
        #my_s = my_s + "MASTER        0    0    0    0    0    0    0    0    1    0    1    0\nEND\n"
    except KeyError:
        print atom
        return ""
    return my_s


def find_dist(mol_1, mol_2):
    """Function to find the square distance between two points in 3D
    Takes two len=3 tuples
    Returns a float"""
    return pow((mol_1.x-mol_2.x),2) + pow((mol_1.y-mol_2.y),2) + pow((mol_1.z-mol_2.z),2)


def get_data_for_targ(my_site, target_name, csv_path=""):
    """Function to get the data for a given target
    - depending on what site you are in
    Returns a dict - keys - model_id, chain_id, path_to_pdb, path_to_map, cmpd_id, smiles
    MUST be unique on model_id+chain_id. PDBs returned must be aligned. User can input own functions"""
    if my_site == "SGC":
        from WONKA.datasites import sgc_load
        if csv_path:
            # CSV path is a list of model IDs
            act_path, out_d = sgc_load(target_name, csv_path)
        else:
            act_path, out_d = sgc_load(target_name)
    elif my_site == "PDB":
        from WONKA.datasites import pdb_load
        act_path, out_d = pdb_load(target_name)
    elif my_site == "CSV":
        out_d = {}
        if not csv_path:
            print "You must add a CSV filepath - none added"
            sys.exit()
        if not os.path.isfile(csv_path):
            print "CSV not file"
            sys.exit()
        # Take the info from a CSV file
        in_d = csv.DictReader(open(csv_path))
        # Fields looking for
        all_fields = ["smiles", "path_to_pdb", "cmpd_id", "model_id", "chain", "path_to_map"]
        mandatory_fields = ["smiles", "path_to_pdb", "cmpd_id", "model_id", "chain"]
        # Check to see if fields are missing
        missing_fields = [x for x in mandatory_fields if x not in in_d.fieldnames]
        if len(missing_fields) != 0:
            print " ".join(missing_fields), " fields required"
            sys.exit()
        if len([x for x in all_fields if x not in in_d.fieldnames]) != 0:
            print " ".join([x for x in all_fields if x not in in_d.fieldnames]), " fields missing"
        out_d = {}
        # Loop through the d - populating the totla
        for item in in_d:
            model = item["model_id"]
            chain = item["chain"]
            if model + chain in out_d:
                print "DUPLICATE MODEL ERROR!!!!!"
                continue
            out_d[model + chain] = {}
            out_d[model + chain]["path_to_pdb"] = item["path_to_pdb"]
            out_d[model + chain]["smiles"] = str(item["smiles"]).decode('string-escape')
            out_d[model + chain]["cmpd_id"] = item["cmpd_id"]
            out_d[model + chain]["path_to_map"] = None
            out_d[model + chain]["model_id"] = model
            out_d[model + chain]["chain"] = chain
        # Create the path for the activity data
        act_path = target_name + "ACTS.csv"
    elif my_site == "SIMPLE":
        out_d = csv_path
        act_path = ""#None
    return out_d, act_path


def do_wonka(target_id):
    # Make the MMPs for 3D molecules
    make_mmp_database(target_id, option=None)
    # Now the pharmacophore points
    fdefName = os.path.join(os.path.join(os.path.split(sys.argv[0])[0], 'data/media'), 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    title = Target.objects.get(pk=target_id).title
    mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__contains=title)
    print "Making pharmacophore points..."
    old = -1
    tot = len(mols)
    for i, m in enumerate(mols):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rMaking pharmacophore points %d%% complete..." % old)
            sys.stdout.flush()
        find_pharmacophore_points(m, factory, fdefName)
    print "\nPharmacophores made"
    print "\nNow finding fragment Centre of Mass"
    # Now find the CofM
    find_frag_com(MMPFrag.objects.filter(mol_id__in=mols))
    # Now make the summary clusters
    print "Finding clusters..."
    find_clusts(target_id)


def find_target(target_name, overwrite=None, save_map=None, data_site="SGC", csv_path=None):
    """Function to find / update the data for a given target"""
    from django.contrib.auth.models import User
    from WONKA.views import send_email_user
    import urllib
    # First pull the target
    new_target = Target.objects.get_or_create(title=target_name)
    target = new_target[0]
    if new_target[1]:
        for user in UserData.objects.all():
            user.new_targs.add(target)
            user.save()
    # Now call the method required to get the data in - returns a dict with every
    # Aligned chain and extra information
    tot_d, act_path = get_data_for_targ(data_site, target_name, csv_path)
    # Now register all the information
    prots_made = register_targ_data(tot_d, target, save_map=False, overwrite=False)
    if len(prots_made) == 0:
        print "Nothing updated -  no need to bother with analysis"
    else:
        # Now load up the RMSDs
        out_d, prot_d, clus_d = find_targ_rmsd(target)
        # Now store these residue shifts on a per protein and per target level
        mark_res_shift(out_d, prot_d, clus_d, target)
        # Load the activity data if it exists - in CSV form 
        if os.path.isfile(act_path):
            load_activity_data(target, act_path)
        # Do the WONKA processing
        do_wonka(target.pk)
        # Title
        subject = "New data available for " + target.title
        # Now e-mail everyone to say that the models have been updated
        users = User.objects.filter(pk__in=UserData.objects.filter(cool_targs=target).values_list("user_id", flat=True))
        # Now  send the e-mails to the users
        out_m = "There is new data for a target you are watching: " + target.title
        out_m = "<br>" + str(len(prots_made)) + " new proteins were deposited"
        mols = []
        for prot in prots_made:
            mols.extend(prot.molecule_set.all())
        out_m += "<br> " + str(len(mols)) + " new molecules were found <br> They are the following:<br>"
        for mol in mols:
            internal_id = InternalIDLink.objects.filter(mol_id=mol)
            if internal_id:
                my_int_id = internal_id[0].internal_id
            else:
                my_int_id = "ID NOT KNOWN"
            out_m += '<figure><img src="http://oommppaa.sgc.ox.ac.uk/Viewer/loader/?function=2DMOL&choice=' + urllib.quote(str(mol.smiles))+'"></img>  <figcaption>'+my_int_id+'</figcaption></figure>'
        for us in users.filter():
            send_email_user(out_m, subject, us.email)


def remove_hetatm(lines):
    """Function to remove HETATM and CONECT record"""
    out_lines = []
    for line in lines:
        if line[:6] in ["HETATM", "CONECT", "ANISOU"]:
            if  line[:6] == "HETATM" and len(line[17:20].strip()) == 3:
                continue
            elif  line[:6] == "HETATM":
                pass
            else:
                continue
        out_lines.append(line)
    return "".join(out_lines)


def assign_temp(block, smiles, model):
    """Helper function to create an RDMol from PDB information"""
    tmp = Chem.MolFromPDBBlock(block)
    template = Chem.MolFromSmiles(smiles)
    if not template:
        print smiles, " not recognised"
        return None
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, tmp)
    except ValueError:
        print "DOESN'T FIT", smiles
        mol = None
    except:
        print "UNSPECIFED ERRROR"
        mol = None
    return mol


def get_waters(lines):
    """Helper function to extract waters from a PDB file"""
    return "".join([line for line in lines if line[17:20] == "HOH"])


def get_ligs(lines):
    """Helper function to extract LIGAND data from a PDB file"""
    old_id = "NOTME"
    old_chain = ""
    out_d = []
    for line in lines:
        if line[:6] != "HETATM":
            continue
        res = line[17:20]
        if res in ["HOH", "EDO", "ACT"]:
                continue
        if old_id == "NOTME":
            line_list = [line]
            old_id = res
        elif res != old_id or int(line[22:26]) != old_chain:
            out_d.append("".join(line_list))
            line_list = [line]
            old_id = res
        else:
            line_list.append(line)
        old_chain = int(line[22:26])
    if old_id == "NOTME":
        return []
    out_d.append("".join(line_list))
    return out_d


def register_targ_data(tot_d, target, save_map=None, overwrite=None):
    """Function to register a targets data"""
    import gzip
    old = -1
    tot = len(tot_d)
    prots_made = []
    for tot_c, chain in enumerate(tot_d):
        if tot_c * 100 / tot != old:
            old = tot_c * 100 / tot
            sys.stdout.write("\rRegistering proteins %d%% complete..." % old)
            sys.stdout.flush()
        smiles = tot_d[chain]["smiles"]
        model = tot_d[chain]["model_id"]
        # Get the PDB file
        if tot_d[chain]["path_to_pdb"][-3:] == ".gz":
            file_lines = gzip.open(tot_d[chain]["path_to_pdb"]).readlines()
        else:
            file_lines = open(tot_d[chain]["path_to_pdb"]).readlines()
        #If we're saving maps and the map exists
        if save_map and tot_d[chain]["path_ to_map"]:
        # Get the map file
            map_lines = open(tot_d[chain]["path_to_map"]).readlines()
        # Get the mols from this
        mols = [assign_temp(block, smiles, model) for block in get_ligs(file_lines) if assign_temp(block, smiles, model) is not None]
        # So now we check that the model hasn't been updated OR exists
        # Check that everything's ok with the protein
        apo_prot = Chem.MolFromPDBBlock(remove_hetatm(file_lines), sanitize=False)
        if not apo_prot:
            print "NONE PROTEIN: ", model + "_" + tot_d[chain]["chain"] + "_" + str(target.title)
            continue
        # Check the protein exists
        prot_code = model + "_" + tot_d[chain]["chain"] + "_" + str(target.pk)
        prot_me = Protein.objects.get_or_create(target_id=target, code=prot_code)
        # Only proceed if this has been created OR overwrite / refresh is set
        my_prot = prot_me[0]
        # If this is a newly created prot - or overwrite is on carry on
        if prot_me[1] or overwrite:
            prots_made.append(my_prot)
            # Delete the molecules for this protein
            Molecule.objects.filter(prot_id=my_prot).delete()
            # Loop through them
            for i, mol in enumerate(mols):
                # Find the reference compounds
                comp_ref = add_new_comp(mol)
                # Find the molecules
                new_mol = Molecule.objects.get_or_create(smiles=smiles, sdf_info=Chem.MolToMolBlock(mol, includeStereo=True), lig_id=str(i), chain_id="A", cmpd_id=comp_ref, occupancy=0.0, prot_id=my_prot)
                my_mol = new_mol[0]
                ### Now add this to the user data
                if new_mol[1]:
                    for user in UserData.objects.all():
                        user.new_mols.add(my_mol)
                        user.save()
                # Calculate the RMSD between ligands
                if i > 0:
                    my_mol.rmsd = AllChem.GetBestRMS(mols[0], mol)
                    my_mol.save()
                else:
                    my_mol.rmsd = 0.0
                # Give the internal ID
                iidl = InternalIDLink()
                iidl.mol_id = my_mol
                iidl.internal_id = tot_d[chain]["cmpd_id"]
                iidl.save()
            # DEFINE THIS CLUSTER 
            # Apo protein - within 5A of the ligand for this protein
#            if not Chem.MolToPDBBlock(apo_prot):
#                my_prot.delete()
#                print "APO PROTEIN NOT REAL"
#                print model
#                sys.exit()
            if prot_me[1] or overwrite:
                my_prot.pdb_info = remove_hetatm(file_lines)
                # If we have a map - add it
                if save_map:
                    my_prot.cif_info = "".join(map_lines)
                my_prot.save()
            # Waters - within 5A of the main ligand cluster
            waters = Chem.MolFromPDBBlock(get_waters(file_lines))
            if waters is not None:
                # Check the waters exist
                if len(Water.objects.filter(prot_id=my_prot)) != 0 and not overwrite:
                    pass
                else:
                    conf = waters.GetConformer()
                    # Delete them for this protein
                    for w in Water.objects.filter(prot_id=my_prot):
                        w.delete()
                    for i in range(waters.GetNumAtoms()):
                        cp = conf.GetAtomPosition(i)
                        if waters.GetAtomWithIdx(i).GetSmarts() != "O":
                            continue
                        Water.objects.get_or_create(water_num=i + 1, prot_id=my_prot, target_id=target,x_com=cp.x,y_com=cp.y,z_com=cp.z)
    print "\nRegistered proteins"
    return prots_made


def mark_res_shift(out_d, prot_d, clus_d, target, lam=2.5):
    """Function to make the residue shifts for a target"""
    # Anaylsing reisude shifts
    print "Analysing residue shifts"
    tot = len(out_d)
    old = -1
    for i, res in enumerate(out_d):
        old = i * 100 / tot
        sys.stdout.write("\rSummarising residues %d%% complete..." % old)
        sys.stdout.flush()
        rp = ResShift()
        rp.target_id = target
        rp.res_num = int(res[3:])
        rp.res_name = res[:3]
        rp.max_shift = max(out_d[res])
        if [x for x in out_d[res] if x != 0.0]:
            rp.min_shift = min([x for x in out_d[res] if x != 0.0])
        else:
            rp.min_shift = 0.0
        rp.avg_shift = float(sum([x for x in out_d[res] if x])) / float(len(out_d[res]))
        try:
            rp.validate_unique()
            rp.save()
        except ValidationError:
            # Otherwise updateeverything
            rp = ResShift.objects.get(target_id=target, res_num=int(res[3:]), res_name=res[:3])
            rp.max_shift = max(out_d[res])
            if [x for x in out_d[res] if x != 0.0]:
                rp.min_shift = min([x for x in out_d[res] if x != 0.0])
            else:
                rp.min_shift = 0.0
            rp.avg_shift = float(sum([x for x in out_d[res] if x])) / float(len(out_d[res]))
            rp.save()
    # Now lets relate these to the residue objects
    print "General residues summarised"
    tot = len(prot_d)
    old = -1
    for i, prot in enumerate(prot_d):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rSummarising individual residues %d%% complete..." % old)
            sys.stdout.flush()
        for res in prot_d[prot]:
            res_num = int(res[3:])
            res_name = res[:3]
            # Pull it down and storethis info
            my_res = Residue.objects.filter(prot_id=prot, res_num=res_num, res_name=res_name)
            if my_res:
                my_res = my_res[0]
            else:
                continue
            my_res.max_shift = max(prot_d[prot][res])
            if [x for x in prot_d[prot][res] if x != 0.0]:
                my_res.min_shift = min([x for x in prot_d[prot][res] if x != 0.0])
            else:
                my_res.min_shift = 0.0
            my_res.avg_shift = float(sum([x for x in prot_d[prot][res] if x])) / float(len(prot_d[prot][res]))
            my_res.save()
    print "Individual residues summarised"
    # Now save the residues into clusters
    tot = len(clus_d)
    old = -1
    for i, res in enumerate(clus_d):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rClustering residues %d%% complete..." % old)
            sys.stdout.flush()
        # Now cluster this residue
        dp = dpmeans([clus_d[res][prot] for prot in clus_d[res]], lam, 0, False)
        dp.run()
        # Now save these clusters
        for j, prot in enumerate(clus_d[res]):
            my_res = Residue.objects.filter(prot_id__pk=prot, res_num=int(res[3:]), res_name=res[:3])
            if my_res:
                my_res = my_res[0]
                # Update the cluster id
                my_res.clust_id = dp.dataClusterId[j]
                my_res.save()
    print "Residues clustered"


def return_numbered_grid_models(model_list, target_title):
    """Fucntion to return a numbered grid of models. Model list is a list like m023 etc"""
    from rdkit.Chem import Draw
    t = Target.objects.get(title=target_title)
    # First get all the models from this list
    code_list = [x + "_a_" + str(t.pk) for x in model_list]
    cmps = Protein.objects.filter(code__in=code_list).values_list("molecule__cmpd_id__pk", flat=True)
    my_cmpds = Compound.objects.filter(pk__in=cmps)
    my_mols = [Chem.MolFromSmiles(str(x.smiles)) for x in my_cmpds]
    img = Draw.MolsToGridImage(my_mols, molsPerRow=6,subImgSize=(200,200),legends=["COMPOUND "+str(i) for i, cmp in enumerate(my_mols)])
    return img


def score_point(target_point, score_points):
    """Function to score a given ph4 point.
Takes a target point as input and a list of score points"""
    out_score = 0
    for score_p in score_points:
        out_score += exp(-eucl_dist(target_point, score_p))
    return out_score

from rdkit import DataStructs


def score_frag(target_frag, score_frags, fp_dict):
    """Function to score a given fragment
    Takes a target point as input and a list of score points,
    and a dictionary containing morgan fingerprints"""
    out_score = 0
    for score_p in score_frags:
        # Find the morgan fingerprint similarity
        morg_sim = DataStructs.TanimotoSimilarity(fp_dict[target_frag.pk], fp_dict[score_p.pk])
        out_score += exp(-eucl_dist(target_frag, score_p)) * morg_sim
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


def create_frag_scores(target_id, opt="XTAL"):
    """Function to create FragScores"""
    target = Target.objects.get(pk=target_id)
    tot_mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__contains=target.title)
    tot_set = MMPFrag.objects.filter(mol_id__in=tot_mols)
    # Evaluate the set
    len(tot_set)
    # Now make a dict of the Morgan fingerprints
    fp_dict = {}
    for mmp_f in tot_set:
        rdm = Chem.MolFromSmiles(str(mmp_f.smistore))
        fp_dict[mmp_f.pk] = AllChem.GetMorganFingerprint(rdm, 2)
    # Now score the mols
    len_mols = len(tot_mols)
    for i, m in enumerate(tot_mols):
        print (i + 1) * 100 / len_mols
        # Set up the score set (set of ph4s to score against)
        score_frags = tot_set.exclude(mol_id=m)
        # And the test set (set of ph4s to test)
        test_frags = tot_set.filter(mol_id=m)
        # Now test these
        for frag in test_frags:
            frag_score = FragScore()
            frag_score.type = "XTAL"
            frag_score.frag_id = frag
            try:
                frag_score.validate_unique()
            except ValidationError:
                frag_score = FragScore.objects.get(frag_id=frag, type="XTAL")
            # Now get the score for this one
            frag_score.score = score_frag(frag, score_frags, fp_dict)
            frag_score.save()

from rdkit import Chem
from rdkit.Chem import AllChem
from math import sqrt,pow
import sys



def eucl_dist(mol_x_com, mol_y_com, mol_z_com, x_com, y_com, z_com):
    """Function to find the Euclidean distance between two points
    Takes 6 floats (coords)
    Returns a Euclidean distance"""
    return sqrt(pow((mol_x_com - x_com), 2) + pow((mol_y_com - y_com), 2) + pow((mol_z_com - z_com), 2))


def prepare_for_confs(mol, core, coreConfId=-1):
    """Function """
    import random
    # Step to help failures - removes hydrogens from molecule
    mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
    # Get the matches - and don't uniquify
    matches = mol.GetSubstructMatches(core, uniquify=False)
    if not matches:
        print "CORE ", Chem.MolToSmiles(core, isomericSmiles=True)
        print "MOL ", Chem.MolToSmiles(mol, isomericSmiles=True)
        raise ValueError, "molecule doesn't match the core"
    match = matches[random.randint(0, len(matches) - 1)]
    coordMap = {}
    coreConf = core.GetConformer(coreConfId)
    for i, idxI in enumerate(match):
        corePtI = coreConf.GetAtomPosition(i)
        coordMap[idxI] = corePtI
    return mol, match, coordMap


def find_conformations(mol, core, match, coordMap, useTethers=True, coreConfId=-1, randomseed=2342, max_iters=200, opt=None):
    """Function to generate conformations. Heavily based on ConstrainedEmbed in the RDKit
    Uses a forcefield (default MMFF) to generate conformations constrained to a
    core smiles. Does energy minimisation. Calculates the RMSD
    Takes an RDKit molecule and a core. Options are to useTethers,
    coreConfId - the conformer ID to use, randomseed - the randomseed to use,
    maxIts - the maximum number of iterations for the minimisation,
    opt -  the forcefield to use.
    Returns an RDKit molecule
    """
    ci = AllChem.EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, useRandomCoords=True)
    if ci < 0:
        print Chem.MolToMolBlock(mol)
        print Chem.MolToMolBlock(core)
        #raise ValueError, 'Could not embed molecule.'
        print "COULD NOT  EMBED"
        return None
    # Now make a map of the points to tether
    algMap = [(j, i) for i, j in enumerate(match)]
    if not useTethers:
        # clean up the conformation
        if opt is "MMFF":
            try:
                mmff_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), sanitize=False, removeHs=False)
                myff = Chem.rdForceFieldHelpers.SetupMMFFForceField(mmff_mol, mmffVerbosity=0)
                ff = AllChem.MMFFGetMoleculeForceField(mol, myff, confId=0)
            # Because the newer version of RDKit has this difference
            except AttributeError:
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
        for i, idxI in enumerate(match):
            for j in range(i + 1, len(match)):
                idxJ = match[j]
                d = coordMap[idxI].Distance(coordMap[idxJ])
                ff.AddDistanceConstraint(idxI, idxJ, d, d, 300.)
        ff.Initialize()
        n = 4
        more = ff.Minimize()
        while more and n:
            more = ff.Minimize()
            n -= 1
        # rotate the embedded conformation onto the core:
        rms = AllChem.AlignMol(mol, core, atomMap=algMap)
    else:
        # rotate the embedded conformation onto the core:
        rms = AllChem.AlignMol(mol, core, atomMap=algMap)
        if opt is "MMFF":
            try:
                mmff_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), sanitize=False, removeHs=False)
                myff = Chem.rdForceFieldHelpers.SetupMMFFForceField(mmff_mol, mmffVerbosity=0)
                ff = AllChem.MMFFGetMoleculeForceField(mol, myff, confId=0)
            # Because the newer version of RDKit has this difference
            except AttributeError:
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
        conf = core.GetConformer()
        if ff is None:
            sys.stderr.write("FORCEFIELD IS NONE\n" + Chem.MolToSmiles(mol))
            return None

        for i in range(core.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
            ff.AddDistanceConstraint(pIdx, match[i], 0, 0.0, 300.)
        ff.Initialize()
        # Do an energy minimisation
        # Forcefield parameters taken from Greg Landrum
        more = ff.Minimize(maxIts=max_iters, energyTol=1e-4, forceTol=1e-3)
        # Four extra steps of minimisation -> as prescribed in Greg's method
        n = 4
        while more and n:
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            n -= 1
        # Realign
        rms = AllChem.AlignMol(mol, core, atomMap=algMap)
    mol.SetProp('EmbedRMS', str(rms))
    return (mol, ff.CalcEnergy())


def get_vdw_rad(atomic_num):
    """Function to get the user defined atomic radius"""
    atomic_rad_dict = {6: 1.7, 7: 1.55, 8: 1.52, 9: 1.47}
    if atomic_num in atomic_rad_dict:
        return atomic_rad_dict[atomic_num]
    else:
        return float(Chem.GetPeriodicTable().GetRvdw(atomic_num))


def check_for_clashes(rdmol, protein_coords=None, threshold=None):
    """Function to look for clashes between the protein and the RDKit molecule
    rdmol - the RDKit molecule
    protein - a list of geometry points for the protein
    threshold - the min distance between points - can be negative if clashes are allowed"""
    # If the protein restrictions do not exist- 
    tot_dists = []
    if not protein_coords:
        return True
    # Find the molecule coords
    conf = rdmol.GetConformer()
    # Now get the radius
    mol_pos = [[conf.GetAtomPosition(atm.GetIdx()), get_vdw_rad(atm.GetAtomicNum())] for atm in rdmol.GetAtoms()]
    # If it is a list of list
    # Go throough the molecule coordinates
    # Go through the proteisn
    for protein in protein_coords:
        this_prot_list = []
        # Go through the individidual coord
        for prot_pos in protein:
            for pos, rad in mol_pos:
                this_prot_list.append(eucl_dist(pos.x, pos.y, pos.z, prot_pos.x, prot_pos.y, prot_pos.z) - (rad + get_vdw_rad(prot_pos.atom_num)))
        tot_dists.append(min(this_prot_list))
    #Find the closest approach in the most distant protein
    return max(tot_dists)


def eucl_dist_2(mol_1, mol_2):
    """Function to find the Euclidean distance between two molecules in 3D
    Takes two len=3 tuples
    Returns a float"""
    from math import sqrt
    return sqrt(pow((mol_1.x - mol_2.x), 2) + pow((mol_1.y - mol_2.y),2) + pow((mol_1.z - mol_2.z), 2))


def find_dist(mol_1, mol_2):
    """Function to find the square difference between two vectors"""
    return pow((mol_1.x-mol_2.x),2) + pow((mol_1.y-mol_2.y),2) + pow((mol_1.z-mol_2.z),2)


def find_rmsd(mol1, mol2):
    """Calculate the differences between the atom coordinates of two identical structures"""
    import math
    if (not mol1) or (not mol2):
        print "NONE MOL"
        return None
    # Gets atoms in mol1 (e.g. 14,5,3...) that match mol2 (1,2,3...)
    matchpatterns = mol1.GetSubstructMatches(mol2, uniquify=False)
    # Check to see if the molecules actually DO contain common substructures
    if not matchpatterns:
        # In this instance it may be only partial occupancy
        matchpatterns = mol2.GetSubstructMatches(mol1, uniquify=False)
        if not  matchpatterns:
            matchpatterns = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol1)).GetSubstructMatches(Chem.MolFromMolBlock(Chem.MolToMolBlock(mol2)), uniquify=False)
            if not matchpatterns:
                print "NO MATCH"
                return 0.0
            else:
                mol1 = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol1))
                mol2 = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol2))
    differences = []
    # Get the conformers to access the coords
    conf1 = mol1.GetConformer(0)
    conf2 = mol2.GetConformer(0)
    # May be more than one matching pattern. Calculate all of them.
    for matchlist in matchpatterns:
        # The total distance => set to zero
        tot_dist = 0
        # Loop through the matches
        for idx2, idx1 in enumerate(matchlist):
            # Get the atom coords
            try:
                atm1 = conf1.GetAtomPosition(idx1)
                atm2 = conf2.GetAtomPosition(idx2)
            except:
                return None
            # Find the distance
            dist = find_dist(atm1, atm2)
            # Add this distance to the sum
            tot_dist += dist
        # Find the mean dists
        mean_dist = float(tot_dist) / float(len(matchlist))
        # Now add the route of this to the possible differences (there may be more than one substructure match
        differences.append(math.sqrt(mean_dist))
    # Return the differences corresponding to all of the ways of matching the molecules -> we want the minimum
    return min(differences)


def gen_confs_lloommppaa(rdmol, core_smi, num_confs=100, ff="MMFF", prot_pos=None, max_iters=200):
    """Function to make all the conformations - we then store everything against a LLConf object and filter from that"""
    import random
    try:
        Chem.SanitizeMol(core_smi)
    except ValueError:
        print "ERROR LOADING MOL"
        return []
    # Now add Hydrogens
    rdmol = AllChem.AddHs(rdmol, addCoords=True)
    # If this is None too - return None
    if not rdmol:
        return None
    conf_list = []
    conf_counter = 0
    while True:
        if conf_counter == num_confs:
            break
        rand_val = random.randint(1, 10000)
        # Find a random conformation
        mol, match, coordMap = prepare_for_confs(rdmol, core_smi, coreConfId=-1)
        # Add hydrogens to the returned conforation
        try:
            mol = AllChem.AddHs(mol, addCoords=True)
            retval = find_conformations(mol, core_smi, match, coordMap, useTethers=True,
                                        randomseed=rand_val, max_iters=max_iters,
                                        opt=ff)
        except:
            print "ADDING HYDROGENS FAILS"
            retval = find_conformations(mol, core_smi, match, coordMap, useTethers=True,
                                        randomseed=rand_val, max_iters=max_iters,
                                        opt=ff)
        if retval is None:
            return []
        # Now take the strain and retval values
        strain = retval[1]
        retval = retval[0]
        if retval == None:
            # If there is an error in the conformation generatin - carry on
            conf_counter += 1
            continue
        # When we're off the ground
        if conf_counter > 0:
            my_min_diff = min([find_rmsd(Chem.MolFromMolBlock(item[0]), retval) for item in conf_list])
        else:
            my_min_diff = 100.0
        my_clash = check_for_clashes(retval, prot_pos)
        try:
            conf_list.append([Chem.MolToMolBlock(retval), my_clash, strain, my_min_diff])
        except:
            print "ERROR IN MOL"
            return []
        conf_counter += 1
    return conf_list


def generate_conformations(rdmol, core_smi, num_confs=100, num_fails=150, max_iters=200, ff="MMFF", prot_pos=None, min_diff=0.35, threshold=-0.5, clash_fail=True):
    """Function to generate conformations for the activity point molecule.
    Removes overly similar conformations.
    Input RDKit molecule of the ActivityPoint molecule (already given 3D coords),
    core_smi - a smiles string of the shared core or a molecule of it ,
    num_confs - the number of conformations required,
    num_fails - the number of potential fails where the new conformation is too similar,
    max_iters - the maximum number of iterations,
    ff - the forcefield used,
    prot_pos - the coords of steric hindrance regions,
    min_diff - the min RMSD between conformers,
    threshold - the threshold for the clash detection,
    clash_fail - whether to recognise a clash as a fail.
    Output is a list of RDKit molecule in different conformations"""
    import random
    try:
        Chem.SanitizeMol(core_smi)
    except ValueError:
        print "ERROR LOADING MOL"
        return []
    # Now add Hydrogens
    rdmol = AllChem.AddHs(rdmol, addCoords=True)
    fail_counter = 0
    conf_list = []
    conf_counter = 0
    while True:
        if fail_counter == num_fails:
            break
        if conf_counter == num_confs:
            break
        rand_val = random.randint(1, 10000)
        # Find a random conformation
        mol, match, coordMap = prepare_for_confs(rdmol, core_smi, coreConfId=-1)
        # Add hydrogens to the returned conforation
        try:
            mol = AllChem.AddHs(mol, addCoords=True)
            retval = find_conformations(mol, core_smi, match, coordMap, useTethers=True,
                                        randomseed=rand_val, max_iters=max_iters,
                                        opt=ff)
        except:
            print "ADDING HYDROGENS FAILS"
            retval = find_conformations(mol, core_smi, match, coordMap, useTethers=True,
                                        randomseed=rand_val, max_iters=max_iters,
                                        opt=ff)
        if retval is None:
            # If there is an error in the conformation generation
            # This means all will fail - so leave the function
            return None
        else:
            strain = retval[1]
            retval = retval[0]
        # When we're off the ground
        if conf_counter > 0:
          #  Check for overly similar conformations
            if min([find_rmsd(Chem.MolFromMolBlock(item[0]), retval) for item in conf_list]) < min_diff:
                # Add to the fail counter
                fail_counter += 1
            else:
                # Add to the conf_counter and add the conf
                # Check for steric clash with the protein (currently None)
                my_clash = check_for_clashes(retval, prot_pos, threshold)
                if my_clash:
                  # Add this molecule
                    conf_counter += 1
                    conf_list.append([Chem.MolToMolBlock(retval), my_clash, strain])
                else:
                    if clash_fail:
                        fail_counter += 1
                    else:
                        conf_counter += 1
                        conf_list.append([Chem.MolToMolBlock(retval), False, strain])
        # To begin with
        else:
            # Or if this is a new one check for clashes only
            my_clash = check_for_clashes(retval, prot_pos, threshold)
            if my_clash:
                conf_counter += 1
                # Add a new molecule
                conf_list.append([Chem.MolToMolBlock(retval), my_clash, strain])
            else:
                # Otherwise - score the molecule as a clash
                if clash_fail:
                    fail_counter += 1
                else:
                    conf_counter += 1
                    conf_list.append([Chem.MolToMolBlock(retval), False, strain])
    # Return the conformations generated
    return conf_list


##### SANIFIX 4 ##################
### This code is for fixing aromatic nitrogens (pyrrolic) 
### Part of the RDKit (not my code)
def fragment_inds_to_mol(oMol, indices):
    em = Chem.EditableMol(Chem.Mol())
    newIndices = {}
    for i, idx in enumerate(indices):
        em.AddAtom(oMol.GetAtomWithIdx(idx))
        newIndices[idx] = i

    for i, idx in enumerate(indices):
        at = oMol.GetAtomWithIdx(idx)
        for bond in at.GetBonds():
            if bond.GetBeginAtomIdx() == idx:
                oidx = bond.GetEndAtomIdx()
            else:
                oidx = bond.GetBeginAtomIdx()
            # make sure every bond only gets added once:
            if oidx < idx:
                continue
            em.AddBond(newIndices[idx], newIndices[oidx], bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    Chem.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap = newIndices
    return res


def _recursivelyModifyNs(mol, matches, indices=None):
    if indices is None:
        indices = []
    res = None
    while len(matches) and res is None:
        tIndices = indices[:]
        nextIdx = matches.pop(0)
        tIndices.append(nextIdx)
        nm = Chem.Mol(mol.ToBinary())
        nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
        nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(1)
        cp = Chem.Mol(nm.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            res, indices = _recursivelyModifyNs(nm, matches, indices=tIndices)
        else:
            indices = tIndices
            res = cp
    return res, indices


def adjust_arom_Ns(m, nitrogenPattern='[n&D2&H0;r5,r6]'):
    """
       default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
       to fix: O=c1ccncc1
    """
    Chem.GetSymmSSSR(m)
    m.UpdatePropertyCache(False)

    # break non-ring bonds linking rings:
    em = Chem.EditableMol(m)
    linkers = m.GetSubstructMatches(Chem.MolFromSmarts('[r]!@[r]'))
    plsFix = set()
    for a, b in linkers:
        em.RemoveBond(a, b)
        plsFix.add(a)
        plsFix.add(b)
    nm = em.GetMol()
    for at in plsFix:
        at = nm.GetAtomWithIdx(at)
        if at.GetIsAromatic() and at.GetAtomicNum() == 7:
            at.SetNumExplicitHs(1)
            at.SetNoImplicit(True)

    # build molecules from the fragments:
    fragLists = Chem.GetMolFrags(nm)
    frags = [fragment_inds_to_mol(nm, x) for x in fragLists]

    # loop through the fragments in turn and try to aromatize them:
    ok = True
    for i, frag in enumerate(frags):
        cp = Chem.Mol(frag.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
            lres, indices = _recursivelyModifyNs(frag, matches)
            if not lres:
                #print 'frag %d failed (%s)'%(i,str(fragLists[i]))
                ok = False
                break
            else:
                revMap = {}
                for k, v in frag._idxMap.iteritems():
                    revMap[v] = k
                for idx in indices:
                    oatom = m.GetAtomWithIdx(revMap[idx])
                    oatom.SetNoImplicit(True)
                    oatom.SetNumExplicitHs(1)
    if not ok:
        return None
    return m

######################### OOMMPPAA FUNCTIONS ###################################

# The following  were developed to deal with a bug in the RDKit that would not 
# allow constrained embed using aromatic sulphurs


def clear_aromatic_S(my_mol):
    """Function to replace an aromatic sulphur with an aromatic oxygen with an unusual isotope
    Takes an RDKit molecule
    Returns the updated molecule"""
    for atom in my_mol.GetAtoms():
        if atom.GetIsAromatic()is True and atom.GetAtomicNum() == 16:
            atom.SetAtomicNum(8)
            atom.SetIsotope(17)
    return my_mol


def clear_isotopic_O(my_mol):
    """Function to replace an aromatic oxygen with an unusual isotope with an aromatic sulphur
    Takes an RDKit molecule
    Returns the updated molecule"""
    for atom in my_mol.GetAtoms():
        if atom.GetIsotope() == 17 and atom.GetAtomicNum() == 8:
            atom.SetAtomicNum(16)
            atom.SetIsotope(0)
    return my_mol


def remove_breaks(my_mol):
    """Function to remove the break point
    Takes an RDKit molecule
    Returns the upadte molecule"""
    for atom in my_mol.GetAtoms():
        if "*" in atom.GetSmarts():
            atom.SetAtomicNum(1)
            atom.SetIsotope(0)
    return my_mol


def canonicalise_context(my_mol):
    """Function to canonicalise the context
    Takes an RDKit molecule
    Returns the upadte molecule"""
    return Chem.MolFromSmiles(Chem.MolToSmiles(my_mol,isomericSmiles=True))


def clean_up_frags(frag_out1, frag_out2=None):
    """Function to clean up fragments that are broken becasue of nH aromaticity issues
    Takes two fragments
    Runs sanfix4 programs
    Returns the fixed fragments"""
    try:
        Chem.SanitizeMol(frag_out1)
    except ValueError:
        nm = adjust_arom_Ns(frag_out1)
        if nm is not None:
            Chem.SanitizeMol(nm)
            sys.stderr.write('Fixed aromaticity:' + Chem.MolToSmiles(nm))
            frag_out1 = nm
        else:
            sys.stderr.write('Aromaticity still broken')
            return None,None
    if frag_out2:
        try:
            Chem.SanitizeMol(frag_out2)
        except ValueError:
            nm = adjust_arom_Ns(frag_out2)
            if nm is not None:
                Chem.SanitizeMol(nm)
                sys.stderr.write('Fixed aromaticity:' + Chem.MolToSmiles(nm))
                frag_out2 = nm
            else:
                sys.stderr.write('Aromaticity still broken')
                return None,None
    return frag_out1, frag_out2


def find_attachment_point(match1, match2, mol1, mol2):
    """Function to find the attachment point for two molecules
    Takes two lists of matches and the two molecules
    Returns the updated matche lists."""
# Now align the fragments as best you can
    listmatch1 = list(match1)
    listmatch2 = list(match2)
    # Find the linking point -> an atom in match but bonded to an atom not in match
    for atm in  mol1.GetAtoms():
        if atm.GetIdx() in match1:
            if len([at for at in atm.GetNeighbors() if at not in match1]) != 0:
                listmatch1.append(atm.GetIdx())
                #Now find the substituent closest to this guy
                if len([at for at in atm.GetNeighbors() if at not in match1]) == 1:
                    listmatch1.append(atm.GetNeighbors()[0].GetIdx())
                else:
                    listmatch2.append(atm.GetNeighbors()[0].GetIdx())
    for atm in  mol2.GetAtoms():
        if atm.GetIdx() in match2:
            if len([at for at in atm.GetNeighbors() if at not in match2]) != 0:
                listmatch2.append(atm.GetIdx())
                if len([at for at in atm.GetNeighbors()if at not in match2]) == 1:
                    listmatch2.append([at for at in atm.GetNeighbors()if at not in match2][0].GetIdx())
                else:
                # Pick one at random
                    listmatch2.append([at for at in atm.GetNeighbors()if at not in match2][0].GetIdx())
    return tuple(listmatch1), tuple(listmatch2)

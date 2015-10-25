import os
import math
import sys
import re
import cStringIO
import ast
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import Lipinski
from rdkit import Chem
from django.core.exceptions import ValidationError
from Pharmacophore.models import PharmaPoint,find_pharmacophore_points,UniquePharma
from Pharmacophore.functions import eucl_dist
from IOhandle.models import Target,ActivityPoint,Molecule,Protein
from IOhandle.functions import centre_of_mass_from_SD_block,set_pH_sdf,set_pH,make_smarts_from_frag
from MMPMaker.models import ActPharmaPoint,make_new_mmp,MMP,MMPFrag,MMPDiffMap,MMPComparison,Act3DFrag,make_new_3d_frag,ActMapPoint
from MMPMaker.helpers import clear_aromatic_S,clear_isotopic_O,canonicalise_context,remove_breaks,clean_up_frags,generate_conformations


def find_core(mol1, mol2, context, option=None):
    """Function to find the core, with attachement point, of a pair of like
    molecules.
    Takes two RDKit molecules and a matching part. Option "SIMPLE" means
    attachment point isn't used.
    Returns an RDKit molecule of the core"""
    match2 = mol2.GetSubstructMatch(context)
    myatm = None
    for idx in range(mol2.GetNumAtoms() - 1, - 1, - 1):
        if idx not in match2:
            # Find the atom adjoining this core
            if len([atm.GetIdx() for atm in  mol2.GetAtomWithIdx(idx).GetNeighbors() if atm.GetIdx() in match2]) == 1:
                # Replace it with the atom of the joining position of the new molecule
                myatm = mol2.GetAtomWithIdx(idx).GetAtomicNum()
            else:
                pass
    match = mol1.GetSubstructMatch(context)
    em = Chem.EditableMol(mol1)
    omit = []
    non_arom = []
    if option == "SIMPLE":
        # We don't do any replacement if option is "SIMPLE"
        myatm = False

    # Now loop through and find the match
    for idx in range(mol1.GetNumAtoms() - 1, - 1, - 1):
        if idx not in match:
            # I fwe're doing this extra step
            if myatm:
                if len([atm.GetIdx() for atm in  mol1.GetAtomWithIdx(idx).GetNeighbors() if atm.GetIdx() in match]) == 1:
            # Replace it with the atom of the joining position of the new molecule
                    em.ReplaceAtom(idx, Chem.Atom(myatm))
                    omit.append(idx)
                    non_arom.append(idx)
                else:
                    em.RemoveAtom(idx)
                # Find one with the correct atomic number and add that to the list
            else:
                em.RemoveAtom(idx)
    coreMol = em.GetMol()
    return coreMol


def find_correct(fragmented_smi, fragmented_smi_noIsotopes):
    """Function to find the correct co-ordinates of the 3D molecule.
    Takes the string of the fragments seperated by "."
    Returns the strings of the core and sidechains"""
    core = ""
    side_chains = ""
    f_array = fragmented_smi.split(".")
    # Now loop through the fragments
    for f in f_array:
        attachments = f.count("*")
        # If the guy only has one attachemnt it is the core
        if (attachments == 1):
            side_chains = "%s.%s" % (side_chains, f)
        else:
            core = f
    if core == "":
        # Error check for null core
        sys.stderr.write("ERROR NULL CORE\n" + f_array + "\n" + \
                         fragmented_smi_noIsotopes)
    side_chains = side_chains.lstrip('.')

    #cansmi the side chains
    temp = Chem.MolFromSmiles(side_chains)
    side_chains = Chem.MolToSmiles(temp, isomericSmiles=True)

    #and cansmi the core
    temp = Chem.MolFromSmiles(core)
    if temp is None:
        sys.stderr.write("ERROR NONE CORE: " + core)
    core = Chem.MolToSmiles(temp, isomericSmiles=True)
    return core, side_chains


def delete_bonds(bonds, mol, smi, outlines, my_id, option=None):
    """Function to delete bonds (specified in a list) of a given molecule and
    then deal with the fragments.
    Updates the outlines dictionary
    Returns the outlines dictionary"""
    #use the same parent mol object and create editable mol
    em = Chem.EditableMol(mol)
    #loop through the bonds to delete
    isotope = 0
    isotope_track = {}
    # Everything maps onto everything
    mymap = [(i, i) for i in range(len(mol.GetAtoms()))]
    for i in bonds:
        isotope += 1
        #remove the bond
        em.RemoveBond(i[0], i[1])
        # Take and store either side of this
        #now add attachement points
        newAtomA = em.AddAtom(Chem.Atom(0))
        em.AddBond(i[0], newAtomA, Chem.BondType.SINGLE)
        mymap.append((newAtomA, i[0]))

        newAtomB = em.AddAtom(Chem.Atom(0))
        em.AddBond(i[1], newAtomB, Chem.BondType.SINGLE)
        mymap.append((newAtomB, i[1]))
        #keep track of where to put isotopes
        isotope_track[newAtomA] = isotope
        isotope_track[newAtomB] = isotope
    #should be able to get away without sanitising mol
    #as the existing valencies/atoms not changed
    modifiedMol = em.GetMol()
    # Now align this mol to the old molecule so that R groups are
    # in the correct positions
    Chem.SanitizeMol(modifiedMol)
    for key in isotope_track:
            #to add isotope label
        modifiedMol.GetAtomWithIdx(key).SetIsotope(isotope_track[key])

    if option != None:
        return modifiedMol, isotope

    #canonical smiles can be different with and without the isotopes
    #hence to keep track of duplicates use fragmented_smi_noIsotopes
    fragmented_smi_noIsotopes = Chem.MolToSmiles(modifiedMol, isomericSmiles=True)
    valid = True
    fragments = fragmented_smi_noIsotopes.split(".")
    #Now take these fragments and map them onto those in
    #fragmented_smi_no_istopes
    #check if its a valid triple cut
    if(isotope == 3):
        valid = False
        for f in fragments:
            matchObj = re.search('\*.*\*.*\*', f)
            if matchObj:
                valid = True
                break

    if valid:
        if(isotope == 1):
            fragmented_smi_noIsotopes = re.sub('\[\*\]', '[*:1]', fragmented_smi_noIsotopes)
            fragments = fragmented_smi_noIsotopes.split(".")

            #print fragmented_smi_noIsotopes
            s1 = Chem.MolFromSmiles(fragments[0])
            s2 = Chem.MolFromSmiles(fragments[1])
            #need to cansmi again as smiles can be different
            core = Chem.MolToSmiles(s1, isomericSmiles=True)
            side_chains = Chem.MolToSmiles(s2, isomericSmiles=True)
            output = '%s,%s,,%s.%s' % (smi, my_id, core,side_chains)
            my_dict = {}
            # Get our fragments
            frgs = Chem.GetMolFrags(modifiedMol, asMols=True)
            # Now we find the core and the sidechains of the 3D frag
            # Do this by matching the canonical smiles of the fragment
            # to the canonical smiles of the original smiles fragments
            # ensures consistency with Hussain algorithm
            frag_zero = Chem.CanonSmiles(Chem.MolToSmiles(frgs[0], isomericSmiles=True).replace("[*]", "[H]"))
            frag_one = Chem.CanonSmiles(Chem.MolToSmiles(frgs[1], isomericSmiles=True).replace("[*]", "[H]"))
            smiles_one = Chem.CanonSmiles(Chem.MolToSmiles(s1, isomericSmiles=True).replace("[*:1]", "[H]"))
            if  frag_zero == smiles_one:
                my_dict["core"] = None
                my_dict["sidechains"] = []
                my_dict["sidechains"].append((frgs[0]))
                my_dict["sidechains"].append((frgs[1]))
            elif frag_one == smiles_one:
                my_dict["core"] = None
                my_dict["sidechains"] = []
                my_dict["sidechains"].append((frgs[1]))
                my_dict["sidechains"].append((frgs[0]))
            else:
                sys.stderr.write(Chem.MolToInchi(mol) + "\n" + [x for x  in fragments])
                sys.stderr.write(frag_zero + frag_one + smiles_one + "ERROR")
                sys.stderr.write("ERROR CAN'T FIND A FRAGMENT")
                sys.exit()

            if output not in outlines:
                outlines[output] = my_dict

        elif (isotope >= 2):
            for key in isotope_track:
                #to add isotope labels
                modifiedMol.GetAtomWithIdx(key).SetIsotope(isotope_track[key])
            fragmented_smi = Chem.MolToSmiles(modifiedMol, isomericSmiles=True)
            #change the isotopes into labels - currently can't add SMARTS or 
            #labels to mol
            fragmented_smi = re.sub('\[1\*\]', '[*:1]', fragmented_smi)
            fragmented_smi = re.sub('\[2\*\]', '[*:2]', fragmented_smi)
            fragmented_smi = re.sub('\[3\*\]', '[*:3]', fragmented_smi)

            fragments = fragmented_smi.split(".")

            #identify core/side chains and cansmi them
            core, side_chains = find_correct(fragmented_smi, fragmented_smi_noIsotopes.split("."))

            #now change the labels on sidechains and core
            #to get the new labels, cansmi the dot-disconnected side chains
            #the first fragment in the side chains has attachment label 1, 2nd:
            # 2, 3rd: 3
            #then change the labels accordingly in the core

            #this is required by the indexing script, as the side-chains are 
            #"keys" in the index
            #this ensures the side-chains always have the same numbering
            isotope_track = {}
            side_chain_fragments = side_chains.split(".")

            for s in xrange(len(side_chain_fragments)):
                matchObj = re.search('\[\*\:([123])\]', side_chain_fragments[s])
                if matchObj:
                    #add to isotope_track with key: old_isotope, value:
                    isotope_track[matchObj.group(1)] = str(s + 1)
            #change the labels if required, find all the old ones and how they 
            #map onto new ones
            for atom in modifiedMol.GetAtoms():
                #to add isotope labels
                if atom.GetIsotope() == (int(key)):
                    atom.SetIsotope(int(isotope_track[key]))
            # Do the same for the other guy
            if(isotope_track['1'] != '1'):
                core = re.sub('\[\*\:1\]', '[*:XX' + isotope_track['1'] + 'XX]', core)
                side_chains = re.sub('\[\*\:1\]', '[*:XX' + isotope_track['1'] + 'XX]', side_chains)
            if(isotope_track['2'] != '2'):
                core = re.sub('\[\*\:2\]', '[*:XX' + isotope_track['2'] + 'XX]' , core)
                side_chains = re.sub('\[\*\:2\]', '[*:XX' + isotope_track['2'] + 'XX]', side_chains)

            if(isotope == 3):
                if(isotope_track['3'] != '3'):
                    core = re.sub('\[\*\:3\]', '[*:XX' + isotope_track['3'] + 'XX]', core)
                    side_chains = re.sub('\[\*\:3\]', '[*:XX' + isotope_track['3'] + 'XX]', side_chains)
            #now remove the XX
            core = re.sub('XX', '', core)
            side_chains = re.sub('XX', '', side_chains)
            # Here we want to link core and sidechain to eachother and the molecule
            output = '%s,%s,%s,%s' % (smi, my_id, core, side_chains)
            # Now make it go into a dictionary
            if output not in outlines:
            # Now make a dict of the fragments as rdkit objects in a dict core: obj, side_chains: [obj,obj,obj,etc]
                my_dict = {}
            # Firstly lets find the core
                core_rd = [frag for frag in Chem.GetMolFrags(modifiedMol, asMols=True) if len(Chem.MolToSmiles(frag).split("*"))>2]
                if len(core_rd) >1:
                    sys.stderr.write("Found two options")
                    sys.exit()
                elif len(core_rd) == 0:
                    sys.stderr.write("Found no options" + str([Chem.MolToSmiles(frag) for frag in Chem.GetMolFrags(modifiedMol, asMols=True)]))
                    sys.exit()
                else:
                    my_dict["core"] = core_rd[0]
                my_dict["sidechains"] = []
                for side_chain in side_chains.split("."):
                    core_rd = [x for x in Chem.GetMolFrags(modifiedMol, asMols=True) if Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(x).replace("*", "H"))))) == Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(core.replace("*:1", "H").replace("*:2", "H").replace("*:3", "H")))))]
                    if len(core_rd) == 0:
                        print "CAN'T FIND IT TWO"
                        sys.stderr.write("ERROR CAN'T FIND FRAG" + Chem.MolToSmiles(Chem.MolFromSmiles(side_chain.replace("*:1", "H").replace("*:2", "H").replace("*:3", "H"))) + str([Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(x).replace("*", "H"))) for x in Chem.GetMolFrags(modifiedMol, asMols=True)]))
                    elif len(core_rd) == 1:
                        my_dict["sidechains"].append((core_rd[0]))
                    else:
                    # So now we need to pick from the possible options based on connections
                    # If it's a SIDECHAIN then it's the one with correct attachment points
                        core_rd = [x for x in Chem.GetMolFrags(modifiedMol,asMols=True) if re.findall(r'(\d+)\*', Chem.MolToSmiles(x, isomericSmiles=True)) == re.findall(r'\*:(\d+)', side_chain)]
                        if len(core_rd) != 1:
                            sys.stderr.write("ERROR IN FINDING BAD ONES IN SIDE CHAINS")
                            sys.stderr.write(re.findall(r'\*:(\d+)', side_chain))
                            sys.stderr.write(str([Chem.MolToSmiles(x, isomericSmiles=True) for x in Chem.GetMolFrags(modifiedMol, asMols=True)]))
                            sys.stderr.write(str([re.findall(r'(\d+)\*', Chem.MolToSmiles(x, isomericSmiles=True)) for x in core_rd]))
                            sys.exit()
                        else:
                            my_dict["sidechains"].append((core_rd[0]))
                    # Now find where there are two options
                    if my_dict is None:
                        sys.stderr.write("ERROR DICTIONARY EMPTY")
                        sys.exit()
                    outlines[output] = my_dict
    return outlines


def store_mmp(outlines, context_dict, i):
    """Function to store an MMP and MMPFrag in the database.
    Takes a list of lines in outlines.
    Returns None"""
    for line in outlines:
        my_dict = outlines[line]
        smi, my_id, core, context = line.split(',')
        if(len(core) == 0):
            side_chains = context.split('.')
            context = side_chains[0]
            core = side_chains[1]
            # FOR EACH OF THESE MAKE A NEW MMP OBJECT
            #print "CHAIN",side_chains[1],side_chains[0]
            context, mmp = make_new_mmp(context, my_id, core, my_dict["sidechains"][1],1,core, i)
            if context in context_dict:
                context_dict[context].append(mmp)
            else:
                context_dict[context] = [mmp]
            context = side_chains[1]
            core = side_chains[0]
            #print "CHAIN",side_chains[1],side_chains[0]
            context, mmp = make_new_mmp(context, my_id, core,my_dict["sidechains"][0],1,core, i)
            if context in context_dict:
                context_dict[context].append(mmp)
            else:
                context_dict[context] = [mmp]
    #double or triple cut
        else:
            attachments = core.count('*')
            #add the context with id to index
            context, mmp = make_new_mmp(context, my_id, core,my_dict["core"],attachments,core, i)
            if context in context_dict:
                context_dict[context].append(mmp)
            else:
                context_dict[context] = [mmp]
    return context_dict


def make_mol_mmp(mol, context_dict, myid=None, target_id=None):
    """Function to derive the MMPs for a given molecule.
    Takes in the molecule and fragments it up.
    Returns None."""
    if mol is None:
        sys.stderr.write("ERROR MOL IS NONE")
        return
    # Now get the smiles
    smi = Chem.MolToSmiles(mol)
    #different cuts can give the same fragments
    #to use outlines to remove them
    outlines = {}
    #SMARTS for "acyclic and not in a functional group"
    smarts = Chem.MolFromSmarts("[#6+0;!$(*=,#[!#6])]!@!=!#[*]")
    #finds the relevant bonds to break
    #find the atoms maches
    matching_atoms = mol.GetSubstructMatches(smarts)
    total = len(matching_atoms)
    #catch case where there are no bonds to fragment
    if(total == 0):
        output = '%s,%s,,' % (smi, myid)
        # Fix to prevent the dictionary being filled
        if (output in outlines) == False:
            # Sometimes there will be nothing to fragment
            outlines[output] = None
            return context_dict
    bonds_selected = []
    #loop to generate every single, double and triple cut in the molecule
    for x in xrange(total):
        #print matches[x]
        bonds_selected.append(matching_atoms[x])
        outlines = delete_bonds(bonds_selected, mol, smi, outlines, myid)
        # Now store this as the mols
        context_dict = store_mmp(outlines, context_dict, x)
        outlines = {}
        bonds_selected = []
##### Only do single cuts in this implementation. Uncomment to do double  and triple too ###
        for y in xrange(x+1,total):
        #print matching_atoms[x],matching_atoms[y]
            bonds_selected.append(matching_atoms[x])
            bonds_selected.append(matching_atoms[y])
            outlines = delete_bonds(bonds_selected,mol,smi,outlines,myid)
            bonds_selected = []
      # Now store this as the mols
            context_dict = store_mmp(outlines, context_dict, str(x)+"_"+str(y))
            outlines = {}
            for z in xrange(y+1, total):
        #print matching_atoms[x],matching_atoms[y],matching_atoms[z]
                bonds_selected.append(matching_atoms[x])
                bonds_selected.append(matching_atoms[y])
                bonds_selected.append(matching_atoms[z])
                outlines = delete_bonds(bonds_selected,mol,smi,outlines,myid)
                context_dict = store_mmp(outlines, context_dict, str(x)+"_"+str(y)+"_"+str(z))
                # Now store this as the mols
                outlines = {}
                bonds_selected = []
    return context_dict


def bulk_commit(context_dict):
    """Function to bulk create a series of MMPFrags"""
    tot_list = []
    for context in context_dict:
        # First get or create this guy and make him
        mmp = MMP.objects.get_or_create(context=context)[0]
        for mmpfrag_obj in context_dict[context]:
            mmpfrag_obj.mmp_link = mmp
            try:
                mmpfrag_obj.validate_unique()
                tot_list.append(mmpfrag_obj)
            except ValidationError:
                continue
    MMPFrag.objects.bulk_create(tot_list)


def find_mmps(rdmols, target_id):
    """Function to make all the MMPs from a list of molecules.
    Takes a dict of RDKit molecules with their corresponding ID as the value. And a target id
    Returns None"""
    tot = len(rdmols)
    old = -1
    context_dict = {}
    for i, mol in enumerate(rdmols):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
        val = rdmols[mol]
        if mol == "":
            continue
        if rdmols[mol][0:3] == "mol":
            mymol = Chem.MolFromMolBlock(mol)
            if mol is None:
                print mol
        else:
            mymol = Chem.MolFromSmiles(mol)
        if mol is None:
            print mol
        elif len(Chem.MolToSmiles(mymol).split(".")) > 1:
            # If the molecule is fragmented then it wont work
            print "ERROR FRAGMENTED MOL -> SKIPPING"
            continue
        context_dict = make_mol_mmp(mymol, context_dict, val, target_id)
    # Every 100 mols commit this stuff to the database
        if i % 100 == 0 and i != 0:
            bulk_commit(context_dict)
            context_dict = {}
    bulk_commit(context_dict)
    old = 100
    sys.stdout.write("\r%d%% complete..." % old)
    sys.stdout.flush()
    print "Completed"
    return None


def make_list_mmps(target_id, option=None, max_size=10, ratio=3, use_ratio=False):
    """Function to find MMPs from the database for agiven target.
    Takes a target ID and whether to do ActivityPoint to Molecule comparison (ACT) or Molecule to Molecule comparison (None)
    max_size is the max size of the non-matching part. ratio is the maximum ratio of the sizes of the matched part and non-matching part
    Returns the list of matched pairs"""
    print "Finding matched pairs"
    # Find all the context for this target
    myKeys = MMP.objects.filter(mmpfrag__mol_id__prot_id__target_id=target_id).distinct()
    out_mmps = []
    pkset = {}
    pkout = {}
    for key in myKeys:
        cont_mol = Chem.MolFromSmiles(str(key.context))
        if cont_mol.GetNumHeavyAtoms() < 3:
            continue
        # If we're dealing with 3D molecules
        if option is None:
            # Pull down all the matches
            myfrags = MMPFrag.objects.filter(mmp_link=key, mol_id__prot_id__target_id=target_id, core_size__lte=max_size)
            # If there is only one match
            if len(myfrags) == 1:
                continue
            # Now find all combinations of these two lists
            for i, frag1 in enumerate(myfrags):
                for frag2 in myfrags[i + 1:]:
                #now generate the pairs
                    id_a, core_a = frag1.mol_id, str(frag1.smistore)
                    id_b, core_b = frag2.mol_id, str(frag2.smistore)
                    #make sure pairs are not same molecule
                    if(id_a.pk != id_b.pk):
                        #make sure LHS and RHS of SMIRKS are not the same
                        if (core_a != core_b):
                            smirks, context = cansmirk(core_a, core_b, key.context)
                            out_mmps.append((id_a, id_b, frag1, frag2))
        # Or wif we're dealing with 2D activity molecules
        elif option is "ACT":
            actfrags = MMPFrag.objects.filter(mmp_link=key, act_id__target_id=target_id, core_size__lte=max_size)
            myfrags = MMPFrag.objects.filter(mmp_link=key, mol_id__prot_id__target_id=target_id, core_size__lte=max_size)
            if len(myfrags) == 0 or len(actfrags) == 0:
                continue
            for i, frag1 in enumerate(myfrags):
                for j, frag2 in enumerate(actfrags):
                    #now generate the pairs
                    id_a, core_a = frag1.mol_id.cmpd_id.pk, str(frag1.smistore)
                    id_b, core_b = frag2.act_id.cmpd_id.pk, str(frag2.smistore)
                    # Make sure the other compound does not have an xtal structure
                    target = Target.objects.get(pk=target_id)
                    if len(Molecule.objects.filter(cmpd_id=frag2.act_id.cmpd_id, prot_id__target_id=target_id).exclude(prot_id__code__contains=target.title)) > 0:
                        continue
                    #make sure pairs are not the same compound OR the same GSK compound
                    if id_a != id_b:
                        #make sure LHS and RHS of SMIRKS are not the same
                        if (core_a != core_b):
                            # At this point we should only add in MMPs if they are a new pair and contain the biggest context
                            # Only find one interaction per pair of compounds.
                            if (frag1.mol_id.cmpd_id, frag2.act_id.cmpd_id) in pkset:
                                if pkset[(frag1.mol_id.cmpd_id, frag2.act_id.cmpd_id)] < Lipinski.HeavyAtomCount(Chem.MolFromSmiles(str(key.context))):
                                    pkset[(frag1.mol_id.cmpd_id, frag2.act_id.cmpd_id)] = Lipinski.HeavyAtomCount(Chem.MolFromSmiles(str(key.context)))
                                    smirks, context = cansmirk(core_a, core_b, key.context)
                                    pkout[(frag1.mol_id.cmpd_id, frag2.act_id.cmpd_id)] = (frag1.mol_id, frag2.act_id, frag1, frag2, context, smirks)
                            else:
                                pkset[(frag1.mol_id.cmpd_id, frag2.act_id.cmpd_id)] = Lipinski.HeavyAtomCount(Chem.MolFromSmiles(str(key.context)))
                                smirks, context = cansmirk(core_a, core_b, key.context)
                                pkout[(frag1.mol_id.cmpd_id, frag2.act_id.cmpd_id)] = (frag1.mol_id, frag2.act_id, frag1, frag2, context,  smirks)
    for val in pkout:
        # Now populate the list
        out_mmps.append(pkout[val])
    # Give this warning if we haven't found any
    if len(out_mmps) == 0:
        print "WARNING!!! No matched pairs found!!!"
    else:
        print len(out_mmps), " matched pairs found in database"
    return out_mmps


def heavy_atom_count(smi):
    """Function compute the heavy atom count of a smiles string.
    Takes a smiles string
    Returns the heavy atom count as an int"""
    m = Chem.MolFromSmiles(smi)
    return m.GetNumAtoms()

### Functions of Jammeed Hussain from RDKit MMP contrib


def get_symmetry_class(smi):
    symmetry = []
    m = Chem.MolFromSmiles(smi)
    #determine the symmetry class
    #see: http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01894.html
    #A thanks to Greg (and Alan)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    #get the symmetry class of the attachements points
    #Note: 1st star is the zero index,
    #2nd star is first index, etc
    for atom in m.GetAtoms():
        if(atom.GetMass() == 0):
            symmetry.append(atom.GetProp('_CIPRank'))
    return symmetry


def cansmirk(lhs, rhs, context):
        #cansmirk algorithm
        #1) cansmi the LHS.
        #2) For the LHS the 1st star will have label 1, 2nd star will have label 2 and so on
        #3) Do a symmetry check of lhs and rhs and use that to decide if the labels on
        #   RHS or/and context need to change.
        #4) For the rhs, if you have a choice (ie. two attachement points are symmetrically
        #   equivalent), always put the label with lower numerical value on the earlier
        #   attachement point on the cansmi-ed smiles

        #print "in: %s,%s" % (lhs,rhs)
    isotope_track={}
    #if the star count of lhs/context/rhs is 1, single cut
    stars = lhs.count("*")

    if(stars > 1):
                #get the symmetry class of stars of lhs and rhs
        lhs_sym = get_symmetry_class(lhs)
        rhs_sym = get_symmetry_class(rhs)
        try:
            (lhs_sym[0] != lhs_sym[1]) and (rhs_sym[0] != rhs_sym[1])
        except IndexError:
            sys.stderr.write("ERROR IN CANSMI")
            sys.stderr.write(lhs+" "+rhs+" "+context)
            smirk = "%s>>%s" % (lhs,rhs)
            return smirk,context

    #deal with double cuts
    if(stars == 2):
        #simple cases
        #unsymmetric lhs and unsymmetric rhs
        if( (lhs_sym[0] != lhs_sym[1]) and (rhs_sym[0] != rhs_sym[1]) ):
                        #get 1st and 2nd labels and store the new label for it in isotope_track
                        #structure: isotope_track[old_label]=new_label (as strings)
            isotope_track = build_track_dictionary(lhs,stars)

            #switch labels using isotope track
            lhs = switch_labels_on_position(lhs)
            rhs = switch_labels(isotope_track,stars,rhs)
            context = switch_labels(isotope_track,stars,context)

        #symmetric lhs and symmetric rhs
        elif( (lhs_sym[0] == lhs_sym[1]) and (rhs_sym[0] == rhs_sym[1]) ):
            #the points are all equivalent so change labels on lhs and rhs based on position
            #labels on context don't need to change
            lhs = switch_labels_on_position(lhs)
            rhs = switch_labels_on_position(rhs)

        #more difficult cases..
        #symmetric lhs and unsymmetric rhs
        elif( (lhs_sym[0] == lhs_sym[1]) and (rhs_sym[0] != rhs_sym[1]) ):
            #switch labels lhs based on position
            lhs = switch_labels_on_position(lhs)
            #change labels on rhs based on position but need to record
            #the changes as need to appy them to the context
            isotope_track = build_track_dictionary(rhs,stars)
            rhs = switch_labels_on_position(rhs)
            context = switch_labels(isotope_track,stars,context)

        #unsymmetric lhs and symmetric rhs
        elif( (lhs_sym[0] != lhs_sym[1]) and (rhs_sym[0] == rhs_sym[1]) ):
            #change labels on lhs based on position but need to record
            #the changes as need to appy them to the context
            isotope_track = build_track_dictionary(lhs,stars)
            lhs = switch_labels_on_position(lhs)
            context = switch_labels(isotope_track,stars,context)
            #as rhs is symmetric, positions are equivalent so change labels on position
            rhs = switch_labels_on_position(rhs)

    #deal with triple cut
    #unwieldy code but most readable I can make it
    elif(stars == 3):
        #simple cases
        #completely symmetric lhs and completely symmetric rhs
        if( ( (lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2]) ) and
        ( (rhs_sym[0] == rhs_sym[1]) and (rhs_sym[1] == rhs_sym[2]) and (rhs_sym[0] == rhs_sym[2]) ) ):
            #the points are all equivalent so change labels on lhs and rhs based on position
            #labels on context don't need to change
            lhs = switch_labels_on_position(lhs)
            rhs = switch_labels_on_position(rhs)

        #completely symmetric lhs and completely unsymmetric rhs
        elif( ( (lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2]) ) and
        ( (rhs_sym[0] != rhs_sym[1]) and (rhs_sym[1] != rhs_sym[2]) and (rhs_sym[0] != rhs_sym[2]) ) ):
            #alter lhs in usual way
            lhs = switch_labels_on_position(lhs)
            #change labels on rhs based on position but need to record
            #the changes as need to appy them to the context
            isotope_track = build_track_dictionary(rhs,stars)
            rhs = switch_labels_on_position(rhs)
            context = switch_labels(isotope_track,stars,context)

        #completely unsymmetric lhs and completely unsymmetric rhs
        elif( ( (lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2]) ) and
        ( (rhs_sym[0] != rhs_sym[1]) and (rhs_sym[1] != rhs_sym[2]) and (rhs_sym[0] != rhs_sym[2]) ) ):
            #build the isotope track
            isotope_track = build_track_dictionary(lhs,stars)
            #alter lhs in usual way
            lhs = switch_labels_on_position(lhs)
            #change rhs and context based on isotope_track
            rhs = switch_labels(isotope_track,stars,rhs)
            context = switch_labels(isotope_track,stars,context)

        #completely unsymmetric lhs and completely symmetric rhs
        elif( ( (lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2]) ) and
        ( (rhs_sym[0] == rhs_sym[1]) and (rhs_sym[1] == rhs_sym[2]) and (rhs_sym[0] == rhs_sym[2]) ) ):
            #build isotope trach on lhs
            isotope_track = build_track_dictionary(lhs,stars)
            #alter lhs in usual way
            lhs = switch_labels_on_position(lhs)
            #change labels on context
            context = switch_labels(isotope_track,stars,context)
            #all positions on rhs equivalent so add labels on position
            rhs = switch_labels_on_position(rhs)

        #more difficult cases, partial symmetry
        #completely unsymmetric on lhs and partial symmetry on rhs
        elif( (lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2]) ):
            #build the isotope track
            isotope_track = build_track_dictionary(lhs,stars)
            #alter lhs in usual way
            lhs = switch_labels_on_position(lhs)
            #change rhs and context based on isotope_track
            rhs = switch_labels(isotope_track,stars,rhs)
            context = switch_labels(isotope_track,stars,context)

            #tweak positions on rhs based on symmetry
            #rhs 1,2 equivalent
            if(rhs_sym[0] == rhs_sym[1]):
                                #tweak rhs position 1 and 2 as they are symmetric
                rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,1,2)

            #rhs 2,3 equivalent
            elif(rhs_sym[1] == rhs_sym[2]):
                #tweak rhs position 1 and 2 as they are symmetric
                rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,2,3)

            #rhs 1,3 equivalent - try for larger set in future
            elif(rhs_sym[0] == rhs_sym[2]):
                #tweak rhs position 1 and 2 as they are symmetric
                rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,1,3)

        #now we are left with things with partial symmetry on lhs and not completely symmetric or unsymmetric on rhs
        else:
            #lhs 1,2,3 equivalent and any sort of partial symmetry on rhs
            if( (lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2]) ):

                #alter lhs in usual way
                lhs = switch_labels_on_position(lhs)
                #change labels on rhs based on position but need to record
                #the changes as need to appy them to the context
                isotope_track = build_track_dictionary(rhs,stars)
                rhs = switch_labels_on_position(rhs)
                context = switch_labels(isotope_track,stars,context)


            #now deal partial symmetry on lhs or rhs.
            #Cases where:
            #lhs 1,2 equivalent
            #lhs 2,3 equivalent
            #lhs 1,3 equivalent
            else:
                #build isotope track on lhs
                isotope_track = build_track_dictionary(lhs,stars)
                #alter lhs in usual way
                lhs = switch_labels_on_position(lhs)
                #change rhs and context based on isotope_track
                rhs = switch_labels(isotope_track,stars,rhs)
                context = switch_labels(isotope_track,stars,context)

                #tweak positions on rhs based on symmetry

                #lhs 1,2 equivalent
                if(lhs_sym[0] == lhs_sym[1]):
                    #tweak rhs position 1 and 2 as they are symmetric on lhs
                    rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,1,2)

                #lhs 2,3 equivalent
                elif(lhs_sym[1] == lhs_sym[2]):
                    #tweak rhs position 1 and 2 as they are symmetric on lhs
                    rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,2,3)

                #lhs 1,3 equivalent - try for larger set in future
                elif(lhs_sym[0] == lhs_sym[2]):
                    #tweak rhs position 1 and 2 as they are symmetric on lhs
                    rhs = switch_specific_labels_on_symmetry(rhs,rhs_sym,1,3)

    smirk = "%s>>%s" % (lhs,rhs)
    return smirk,context


def switch_specific_labels_on_symmetry(smi, symmetry_class, a, b):

    #check if a and b positions are symmetrically equivalent
    #if equivalent, swap labels if the lower numerical label is not on the
    #1st symmetrically equivalent attachment points in the smi

    if(symmetry_class[a-1] == symmetry_class[b-1]):
        #what are the labels on a and b

        matchObj = re.search( r'\[\*\:([123])\].*\[\*\:([123])\].*\[\*\:([123])\]', smi )
        if matchObj:
            #if the higher label comes first, fix
            if(int(matchObj.group(a)) > int(matchObj.group(b))):
            #if(int(matchObj.group(1)) > int(matchObj.group(2))):
                smi = re.sub(r'\[\*\:'+matchObj.group(a)+'\]', '[*:XX' + matchObj.group(b) + 'XX]' , smi)
                smi = re.sub(r'\[\*\:'+matchObj.group(b)+'\]', '[*:XX' + matchObj.group(a) + 'XX]' , smi)
                smi = re.sub('XX', '' , smi)

    return smi


def switch_labels_on_position(smi):

    #move the labels in order of position
    smi = re.sub(r'\[\*\:[123]\]', '[*:XX1XX]' , smi, 1)
    smi = re.sub(r'\[\*\:[123]\]', '[*:XX2XX]' , smi, 1)
    smi = re.sub(r'\[\*\:[123]\]', '[*:XX3XX]' , smi, 1)
    smi = re.sub('XX', '' , smi)

    return smi


def switch_labels(track, stars, smi):

    #switch labels based on the input dictionary track
    if(stars > 1):
        #for k in track:
        #        print "old: %s, new: %s" % (k,track[k])

        if(track['1'] != '1'):
            smi = re.sub(r'\[\*\:1\]', '[*:XX' + track['1'] + 'XX]' , smi)

        if(track['2'] != '2'):
            smi = re.sub(r'\[\*\:2\]', '[*:XX' + track['2'] + 'XX]' , smi)

        if(stars == 3):
            if(track['3'] != '3'):
                smi = re.sub(r'\[\*\:3\]', '[*:XX' + track['3'] + 'XX]' , smi)

        #now remove the XX
        smi = re.sub('XX', '' , smi)

    return smi


def build_track_dictionary(smi, stars):

    isotope_track = {}

    #find 1st label, record it in isotope_track as key, with value being the
    #new label based on its position (1st star is 1, 2nd star 2 etc.)
    if(stars == 2):
        matchObj = re.search( r'\[\*\:([123])\].*\[\*\:([123])\]', smi )
        if matchObj:
            isotope_track[matchObj.group(1)] = '1'
            isotope_track[matchObj.group(2)] = '2'

    elif(stars == 3):
        matchObj = re.search( r'\[\*\:([123])\].*\[\*\:([123])\].*\[\*\:([123])\]', smi )
        if matchObj:
            isotope_track[matchObj.group(1)] = '1'
            isotope_track[matchObj.group(2)] = '2'
            isotope_track[matchObj.group(3)] = '3'

    return isotope_track


### OOMMPPAA functions
def make_h_frag(link, mmp, my_smi):
    """Function to make a H MMPFrag object
    Takes a Molecule or ActivityPoint object, an MMP and a smiles string.
    Returns None"""
    new_mmp = MMPFrag()
    # Assign the appropriate point
    if type(link) == Molecule:
        new_mmp.mol_id = link
    elif type(link) == ActivityPoint:
        new_mmp.ActlID = link
    # The smiles that is stored
    new_mmp.smistore = "*[H]"
    # The 3D information for the fragment if it exists
    new_mmp.sdf_info = ""
    # The molecule it is associated to
    # Link back to the appropriate context
    new_mmp.mmp_link = mmp
    # To give the core size
    new_mmp.core_size = 1.0
    # To give the core to MMP ratio
    new_mmp.core_ratio = 1.0 / Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(my_smi))
    try:
        new_mmp.validate_unique()
        new_mmp.save()
        return new_mmp
    except ValidationError:
        return MMPFrag.objects.get(smistore="*[H]", sdf_info="", mmp_link=mmp)
    return None


def index_hydrogen_change(target_id):
    """Function to implement index h change using the Django data structure
    Takes a target id
    Returns None
    Algorithm ->
    1) Get the smiles string of all ActivityPoints and Molecules
    2) Get all the MMP objects and add Hydrogen
    3) Canonicalise and see if this matches the smiles string of the ActivityPoint or Molecules
    4) If yes add an MMPFrag object showing this H-addition.
    """
    print "Indexing H change"
    # Pull out all my activity pints
    acts = ActivityPoint.objects.filter(target_id=target_id)
    # I need to set_pH to charged -> as this is how the contexts will be
    act_mols = {set_pH(str(act.cmpd_id.smiles)): act.pk for act in acts }
    # Canonicalises the smiles
    act_smis = {Chem.MolToSmiles(Chem.MolFromSmiles(str(act_mol)), isomericSmiles=True): act_mols[act_mol] for act_mol in act_mols if Chem.MolFromSmiles(str(act_mol)) is not None }
    # Do the same for the rdkit molecule
    mol_smis = {Chem.MolToSmiles(Chem.MolFromMolBlock(str(m.sdf_info)), isomericSmiles=True): m.pk for m in Molecule.objects.filter(prot_id__target_id=target_id)}
    # Now loop through the contexts, excluding any double or triple cutes
    mmps = MMP.objects.exclude(context__contains="[*:2]")
    tot = len(mmps)
    old = -1
    for i, m in enumerate(mmps):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." %old)
            sys.stdout.flush()
        # Canonicalise the smiles of the context
        my_smi = Chem.MolToSmiles(Chem.MolFromSmiles(str(m.context).replace("*:1", "H")), isomericSmiles=True)
        if my_smi in act_smis:
            # If I can find it in the activity list
            # Pullout this activity point
            act_id = ActivityPoint.objects.get(pk=act_smis[my_smi])
            # Check this is still correct
            try:
                assert act_id.target_id.pk == target_id
            except AssertionError:
                print "Error: activity point no longer corresponds to original target..."
                sys.exit()
            # Now make the fragment
            make_h_frag(act_id, m, my_smi)
        if my_smi in mol_smis:
            mol_id = Molecule.objects.get(pk=mol_smis[my_smi])
            # Check this is still correct
            try:
                assert mol_id.prot_id.target_id.pk == target_id
            except AssertionError:
                print "Error: molecule no longer corresponds to original target..."
                sys.exit()
            # Now make the fragment
            make_h_frag(mol_id, m, my_smi)
    old = 100
    sys.stdout.write("\r%d%% complete..." % old)
    sys.stdout.flush()
    print "Indexing H change complete"
    return


def make_space_point(myfrag, actchange, target, my_type):
    """Function to find the centre of mass of the attached fragment for each molecule. It then makes a pharmacophore point object.
    Takes an Act3DFrag an activity change a Target and a type (indicating the type of activity used)
    Returns None."""
    # First of all make the PharmaPoint object
    ph_pt1 = PharmaPoint()
    ph_pt1.smiles = "SHAPE"
    ph_pt1.method_used = "PROTRUDEDIST"
    ph_pt1.git_version = "DISTRIB"
# Links back to the molecule found
    ph_pt1.mol_id = myfrag.mol_id
    ph_pt1.uniq_id = UniquePharma.objects.get_or_create(smiles="SHAPE")[0]
    ph_pt1.pharma_id = 10000
    # Make the 3D Fragments between the two molecules
    mol1 = Chem.MolFromMolBlock(str(myfrag.mol_id.sdf_info))
    # Define weight as the size differences
    ph_pt1.Weight = myfrag.mmp_frag_id.core_size
    # Now make the fragments
    # Make the central substructure fragment
    context = Chem.MolFromSmarts(make_smarts_from_frag(str(myfrag.mmp_frag_id.mmp_link.context)))
    # Produce the 3D Fragment that will match what we want
    frag_1 = AllChem.DeleteSubstructs(mol1, context)
    frag_1, fr_none = clean_up_frags(frag_1)
    if frag_1 is None:
        sys.stderr.write("ERROR NOT VALID FRAG")
        sys.stderr.write("ERROR "+Chem.MolToSmiles(mol1))
        sys.stderr.write("ERROR "+myfrag.mmp_frag_id.mmp_link.context)
        return None
    # Now register the COM of this fragment -> THIS IS WHERE THE POINT SITS
    ph_pt1.x_com, ph_pt1.y_com, ph_pt1.z_com = centre_of_mass_from_SD_block(Chem.MolToMolBlock(frag_1))
    # Catch for their being no info
    if ph_pt1.x_com is None:
        return
    # Now lets save the Ph4 point
    try:
        ph_pt1.validate_unique()
        ph_pt1.save()
    except ValidationError:
        ph_pt1 = PharmaPoint.objects.get(pharma_id=10000, mol_id=myfrag.mol_id, method_used="PROTRUDEDIST")
    # Add this point to the list -> It shows the addition this molecule.
    newActPha1 = ActPharmaPoint()
    newActPha1.act_3d_frag_id = myfrag
    newActPha1.from_xtal = True
    newActPha1.target_id = target
    newActPha1.act_change = actchange
    newActPha1.type = my_type
    newActPha1.pharma_id = ph_pt1
    try:
        newActPha1.validate_unique()
        newActPha1.save()
    except ValidationError:
        pass


def find_points(target_id, act_type):
    """Function to return all the Activity Inactivity and Shape points for the OOMMPPAA viewer
    Takes a Target and a Type of activity change.
    Returns three lists of points with associated data."""
    # Find the list of activity points
    act_list = [(act.act_3d_frag_id.over_mol_id.pk, act.act_3d_frag_id.mol_id.pk,
    act.pharma_id.pk, act.pharma_id.smiles, act.pharma_id.x_com, act.pharma_id.y_com,
    act.pharma_id.z_com, act.act_change, act.num_diff_15) for act in
    ActPharmaPoint.objects.filter(target_id=target_id, act_change__gt=0.0, in_diff_15=True,
    num_diff_15__lte=6, type=act_type).order_by("act_3d_frag_id")]
    # Find the list of inactivity points
    inact_list = [(act.act_3d_frag_id.over_mol_id.pk, act.act_3d_frag_id.mol_id.pk,
    act.pharma_id.pk, act.pharma_id.smiles, act.pharma_id.x_com, act.pharma_id.y_com,
    act.pharma_id.z_com, act.act_change, act.num_diff_15) for act in 
    ActPharmaPoint.objects.filter(target_id=target_id, act_change__lt=0.0, in_diff_15=True,
    num_diff_15__lte=6, type=act_type).order_by("act_3d_frag_id")]
    # Find all the points
    shape_list = [(act.act_3d_frag_id.over_mol_id.pk, act.act_3d_frag_id.mol_id.pk,
    act.pharma_id.pk, act.pharma_id.smiles, act.pharma_id.x_com, act.pharma_id.y_com,
    act.pharma_id.z_com, act.act_change * 3.0 / float(act.pharma_id.Weight),
    act.pharma_id.Weight) for act in ActPharmaPoint.objects.filter(target_id=target_id,
    pharma_id__smiles="SHAPE", type=act_type).order_by("act_3d_frag_id")]
    # Return them to the function
    return act_list, inact_list, shape_list


def make_mmp_database(target_id, option=None):
    """Function to find all the mmps for a given target.
    Takes a Target.
    Returns None"""
    # Now lets pull all the activity information for my molecules in 3D
    target = Target.objects.get(pk=target_id)
    # Define the possible activity types
    poss_list = ["IC50", "Ki", "TM"]
    # First of all get all compound ids that have ligand coordinates AND that have activity
    # Now pull all the molecules relating to this
    rdmols = {}

    # Finding activity molecuels
    if option is "TWOD":
        old = -1
        # Get all the activity points that haven't been done
        act_list = ActivityPoint.objects.filter(cmpd_id__mol_wt__lte=5000000,
        target_id=target_id, source__in=poss_list).filter(mmpfrag__isnull=True)
        # Find the length of the lists
        tot = len(act_list)
        print "CONVERTING MOLECULES TO APPROPRIATE pH..."
        for i, act in enumerate(act_list):
            if i * 100 / tot != old:
                old = i * 100 / tot
                sys.stdout.write("\r%d%% complete..." % old)
                sys.stdout.flush()
            # Now set the pH
            new_s = set_pH(str(act.cmpd_id.smiles))
            new_m = Chem.MolFromSmiles(new_s)
            if new_m is None:
                sys.stderr.write("NONE MOLECULE: " + new_s)
                continue
            rdmols[new_s] = "act" + str(act.pk)
        old = 100
        sys.stdout.write("\r%d%% complete..." % old)
        sys.stdout.flush()
        print "\npH SET"
        print "Now making MMPDB for activty molecules..."
    elif option is None:
        mol_list = Molecule.objects.exclude(prot_id__code__contains=target.title).filter(prot_id__target_id=target_id).order_by("pk")
        for i, act in enumerate(mol_list):
            # Add the 3D Molecules - assuming user has input at appropriate pH
            rdmols[str(act.sdf_info)] = "mol" + str(act.pk)
        if len(rdmols) != 0:
            print "Now making MMPDB for 3D molecules..."
    if len(rdmols) == 0:
        print "NO AVAILABLE DATA"
        return [], []
    find_mmps(rdmols, target_id)
    return


def act_mmp_3d(out_mmps, target_id):
    """Function to attribute 3D coordinates to the activity only molecule. Uses MMFF and maximises shape overlap with existing fragment.
    Takes a list of mmps and a Target
    Returns None"""
    print "Generating 3D conformations"
    # First of all make sure the protein exists to put all these artificial coordinates on
    my_target = Target.objects.get(pk=target_id)
    new_protein = Protein()
    new_protein.code = my_target.title + "ChEMBL"
    new_protein.target_id = my_target
    try:
        new_protein.validate_unique()
        new_protein.save()
    except ValidationError:
        new_protein = Protein.objects.get(code=my_target.title + "ChEMBL")
    # Now let's loop through the MMPs and overlay the act mol onto the other one
    tot = len(out_mmps)
    old = -1
    for i, ans in enumerate(out_mmps):
        # Print the progress
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
    # Pull out the molecules and the fragments mol1 is in  3D mol2 is an activity molecule
        mol1 = Chem.MolFromMolBlock(str(ans[0].sdf_info))
        # Set the protonation state here
        new_s = set_pH(str(ans[1].cmpd_id.smiles))
        if Chem.MolFromSmiles(new_s) is None:
            continue
        mol2 = Chem.MolFromSmiles(new_s)
        frag2 = ans[3]
        activity = ans[1]
        s_context = str(ans[4])
    # Split the context into it's fragments
        my_frags = s_context.split(".")
    # If its a single cut
        if len(my_frags) == 1:
            context = canonicalise_context((Chem.MolFromSmiles(s_context)))
        elif len(my_frags) == 2 or len(my_frags) == 3:
            # If it is a double or triple split ignore
            continue
        else:
            sys.stderr.write("ERROR MORE THAN THREE FRAGMENTS...")
            sys.exit()
        # Now filter out if it's too small
        if context.GetNumHeavyAtoms() < 3:
            continue
        # Now find the best core, using make_smarts_from_frag to make a 
        # relevant smarts pattern to do substructure subs
        smarts_mol = Chem.MolFromSmarts(make_smarts_from_frag(Chem.MolToSmiles(context, isomericSmiles=True)))
        core_mol = find_core(mol1, mol2, smarts_mol)
        # Check there is only one option. If there are two then we have a
        # problem. With the current implementation we can only
        if len(mol2.GetSubstructMatches(context)) > 1:
            sys.stderr.write("ERROR TWO POSSIBLE MATCHES")
            sys.stderr.write("ERROR " + Chem.MolToSmiles(context))
            sys.stderr.write("ERROR " + Chem.MolToSmiles(mol2))
            continue
        #Now the function to do the constraining
        try:
            mol2 = AllChem.ConstrainedEmbed(mol2, core_mol)
        except ValueError:
            # If it is a nitrile it needs extra sanitisation before embedding.
            Chem.SanitizeMol(core_mol)
            try:
                mol2 = AllChem.ConstrainedEmbed(mol2, core_mol)
            except ValueError:
                try:
                    # Setting the mol as a molblock and reading back in 
                    # can fix this issue
                    mol2 = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol2))
                    core_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(core_mol))
                    mol2 = AllChem.ConstrainedEmbed(mol2, core_mol)
                except ValueError:
                    # Except when it doesn't so print out this
                    sys.stderr.write(Chem.MolToMolBlock(mol2))
                    sys.stderr.write(Chem.MolToMolBlock(core_mol))
                    sys.stderr.write("ERROR MOL ONE:\n" + Chem.MolToSmiles(mol1))
                    sys.stderr.write("ERROR MOL TWO:\n" + Chem.MolToSmiles(mol2))
                    sys.stderr.write("ERROR CORE MOL:\n" + Chem.MolToSmiles(core_mol))
                    sys.stderr.write("ERROR CONTEXT:\n" + Chem.MolToSmiles(context))
                    sys.stderr.write("ERROR CANNOT CONSTRAINED EMBED MOLECULE")
                    continue
            # If the molecule is the same as the core then why bother
            if Chem.MolToSmiles(mol2, isomericSmiles=True) == Chem.MolToSmiles(core_mol, isomericSmiles=True):
                continue
        #NOW FILTER ON SHAPE FOR THE BEST ONE
        out_confs = generate_conformations(mol2, core_mol, num_confs=30, num_fails=30, max_iters=100, ff="MMFF")
        if out_confs is None:
            continue
        my_mols = [Chem.MolFromMolBlock(x[0]) for x in out_confs]
        mols = []
        for conf in my_mols:
            # Now calculate the shape distance
            mols.append((AllChem.ShapeTanimotoDist(conf, mol1, ignoreHs=True), conf))
        mols = sorted(mols, key=lambda x: x[0])
        # Now pick the best (lowest) one
        mol2 = mols[0][1]
#        # code to find the most energetically favourable one - if several have
#        if len(mols) > 1:
#            # If the molecules are identical in terms of shapeprotrudedist
#            # then use energy to get the value
#            if mols[0][0] == mols[1][0]:
#                # Check there aren't more like this
#                samemols = [x for x in mols if x[0] == mols[0][0]]
#                mineng = 100000
#                # Check the energy
#                for smol in samemols:
#                    mysmol = smol[1]
#                    # Calculate the energy
#                    try:
#                        mmff_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mysmol), sanitize=False)
#                        myff = Chem.rdForceFieldHelpers.SetupMMFFForceField(mmff_mol, mmffVerbosity=0)
#                        ff = AllChem.MMFFGetMoleculeForceField(mysmol, myff, confId=0)
#                    # Because the neweer version of RDKit has this difference
#                    except AttributeError:
#                        ff = AllChem.MMFFGetMoleculeForceField(mysmol, AllChem.MMFFGetMoleculeProperties(mysmol))
#                    if ff.CalcEnergy() < mineng:
#                        mineng = ff.CalcEnergy()
#                        mol2 = mysmol
#                # Derive the fragments
        me = Chem.MolFromSmarts(make_smarts_from_frag(s_context))
        frag_out2 = AllChem.DeleteSubstructs(mol2,me)
        # Clean up the fragments
        frag_out2, fr_none = clean_up_frags(frag_out2)
        # Now add the coordinates for the fragment to the fragment and add the molecule linker -> the fact that the molecule is linked to an activity protein singles it out as being Activity firs
        make_new_3d_frag(frag2, activity, ans[0], frag_out2, new_protein, mol2)
    old = 100
    sys.stdout.write("\r%d%% complete..." % old)
    sys.stdout.flush()
    print "Completed generating all 3D coordinates"


def make_ph4_difference_points(target_id, opt_type=["SPR"]):
    """Function to find pharmacophore differences between compounds.
    Takes a Target and activity type to be used.
    Returns None
    """
    # Define the allowed PH4 objects
    allowed_list = ["SingleAtomDonor", "SingleAtomAcceptor", "RH5_5", "RH6_6",
                    "Arom6", "Arom5", "ThreeWayAttach"]
    # Get the target object
    my_target = Target.objects.get(pk=target_id)
    # Set up the PH4 factory to make the Ph4 points
    fdefName = os.path.join(os.path.join(os.path.split(sys.argv[0])[0], 'data/media'), 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    # Now pull out the activity point for the original compound
    new_frags = []
    comb_list = []
    # Pull down all the Act3D frags -> we want the smallest possible one
    myfrags = Act3DFrag.objects.filter(act_id__target_id=my_target).distinct()
    # The below is to select unique of act_id and Overmol_id combination
    # Make a list of unique pairs
    print "Finding pharmacophore differences"
    for f in myfrags:
        comb_vals = (f.act_id, f.over_mol_id)
        # If it's already in the comb list then continue
        if comb_vals in comb_list:
            continue
        # Otherwwise add it to the list
        comb_list.append(comb_vals)
        # And add this fragment to the list to iterate through.
        new_frags.append(f)
    # Now pull out the activity point for the original compound
    tot = len(new_frags)
    old = -1
    for i, myfrag in enumerate(new_frags):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
        # Pick out the smallest one
        myfrag = Act3DFrag.objects.filter(act_id=myfrag.act_id, over_mol_id=myfrag.over_mol_id).order_by("mol_wt")
        if len(myfrag) == 0:
            print "No frags"
            continue
        else:
            myfrag = myfrag[0]
        # Now lets pull out the activity point
        over_act_point = ActivityPoint.objects.filter(cmpd_id=myfrag.over_mol_id.cmpd_id, target_id=target_id, source__in=opt_type).order_by("activity")
        if len(over_act_point) == 0:
            # Here a compound does not have activity data
            continue
        elif len(over_act_point) > 1:
            # If there is multiple activity points we simply take the biggest one.
            my_act_point = over_act_point[len(over_act_point) - 1].activity
            my_operator = over_act_point[len(over_act_point) - 1].operator
            my_source = over_act_point[len(over_act_point) - 1].source
        else:
            my_act_point = over_act_point[0].activity
            my_operator = over_act_point[0].operator
            my_source = over_act_point[0].source
        # Now pull out the SPR data for the given compound
        act_dat_point = ActivityPoint.objects.filter(cmpd_id=myfrag.act_id.cmpd_id, target_id=target_id, source__in=opt_type).order_by("activity")
        if len(act_dat_point) == 0:
            print "No activity points"
            continue
        elif len(act_dat_point) > 1:
            # If there is multiple activity points we simply take the biggests one.
            fr_act_point = act_dat_point[len(act_dat_point) - 1].activity
            fr_operator = act_dat_point[len(act_dat_point) - 1].operator
            fr_source = act_dat_point[len(act_dat_point) - 1].source
        else:
            fr_act_point = act_dat_point[0].activity
            fr_operator = act_dat_point[0].operator
            fr_source = act_dat_point[0].source
        # Generate a point that is the centre of mass of the attachment frag
        # that shows the shape difference
        # So now define the space differences. This also gives a general picture
        # of where things have been
        # If they are not the same type -> ignore
        if  my_operator == "<" and fr_operator == "<":
            # Both compoounds are inactive so ignore this one
            continue
        elif  my_operator == "<":
            # Check if the inactive compound has a higher measured activity
            if my_act_point > fr_act_point:
                # This is now invalid
                continue
            else:
                diff_act = my_act_point - fr_act_point
        elif fr_operator == "<":
            # Check if the inactivity compound has a higher measured activity
            if fr_act_point > my_act_point:
                # This is also invalid
                continue
            else:
                diff_act = my_act_point - fr_act_point
        else:
            # If its just a normal one
            diff_act = my_act_point - fr_act_point
        make_space_point(myfrag, diff_act, my_target, opt_type)
        # Make the PH4 points in space for the two mols and then define their 
        # differences
        find_pharmacophore_points(myfrag.over_mol_id, factory, fdefName)
        find_pharmacophore_points(myfrag.mol_id, factory, fdefName)
        over_mol_ph4 = PharmaPoint.objects.filter(smiles__in=allowed_list, mol_id=myfrag.over_mol_id)
        mol_ph4 = PharmaPoint.objects.filter(smiles__in=allowed_list, mol_id=myfrag.mol_id)
        # So now contrast these pharmapoints
        # Now we just loop through both sets and add them in
        match_set15 = set()
        match_set10 = set()
        match_set25 = set()
        match_set30 = set()
        match_set35 = set()
        dist_dict = {}
        # Take all the ph4 points and the one that's closest
        for ph_pt in over_mol_ph4:
            # Loop over the points from the other molecule
            for ph_pt2 in mol_ph4:
                #  Check that they are of the same type -> this is uniqueID so
                # Arom5 != Arom6
                if ph_pt.smiles == ph_pt2.smiles:
                    # Add the pk to a dictionrary. Add this small distance to
                    # the dict -> i.e. keep track of the smallest distance
                    if ph_pt.pk not in dist_dict:
                        dist_dict[ph_pt.pk] = eucl_dist(ph_pt, ph_pt2)
                    elif eucl_dist(ph_pt, ph_pt2) < dist_dict[ph_pt.pk]:
                        dist_dict[ph_pt.pk] = eucl_dist(ph_pt, ph_pt2)
                    # Do the same for the other point. This could certainly be done more efficiently
                    if ph_pt2.pk not in dist_dict:
                        dist_dict[ph_pt2.pk] = eucl_dist(ph_pt, ph_pt2)
                    elif eucl_dist(ph_pt, ph_pt2) < dist_dict[ph_pt2.pk]:
                        dist_dict[ph_pt2.pk] = eucl_dist(ph_pt, ph_pt2)

                    # Now if the points fall below certain thresholds add them into the appropriate sets
                    if eucl_dist(ph_pt, ph_pt2) < 3.5:
                        match_set35.add(ph_pt2.pk)
                        match_set35.add(ph_pt.pk)
                    if eucl_dist(ph_pt, ph_pt2) < 3.0:
                        match_set35.add(ph_pt2.pk)
                        match_set35.add(ph_pt.pk)
                        match_set30.add(ph_pt2.pk)
                        match_set30.add(ph_pt.pk)
                    if eucl_dist(ph_pt, ph_pt2) < 2.5:
                        match_set25.add(ph_pt2.pk)
                        match_set25.add(ph_pt.pk)
                        match_set30.add(ph_pt2.pk)
                        match_set30.add(ph_pt.pk)
                        match_set35.add(ph_pt2.pk)
                        match_set35.add(ph_pt.pk)
                    if eucl_dist(ph_pt, ph_pt2) < 1.5:
                        match_set15.add(ph_pt2.pk)
                        match_set15.add(ph_pt.pk)
                        match_set25.add(ph_pt2.pk)
                        match_set25.add(ph_pt.pk)
                        match_set30.add(ph_pt2.pk)
                        match_set30.add(ph_pt.pk)
                        match_set35.add(ph_pt2.pk)
                        match_set35.add(ph_pt.pk)
                    if eucl_dist(ph_pt, ph_pt2) < 1.0:
                        match_set10.add(ph_pt2.pk)
                        match_set10.add(ph_pt.pk)
                        match_set15.add(ph_pt2.pk)
                        match_set15.add(ph_pt.pk)
                        match_set25.add(ph_pt2.pk)
                        match_set25.add(ph_pt.pk)
                        match_set30.add(ph_pt2.pk)
                        match_set30.add(ph_pt.pk)
                        match_set35.add(ph_pt2.pk)
                        match_set35.add(ph_pt.pk)

        #  Find the number of points in each set
        mt_set10 = len([x for x in PharmaPoint.objects.filter(smiles__in=allowed_list, mol_id__in=[myfrag.over_mol_id, myfrag.mol_id]) if x.pk not in match_set10])
        mt_set15 = len([x for x in PharmaPoint.objects.filter(smiles__in=allowed_list, mol_id__in=[myfrag.over_mol_id ,myfrag.mol_id]) if x.pk not in match_set15])
        mt_set25 = len([x for x in PharmaPoint.objects.filter(smiles__in=allowed_list, mol_id__in=[myfrag.over_mol_id, myfrag.mol_id]) if x.pk not in match_set25])
        mt_set30 = len([x for x in PharmaPoint.objects.filter(smiles__in=allowed_list, mol_id__in=[myfrag.over_mol_id, myfrag.mol_id]) if x.pk not in match_set30])
        mt_set35 = len([x for x in PharmaPoint.objects.filter(smiles__in=allowed_list, mol_id__in=[myfrag.over_mol_id, myfrag.mol_id]) if x.pk not in match_set35])

        # Now loop throigh the pharmacophore points on each molecule and make this ActPharmaPoint
        for i, ph_pt in enumerate(over_mol_ph4):
        # Add this point in and then
            new_act_pharma = ActPharmaPoint()
            new_act_pharma.pharma_id = ph_pt
            try:
                new_act_pharma.dist_same = dist_dict[ph_pt.pk]
            except KeyError:
                new_act_pharma.dist_same = 100000000.0
            new_act_pharma.act_3d_frag_id = myfrag
            new_act_pharma.from_xtal = True
            new_act_pharma.num_diff_10 = mt_set10
            new_act_pharma.in_diff_10 = ph_pt.pk not in match_set10
            new_act_pharma.num_diff_15 = mt_set15
            new_act_pharma.in_diff_15 = ph_pt.pk not in match_set15
            new_act_pharma.num_diff_25 = mt_set25
            new_act_pharma.in_diff_25 = ph_pt.pk not in match_set25
            new_act_pharma.num_diff_30 = mt_set30
            new_act_pharma.in_diff_30 = ph_pt.pk not in match_set30
            new_act_pharma.num_diff_35 = mt_set35
            new_act_pharma.in_diff_35 = ph_pt.pk not in match_set35
            new_act_pharma.target_id = my_target
            new_act_pharma.act_change = diff_act
            new_act_pharma.type = opt_type
            try:
                new_act_pharma.validate_unique()
                new_act_pharma.save()
            except ValidationError:
                pass

        # And the other molecule
        for i, ph_pt in enumerate(mol_ph4):
            new_act_pharma = ActPharmaPoint()
            new_act_pharma.pharma_id = ph_pt
            try:
                new_act_pharma.dist_same = dist_dict[ph_pt.pk]
            except KeyError:
                new_act_pharma.dist_same = 100000000.0
            new_act_pharma.act_3d_frag_id = myfrag
            new_act_pharma.from_xtal = True
            new_act_pharma.num_diff_10 = mt_set10
            new_act_pharma.in_diff_10 = ph_pt.pk not in match_set10
            new_act_pharma.num_diff_15 = mt_set15
            new_act_pharma.in_diff_15 = ph_pt.pk not in match_set15
            new_act_pharma.num_diff_25 = mt_set25
            new_act_pharma.in_diff_25 = ph_pt.pk not in match_set25
            new_act_pharma.num_diff_30 = mt_set30
            new_act_pharma.in_diff_30 = ph_pt.pk not in match_set30
            new_act_pharma.num_diff_35 = mt_set35
            new_act_pharma.in_diff_35 = ph_pt.pk not in match_set35
            new_act_pharma.act_change = -diff_act
            new_act_pharma.target_id = my_target
            new_act_pharma.type = opt_type
            try:
                new_act_pharma.validate_unique()
                new_act_pharma.save()
            except ValidationError:
                pass
    old = 100
    sys.stdout.write("\r%d%% complete..." % old)
    sys.stdout.flush()
    print "\nCompleted finding pharmacophore points"
    return None


def find_frag_com(frags=None):
    """Function to find the centre of mass of all MMPFrag objects
    Takes no args
    Returns None"""
    # Get all single cut MMPFrag objects related to an xtal structure with a reasonable core size
    if frags is None:
        frags = MMPFrag.objects.exclude(smistore__contains="[*:2]").exclude(mol_id__smiles="").filter(core_size__lte=10)
    tot = len(frags)
    old = -1
    for i, frag in enumerate(frags):
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
        try:
            frag.x_com, frag.y_com, frag.z_com = centre_of_mass_from_SD_block(frag.sdf_info)
        except:
            if frag.sdf_info == "":
                continue
        frag.save()


def find_context_compound():
    """Function to add foreign key link to MMP object for Compound object
    Takes no args
    Returns None"""
    frags = MMPFrag.objects.exclude(smistore__contains="[*:2]").filter(core_size__lte=10)
    for i, frag in enumerate(frags):
        print i
        mmp = frag.mmp_link
        if frag.mol_id is not None:
            cmpd = frag.mol_id.cmpd_id
        elif frag.act_id is not None:
            cmpd = frag.act_id.cmpd_id
        mmp.compound_id.add(cmpd)
        mmp.save()


def find_pharma_changes(target_id, my_type):
    """Function to find the relevant pharamcophore changes and store in the DB
    Takes a Target and an activity change
    Returns None"""
    print "FINDING RELEVANT SHIFTS"
    act_list, inact_list, shape_list = find_points(target_id, my_type)
    print "MAKING ACTIVITY MAPS"
    targ = Target.objects.get(pk=target_id)
    commit_pharma_changes(act_list, inact_list, shape_list, my_type, targ)


def render_pdb_data(type_dict, item, index, option=None):
    """Function to render atomic data into a HETATM record
    Takes a pharma type, an occupancy, a B-factor, x,y and z coords, a
    dictionary of atom types to convert to and an index of the atom
    Returns a string"""
    if option == "SHAPE":
        my_s = "HETATM"+" "*(5-len(str(index+1)))+str(index+1)+" "+type_dict[item[3]]+" "*(5-len(type_dict[item[3]]))+type_dict[item[3]]+" "*(5-len(type_dict[item[3]]))+" "*(4-len(str(index+1)))+str(index+1)+" "*(12-len("{0:.3f}".format(item[4])))+"{0:.3f}".format(item[4])+" "*(8-len("{0:.3f}".format(item[5])))+"{0:.3f}".format(item[5])+" "*(8-len("{0:.3f}".format(item[6])))+"{0:.3f}".format(item[6])+" "*(6-len("{0:.2f}".format(min(9.0,math.fabs(item[7]))*0.1)))+"{0:.2f}".format(min(9.0,math.fabs(item[7]))*0.1)+" "*(6-len("{0:.2f}".format(item[8])))+"{0:.2f}".format(item[8])+" "*(12-len(str(type_dict[item[3]])))+type_dict[item[3]]+"1+\n"
    elif option is None:
        my_s = "HETATM"+" "*(5-len(str(index+1)))+str(index+1)+" "+type_dict[item[3]]+" "*(5-len(type_dict[item[3]]))+type_dict[item[3]]+" "*(5-len(type_dict[item[3]]))+" "*(4-len(str(index+1)))+str(index+1)+" "*(12-len("{0:.3f}".format(item[4])))+"{0:.3f}".format(item[4])+" "*(8-len("{0:.3f}".format(item[5])))+"{0:.3f}".format(item[5])+" "*(8-len("{0:.3f}".format(item[6])))+"{0:.3f}".format(item[6])+" "*(6-len("{0:.2f}".format(min(1.0,math.fabs(item[7])))))+"{0:.2f}".format(min(1.0,math.fabs(item[7])))+" "*(6-len("{0:.2f}".format(item[8])))+"{0:.2f}".format(item[8])+" "*(12-len(str(type_dict[item[3]])))+type_dict[item[3]]+"1+\n"
    my_s = "COMPND    UNNAMED\nAUTHOR    GENERATED BY OPEN BABEL 2.3.2\n" + my_s
    my_s = my_s + "MASTER        0    0    0    0    0    0    0    0    1    0    1    0\nEND\n"
    return my_s


def commit_pharma_changes(act_list, inact_list, shape_list, activity_type, target_id):
    """Function to store information about activity,inactivity and shape points in the DB
    PDB file is written of all the points and ActMapPoint objects are stored.
    Takes a list of activity, incactivty and shape points, activity type and Target"""
    # This should be replaced by making them into PharmaPoints in the DB with their assoicated meta-data
    lists = [act_list, inact_list, shape_list]
    # Define the points to draw
    my_fs = ["act", "inact", "shape"]
    old_list =[]
    comps = []
    i = 0
    # Dictionary to relate Atom types to feature types
    type_dict = {u'SHAPE': 'Na', u'SingleAtomAcceptor': 'Li',
                 u'SingleAtomDonor': 'Na', u'RH5_5': 'K',
                 u'RH6_6': 'K', u'Arom6': 'Zn', u'Arom5': 'Zn',
                 u'ThreeWayAttach': 'K'}
    # Loop through act,inact and non
    for j, my_list in enumerate(lists):
        print "MAKING MAP ", j + 1, " of ", len(lists)
        # Make a new map object and initialise its values
        new_map = MMPDiffMap()
        new_map.type = my_fs[j]
        new_map.num_change = 10
        new_map.activity_change = 9.0
        new_map.activity_type = activity_type
        new_map.target_id = target_id
        new_map.pdb_info = ""
        try:
            new_map.validate_unique()
            new_map.save()
        except ValidationError:
            continue
        # Count the individual points in the map
        counter =0
        # Loop through the points
        old = -1
        tot = len(my_list)
        for k, item in enumerate(my_list):
            if k * 100 / tot != old:
                old = k * 100 / tot
                sys.stdout.write("\r%d%% complete..." % old)
                sys.stdout.flush()
            # If the pharmacophore point is not represeented by a metal then skip it
            try:
                type_dict[item[3]]
            except KeyError:
                continue
            # This is a point to be stored in the map
            new_point = ActMapPoint()
            # Associate it to it's underlying PharmaPoint
            new_point.pharma_id = PharmaPoint.objects.get(pk=item[2])
            counter +=1
            # if it's a new combination make a new MMP combo pair
            if (item[0], item[1]) != old_list:
                new_com = MMPComparison()
                # Now inititalise the unique identifiers
                new_com.num_change = int(item[8])
                new_com.activity_change = abs(float(item[7]))
                new_com.activity_type = activity_type
                new_com.target_id = target_id
                new_com.type = new_map.type
                comps.append(new_com)
                # A list to hold the indices in
                new_com.my_list = []
                # Add this combination to the checking list
                old_list = (item[0], item[1])
                # Make a new string IO object to write the molecular information
                new_com.out_s = cStringIO.StringIO()
                # Make an RDKIT SD Writer object for the visualisations
                new_com.c_f = Chem.SDWriter(new_com.out_s)
                # Add in the foreign key links to the molecules
                new_com.xtal_mol = Molecule.objects.get(pk=item[0])
                new_com.chembl_mol = Molecule.objects.get(pk=item[1])
                # Add in the foreign key links to the activity points
                act_points = ActivityPoint.objects.filter(cmpd_id=new_com.xtal_mol.cmpd_id, target_id=target_id, source__in=ast.literal_eval(activity_type)).order_by("activity")
                if len(act_points) > 1:
                    # If the user has added multiple points. The the biggest one is taken
                    new_com.xtal_act = act_points[len(act_points) - 1]
                elif len(act_points) == 0:
                    sys.stderr.write("No activity points in list")
                else:
                    new_com.xtal_act = act_points[0]
                # Now do this for the other compounds
                act_points = ActivityPoint.objects.filter(cmpd_id=new_com.chembl_mol.cmpd_id, target_id=target_id, source__in=ast.literal_eval(activity_type)).order_by("activity")
                if len(act_points) > 1:
                    # If the user has added multiple points. The the biggest one is taken
                    new_com.chembl_act = act_points[len(act_points) - 1]
                elif len(act_points) == 0:
                    sys.stderr.write("No activity points in list")
                else:
                    new_com.chembl_act = act_points[0]
                # Now write the SDInfo to this st
                new_com.c_f.write(Chem.MolFromMolBlock(str(Molecule.objects.get(pk=item[0]).sdf_info)))
                new_com.c_f.write(Chem.MolFromMolBlock(str(Molecule.objects.get(pk=item[1]).sdf_info)))
                new_com.sdf_info = ""
                try:
                    new_com.validate_unique()
                    new_com.save()
                    new_com_pk = new_com.pk
                except ValidationError:
                  # Find the comparison PK
                    new_com_pk = MMPComparison.objects.get(target_id=new_com.target_id,
                    num_change=new_com.num_change, activity_change=new_com.activity_change,
                    activity_type=new_com.activity_type, xtal_mol=new_com.xtal_mol,
                    chembl_mol=new_com.chembl_mol, type=new_com.type).pk

            new_point.mmp_comparsion_id = MMPComparison.objects.get(pk=new_com_pk)
            # Now add the act change and numchange information to this
            if new_map.type != "shape":
              # Two strings to display the HETATM records
                my_s = render_pdb_data(type_dict, item, i, option="SHAPE")
            else:
                my_s = render_pdb_data(type_dict, item, i)
            # Now add that this is linked to that point in the map
            new_com.my_list.append(my_fs[j] + str(counter))
            new_map.pdb_info += my_s
            new_map.save()
            new_point.point_id = counter
            new_point.map_id = new_map
            new_point.validate_unique()
            new_point.save()
        # Save the map to the DB
        new_map.save()
        old = 100
        sys.stdout.write("\r%d%% complete..." % old)
        sys.stdout.flush()
        print "\nMADE MAP ", j + 1, " of ", len(lists)
    # Now write the completed SD file to the comp objects
    for cm in comps:
        cm.c_f.flush()
        cm.sdf_info = cm.out_s.getvalue()
        cm.map_inds = str(cm.my_list)
        try:
            cm.validate_unique()
            cm.save()
        except ValidationError:
            pass
    print "COMPLETED MAKING MAPS"
    return None

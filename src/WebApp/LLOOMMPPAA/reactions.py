### List of functions to identify and peform reactions
from IOhandle.models import Target, Molecule
from IOhandle.functions import add_new_comp
from rdkit import Chem
from rdkit.Chem import AllChem
from LLOOMMPPAA.models import Reaction, PossReact, Reactant, ReactantLib, ProductLib, Product, Process, ReactionProcess
from django.core.exceptions import ValidationError
from LLOOMMPPAA.functions import make_synth_points, do_diff_analysis
from PLIFS.functions import restore_probs
import uuid
from djutils.decorators import async
import sys


def define_reacts():
    """Function to define all possible reactions"""
    # Take the hits and find all AMIDES, THIOUREAS, SONAGASHIRA couplings
    funct_g_dict = {"AMIDES": {"RETRO_SMARTS": "[C:1](=[O:2])[NX3;R0:3]>>[C:1](=[O:2])-[ClA1].[N1:3]", "REACT_SMARTS": "[C:1](=[O:2])-[ClD1,OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]", "CONT": '["[N1:1]>>[N1:1][CA1](=[OA2])","[1ClD1:1]>>[NA1:1]"]'},
    "THIOUREA": {"RETRO_SMARTS": "[NX3;R0:1][CX3:2](=[SX1:3])[NX3:4]>>[N1:1].[ClA1][CX3:2](=[SX1:3])[NX3:4]", "REACT_SMARTS": "[N!H0:1].[ClD1][CX3:2](=[SX1:3])[NX3:4]>>[NX3:1][CX3:2](=[SX1:3])[NX3:4]", "SMARTS": "[NX3;R0][CX3](=[SX1])[NX3]", "COORDS": [0, 1], "CONT": '["[N1:1]>>[N1:1][CX3A1](=[SX1A2])[NX3A3]","[1ClD1:1]>>[NA1:1]"]'},
    # This one needs sorting - it will do the retro step - but the forward step will not be at all selective
    "ARYL COUPLING": {"RETRO_SMARTS": "[c:1]!@[c:2]>>[c:1].[c:2]", "REACT_SMARTS": "[c:1].[c:2]>>[c:1]-[c:2]", "SMARTS": "[c]!@[c]", "CONT": "[]", "COORDS": [0, 1]}
    }
    for item in funct_g_dict:
        my_reaction = Reaction()
        my_reaction.name = item
        try:
            my_reaction.validate_unique()
        except ValidationError:
            my_reaction = Reaction.objects.get(name=item)
        my_reaction.react_smarts = funct_g_dict[item]["REACT_SMARTS"]
        my_reaction.retro_smarts = funct_g_dict[item]["RETRO_SMARTS"]
        my_reaction.cont_smarts = funct_g_dict[item]["CONT"]
        my_reaction.save()


def find_context(atom_ids, my_mol, my_sub_m, cont_sub):
    """Function to get the context - by adding all the atoms in a substructure and adding them to the fragment"""
    # First combine the tuples
    # Get the atoms in this match
    out_m = []
    for x in my_mol.GetSubstructMatches(cont_sub):
        for y in x:
            # If it's in both substructures
            if y in my_sub_m:
                out_m.append(y)
    # Add the core structure
    out_m.extend(atom_ids)
    # get an editable mol
    em = Chem.EditableMol(my_mol)
    # Loop through the atoms and remove them
    for my_id in my_mol.GetAtoms():
        if my_id.GetIdx() in out_m:
            continue
        else:
            # Replace with a gap
            em.ReplaceAtom(my_id.GetIdx(), Chem.Atom(0))
    # Now remove all these atoms
    star_mol = em.GetMol()
    out_mol = Chem.DeleteSubstructs(star_mol, Chem.MolFromSmarts('[#0]'))
    out_ans = Chem.MolToSmiles(out_mol)
    print "ORIG MOL", Chem.MolToSmiles(my_mol)
    print "FRACT MOL", out_ans
    return out_ans


def remove_zero_atoms(mol_block):
    "Function to take a molecule and remove all atoms with no mass"
    rdmol = Chem.MolFromMolBlock(str(mol_block))
    conf = rdmol.GetConformer()
    em = Chem.EditableMol(rdmol)
    out_lines = []
    [em.ReplaceAtom(x.GetIdx(), Chem.Atom(0)) for x in rdmol.GetAtoms() if [0.0, 0.0, 0.0] == [conf.GetAtomPosition(x.GetIdx()).x,conf.GetAtomPosition(x.GetIdx()).z, conf.GetAtomPosition(x.GetIdx()).z ]]
    star_mol = em.GetMol()
    out_mol = Chem.DeleteSubstructs(star_mol, Chem.MolFromSmarts('[#0]'))
    return Chem.MolToMolBlock(out_mol)


def replace_isotope(smiles):
    """Replace an added atom with an attachment point"""
    rdmol = Chem.MolFromSmiles(str(smiles))
    em = Chem.EditableMol(rdmol)
    for my_id in rdmol.GetAtoms():
        if my_id.GetIsotope() == 0:
            continue
        else:
            # Replace with a gap
            em.ReplaceAtom(my_id.GetIdx(), Chem.Atom(0))
    # Now remove all these atoms
    star_mol = em.GetMol()
    out_mol = Chem.DeleteSubstructs(star_mol, Chem.MolFromSmarts('[#0]'))
    return Chem.MolToSmiles(out_mol, isomericSmiles=True)


# Identify reactions possible for the molecules
def find_follow_ups(target_id=2):
    """Function to find all the fragmentation points for a given targets molecules"""
    # Get the target
    target = Target.objects.get(pk=target_id)
    # Now loop through the compounds in this compound
    out_d = {}
    res_d = {}
    my_mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)
    # Make a unique dict of smiles
    for item in my_mols:
        if item.smiles in out_d:
            out_d[item.smiles].append(item.pk)
        else:
            out_d[item.smiles] = [item.pk]
            res_d[item.smiles] = {}
    # Get the reactions
    my_reacts = Reaction.objects.all()
    # Now loop through this dict and add the matches
    import ast
    for item in out_d:
        my_mol = Chem.MolFromMolBlock(str(Molecule.objects.get(pk=out_d[item][0]).sdf_info))
        for funct_g in my_reacts:
#            my_sub = Chem.MolFromSmarts(funct_g.funct_smarts)
#            # my_
            my_cont = ast.literal_eval(str(funct_g.cont_smarts))
            reacts = []
            for cont in my_cont:
                my_rxn = AllChem.ReactionFromSmarts(cont)
                my_rxn.Initialize()
                reacts.append(my_rxn)
            # Now get the the substructure match#
            if funct_g.retro_smarts:
                retro_rxn = AllChem.ReactionFromSmarts(str(funct_g.retro_smarts))
                these_matches = retro_rxn.RunReactants((my_mol,))
                if not these_matches:
                    continue
            else:
                continue
            # If we have matches -> then produce the objects
            if these_matches:
                for my_match in these_matches:
#                    my_match = mol_and_ids[0]
#                    atom_ids = mol_and_ids[1]
#                    sub_m = mol_and_ids[2]
                    if len(my_match) != 2:
                        print Chem.MolToSmiles(my_match[0], isomericSmiles=True)
                        print Molecule.objects.get(pk=out_d[item][0]).smiles
                    new_poss_react = PossReact()
                    mol_id = Molecule.objects.get(pk=out_d[item][0])
                    retained_frag = Chem.MolToSmiles(my_match[0], isomericSmiles=True)
                    replaced_frag = Chem.MolToSmiles(my_match[1], isomericSmiles=True)
                    new_poss_react.react_id = funct_g
                    new_poss_react.mol_id = mol_id
                    # Now get one half of the retro synthesis
                    new_poss_react.retained_frag = retained_frag
                    # Now get the other half of the retro synthesis
                    new_poss_react.replaced_frag = replaced_frag
                    try:
                        new_poss_react.validate_unique()
                    except ValidationError:
                        new_poss_react = PossReact.objects.get(react_id=funct_g, mol_id=mol_id, retained_frag=retained_frag, replaced_frag=replaced_frag)
                    new_poss_react.ret_frag_sdf_info = remove_zero_atoms(Chem.MolToMolBlock(my_match[0]))
                    new_poss_react.rep_frag_sdf_info = remove_zero_atoms(Chem.MolToMolBlock(my_match[1]))
                    # Now apply all the reactions 
                    in_ret = Chem.MolFromSmiles(str(retained_frag))
                    in_rep = Chem.MolFromSmiles(str(replaced_frag))
                    for rxn in reacts:
                        products_ret = rxn.RunReactants((in_ret,))
                        products_rep = rxn.RunReactants((in_rep,))
                        if products_ret:
                            in_ret = products_ret[0][0]
                        if products_rep:
                            in_rep = products_rep[0][0]
                    # Now set these
                    new_poss_react.retained_frag_frag = str(retained_frag)
                    # Now find the context of this fragment - the compound with the full group on the other side
                    new_poss_react.retained_frag_context = Chem.MolToSmiles(in_ret)
                    new_poss_react.replaced_frag_frag = str(replaced_frag)
                    new_poss_react.replaced_frag_context = Chem.MolToSmiles(in_rep)
                    new_poss_react.save()
    return None


def load_in_reagents(library_title, path_to_file, reaction):
    """Function to load up a library of reagents - e.g. acyl chlorides"""
    # Open the file
    my_mols = Chem.SDMolSupplier(path_to_file)
    # Get the library name this corresponds to
    if len(ReactantLib.objects.filter(library_name=library_title)) == 0:
        out_lib = ReactantLib()
        out_lib.library_name = library_title
        out_lib.save()
    else:
        out_lib = ReactantLib.objects.filter(library_name=library_title)[0]
    # Go through the compounds
    for rdmol in my_mols:
        # Add them as a django thingy
        dj_comp = add_new_comp(rdmol)
        # Link it to the reaction
        my_r = Reactant()
        my_r.cmpd_id = dj_comp
        try:
            my_r.validate_unique()
            my_r.save()
        except ValidationError:
            my_r = Reactant.objects.get(cmpd_id=dj_comp)
        # Now update this
        my_r.is_available = True
        my_r.react_id.add(reaction)
        my_r.save()
        out_lib.reactant_id.add(my_r)
        out_lib.save()
    return out_lib

def load_in_follow_ups(library_title, path_to_file, reaction, mol_id):
    """Function to load in follow ups - ready made"""
    # Open the file
    my_mols = Chem.SDMolSupplier(path_to_file)
    # Get the library name this corresponds to
    if len(ProductLib.objects.filter(lib_name=library_title)) == 0:
        out_lib = ProductLib()
        out_lib.lib_name = library_title
        out_lib.save()
    else:
        out_lib = ProductLib.objects.filter(lib_name=library_title)[0]
    # Get the process
    my_process = Process()
    my_process.mol_id = mol_id
    my_process.is_made_lloommppaa = False
    my_process.reaction_id = reaction
    my_process.save()
    # Go through the compounds
    for rdmol in my_mols:
        # Add them as a django thingy
        dj_comp = add_new_comp(rdmol)
        # Link it to the reaction
        my_p = Product()
        my_p.cmpd_id = dj_comp
        try:
            my_p.validate_unique()
            my_p.save()
        except ValidationError:
            my_p = Product.objects.get(cmpd_id=dj_comp)
        # Now update this
        my_p.process_id.add(my_process)
        my_p.save()
        out_lib.product_id.add(my_p)
        out_lib.save()
    return out_lib

def define_reaction_process(mol_id, context, my_prots, react_id, my_context, react_frag, products_id=None, reactants_id=None):
    """Function to define what we should do next"""
    # Define the reaction process
    react_proc = ReactionProcess()
    # Define the molecule to come off
    react_proc.mol_id = mol_id
    # Which potential reaction to use
    react_proc.react_id = react_id
    # The fragment to use in the reaction
    react_proc.react_frag = react_frag
    # Which context to use - should be stripped of fragment points
    react_proc.context = my_context
    # Now save or get this
    react_proc.stage_completion = 0
    try:
        react_proc.validate_unique()
        react_proc.save()
    except ValidationError:
        react_proc = ReactionProcess.objects.get(mol_id=mol_id, react_id=react_id, context=context)
        print "Warning - Reaction has been run before"
    # Which compounds to use
    if products_id:
        react_proc.is_reaction = False
        react_proc.products_id.add(products_id)
        # Now add these as foreign key links to the queue
        for x in products_id.product_id.all():
            react_proc.product_queue.add(x.cmpd_id)
            react_proc.save()
    if reactants_id:
        react_proc.is_reaction = True
        react_proc.reactants_id.add(reactants_id)
        # Now add these as foreign key links to the queue
        for x in reactants_id.reactant_id.all():
            react_proc.reactant_queue.add(x.cmpd_id)
            react_proc.save()
    # If nothing then get rid of them
    if not products_id and not reactants_id:
        print "You must specifiy a library to use"
        return None
    # Define the proteins to use as clashes
    for prot_id in my_prots:
        react_proc.prot_id.add(prot_id)
    react_proc.save()
    return react_proc


def create_lib(rxn, react_proc, lib_name):
    """Function to  create a library from a reaction"""
    from LLOOMMPPAA.pains_filter import pains_test
    # Make the molecule fit for reaction
    rdmol = Chem.MolFromSmiles(str(react_proc.react_frag))
    p_lib = ProductLib()
    p_lib.lib_name = str(uuid.uuid4())
    p_lib.save()
    # Get the process
    my_process = Process()
    my_process.mol_id = react_proc.mol_id
    my_process.is_made_lloommppaa = False
    my_process.reaction_id = react_proc.react_id
    my_process.save()
    # Loop through the library
    # Get the lib
    my_cmpd = react_proc.reactant_queue.all()
    tot = len(my_cmpd)
    old = -1
    for i, cmpd in enumerate(my_cmpd):
        sys.stdout.write("\rCarried out reaction %d of %d..." % (i, tot))
        sys.stdout.flush()
        re_mol = Chem.MolFromSmiles(str(cmpd.smiles))
        out_prods = rxn.RunReactants((re_mol, rdmol))
        if len(out_prods) == 0:
            out_prods = rxn.RunReactants((rdmol, re_mol))
            if len(out_prods) == 0:
                print "NO PRODUCTS"
                print Chem.MolToSmiles(rdmol, isomericSmiles=True)
                print Chem.MolToSmiles(re_mol, isomericSmiles=True)
                continue
        products = out_prods[0]
        if len(products) > 1:
            print "MULTIPLE PRODUCTS"
            print products
            print Chem.MolToSmiles(rdmol)
            print Chem.MolToSmiles(re_mol)
            continue
        if pains_test(products[0]):
            print "PAINS FILTER SKIPPING!!!"
            continue
        # Register the compound
        dj_comp = add_new_comp(products[0])
        # Add it to the list of products
        my_prod = Product.objects.filter(cmpd_id=dj_comp)
        if my_prod:
            my_prod = my_prod[0]
        else:
            my_prod = Product()
            my_prod.cmpd_id = dj_comp
            my_prod.save()
        my_prod.process_id.add(my_process)
        # Add the product to the product library
        p_lib.product_id.add(my_prod)
        p_lib.save()
        # Add it to the product queue
        react_proc.product_queue.add(dj_comp)
        react_proc.reactant_queue.remove(cmpd)
        my_prg = int((float(i) / float(tot)) * 100)
        if my_prg != old:
            react_proc.stage_completion = my_prg
            old = my_prg
            react_proc.save()
    react_proc.products_id.add(p_lib)
    react_proc.save()
    # Now return this
    return react_proc


def run_reaction(react_proc):
    """Function to take a ReactionProcess object - and run a reaction"""
    # Get the reaction smarts
    react_proc.proc_stage = "RUN REACTION"
    react_proc.stage_completion = 0
    react_proc.save()
    try:
        my_rxn = AllChem.ReactionFromSmarts(str(react_proc.react_id.react_smarts))
    except ValueError:
        print "Invalid reaction", react_proc.react_id.name
        print "Invalid snmats", react_proc.react_id.react_smarts
        return None
    # Get the library of compounds
    if not react_proc.reactant_queue.all():
        pass
    else:
        # Lib name is a combination of
        lib_name = react_proc.title + react_proc.react_id.name
        react_proc = create_lib(my_rxn, react_proc, lib_name)
    # Loop through the crystal structures
    react_proc.proc_stage = "GENERATE CONFS"
    react_proc.stage_completion = -1
    react_proc.save()


def gen_confs(react_proc):
    make_synth_points(target_id=react_proc.mol_id.prot_id.target_id.pk, react_proc=react_proc)
    # Now signal this ready to  do analysis for
    react_proc.proc_stage = "MAKING PROTEIN PROBES"
    react_proc.stage_completion = -1
    react_proc.save()


@async
def init_react_proc(react_proc):
    """Function to intialise a reaction process - from a view"""
    print "ADDING COMPOUNDS"
    react_proc.proc_stage = "ADD REACTANTS"
    react_proc.stage_completion = 0
    react_proc.save()
    my_cmpds = react_proc.reactants_id.filter().values_list("reactant_id__cmpd_id", flat=True)
    tot = len(my_cmpds)
    old = -1
    for i, cmpd in enumerate(my_cmpds):
        sys.stdout.write("\rLoaded reactant %d of %d..." % (i, tot))
        sys.stdout.flush()
        react_proc.reactant_queue.add(cmpd)
        my_prg = int((float(i) / float(tot)) * 100)
        if my_prg != old:
            react_proc.stage_completion = my_prg
            old = my_prg
            react_proc.save()
    sys.stdout.write("\rLoaded reactants")
    sys.stdout.flush()    # Now do the same for the products
    my_cmpds = react_proc.products_id.filter().values_list("product_id__cmpd_id", flat=True)
    tot = len(my_cmpds)
    react_proc.proc_stage = "ADD PRODUCTS"
    react_proc.stage_completion = 0
    react_proc.save()
    old = -1
    for i, cmpd in enumerate(my_cmpds):
        sys.stdout.write("\rLoaded product %d of %d..." % (i, tot))
        sys.stdout.flush()
        react_proc.product_queue.add(cmpd)
        my_prg = int((float(i) / float(tot)) * 100)
        if my_prg != old:
            react_proc.stage_completion = my_prg
            old = my_prg
            react_proc.save()
    sys.stdout.write("\rLoaded products")
    sys.stdout.flush()
    react_proc.proc_stage = "RUN REACTION"
    react_proc.stage_completion = -1
    react_proc.save()


def find_procs_to_run():
    """Function to find processes to run"""
    import time
    for i in range(30):
        time.sleep(1)
        print i + 1, "seconds"
        # Get all of the process  that need to run from the beginning
        rp = ReactionProcess.objects.filter(proc_stage = "RUN REACTION", stage_completion=-1)
        # Run the first in the list
        if rp:
            run_react_proc(rp[0])
        else:
            # Now we need to check if any process have stopped - this will be manually decided for now
            rp = ReactionProcess.objects.filter(stage_completion=-1)
            if rp:
                run_react_proc(rp[0])


def find_dead_procs():
    """Find processes still running -> in RUN REACTANTS or GENERATE CONFS
    where stage completion has not changed"""
    import time
    # First find a dict of the ReactionProcess -> hold this dict
    rp = ReactionProcess.objects.filter(proc_stage__in = ["RUN REACTION", "GENERATE CONFS", "MAKING PROTEIN PROBES",
                  "MAKING MOLECULE PROBES", "MAKING LLOOMMPPAA PROBES", "FINDING INTERACTIONS"])
    store_d = {}
    for r in rp:
        store_d[r.pk] = {"PROCESS": r.proc_stage, "COMPLETE": r.stage_completion}
    # Now wait for 20 minutes
    print "WAITING"
    for i in range(1200):
        time.sleep(1)
        print i + 1, "seconds"
    for item in store_d:
        rp = ReactionProcess.objects.filter(pk=item, proc_stage=store_d[item]["PROCESS"], stage_completion=store_d[item]["COMPLETE"])
        if rp:
            for r in rp:
                print "FOUND ONE", r.pk
                r.stage_completion = -1
                r.save()


def restart_and_find_dead_procs():
    """Function to find dead procs and restart them"""
    find_dead_procs()
    find_procs_to_run()


def do_analysis(react_proc):
    """Function to do the analysis for a given  reaction process"""
    # Now do the analysis
    react_proc.proc_stage = "DO ANALYSIS"
    react_proc.stage_completion = 0
    react_proc.save()
    do_diff_analysis(react_proc.mol_id.prot_id.target_id.pk, num_iters=1000, num_trial=20, clash=-1.0, rmsd=0.5, shape_dist=1.0, react_id=react_proc, redo=True)


def run_react_proc(react_proc):
    """Function to just run the reaction  - THIS IS THE MAIN FUNCTION TO DO THIS!!!!"""
    from PLIFS.interaction_defs import init_schemes
    # First add all the things to  the product and reaction queue
    print "STARTING REACTION PROCESS"
    # Now run the reactions
    if react_proc.proc_stage == "RUN_REACTION":
        run_reaction(react_proc)
    if react_proc.proc_stage == "GENERATE CONFS":
    # Now generate the conformations
        gen_confs(react_proc)
    # Now generate the PLIFS
    if react_proc.proc_stage in ["FINDING INTERACTIONS", "MAKING PROTEIN PROBES", "MAKING MOLECULE PROBES", "MAKING LLOOMMPPAA PROBES"]:
        restore_probs(react_proc.mol_id.prot_id.target_id.pk, opt="LLOOMMPPAA", clash=-1.0, rmsd=0.5, shape_dist=1.0, react_id=react_proc)
    # Now initialse the schemes
    if react_proc.proc_stage == "MAKE INT SCHEMES":
        print "MOVING ONTO INTERACTION SCHEMES"
        init_schemes(refresh=True, react_id=react_proc)
    # Now do the analysis
    #do_analysis(react_proc)
    react_proc.proc_stage = "COMPLETE"
    react_proc.stage_completion = 100
    react_proc.save()

from IOhandle.models import Compound,Molecule,Protein,ActivityPoint,Project,Target
from LLOOMMPPAA.models import Reaction
from LLOOMMPPAA.views import get_sdf_file
from MMPMaker.functions import index_hydrogen_change, make_mol_mmp, make_mmp_database, make_list_mmps, act_mmp_3d, make_ph4_difference_points, find_pharma_changes
from MMPMaker.models import MMPDiffMap, MMPComparison, ActPharmaPoint, MMPFrag
from Group.models import Group
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from django.core.exceptions import ValidationError
import sys, os
from IOhandle.functions import add_new_comp
import csv
import uuid, glob
import ast
from WONKA.models import ResShift
import numpy
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import MCS
# Global variable to define what things can be used possibly
__POSS_LIST__ = ["IC50", "Ki", "TM"]


def check_smiles_or_SD(file_path):
    """Function to check if a is a valid smiles or SD file
    Takes a file path
    Returns a lits of Molecules or None"""
    if os.path.isfile(file_path):
        pass
    else:
        print "NO SUCH FILE: ", file_path
        return None
    try:
        mols = Chem.SDMolSupplier(file_path)
    except IOError:
        mols = None
    try:
        mols = Chem.SmilesMolSupplier(file_path)
    except IOError:
        mols = None

    if mols == None:
        print "NOT VALID SMILE OR SD FILE: ", file_path
    else:
        return mols


def read_CSV(file_in):
    """Function to read a csv file.
    Takes a file path
    Returns a csv dict object"""
    return csv.DictReader(file_in)

def get_sdf_file_ll(target_id):
    get_sdf_file(target_id)

def add_new_mol(rdmol, target):
    """Function to add a new bound Molecule object
    Takes an RDKit molecule and a Target
    Returns None"""
    import re
    new_mol = Molecule()
    rdProps = rdmol.GetProp("_Name").split("_")
    comp_ref = add_new_comp(rdmol)
    if comp_ref is None:
        return None
    # To get rid of the .pdb suffix
    pdb_id = rdProps[0].split(".")[0]
    # If it's an SGC model entry -> rename it appropriately
    if rdmol.HasProp("name"):
        if re.match("^m\d\d\d$", rdmol.GetProp("name")):
            pdb_id = rdmol.GetProp("name") + "_" + str(target.pk)
    # Check that the name is unique and more than 3 characters long
    # and doesn't contain the target.title
    if len(pdb_id) < 4 or len(Protein.objects.filter(code=pdb_id)) > 0 or target.title in pdb_id:
        # Make a new uniqid
        # First up check that this molecule has not been added before
        mols = [[Chem.MolFromMolBlock(str(x.sdf_info)), x.pk] for x in Molecule.objects.filter(prot_id__target_id=target, cmpd_id=comp_ref)]
        [x[0].SetProp("_Name", "N") for x in mols]
        rdmol.SetProp("_Name", "N")
        # If this is actually a duplicate then continue
        sd_block = Chem.MolToMolBlock(Chem.MolFromMolBlock(Chem.MolToMolBlock(rdmol)))
        matches = [x for x in mols if sd_block == Chem.MolToMolBlock(x[0])]
        if len(matches) > 0:
            # Just return
            return
        else:
            # We have not put this EXACT mol into the database
            # Now lets make an ID
            molid = uuid.uuid4().hex + "_" + pdb_id
            rdmol.SetProp("_Name", molid)
    else:
        molid = pdb_id
    # Make a protein object by which it is related in the DB
    new_mol.prot_id = Protein.objects.get_or_create(code=molid, target_id=target)[0]
    new_mol.sdf_info = Chem.MolToMolBlock(rdmol)
    new_mol.smiles = Chem.MolToSmiles(rdmol, isomericSmiles=True)
    try:
        new_mol.lig_id = rdProps[1]
        new_mol.chain_id = rdProps[2]
        new_mol.occupancy = float(rdProps[3])
    except IndexError:
        new_mol.lig_id = "UNL"
        new_mol.chain_id = "Z"
        new_mol.occupancy = 0.0
    # Add this to the compound list -> make sure this passes in for the
    # correct molecule. I.e. if it fails where does it go???
    # Now link that compound back
    new_mol.cmpd_id = comp_ref
    try:
        new_mol.validate_unique()
        new_mol.save()
    except ValidationError:
        pass


def initialise_dummys():
    """Function to initialise all the dummy database objects
    Takes no args
    Returns None"""
    # Make the dummy objects required for other objects
    d_targ = Target.objects.get_or_create(title="DUMMY", uniprot_id="DUMMY")[0]
    d_prot = Protein.objects.get_or_create(code="DUMMY", target_id=d_targ)[0]
    d_cmpd = Compound.objects.get_or_create(smiles="DUMMY", inchi="DUMMY",
                                            mol_log_p=0.0,
                                            mol_wt=0.0, heavy_atom_count=0,
                                            heavy_atom_mol_wt=0, nhoh_count=0,
                                            no_count=0, num_h_acceptors=0,
                                            num_h_donors=0, num_het_atoms=0,
                                            num_rot_bonds=0, num_val_electrons=0,
                                            ring_count=0, tpsa=0.0)[0]
    Molecule.objects.get_or_create(prot_id=d_prot, cmpd_id=d_cmpd, rscc=0.0,
                                   lig_id="DUMMY", chain_id="X", occupancy=0.0,
                                   x_com=0.0, y_com=0.0, z_com=0.0, sdf_info="DUMMY")
    #Group.objects.get_or_create(group_number=1, group_text="DUMMY", method_used="DUMMY")
    ActivityPoint.objects.get_or_create(confidence=0, cmpd_id=d_cmpd, target_id=d_targ,
                                        activity=0.0, units="DUMMY", internal_id="DUMMY",
                                        source="DUMMY")


def do_oommppaa_proc(target_id=1):
    """Function to generate the full MMP database and find differences.
    Make the mpps, index h change, find all matched pairs
    make 3D coords and make and find differences
    Takes a target_id
    Returns None"""
    initialise_dummys()
    # Make them for activitv points
    make_mmp_database(target_id, option="TWOD")
    # Make them for 3D molecules
    make_mmp_database(target_id, option=None)
    # Do index H change
    index_hydrogen_change(target_id)
    # Find the mmps
    out_mmps = make_list_mmps(target_id, option="ACT", max_size=10,
                              ratio=3, use_ratio=False)
    # Now generate the  three-d mmps from this
    act_mmp_3d(out_mmps, target_id)
    # Find differences between molecules
    make_ph4_difference_points(target_id, opt_type=__POSS_LIST__)
    # Record them to the database
    find_pharma_changes(target_id, my_type=str(__POSS_LIST__))


def find_3d_mmps(target_id=1):
    """Function to fragment all Molecules for a target
    Takes a target id
    Returns None"""
    make_mmp_database(target_id, option=None)


def loading_index_h_change(target_id=1):
    """Function to index hydrogen change for a target
    Takes a target id
    Returns None"""
    index_hydrogen_change(target_id)


def delete_target(target_id):
    """Function to delete all the molecules, 
    proteins and MMP data for a given target
    Takes a target id
    Returns None
    """
    # We need to delete them on by one for SQLite
    # Get the mols, prots and mmpfragss
    mols = Molecule.objects.filter(prot_id__target_id=target_id)
    prots = Protein.objects.filter(target_id=target_id)
    mmpfrags = MMPFrag.objects.filter(mol_id__prot_id__target_id=target_id)
    MMPDiffMap.objects.filter(target_id=target_id).delete()
    acts = ActivityPoint.objects.filter(target_id=target_id)
    ap = ActPharmaPoint.objects.filter(target_id=target_id)
    mmpcomps = MMPComparison.objects.filter(target_id=target_id)
    # Now delete all of them
    for a in acts:
        a.delete()
    for a in ap:
        a.delete()
    for m in mols:
        m.delete()
    for p in prots:
        p.delete()
    for m in mmpfrags:
        m.delete()
    for m in mmpcomps:
        m.delete()
    ActPharmaPoint.objects.filter(target_id=target_id).delete()
    # Now clean up by deleting the target
    Target.objects.get(pk=target_id).delete()


def list_targets():
    """Function to list all the targets
    Takes no args
    Returns None"""
    targs = Target.objects.all()
    for t in targs:
        print t.title


def make_3d_confs(target_id=1):
    """Function to make 3D coords for a target
    Takes a target_id
    Returns None"""
    out_mmps = make_list_mmps(target_id, option="ACT", max_size=10, ratio=3, use_ratio=False)
    # Now generate the  three-d mmps from this
    act_mmp_3d(out_mmps, target_id)


def add_new_act(comp_ref, target, activity, units, chembl_id, source, operator):
    """Function to register a new activity point
    Takes a Compound object a Target object, a measured activity, measured units,
    an id for the compound and a source, e.g. IC50
    Returns None"""
    new_act = ActivityPoint()
    try:
        new_act.activity = float(activity)
    except ValueError:
        print "ACTIVITY DATA NOT A NUMBER -> SKIPPING"
        return
    new_act.units = units
    new_act.source = source
    new_act.cmpd_id = comp_ref
    new_act.target_id = target
    new_act.internal_id = chembl_id
    new_act.operator = operator
    try:
        new_act.validate_unique()
        new_act.save()
    except ValidationError:
        pass


def load_mols(target, file_path):
    """Function to load in the 3D molecules
    Takes a Target object and a file path
    Returns None"""
    mols = Chem.SDMolSupplier(file_path)
    if not mols:
        return
    tot = len(mols)
    if tot == 0:
        print "No molecules given"
        return
    old = -1
    print "Adding molecules..."
    for i, m in enumerate(mols):
        # Catch none molecules
        if m is None:
            print "None molecule", sys.exit()
        # Print the progress
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
        Chem.AssignAtomChiralTagsFromStructure(m)
        add_new_mol(m, target)
    old = 100
    sys.stdout.write("\r%d%%" % old)
    sys.stdout.flush()
    print "\nAdding molecules complete"


def load_activity_data(target, file_path):
    """Function to load in a CSV file of activity data
    Takes a Target object and a file path
    Returns None"""
    # Read the file into a CSV dict
    in_d = read_CSV(open(file_path))
    # Fields looking for
    all_fields = ["smiles", "Activity", "ID", "operator"]
    mandatory_fields = ["smiles", "Activity"]
    # Check to see if fields are missing
    missing_fields = [x for x in mandatory_fields if x not in in_d.fieldnames]
    if len(missing_fields) != 0:
        print " ".join(missing_fields), " fields required"
        sys.exit()
    if len([x for x in all_fields if x not in in_d.fieldnames]) != 0:
        print " ".join([x for x in all_fields if x not in in_d.fieldnames]), " fields missing"
    tot = len(open(file_path).readlines()) - 1
    if tot == 0:
        print "No activity data"
        return
    old = -1
    print "Loading activity data"
    for i, l in enumerate(in_d):
        # Do the percent clock
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
        m = Chem.MolFromSmiles(str(l["smiles"]))
        if m is None:
             # Try doing this - in case it needs escaping
             m = Chem.MolFromSmiles(str(l["smiles"])).decode('string-escape')
        if m is None:
            print "Error None molecule", l["smiles"]
            continue
        comp_ref = add_new_comp(m)
        if comp_ref is None:
            continue
        # Now add the required information if no column is entered
        units = l.get("units")
        if units is None:
            units = "pnM"
        chid = l.get("ID")
        if chid is None:
            chid = "NONE"
        source = l.get("Source")
        if source is None:
            source = "IC50"
        operator = l.get("operator")
        if operator is None:
            operator = "NA"
        add_new_act(comp_ref, target, l["Activity"], units, chid, source, operator)
    old = 100
    sys.stdout.write("\r%d%%" % old)
    sys.stdout.flush()
    print "\nAdding activity data complete"
    return None


def load_protein(target, file_path):
    """Function to load in a pdb file
    Takes a Target object and a file path
    Returns None"""
    new_prot = Protein()
    new_prot.code = target.title + "TEMP"
    new_prot.pdb_info = open(file_path).read()
    new_prot.target_id = target
    try:
        new_prot.validate_unique()
        new_prot.save()
    except ValidationError:
        my_prot = Protein.objects.get(code=target.title + "TEMP")
        print "TEMPLATE PROTEIN ALREADY ENTERED"
        print "OVERWRITING PROTEIN"
        my_prot.pdb_info = open(file_path).read()
        my_prot.save()


def load_protein_file(target, file_path, code):
    """Function to load in a pdb file
    Takes a Target object and a file path
    Returns None"""
    new_prot = Protein()
    new_prot.code = code
    new_prot.pdb_info = open(file_path).read()
    new_prot.target_id = target
    try:
        new_prot.validate_unique()
        new_prot.save()
    except ValidationError:
        my_prot = Protein.objects.get(code=code)
        my_prot.pdb_info = open(file_path).read()
        my_prot.save()


def load_proteins(target, dir_path):
    """Function to load all the APO structures that will then relate
    back to the bound ligands.
    Takes a target and a directory path
    Returns None"""
    curr_dir = os.getcwd()
    os.chdir(dir_path)
    # Get all the files
    for my_file in glob.glob("*apo.pdb"):
        my_code = os.path.split(my_file)[1][:4] + "_" + str(target.pk)
        load_protein_file(target, my_file, my_code)
    os.chdir(curr_dir)


def load_res_shifts(target_id, dir_path):
    """Function to register the biggest residue shifts and most shifting residues"""
    in_f = os.path.join(dir_path, "res_split.txt")
    in_d = ast.literal_eval(open(in_f).read())
    # First get and store the proteins
    #for prot in in_d["PROT"]:
    #    open(in_d["PROT"][prot])
    # Now loop through these
    for res in in_d:
        if res == "PROT":
            continue
        if len(in_d[res]) != len(in_d["PROT"]):
            print "RESIDUE "+res+" NOT FULLY ANALYSED - CONTINUE"
            continue
        rp = ResProt()
        rp.target_id = target_id
        rp.res_name = res
        try:
            rp.validate_unique()
        except ValidationError:
            rp = ResProt.objects.get(target_id=target_id, res_name=res)
        rp.avg_shift = numpy.mean(in_d[res])
        rp.max_shift = max(in_d[res])
        rp.min_shift = min(in_d[res])
        rp.save()
        for i, prot in enumerate(in_d["PROT"]):
            rs = ResShift()
            rs.target_id = target_id
            rs.res_id = rp
            # Check 
            try:
                my_prot = Protein.objects.get_or_create(code=prot + "_"  + str(target_id.pk))[0]
            except:
                my_prot = Protein.objects.filter(code__startswith=prot, target_id=target_id)
                if len(my_prot) == 0:
                    print prot
                    print target_id.title
                    continue
                else:
                    my_prot = my_prot[0]
            my_mol = Molecule.objects.filter(prot_id=my_prot)
            rs.prot_id = my_prot
            rs.shift = in_d[res][i]
            try:
                rs.validate_unique()
            except ValidationError:
                rs = ResShift.objects.get(target_id=target_id, prot_id=my_prot, res_id=rp)
                rs.shift = in_d[res][i]
          # If it's APO then we have a none pk
            if len(my_mol) == 0:
                pass
            else:
                rs.mol_id = my_mol[0]
            #Now save it
            rs.save()


def load_waters(target, dir_path):
    """Function to load all the APO structures that will then relate
    back to the bound ligands.
    Takes a target and a directory path
    Returns None"""
    curr_dir = os.getcwd()
    os.chdir(dir_path)
    # Get all the files
    from WONKA.models import Water
    for my_file in glob.glob("*waters.mol2"):
        my_code = os.path.split(my_file)[1][:4] + "_" + str(target.pk)
        my_prot = Protein.objects.filter(code=my_code)
        if not my_prot:
            continue
        my_prot = my_prot[0]
        waters = Chem.MolFromMol2File(my_file)
        conf = waters.GetConformer()
        # Delete them for this protein
        for w in Water.objects.filter(prot_id=my_prot):
            w.delete()
        for i in range(waters.GetNumAtoms()):
            cp = conf.GetAtomPosition(i)
            Water.objects.get_or_create(water_num=i+1, prot_id=my_prot, target_id=target,x_com=cp.x,y_com=cp.y,z_com=cp.z)
    os.chdir(curr_dir)


def load_compounds(file_path):
    """Function to load compounds and make the MMPs
    Takes a file path
    Returns None"""
    mols = Chem.SDMolSupplier(file_path)
    counter = 0
    for m in mols:
        if m is None:
            print "NONE MOL"
            continue
        counter +=1
        print counter
        # add the new compound to the database
        comp_ref = add_new_comp(m)
        if comp_ref is None:
            continue
        new_m = Chem.MolFromSmiles(str(comp_ref.smiles))
        # Filter too big molecules
        if Descriptors.ExactMolWt(new_m) > 560:
            continue
        make_mol_mmp(new_m, id="cmp" + str(comp_ref.pk), target_id=None)


def refresh_mmp_maps(target_id):
    """Function to refresh maps for a target
    Takes a target id
    Returns None"""
    MMPDiffMap.objects.filter(target_id=target_id).delete()
    try:
        MMPComparison.objects.filter(target_id=target_id).delete()
    except:
        # Fix for SQLite
        mols = MMPComparison.objects.filter(target_id=target_id)
        for m in mols:
            m.delete()
    try:
        ActPharmaPoint.objects.filter(target_id=target_id).delete()
    except:
        # Fix for SQLite
        mols = ActPharmaPoint.objects.filter(target_id=target_id)
        for m in mols:
            m.delete()
    print "RECHARGING"
    make_ph4_difference_points(target_id, opt_type=__POSS_LIST__)
    find_pharma_changes(target_id, my_type=str(__POSS_LIST__))


def do_lloommppaa_proc(target_id, pdb_protein, smiles, mol2_protein=None, reactants=None, products=None, context=None):
    """Function to do the processing for LLOOMMPPAA"""
    from LLOOMMPPAA.models import PossReact
    from LLOOMMPPAA.reactions import run_react_proc, define_reacts, load_in_reagents, load_in_follow_ups, find_follow_ups, define_reaction_process
    # Load the data
    load_wonka_prot(target_id, pdb_protein, smiles)
    if reactants:
        # Find the potential synthesis vectors
        define_reacts()
        find_follow_ups(target_id=target_id)
        # Show the sides -> USER MUST SELECT A SIDE
        poss_reacts = PossReact.objects.filter(mol_id__prot_id__target_id=target_id)
        # Loop through all the possible reactions
        for ps_r in poss_reacts:
            print ps_r.react_id.name
            print "1):", ps_r.replaced_frag
            print "2):", ps_r.retained_frag
            print "3):", "SKIP"
            choice = int(raw_input("Select a fragment to replace...")) 
            if choice == 1:
                context = ps_r.retained_frag_context
                react_frag = ps_r.retained_frag
                this_react = ps_r.react_id
                break
            elif choice == 2:
                context = ps_r.replaced_frag_context
                react_frag = ps_r.replaced_frag
                this_react = ps_r.react_id
                break
            else:
                continue
    if not context:
        print "YOU MUST SPECIFY A CONTEXT"
        if products:
            print "AUTO GENERATING FROM PRODUCTS"
            print products
            context = Chem.CanonSmiles(Chem.MolToSmiles(Chem.MolFromSmarts(MCS.FindMCS(Chem.SDMolSupplier(products)).smarts),isomericSmiles=True))
            print context
        else:
            return
    if reactants:
        mol_id = ps_r.mol_id
        prot_id = ps_r.mol_id.prot_id
    else:
        mol_id = Molecule.objects.get(prot_id__target_id=target_id,smiles=smiles)
        prot_id = mol_id.prot_id
    my_prots = [prot_id] 
    # Set the mol2 protein for this target -> throw a warning if this doesn't happen
    if mol2_protein:
        from PLIFS.models import PlifProtein
        pp = PlifProtein()
        pp.prot_id = prot_id
        pp.mol2_data = open(mol2_protein).read()
        pp.save()
    # Define the reactants and products
    if reactants:
        react_id = load_in_reagents("RUN_DEF", reactants, ps_r.react_id)
    else:
        react_id = None
    if products:
        my_react = Reaction.objects.get_or_create(name="DUMMY",react_smarts="DUMMY",retro_smarts="DUMMY", cont_smarts="DUMMY")[0]
        this_react = my_react
        react_frag = context
        prod_id = load_in_follow_ups("RUN_PROD", products,my_react, mol_id)
    else:
        prod_id = None
    # Now set up the reaction itself
    react_proc = define_reaction_process(mol_id, context, my_prots, this_react, context, react_frag, products_id=prod_id, reactants_id=react_id)
    if react_id:
        react_proc.proc_stage = "RUN_REACTION"
        react_proc.save()
    if products:
        react_proc.proc_stage = "GENERATE CONFS"
        react_proc.save()
    # Now run the reaction and analysis itself
    run_react_proc(react_proc)


def load_wonka_prot(target_id, protein, smiles):
    """Function to load in data a la WONKA"""
    from WONKA.functions import find_target
    out_d = {"MODEL_ONE": {}}
    out_d["MODEL_ONE"]["path_to_pdb"] = protein
    out_d["MODEL_ONE"]["smiles"] = str(smiles).decode('string-escape')
    out_d["MODEL_ONE"]["cmpd_id"] = "LIG"
    out_d["MODEL_ONE"]["path_to_map"] = None
    out_d["MODEL_ONE"]["model_id"] = "MODEL_ONE"
    out_d["MODEL_ONE"]["chain"] = "A"
    title = Target.objects.get(pk=target_id).title
    find_target(title, overwrite=None, save_map=None, data_site="SIMPLE", csv_path=out_d)


def do_wonka_proc(title, overwrite=None, save_map=None, data_site="SGC", csv_path=None):
    """Function to do the processing for LLOOMMPPAA"""
    from WONKA.functions import find_target
#    # Make the dummys
    initialise_dummys()
    find_target(title, overwrite=None, save_map=None, data_site=data_site, csv_path=csv_path)


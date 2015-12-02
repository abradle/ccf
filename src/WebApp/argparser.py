import sys
import getopt
from loading import *
from IOhandle.models import Target


def myargparse():
    """Function to parse arguments and handle functions
    Takes no args
    Returns None"""
    # Define the paths as none
    my_acts = None
    my_mols = None
    my_prot = None
    my_targ = None
    my_cmps = None
    # Define the flags as none
    make_mmp_flag = None
    lloommppaa_flag = None
   # Defien the new stuff for LLOOMMPPAAA 
    mol2_protein=None 
    reactants=None
    products=None
    context=None
# 
    refresh_map_flag = None
    index_h_flag = None
    find_3d_flag = None
    make_3d_flag = None
    del_targ_flag = None
    list_targ_flag = None
    no_targ = None
    wonka_flag = None
    csv_path = None
    data_site = "SGC"
    preamble = """
 _______  _______  _______  _______  _______  _______  _______  _______
(  ___  )(  ___  )(       )(       )(  ____ )(  ____ )(  ___  )(  ___  )
| (   ) || (   ) || () () || () () || (    )|| (    )|| (   ) || (   ) |
| |   | || |   | || || || || || || || (____)|| (____)|| (___) || (___) |
| |   | || |   | || |(_)| || |(_)| ||  _____)|  _____)|  ___  ||  ___  |
| |   | || |   | || |   | || |   | || (      | (      | (   ) || (   ) |
| (___) || (___) || )   ( || )   ( || )      | )      | )   ( || )   ( |
(_______)(_______)|/     \||/     \||/       |/       |/     \||/     \|
"""

    try:
        d = sys._MEIPASS
        print_option = preamble + "\nmanage.exe"
    # If not use the linux babel or the windows babel in defined locations
    except AttributeError:
        print_option = preamble + "\npython manage.py"
    print_out = print_option + """ --targ <TARGET> --mols <MOLS.sdf> --acts <ACTS.csv> --cmps <COMPOUNDS.sdf> --prot <PDBFILE.pdb> --mmp True

    Core options:
    --targ <TARGET> Name of target which information refers to
    --mols <MOLS.sdf> File path to molecules to add for this target
    --acts <ACTS.csv> File path to activity data for this target
    --mmp <True/False> Flag to run MMP analysis on this data
    --lloommppaa <True/False> Flag to run LLOOMMPPAA analysis for a target
    --wonka <True/False> Flag to run WONKA analysis for a target
    --prot <PDBFILE.pdb> File path for template protein structure for this data

    Refresh options:
    --refreshmaps <True/False> Flag to refresh maps after option changes

    Other options:
    --deleteTarget <True/False> Flag to find delete the target specified
    --listTargets <True/False> Flag to find all the targets and print them out
    --cmps <COMPOUNDS.sdf> SD file for extra compounds to be added (are not attached to a target)
    --pharm <True/False> Flag to generate pharmacophore points for molecules for a target
    --hchange <True/False> Flag to run just the Index H Change option for a target
    --make3d <True/False> Flag to generate 3D conformations for a target
    --find3d <True/False> Flag to find MMPs for all 3D molecules
    --datasite <TYPE> (SGC, GSK, CSV) - the type of data to parse
    --csv_path <PATH/TO/FILE> - the file to pick up
    """
    # List of options
    opt_list = [("target=", "t:"), ("mols=", "m:"), ("acts=", "a:"),
                ("mmp=", "p:"), ("cmps=", "c:"), ("datasite=", "f:"),
                ("refreshmaps=", "e:"), ("hchange=", "q:"), ("make3d=", "w:"),
                ("find3d=", "l:"), ("prot=", "z:"), ("ll_prot=", "j:"),("deleteTarget=", "d:"),("listTargets=", "b:"), ("lloommppaa=", "r:"), ("wonka=", "k:"),("smiles=", "j:"), ("csv_path=", "v"), ("reactants=", "i"), ("products=","h:"), ("mol2_prot=","g:"),("context=","z:")]
    try:
        # Parse the arguments
        opts, args = getopt.getopt(sys.argv[1:],"".join([x[1] for x in opt_list]),[x[0] for x in opt_list])
    except getopt.GetoptError as err:
        print "ERROR IN PARSING"
        print err
        print print_out
        return
    if len(opts) == 0:
        print print_out
        return

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print print_out
            sys.exit()
        elif opt in ("-m", "--mols"):
            my_mols = arg
        elif opt in ("-z", "--prot"):
            my_prot = arg
        elif opt in ("-t", "--target"):
            my_targ = arg.rstrip()
        elif opt in ("-a", "--acts"):
            my_acts = arg
        elif opt in ("-p", "--mmp"):
            make_mmp_flag = arg
        elif opt in ("-c", "--cmps"):
            my_cmps = arg
        elif opt in ("-e", "--refreshmaps"):
            refresh_map_flag = arg
        elif opt in ("-q", "--hchange"):
            index_h_flag = arg
        elif opt in ("-w", "--make3d"):
            make_3d_flag = arg
        elif opt in ("-l", "--find3d"):
            find_3d_flag = arg
        elif opt in ("-d", "--deleteTarget"):
            del_targ_flag = arg
        elif opt in ("-b", "--listTargets"):
            list_targ_flag = arg
        elif opt in ("-r", "--lloommppaa"):
            lloommppaa_flag = arg
        elif opt in ("-k", "--wonka"):
            wonka_flag = arg
        elif opt in ("-f", "--datasite"):
            data_site = arg
        elif opt in ("-v", "--csv_path"):
            csv_path = arg
        elif opt in ("--reactants",):
            reactants = arg
        elif opt in ("--products",):
            products = arg
        elif opt in ("--context",):
            context = arg
        elif opt in ("--mol2_prot",):
            mol2_prot = arg
        elif opt in ("--smiles",):
            smiles = arg
        elif opt in ("--ll_prot",):
            ll_prot = arg
    # Initialise dummy variables
    print preamble
    initialise_dummys()
    if my_cmps or list_targ_flag:
        no_targ = True
    if not my_targ and not no_targ:
        print "You need to specify a target!!!!"
        return
    elif my_cmps and not my_targ:
        load_compounds(my_cmps)
        print "Loaded compounds. No target specified."
        return
    elif list_targ_flag and not my_targ:
        list_targets()
        return
    target = Target.objects.get_or_create(title=my_targ)[0]
    if del_targ_flag:
        print "DELETING TARGET... - ", target.title
        delete_target(target.pk)
        print "DELETED TARGET... - ", target.title
        return
    if refresh_map_flag:
        refresh_mmp_maps(target.pk)
    if index_h_flag:
        loading_index_h_change(target.pk)
    if find_3d_flag:
        find_3d_mmps(target.pk)
    if make_3d_flag:
        make_3d_confs(target.pk)
    if my_mols:
        load_mols(target, my_mols)
    if my_prot:
        load_protein(target, my_prot)
    if my_acts:
        load_activity_data(target, my_acts)
    if wonka_flag:
        do_wonka_proc(my_targ, overwrite=None, save_map=None, data_site=data_site, csv_path=csv_path)
    if make_mmp_flag:
        do_oommppaa_proc(target.pk)
    if my_cmps:
        load_compounds(my_cmps)
    if lloommppaa_flag:
        do_lloommppaa_proc(target.pk, ll_prot, smiles, mol2_prot, reactants, products, context)
    if list_targ_flag:
        list_targets()

# Functions to perform analysis on LLOOMMPPAA datasets
import numpy as np
from LLOOMMPPAA.models import LLConf, GPVal, HistBinVal, HistBin, GridPoint, Histogram, DivList
from IOhandle.models import Target
import pickle, os, csv, math, sys
from collections import Counter
from PLIFS.models import Interaction, PlifProbe
import numpy as np
from PLIFS.models import PlifProbe
from rdkit import Chem
from rdkit.Chem import AllChem
from django.db.models import Sum
import scipy.cluster.hierarchy as sch
# do this before importing pylab or pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def delete_hist():
    """Helper function to delete all the stuff in this area"""
    while HistBin.objects.count():
        ids = HistBin.objects.values_list('pk', flat=True)[:100]
        HistBin.objects.filter(pk__in = ids).delete()
    while GridPoint.objects.count():
        ids = GridPoint.objects.values_list('pk', flat=True)[:100]
        GridPoint.objects.filter(pk__in = ids).delete()
    HistBinVal.objects.all().delete()
    GPVal.objects.all().delete()
    Histogram.objects.all().delete()


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


def fill_in_samp_dict(out_d, samp_list, tot_list, target):
    """Function to fill in the counts in the SAMPLE dic"""
    samp_dict = {}
    for num_samps in samp_list:
        print num_samps
        samp_dict[num_samps] = {}
        for test in out_d:
            # Normalise the values or set zero elements
            for item in tot_list:
                if item in out_d[test][num_samps]:
                    out_d[test][num_samps][item] = float(out_d[test][num_samps][item]) / float(num_samps) #### Muliply by the number of iters here
                else:
                    out_d[test][num_samps][item] = 0
            # Now write these out to the new dict 
            for item in out_d[test][num_samps]:
                if item not in samp_dict[num_samps]:
                    samp_dict[num_samps][item] = {test: out_d[test][num_samps][item]}
                else:
                    samp_dict[num_samps][item][test] = out_d[test][num_samps][item]
    return out_d, samp_dict


def write_total_ans_csv(num_samps, target, samp_dict):
    """Function to write out the whole analysis as a CSV file"""
    # File handle for the total answers
    tot_ans = open(str(num_samps) + "__" + target.title + "__" + "ans.csv", "w")
    tot_ans.write("COMPS, ")
    for my_cmp in samp_dict[num_samps]:
        test_list = [test for test in samp_dict[num_samps][my_cmp] if test[:8] != "PLIFSDIV"]
        tot_ans.write(",".join([test for test in test_list]))
        tot_ans.write("\n")
        break
    # Now create the comparison vals for this test
    mat_vals = {}
    chi_vals = {}
    for i in range(len(test_list)):
        for j in range(i, len(test_list)):
            if i == j:
                continue
            mat_vals[test_list[i] + "__" + test_list[j]] = 0
    for i in range(len(test_list)):
        for j in range(i, len(test_list)):
            if i == j:
                continue
            chi_vals[test_list[i] + "__" + test_list[j]] = 0
    for my_cmp in samp_dict[num_samps]:
        tot_ans.write(str(my_cmp) + ",")
        me = ",".join([str(samp_dict[num_samps][my_cmp][test]) for test in test_list])
        tot_ans.write(me)
        tot_ans.write("\n")
        for item in mat_vals:
            test_1 = item.split("__")[0]
            test_2 = item.split("__")[1]
            diff = math.fabs(samp_dict[num_samps][my_cmp][test_1] - samp_dict[num_samps][my_cmp][test_2])
            my_sum = samp_dict[num_samps][my_cmp][test_1] + samp_dict[num_samps][my_cmp][test_2]
            mat_vals[item] += diff
            # Find the chi squared diff
            if my_sum == 0:
                chi_vals[item] += 0
            else:
                chi_vals[item] += diff / my_sum
    tot_ans.close()
    return mat_vals, chi_vals


def write_comp_mats(num_samps, target, comp_mats):
    """Function to write out the comparison matrices"""
    comp_ans = open(str(num_samps) + "__" + target.title + "__" + "COMPans.csv", "w")
    for mat_vals in comp_mats:
        for item in mat_vals:
            comp_ans.write(item)
            comp_ans.write(",")
            comp_ans.write(str(mat_vals[item]))
            comp_ans.write("\n")
        comp_ans.write("\n\n\n")


def write_out_ans(samp_dict, target):
    """Function to write out the full analysis for this set"""
    for num_samps in samp_dict:
        mat_vals, chi_vals = write_total_ans_csv(num_samps, target, samp_dict)
        write_comp_mats(num_samps, target, [mat_vals, chi_vals])
    return None


def init_data(target_id, file_path):
    """Function to init the data at the start into useable dicts"""
    # Get the information
    target, my_d = parse_pickle(target_id, file_path)
    # Generate histograms - compound = bin, for each method / number of selections
    out_d, samp_list, tot_list = get_counts_for_each(my_d)
    # Now fill in all the zero elements
    out_d, samp_dict = fill_in_samp_dict(out_d, samp_list, tot_list, target)
    return out_d, samp_dict, target, my_d, tot_list


def parse_pickle(target_id, file_path):
    """Function to parse the out answer"""
      # Get the target=
    if target_id:
        target = Target.objects.get(pk=target_id)
    else:
        target = None
#    if not file_path:
#        file_path = os.path.join(r"W:\Informatics\Pharmacophore\anthony\DPhil\CODE\CHOC\src\WebApp", target.title + "20___" + str(rp_id) + "_out.ans")
    # Open the file
    my_d = pickle.load(open(file_path))
    return target, my_d


def render_pdb_data(target_id, test, num_samps, this_item,
                    this_samp, file_name, react_anal_id,
                    method_id, type_id, opt=None):
    """Function to take a dict of PlifProb PKs with corresponding values
    Use the value to fill the BFactor column"""
    # Define this histogram
    new_hist = Histogram.objects.get_or_create(test_made=test, num_samps=num_samps, hist_title= this_item)[0]
    mol = Chem.MolFromSmiles("C")
    mol_txt = ""
    for item in this_samp:
        # First get the PlifProbe if we need
        if opt:
            pp = PlifProbe()
            pp.x_com = item[0]
            pp.y_com = item[1]
            pp.z_com = item[2]
            pp.intensity = item[3]
        else:
            pp = PlifProbe.objects.get(pk=item)
        # Make an editable molecule
        em = AllChem.EditableMol(mol)
        em.RemoveAtom(0)
        sv = em.AddAtom(Chem.Atom(11))
        gm = em.GetMol()
        Chem.SanitizeMol(gm)
        AllChem.EmbedMolecule(gm)
        cnf = gm.GetConformer()
        sp = AllChem.rdGeometry.Point3D()
        # Now add the coords
        sp.x = pp.x_com
        x_char = str(pp.x_com)[:9]
        sp.y = pp.y_com
        y_char = str(pp.y_com)[:9]
        sp.z = pp.z_com
        z_char = str(pp.z_com)[:9]
        cnf.SetAtomPosition(0, sp)
        out_mol = cnf.GetOwningMol()
        out_mol = Chem.MolFromPDBBlock(Chem.MolToPDBBlock(out_mol))
        atm = out_mol.GetAtomWithIdx(0)
        if opt:
            atm.GetPDBResidueInfo().SetTempFactor(pp.intensity)
            val_out = pp.intensity
        else:
            atm.GetPDBResidueInfo().SetTempFactor(this_samp[item])
            val_out = this_samp[item]
        # So now create this data
        my_gp = GridPoint.objects.get_or_create(x_char=x_char, y_char=y_char, z_char=z_char, target_id=target_id)[0]
        my_grid_val = GPVal.objects.get_or_create(my_anal_id=react_anal_id, num_samps=num_samps,
                                                  type_id=type_id, method_id=method_id,
                                                  gp_id=my_gp, target_id=target_id)[0]
        # Now add this value
        my_grid_val.value = val_out
        my_grid_val.pdb_info = Chem.MolToPDBBlock(out_mol)
        my_grid_val.sdf_info = Chem.MolToMolBlock(out_mol)
        my_grid_val.save()
        mol_txt += Chem.MolToPDBBlock(out_mol)+"\n"
    out_f = open(file_name, "w")
    out_f.write(mol_txt)
    out_f.close()


def get_grid_map(target_id, cmp_list):
    """Function to get the grid map for a given target"""
    # First find the x_max
    my_vals = PlifProbe.objects.filter(target_id=target_id, mol_id__cmpd_id__in=cmp_list).values_list("x_com", "y_com", "z_com")
    x_coords = [x[0] for x in my_vals]
    y_coords = [x[1] for x in my_vals]
    z_coords = [x[2] for x in my_vals]
    # Get the min and max coords
    return [max(x_coords), min(x_coords), max(y_coords), min(y_coords), max(z_coords), min(z_coords)]


def get_samp_bins(samp_ints, grid_map, target_id, cmp_list, spacing=1.5, boundary=6.0):
    """Function to find the interactions for a given sample - given this target's
    grid map"""
    # Get all of them for this set
    in_pps = PlifProbe.objects.filter(target_id=target_id, mol_id__cmpd_id__in=cmp_list).filter(pk__in=samp_ints)
    # Make sure we've got the list - then normalise the contributions by this value too - so each set adds up to the same value
    # Normalise by this answer
    len(in_pps)
    # Work out the grid spacings
    x_range = np.arange(grid_map[1] - boundary,
                        grid_map[0] + boundary,
                        spacing)
    y_range = np.arange(grid_map[3] - boundary,
                        grid_map[2] + boundary,
                        spacing)
    z_range = np.arange(grid_map[5] - boundary,
                        grid_map[4] + boundary,
                        spacing)
    # List to hold the points
    out_dict = {}
    old_x = None
    old_y = None
    old_z = None
    # Work out the different types of PLIF types
    type_list = list(set(in_pps.values_list("type", flat=True)))
    type_list.append("ALL")
    old = 0
    # Do a % counter
    step = float(float(1) / float(len(x_range))) * 100
    # Set up empty lists
    my_volume = (spacing * spacing * spacing)
    tot_prop = 0.0
    for my_ph4_type in type_list:
        out_dict[my_ph4_type] = []
    for x in x_range:
        sys.stdout.write("\r%d%% complete..." % old)
        sys.stdout.flush()
        old += step
        # Now calculate the bonds of this square
        x_max = x + spacing
        x_min = x - spacing
        for y in y_range:
            y_max = y + spacing
            y_min = y - spacing
            for z in z_range:
                z_max = z + spacing
                z_min = z - spacing
                # Set the coords
#                pp = PlifProbe()
#                pp.x_com = x
#                pp.y_com = y
#                pp.z_com = z
                these_ps = in_pps.filter(x_com__lte=x_max).filter(y_com__lte=y_max).filter(z_com__lte=z_max)
                these_ps = these_ps.filter(x_com__gt=x_min).filter(y_com__gt=y_min).filter(z_com__gt=z_min)
                # Now loop through the list
                for my_ph4_type in type_list:
                    # Now we've got these points - multiply them by the number of cmps in this set
                    if my_ph4_type == "ALL":
                        pass
                    else:
                        # Now filter them on pharma type
                        my_ps = these_ps.filter(type=my_ph4_type)
                    intensity = 0.0
                    for this_p in my_ps:
                        # First work out how much this one goes into the box
                        # Get the absolute distance
                        x_dist = math.fabs(x - this_p.x_com)
                        y_dist = math.fabs(y - this_p.y_com)
                        z_dist = math.fabs(z - this_p.z_com)
                        if x_dist > spacing or y_dist > spacing or z_dist > spacing:
                            print "DISTANCE ERROR"
                        # Then interpolate this point - normalised by the contribution
                        this_prop = ((spacing - x_dist) * (spacing - y_dist) * (spacing - z_dist)) / my_volume
                        # Now multiply this by the number of compounds with this point
                        intensity += samp_ints[this_p.pk] * this_prop
                    out_dict[my_ph4_type].append((x, y, z, intensity,))
        # Append this list appropriately
    # Now print this out
    sys.stdout.write("\r%d%% complete..." % old)
    sys.stdout.flush()
    print out_dict.keys()
    print "MADE COMP MAP"
    return out_dict


def get_data_from_db(target_id, react_anal_id):
    """Function to take in the dictionary and return counts for each compound set"""
    ele_set = set()
    samp_set = set()
    out_d = {}
    db_lists = DivList.objects.filter(run_id=react_anal_id)
    my_tests = list(set(db_lists.values_list("method_id", flat=True)))
    my_samps = list(set(db_lists.values_list("num_samps", flat=True)))
    for test in my_tests:
        out_d[test] = {}
        for num_samps in my_samps:
            # Now count the number of each ele
            samp_set.add(num_samps)
            out_d[test][num_samps] = {}
            # Tot_ans is a dict key = comp, val = count
            #tot_ans = dict(Counter(tot_l))
            tot_ans = db_lists.get(method_id=test, num_samps=num_samps).score_id.all()
            # Now make the dict for this
            for ans in tot_ans:
                out_d[test][num_samps][ans.cmpd_id.pk] = ans.score
            # Go through and add these to the total set
            for my_ans in tot_ans:
                ele_set.add(my_ans)
    #  Now return this dict
    return out_d, list(ele_set)


def find_populated_interactions(target_id, react_anal_id):
    """For each method / number of picks - find the position and number
    of interactions.
    Writes out a PDB file - BFactor = scale of interactions
    Returns None"""
    from LLOOMMPPAA.functions import calc_cmpd_complex
    # Get all the data in dicts. Out_d = test: num_samps: compd: NUMBER
    # orig_d = test: num_samps: samp_num: CMPD_LIST
    # samp_dict = num_samps: compd: test: NUMBER
    #from IOhandle.models import Compound
    out_d, tot_cmp_list = get_data_from_db(target_id, react_anal_id)
    # Now write out this comparison
    #write_out_ans(samp_dict, target)
    # Get the grid map for this target - based on interaction points
#    grid_map = get_grid_map(target_id, tot_cmp_list)
    grid_map = get_grid_map(target_id, [react_anal_id.react_proc_id.mol_id.cmpd_id.pk])
    base_path = ""
    comp_list = []
    for test in out_d:
        for num_samps in out_d[test]:
            this_samp_ints = {}
            this_samp_plifs = {}
            # Find the compounds for this set
            for comp in out_d[test][num_samps]:
                comp_list.append(comp)
                # Get the normalised number of times this compound is selected
                mult = out_d[test][num_samps][comp]
                # If it wasn't selected
                if not mult:
                    continue
                # Get the interactions associated to this compound
                my_ints = Interaction.objects.filter(mol_id__cmpd_id__pk=comp, target_id=target_id, scheme_id__scheme_name="broad").values_list("interaction_id__probe_dest_id__pk", flat=True)
                # Get its PLIF probes associated to this compounds
                my_plifs = PlifProbe.objects.filter(target_id=target_id, mol_id__cmpd_id=comp)
                # Multiply them by the number of compounds in this group - weighted by the amount it is in this grouping
                # Normalise so the contribution of each compound is mult - so should sum to 100
                for my_int in my_ints:
                    if my_int in this_samp_ints:
                        this_samp_ints[my_int] += float(mult)
                    else:
                        this_samp_ints[my_int] = float(mult)
                # Make a dict for the PLIF probes
                for my_int in my_plifs:
                    if my_int.pk in this_samp_plifs:
                        this_samp_plifs[my_int.pk] += float(mult)
                    else:
                        this_samp_plifs[my_int.pk] = float(mult)
            # Now normalise this_samp_ints
#            for my_int in this_samp_ints:
#                this_samp_ints[my_int] = this_samp_ints[my_int]
#            for my_plif in this_samp_plifs:
#                this_samp_plifs[my_plif] = this_samp_plifs[my_plif]
            # Now make the PDB data
            # Bin samps
#            file_path = os.path.join(base_path, test + str(num_samps) + ".INTS.pdb")
#            render_pdb_data(Target.objects.get(pk=target_id), test, num_samps,
#                            "INTS", this_samp_ints, file_path, react_anal_id,
#                            test, "ALL")
            # Now work out the compound places
            this_samp_bins = get_samp_bins(this_samp_plifs, grid_map, target_id, list(set(comp_list)))
            ### We want to do it for ALL and then per pharmacophore point ###
            for ph4_type in this_samp_bins:
                render_pdb_data(Target.objects.get(pk=target_id), test, num_samps,
                            "INTS", this_samp_bins[ph4_type], "me.pdb", react_anal_id,
                            test, ph4_type, opt="HELLO")
            # Send these to be saved in the database
            calc_cmpd_complex()
            find_complexity(Target.objects.get(pk=target_id), out_d[test][num_samps], num_samps, test, react_anal_id)


def find_my_hists(target_id, react_anal_id):
    """"""
    print "GETTING DATA"
    from LLOOMMPPAA.functions import calc_cmpd_complex
    calc_cmpd_complex()
    out_d, tot_cmp_list = get_data_from_db(target_id, react_anal_id)
    for test in out_d:
        print "DOING ANALYSIS FOR ",test
        for num_samps in out_d[test]:
            find_complexity(Target.objects.get(pk=target_id), out_d[test][num_samps], num_samps, test, react_anal_id)


def get_dend_data(react_anal_id):
    """Function to get the data in the format for a dendrogram"""
    ele_set = set()
    samp_set = set()
    out_d = {}
    db_lists = DivList.objects.filter(run_id=react_anal_id).filter(method_id__in=["MACCS", "USRCAT", "MORGAN", "RANDOM", "PLIFS_default", "PLIFS_broad"])
    my_tests = list(set(db_lists.values_list("method_id", flat=True)))
    my_samps = list(set(db_lists.values_list("num_samps", flat=True)))
    if len(my_samps) > 1:
        print "TOO MANY SAMPLES"
        return None
    for test in my_tests:
        # Tot_ans is a dict key = comp, val = count
        out_d[test] = {}
        #tot_ans = dict(Counter(tot_l))
        tot_ans = db_lists.filter(method_id=test, num_samps=my_samps[0])[0].score_id.all()
        # Now make the dict for this
        for ans in tot_ans:
            out_d[test][ans.cmpd_id.pk] = ans.score
            ele_set.add(ans.cmpd_id.pk)
        # Go through and add these to the total set
    #  Now return this dict
    # Now make a list of tests
    out_ans = []
    for i, test in enumerate(my_tests):
        out_ans.append([])
        print out_d[test]
        for cmpd in ele_set:
            if cmpd in out_d[test]:
                out_ans[i].append(out_d[test][cmpd])
            else:
                out_ans[i].append(0.0)
    return out_ans, my_tests


def make_dendrogram(react_anal_id):
    """Function to make a dendrogram from the data available"""
    ### HERE GENERATE THE APPROPRIATE DATA
    import StringIO
    # Generate histograms - compound = bin, for each method / number of selections
    # Get the data
    output = StringIO.StringIO()
    fig = plt.figure()
    X, labels = get_dend_data(react_anal_id)
    d = sch.distance.pdist(X)
    Z = sch.linkage(d, method='complete')
    P = sch.dendrogram(Z, labels=labels, orientation='right')
    fig.savefig(output, format="PNG")
    contents = output.getvalue()
    return contents


def find_complexity(target_id, comp_d, num_samps, test, react_anal_id):
    """"For each method / number of picks - find the complexity distributions for each metric
    Writes out a CSV file for each one, with predefined binned histograms
    Returns None"""
    from LLOOMMPPAA.models import Cmpd_Complexity

    from django.db.models import Max, Min
    # Loop through the comps
    # Get all the compounds
    all_cmps = Cmpd_Complexity.objects.all()
    all_bins = {}
    my_tests = ['num_hba', 'num_hbd', 'num_rot_bonds','num_rings','num_chir','num_maccs','num_morg']
    for item in my_tests:
        all_bins[item] = {}
        for x in range(all_cmps.aggregate(Max(item))[item + "__max"] + 1):
            all_bins[item][x] = 0.0
    # Now for the continous tests
    cont_tests = ["sa_score", "size"]
    for item in cont_tests:
        all_bins[item] = {}
        agg = all_cmps.aggregate(Min(item), Max(item))
        my_bins = np.linspace(agg[item + "__min"], agg[item + "__max"], 30)
        for i, x in enumerate(my_bins):
            if i == len(my_bins) -1:
                continue
            all_bins[item][i] = {"min": x, "max": my_bins[i+1], "val": 0.0}
    for comp in comp_d:
        dj_comp = Cmpd_Complexity.objects.get(cmpd_id__pk=comp)
        # Get the data for this comp - extra 0.5 to compensate for slight FP calc errors
        dj_comp_num = comp_d[comp]
        if dj_comp_num == 0.0:
            continue
        ### Histogram for each property for this experiment
        for item in my_tests:
            num_val = getattr(dj_comp, item)
            all_bins[item][num_val] += dj_comp_num
        ## Now for sa_score and size
        for item in cont_tests:
            num_val = getattr(dj_comp, item)
            for my_bin in all_bins[item]:
                if num_val >= all_bins[item][my_bin]["min"] and num_val < all_bins[item][my_bin]["max"]:
                    all_bins[item][my_bin]["val"] += dj_comp_num
    # Now write this out to a file
    for item in my_tests:
        # build a new hist
        new_hist = Histogram.objects.get_or_create(test_made=test, num_samps=num_samps, hist_title= item + "__complex")[0]
        for my_bin in all_bins[item]:
            min_char = str(my_bin)
            max_char = str(my_bin + 1)
            hist_bin = HistBin.objects.get_or_create(hist_id=new_hist, min_char=min_char, max_char=max_char)[0]
            bin_val = HistBinVal.objects.get_or_create(hist_id=new_hist, bin_id=hist_bin, target_id=target_id,
                                                       type_id=item, my_anal_id=react_anal_id, method_id=test,
                                                       num_samps=num_samps)[0]
            bin_val.value = all_bins[item][my_bin]
            bin_val.save()

    ## Now do it for the continuous things
    for item in cont_tests:
        new_hist = Histogram.objects.get_or_create(test_made=test, num_samps=num_samps, hist_title= item + "__complex")[0]
        for my_bin in all_bins[item]:
            min_char = str(all_bins[item][my_bin]["min"])[:7]
            max_char = str(all_bins[item][my_bin]["max"])[:7]
            hist_bin = HistBin.objects.get_or_create(hist_id=new_hist, min_char=min_char, max_char=max_char)[0]
            bin_val = HistBinVal.objects.get_or_create(hist_id=new_hist, bin_id=hist_bin, target_id=target_id,
                                                   type_id="ALL", my_anal_id=react_anal_id, method_id=test,
                                                   num_samps=num_samps)[0]
            bin_val.value = all_bins[item][my_bin]["val"]
            bin_val.save()


def plot_confs_against_vals(target_id):
    """Function to plot 1) Number of compounds modelled
    2) Number of conformations per compound, against
    a) strain
    b) closest approach
    c) rmsd from other structures
    Pop these into PLOTLY and then plot for each target - gives the cutoff value for conformations. 
    """
    # Define the Strain bins
    strain_bins = np.linspace(0, 100, 20)
    # Define the Clash bins
    clash_bins = np.linspace(-0.1, 0.1, 20)
    # Define the RMSD bins
    rmsd_bins = np.linspace(0.05, 1.0, 20)
    # Get the LLConfs down
    my_confs = LLConf.objects.filter(target_id=target_id)
    # First register everything in a CSV file
    target = Target.objects.get(pk=target_id)
    # Find the total number of modelled compounds
    tot_cmps = len(set(my_confs.values_list("mol_id__cmpd_id", flat=True)))
    # Now get the strain vals
    strain_ans = []
    out_f = open(target.title + ".STRAIN.csv", "w")
    out_f.write("strain, num_comps, completeness, num_confs\n")
    for i, val in enumerate(strain_bins):
        # Now get all those in this bin - we're going for appreciating
        strain_confs = my_confs.filter(strain__lte=val)
        # Calculate the num comps
        num_comps = len(set(strain_confs.values_list("mol_id__cmpd_id", flat=True)))
        # Calculate the % of total
        completeness = float(num_comps) / float(tot_cmps)
        # Calculate the num confs / comp
        if num_comps == 0:
            num_confs = "NA"
        else:
            num_confs = float(len(strain_confs)) / float(num_comps)
        # Append to the list
        strain_ans.append({"strain": val, "number of compounds modelled": num_comps, "completeness": completeness, "number of conformations per compound": num_confs})
        out_f.write(str(val) + ", " + str(num_comps) + ", " + str(completeness) + ", " + str(num_confs)+"\n")
    out_f.close()
    # Now do the clash
    clash_ans = []
    out_f = open(target.title + ".CLASH.csv", "w")
    out_f.write("clash, num_comps, completeness, num_confs\n")
    for i, val in enumerate(clash_bins):
        # Now get all those in this bin - we're going for appreciating
        strain_confs = my_confs.filter(clash__gte=-val)
        # Calculate the num comps
        num_comps = len(set(strain_confs.values_list("mol_id__cmpd_id", flat=True)))
        # Calculate the % of total
        completeness = float(num_comps) / float(tot_cmps)
        # Calculate the num confs / comp
        if num_comps == 0:
            num_confs = "NA"
        else:
            num_confs = float(len(strain_confs)) / float(num_comps)
        # Append to the list
        clash_ans.append({"clash": val, "number of compounds modelled": num_comps, "completeness": completeness, "number of conformations per compound": num_confs})
        out_f.write(str(-val) + ", " + str(num_comps) + ", " + str(completeness) + ", " + str(num_confs) + "\n")
    out_f.close()
    # Now do the RMSD
    rmsd_ans = []
    out_f = open(target.title + ".RMSD.csv", "w")
    out_f.write("clash, num_comps, completeness, num_confs\n")
    for i, val in enumerate(rmsd_bins):
        # Now get all those in this bin - we're going for appreciating
        strain_confs = my_confs.filter(rmsd__gte=val)
        # Calculate the num comps
        num_comps = len(set(strain_confs.values_list("mol_id__cmpd_id", flat=True)))
        # Calculate the % of total
        completeness = float(num_comps) / float(tot_cmps)
        # Calculate the num confs / comp
        if num_comps == 0:
            num_confs = "NA"
        else:
            num_confs = float(len(strain_confs)) / float(num_comps)
        # Append to the list
        rmsd_ans.append({"rmsd": val, "number of compounds modelled": num_comps, "completeness": completeness, "number of conformations per compound": num_confs})
        out_f.write(str(val) + ", " + str(num_comps) + ", " + str(completeness) + ", " + str(num_confs) + "\n")
    out_f.close()


def find_overlap_methods(target_id=None, file_path=None):
    """Function to analyse the overlap between different methods"""
    # Get the information
    out_d, samp_dict, target, orig_d, tot_list = init_data(target_id, file_path)
    # Now do the accumulative analysis for the different sample sizes
    write_out_ans(samp_dict, target)
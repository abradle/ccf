from IOhandle.models import Molecule, Target
from WONKA.models import InternalIDLink
from PLIFS.models import PlifProbe, PlifProbeBit
import json
from SUPERSTAR.models import PlifVis, PlifProbeScore, PlifVisGrid
from SUPERSTAR.config import SuperStarConfig
from django.core.exceptions import ValidationError


def make_json(targ, grid_space=0.5):
    """Function to make a JSON fingerprint"""
    # Get the superstar config
    ss_config = SuperStarConfig()
    # A dictionary to map interactions
    int_d = {
 #u'Arom6_Carbonyl': "Aromatic",
 #u'SingleAtomAcceptor_WeakDonor': "Weak_H_bonding",
 #u'Arom5_Carbonyl': "Aromatic",
 u'SingleAtomAcceptor_SingleAtomDonor': "H_bonding",
 u'Arom5_Arom6': "Aromatic",
 u'AcidicGroup_BasicGroup': "Acid_Base",
 u'Hydrophobe_Hydrophobe': "Hydrophobic",
 u'Arom6_Arom6': "Aromatic"}
    res_two_three_dict = {'ALA': "A",
                      'ARG': "R",
                      'ASN': "N",
                      'ASP': "D",
                      'CYS': "C", 
                      'GLU': "E",
                      'GLN': "Q",
                      'GLY': "G",
                      'HIS': "H",
                      'ILE': "I",
                      'LEU': "L",
                      'LYS': "K",
                      'MET': "M",
                      'PHE': "F",
                      'PRO': "P",
                      'SER': "S",
                      'THR': "T",
                      'TRP': "W",
                      'TYR': "Y",
                      'VAL': "V"}
    # List of allowed PH4s
    ph4_l = ss_config.allowed_ph4
    # Get the target
      # The return dicts
    out_list = {}
    # Get the rows - molecules
    # Get the molecules
    mol_refs = None
    # Get the mols we are considering
    mols = Molecule.objects.filter(prot_id__target_id=targ).exclude(prot_id__code__startswith=targ.title)
    # Ensure they are in the same order
    len(mols)
    out_list["rows"] = []
    out_list["mol_freqs"] = {}
    out_list["row_keys"] = {}
    mol_d = {}
    print "MAKING MOLS"
    # Order them by id
    for mol in mols:
        if len(InternalIDLink.objects.filter(mol_id=mol)) == 0:
            out_list["rows"].append(str(mol.pk))
            out_list["row_keys"][str(mol.pk)] = mol.smiles
            mol_d[mol.pk] = str(mol.pk)
            out_list["mol_freqs"][str(mol.pk)] = 0
        else:
            for my_id in InternalIDLink.objects.filter(mol_id=mol):
                out_list["rows"].append(my_id.internal_id + "_" + str(mol.pk))
                out_list["row_keys"][my_id.internal_id + "_" + str(mol.pk)] = my_id.internal_id
                mol_d[mol.pk] = my_id.internal_id + "_" + str(mol.pk)
                out_list["mol_freqs"][my_id.internal_id + "_" + str(mol.pk)] = 0
            # Order by 
    out_list["order_mols_on"] = sorted(out_list["rows"])
    print "DONE MOLS"
    # Get the columns - all possible interactions....
    # List of  residues
    print "DOING PROBES"
    ## Get only filled probes
    res_list = set(PlifProbe.objects.exclude(plifs_plifprobebit_source__isnull=True, plifs_plifprobebit_dest__isnull=True).filter(res_id__isnull=False, target_id=targ).filter(type__in=ph4_l).values_list("res_id__res_name","res_id__res_num","type"))
    print res_list
    # For each residue - get all the probes possible
    print "DONE PROBES"
    out_list["columns"] = []
    out_list["my_data"] = []
    print "DOING RES"
    my_tot_probes = PlifProbe.objects.filter(res_id__res_num__in=[x[1] for x in res_list], type__in=ph4_l)
    print len(my_tot_probes), "TOT PROBES"
    my_tot_probebits = PlifProbeBit.objects.filter(mol_id__prot_id__target_id=targ, dist__lte=5.5)
    print len(my_tot_probebits), "PROBE BITS"
    out_list["col_keys"] = {}
    counter = 0
    unsorted_cols = []
    print "Making columns"
    for res in res_list:
        #my_probes = my_tot_probes.filter(res_id__res_name=res[0], res_id__res_num=res[1], target_id=targ)
        #for probe in my_probes:
        # Assign the name
        my_name = res_two_three_dict[str(res[0])] + str(res[1]) + "_" + str(res[2])
        if my_name not in out_list["columns"]:
            counter += 1
            out_list["col_keys"][my_name] = str(res[0]) + str(res[1])
            #out_list["columns"].append(my_name)
            unsorted_cols.append((my_name, res[0], res[1], res[2]))
    print "Made columns"
    # Sort on res name
    unsorted_cols.sort(key=lambda x: x[1])
    out_list["ressort"] = [x[0] for x in unsorted_cols]
    # Sort on pharmacophore type
    unsorted_cols.sort(key=lambda x: x[2])
    out_list["columns"] = [x[0] for x in unsorted_cols]
    unsorted_cols.sort(key=lambda x: x[3])
    out_list["ph4sort"] = [x[0] for x in unsorted_cols]
    print "DOING PROBE BITS"
    done_list = []
    counter = 0
    tot = 0
    for mol_counter, m in enumerate(mols):
        print m.pk
        my_mol_pbs = my_tot_probebits.filter(mol_id=m) ###### RESTORE THIS - WE WANT ON NATIVE CONTACTS, prot_id=mol.prot_id)
        for p in my_mol_pbs:
            if p.type == "SingleAtomAcceptor_SingleAtomDonor":
                if p.dist > 3.0:
                    continue
            # Now get the score
            print "PROBEBIT", p.pk
            if p.probe_dest_id.res_id:
                res = res_two_three_dict[p.probe_dest_id.res_id.res_name] + str(p.probe_dest_id.res_id.res_num) + "_" + p.probe_dest_id.type
                # Now get the vals
                my_vals = p.probe_source_id.plifprobegridscorenew_set.filter(grid_space=grid_space)
                if my_vals:
                    my_score = my_vals[0].score
                else:
                    my_score = None
            elif p.probe_source_id.res_id:
                res = res_two_three_dict[p.probe_source_id.res_id.res_name] + str(p.probe_source_id.res_id.res_num) + "_" + p.probe_source_id.type
                my_vals = p.probe_dest_id.plifprobegridscorenew_set.filter(grid_space=grid_space)
                if my_vals:
                    my_score = my_vals[0].score
                else:
                    my_score = None
            if my_score < 0.4:
                my_type = "v_bad"
            elif my_score < 1.0:
                my_type = "bad"
            elif my_score < 1.2:
                my_type = "avg"
            elif my_score < 3.0:
                my_type = "good"
            else:
                my_type = "v_good"
            counter += 1
            if p.type not in int_d:
                continue
            #  Add this to the out list
            # -> OLD WAY OF MAKING INTERACTIONS int_d[p.type]
            out_list["my_data"].append([mol_d[m.pk], res, [my_type], str(p.dist), p.type, ["XTAL"], my_score]) 
            # To capture what has already been done
            done_list.append(str(mol_d[m.pk]) + str(res) + str(int_d[p.type]))
        # Now get the associated LLOOMMPPAA bits
    if tot:
        for key in out_list["mol_freqs"]:
            out_list["mol_freqs"][key] = float(out_list["mol_freqs"][key]) / float(tot)
        # Now order the molecules like this
    out_list["order_mols_two"] = sorted(out_list["mol_freqs"], key=out_list["mol_freqs"].get, reverse=True)
    # Order the 
    print "DONE RES"
    # Get the data for each molecule
    # Return a JSON
    tot_list = json.dumps(out_list)
    my_fps = PlifVisGrid()
    my_fps.target_id = targ
    my_fps.grid_space = grid_space
    try:
        my_fps.validate_unique()
    except ValidationError:
        my_fps = PlifVisGrid.objects.get(target_id=targ, grid_space=grid_space)
    my_fps.json_text = tot_list
    my_fps.save()
    return tot_list


def make_grids(target_id, grid_list=[0.3,0.5,0.7,1.0]):
    """Function to make the JSON grids for different grid spacings"""
    targ = Target.objects.get(pk=target_id)
    # Now loop over them and make them
    for grid_space in grid_list:
        make_json(targ, grid_space)
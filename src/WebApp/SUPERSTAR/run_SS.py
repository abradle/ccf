import subprocess
import os
from multiprocessing import Pool
from SUPERSTAR.config import SuperStarConfig
from SUPERSTAR.models import PlifProbeScore, PlifProbeGridScoreNew
from rdkit import Chem
from rdkit.Chem import MCS
from rdkit.Chem import AllChem
import math
from IOhandle.models import Protein, Molecule, Target
from PLIFS.models import PlifProbe
from PLIFS.functions import restore_probs
from SUPERSTAR.parse_grid import parse_grid_file, interpolate
from SUPERSTAR.functions import make_json
import time, ast
from django.core.exceptions import ValidationError
import numpy as np
import threading

num_workers = 0
p_dict = {}


def get_centroid(rdmol):
    """Function to return the centroid of a molecule"""
    try:
        for conf in rdmol.GetConformers():
            centroid = Chem.rdMolTransforms.ComputeCentroid(conf)
            x = centroid[0]
            y = centroid[1]
            z = centroid[2]
        return "%s %s %s" % (x, y, z)
    except AttributeError:
        return False


def icm_add_hyds(pdb_code):
    """ """
    out_file = 'W:/Informatics/Pharmacophore/anthony/DPhil/CODE/CHOC/src/WebApp/SUPERSTAR/PROTEINS/' + pdb_code + "_OUT.pdb" 
    pdb_file = 'W:/Informatics/Pharmacophore/anthony/DPhil/CODE/CHOC/src/WebApp/SUPERSTAR/PROTEINS/' + pdb_code + ".pdb" 
    out_icm = """call _macro
# Load the file_in
openFile '""" + pdb_file + """'
# Convert 
convertObject a_1.* yes yes yes yes yes yes no
#Save it 
my_s = '""" + out_file + """'
write pdb a_1.* my_s
q
    """
    out_f = open("out_icm", "w")
    out_f.write(out_icm)
    out_f.close()
    #subprocess.call(["C:\Program Files (x86)\Molsoft LLC\ICM-Pro\icmconsole.exe", "out_icm"])
    return out_file


def run_worker(jobname, ss_c, probe, centroid, mol, rd_conf, h_path, my_probes, grid_space):
    """Function to run each thread as different thread"""
    # NOW DEFINE THE JOBNAME
    global num_workers
    global p_dict
    num_workers += 1
    print jobname, "WORKER", num_workers
    ins_file = '''JOBNAME ''' + jobname + '''
PROBENAME ''' + ss_c.probe_dict[probe] + '''
GRID_SPACING ''' + str(grid_space) + '''
CALC_CONTOUR_SURFACE 0
RUN_TYPE SUPERSTAR
DATA_SOURCE CSD
SHELL_VALUES 0.5 0
SIGCHECK_METHOD POISSON SHELL
PROPENSITY_CORRECTION LOGP DEFAULT DEFAULT
MIN_CAVITY_VOLUME 1
CAVITY_RADIUS 3
DATA_TYPE ISOSTAR
MOLECULE_FILE ''' + h_path + '''
CAVITY_DETECTION 1
MIN_PSP 4
MIN_PROPENSITY 0.0
SAVE_CAVITY MESH
USE_TORSIONS 1
MAP_FORMAT INSIGHT_ASCII_GRID

  '''
    # Now check if there is a centroid
    if centroid:
        ins_file = ins_file + "CAVITY_ORIGIN " + str(centroid) + "\n"
    # Now define the ins file
    ins_file_path = os.path.join(ss_c.ins_path, jobname + ".ins")
    input_file = open(ins_file_path, "w")
    input_file.write(ins_file)
    input_file.close()
    s_path = r"C:\Program Files (x86)\CCDC\SuperStar 2.1.2\bin\superstar_app.exe"
    cmd = [s_path, ins_file_path]
    # Now run it
    in_path = os.path.join(ss_c.out_path, jobname + ".grd")
    print "SPAWNING SUPERSTAR"
    if not os.path.isfile(in_path):
        subprocess.Popen(cmd, shell=True)
    # Wait for the files to write
    while 1 == 1:
        if not os.path.isfile(in_path):
            time.sleep(1)
            continue
        break
    # Now get the grid back
    print "PARSING GRID"
    grid_space, grid_dict = parse_grid_file(in_path)
    print "GRID PARSED"
    # Now read in the grid
    # Now go through the relevant PLIF probes
    # Now loop through these and score them
    for p in my_probes:
        # Get the atom ids
        propensities = []
        atom_ids = ast.literal_eval(p.atom_ids)
        for atom_id in atom_ids:
            # Get the atom
            my_atom = rd_conf.GetAtomPosition(atom_id)
            # Now loop through the atom ids
            propensities.append(interpolate(grid_space, grid_dict, my_atom.x, my_atom.y, my_atom.z))
        # Now get the geometric mean of this
        propensity = np.exp(np.sum(np.log(propensities))) / len(atom_ids)
        p_dict[p.pk] = propensity
    num_workers -= 1
    return None


class myThread (threading.Thread):
    """Class to start new threads"""
    def __init__(self, threadID, jobname, ss_c, probe, centroid, mol, rd_conf, h_path, my_probes, grid_space):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = jobname
        #self.counter = counter
        self.jobname = jobname
        self.ss_c = ss_c
        self.probe = probe
        self.centroid = centroid
        self.mol = mol
        self.rd_conf = rd_conf
        self.h_path = h_path
        self.my_probes = my_probes
        self.grid_space = grid_space

    def run(self):
        print "Starting " + self.name
        run_worker(self.jobname, self.ss_c, self.probe, self.centroid, self.mol, self.rd_conf, self.h_path, self.my_probes, self.grid_space)


def make_ins_file(target_id=2, grid_space=0.5):
    """Function to take a target and do analysis of all pharma points bound to this target
    in xtal structures"""
    global p_dict
    global num_workers
    ss_c = SuperStarConfig()
    target = Target.objects.get(pk=target_id)
    my_mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__contains=target.title)
    # Now change to this directoy
    os.chdir(ss_c.out_path)
    # Loop through the appropriate proteins
    for mol in my_mols:
        # Define the files to clean up
        rd_mol = Chem.MolFromMolBlock(str(mol.sdf_info))
        rd_conf = rd_mol.GetConformer()
        file_list = []
        # Define the pdb_file
        pdb = mol.prot_id.code
        out_path = os.path.join(ss_c.base_path, pdb + ".pdb")
        out_f = open(out_path, "w")
        pdb_mol = Chem.MolFromPDBBlock(str(mol.prot_id.pdb_info))
        out_f.write(Chem.MolToPDBBlock(pdb_mol))
        out_f.close()
        out_file = icm_add_hyds(pdb)
        # Now use ICM to add Hydrogens to this...
        h_path = os.path.join(ss_c.base_path, pdb + "_OUT.pdb")
        centroid = get_centroid(rd_mol)
        my_keys = ss_c.probe_dict.keys()
        counter = -1
        while 1 == 1:
            if num_workers >= 6:
                continue
            elif counter < len(my_keys)-1:
                counter += 1
                probe = my_keys[counter]
                jobname = pdb + "_" + probe + "_" + str(grid_space).replace(".", "_")
                my_probes = PlifProbe.objects.filter(mol_id=mol, type__in=ss_c.type_dict[probe])
                len(my_probes)
                my_t = myThread(counter, jobname, ss_c, probe, centroid, mol, rd_conf, h_path, my_probes, grid_space)
                my_t.start()
            if num_workers == 0 and counter == len(my_keys) -1:
                break
    # Now add the scores - with this multithreaded EG
    for p in p_dict:
        p_score = PlifProbeGridScoreNew()
        plif_p = PlifProbe.objects.get(pk=p)
        p_score.plif_probe = plif_p
        p_score.score = p_dict[p]
        p_score.grid_space = grid_space
        try:
            p_score.validate_unique()
            p_score.save()
        except ValidationError:
            p_score = PlifProbeGridScoreNew.objects.get(plif_probe=plif_p, grid_space=grid_space)
            p_score.score = p_dict[p]
            p_score.save()


def run_anal(target_id=2, grid_spaces=[0.1, 0.2, 0.3, 0.4, 0.5,
                                       0.6, 0.7, 0.8, 0.9, 1.0]):
    """Function to run the analysis on multiple grid spacings"""
    global p_dict
    # Loop over the options
    for grid_s in grid_spaces:
        # Run the function
        # Re-zero the dictionary
        p_dict = {}
        make_ins_file(target_id, grid_space=grid_s)


def run_targs():
    """Function to run the SUPERSTAR analysis for all targets"""
    targs = Target.objects.exclude(title="DUMMY")
    for targ in targs:
        # Now run the analysis for this target
        restore_probs(targ.pk)
        run_anal(targ.pk, grid_spaces=[0.5])
        make_json(targ, grid_space=0.5)
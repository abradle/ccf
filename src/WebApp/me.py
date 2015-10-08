import os,sys
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "WebApp.settings")
from IOhandle.models import Compound,Molecule,Protein,ActivityPoint,Project,Target
from MMPMaker.functions import index_hydrogen_change, make_mol_mmp, make_mmp_database, make_list_mmps, act_mmp_3d, make_ph4_difference_points, find_pharma_changes
from MMPMaker.models import MMPDiffMap, MMPComparison, ActPharmaPoint, MMPFrag
from WONKA.models import Water,InternalIDLink, Residue, ResShift
from Group.models import Group
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from django.core.exceptions import ValidationError
import sys, os
from IOhandle.functions import add_new_comp
import csv
import uuid, glob
from rdkit.Chem import ChemicalFeatures
from IOhandle.functions import add_new_comp
from Pharmacophore.models import find_pharmacophore_points
from LLOOMMPPAA.functions import make_synth_points
from MMPMaker.functions import make_mmp_database,find_frag_com
from loading import load_activity_data
from MMPMaker.functions import index_hydrogen_change, make_mol_mmp,  make_list_mmps, act_mmp_3d, make_ph4_difference_points, find_pharma_changes
from WONKA.functions import find_target


## Main function
from loading import initialise_dummys
initialise_dummys()
error_log = open("ERROR.txt","w")
#####error_log.write("BRD4A")
#####find_target("BRD4A")

#
##error_log.write("JMJD2AA")
#find_target("JMJD2AA")
#error_log.write("BRD1A")
#find_target("BRD1A")
#error_log.write("BAZ2BA")
#find_target("BAZ2BA")
error_log.write("PHIPA")
find_target("PHIPA")
#error_log.write("JMJD2DA")
#find_target("JMJD2DA")
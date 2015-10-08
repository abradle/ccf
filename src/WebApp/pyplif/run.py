from openbabel import OBMol, OBConversion, OBResidueIter
import pyplif as pp
import operator
from IOhandle.models import Protein, Molecule
from LLOOMMPPAA.models import PlifMethod,PlifRes,PlifBit,Plif,SynthPoint,PlifBitInstance
from django.core.exceptions import ValidationError
import os
import ast
import sys
import tempfile
from rdkit import Chem
from auxfuns import *

def django_run(target, opt="XTAL"):
    """Function to take multiple confs of ONE ligand and generate their PLIFS against one template protein"""
    # Set up the OpenBaebel conversion modules
    sdconv = OBConversion()
    ligref = OBMol()
    # Define the residues and the proteisn to analyse
    if os.path.isfile(os.path.join(os.path.split(sys.argv[0])[0], 'data/res_def.py')):
    	res_d = [trans_res(x) for x in ast.literal_eval(open(os.path.join(os.path.split(sys.argv[0])[0], 'data/res_def.py')).read())[target.title].split()]
    print res_d
# Molecules
	  # Now read in the ligand
    plif_method = PlifMethod()
    plif_method.text= "PYPLIF"
    feature_list = ["POLAR","FACE","EDGE","ACCEPTOR","DONOR","NEGATIVE","POSITIVE"]
    try:
        plif_method.validate_unique()
        plif_method.save()
    except ValidationError:
    	  plif_method = PlifMethod.objects.get(text="PYPLIF")
    out_d = {}
    counter = 0
# Create a file for the protein
    t = tempfile.NamedTemporaryFile(suffix=".pdb",delete=False)
    my_prot = Protein.objects.get(code=target.title+"TEMP")
    t.write(my_prot.pdb_info.name)
    t.close()
    protref = read_prot(t.name, res_d)
    t = tempfile.NamedTemporaryFile(suffix=".sdf",delete=False)
    t.close()
    sdconv.SetInFormat("sdf")
    if opt == "XTAL":
        mols = Molecule.objects.exclude(prot_id__code__contains=target.title).filter(prot_id__target_id=target)
    elif opt == "LLOOMMPPAA":
    	  mols = []
    	  sps = SynthPoint.objects.filter(target_id=target)
    	  for s in sps:
    	  	  mols.extend([m for m in s.mol_id.all()])
    else:
    	  print "UNKNOWN OPTION"
    	  return
    for dj_mol in mols:
    	  out_sd = Chem.SDWriter(t.name)
    	  out_sd.write(Chem.MolFromMolBlock(str(dj_mol.sdf_info)))
    	  out_sd.close()
    	  sdconv.ReadFile(ligref, t.name)
    	  # Now make the new plif
    	  new_plif = Plif()
    	  new_plif.mol_id = dj_mol
    	  new_plif.prot_id = my_prot
    	  new_plif.method_id = plif_method
    	  try:
    	  	  new_plif.validate_unique()
    	  	  new_plif.save()
    	  except ValidationError:
    	  	  new_plif = Plif.objects.get(mol_id=dj_mol,prot_id=my_prot,method_id=plif_method)
    	  lig_name = ligref.GetTitle().strip(",")
    	  prot_name = lig_name.split("_")[0]
    	  ligref.AddHydrogens()
    	  counter +=1
    	  refresdict = pp.getresiduedict(protref, res_d)
    	  new_d = get_fp(protref,ligref, res_d)
    	  for res in new_d:
    	  	  new_res =  PlifRes()
    	  	  new_res.res_name = res[:3]
    	  	  new_res.res_num = int(res[3:])
    	  	  new_res.prot_id = my_prot
    	  	  try:
    	  	  	  new_res.validate_unique()
    	  	  	  new_res.save()
    	  	  except ValidationError:
    	  	  	  new_res = PlifRes.objects.get(res_name=res[:3],res_num=int(res[3:]),prot_id=my_prot)
    	  	  new_plif.res_id.add(new_res)
    	  	  for bit_num, bit in enumerate(new_d[res]):
    	  	  	  new_bit = PlifBit()
    	  	  	  new_bit.feature = feature_list[bit_num]
    	  	  	  new_bit.method_id = plif_method
    	  	  	  new_bit.res_id = new_res
    	  	  	  try:
    	  	  	  	  new_bit.validate_unique()
    	  	  	  	  new_bit.save()
    	  	  	  	  my_fun(dj_mol,new_bit,new_plif,bit)
    	  	  	  except ValidationError:
										new_bit = PlifBit.objects.get(feature=feature_list[bit_num],method_id=plif_method,res_id=new_res)
										new_bit.save()
										new_plif.bit_id.add(new_bit)
										my_fun(dj_mol,new_bit,new_plif,bit)
										
    	  ligref = OBMol()
    	  notatend = sdconv.Read(ligref)

def get_fp(protref,ligref,res_d):
	refresdict = pp.getresiduedict(protref, res_d)
	new_d = dict((k, v) for (k, v) in refresdict.iteritems() if k in res_d)
	pp.otherinteractions(new_d, res_d, protref, ligref, [])
	# New_d is now this bit string. Return the output
	return new_d

def write_res(out_d,res_d):
	"""Function to write the results out in a dictionary"""
	# Now make the csv file
	out_f = open("out.csv","w")
	# First the header
	out_f.write("mol,")
	out_l = []
	# Put gaps either side
	for res in res_d:
		out_l.extend([res,"","","","",""])
	# Now write this list
	out_f.write(",".join(out_l))
	out_f.write("\n")
	# Now the val for each mol
	for mol in out_d:
		out_f.write(mol+",")
		for res in res_d:
			out_f.write(",".join([str(int(x)) for x in list(out_d[mol][res])]))
		out_f.write("\n")
	out_f.close()


def read_prot(prot_file, res_d):
	"""Function to read in a protein to an OBMol"""
	conv = OBConversion()
	protref = OBMol()
	conv.SetInFormat("pdb")
	conv.ReadFile(protref,prot_file)
	# Now assign the residue names
	i = 0
	my_res = []
	for residue in OBResidueIter(protref):
		i+=1
		residue.SetName(residue.GetName()+str(residue.GetNum()))
		my_res.append(residue.GetName())
	# Now check that all the residues exist and print out if not
	fail_counter = 0
	fail_list = []
	# Loop through the res and check they are in the list
	for res_me in res_d:
		if res_me not in my_res:
			fail_counter += 1
			fail_list.append(res_me)
	# If it's out of register by one do again
	if fail_counter > 0:
		i = 0
		my_res = []
		for residue in OBResidueIter(protref):
			i+=1
			residue.SetName(residue.GetName()+str(residue.GetNum()))
			my_res.append(residue.GetName())
		# Now check that all the residues exist and print out if not
		fail_counter = 0
		fail_list = []
		# Loop through the res and check they are in the list
		for res_me in res_d:
			if res_me not in my_res:
				fail_counter += 1
				fail_list.append(res_me)		
				out_err.write(prot_file+",")
				out_err.write(str(fail_counter)+"\n")
				out_err.write(str(fail_list))
				out_err.write(str(my_res))
				out_err.write(str(res_d))
	
	protref.AddHydrogens()
	return protref


def trans_res(res):
  """Function to take a RES of form blah./^V64 and output VAL64"""
  d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
  rev_d = dict((v, k) for (k, v) in d.iteritems())
  # Take the res
  out_res = rev_d[res.split("^")[1][0]]
  return out_res+res.split("^")[1][1:]


def get_dict(file):
	"""Function to make a dictionary of residues for the binding site
	of each protein. Based on a text file"""
	res_list = [x.rstrip() for x in open(file).readlines()]
	# Make the output dictionary
	out_d = {}
	# write it
	for res in res_list:
	  if res[2:6] in out_d:
	    out_d[res[2:6]].append(trans_res(res))
	  else:
	    out_d[res[2:6]] = [trans_res(res)]
	# Define the lists to count the aligned residues
	prot_list = []
	tot_list = []
	restrictive_list = []
	tot_dict = {}
	for prot in out_d:
	  if out_d[prot][0] == "ILE10":
	    prot_list.append(prot)
	    tot_list.extend(out_d[prot])
	    for res in  out_d[prot]:
	      if res in tot_dict:
	        tot_dict[res] +=1
	      else:
	        tot_dict[res] =1
	list(set(tot_list))
	sorted_x = sorted(tot_dict.iteritems(), key=operator.itemgetter(1),reverse=True)
	core_res = [x[0] for x in sorted_x[:16]]
	# Now make a restrictive list - only those with all residues
	for prot in out_d:
	  if set(core_res).issubset(set(out_d[prot])):
	  	restrictive_list.append(prot)
	return list(set(core_res)), restrictive_list[:30]



def all_to_all():
	"""Function to compare all to all"""
	# Set up the OpenBaebel conversion modules
	sdconv = OBConversion()
	ligref = OBMol()
	# Define the residues and the proteisn to analyse
	res_d, prot_list = get_dict("myFirstFile.txt")
	# Now read in the ligand
	sdconv.SetInFormat("sdf")
	notatend = sdconv.ReadFile(ligref,"../mols.sdf")
	out_d = {}
	counter = 0
	# Now read the ligand file
	while notatend:
		lig_name = ligref.GetTitle().strip(",")
		prot_name = lig_name.split("_")[0]
		if prot_name not in prot_list:
			ligref = OBMol()
			notatend = sdconv.Read(ligref)
			continue
		ligref.AddHydrogens()
		counter +=1
		print counter
		for j, my_prot in enumerate(prot_list):
			protref = read_prot(r"C:\www\Protoype\media_coninchi\pdb" + "\\" + my_prot + "al.pdb", res_d)
			# Get the reference dictionary
			refresdict = pp.getresiduedict(protref, res_d)
			# Update this dict, to only residues in the binding site
			new_d = get_fp(protref,ligref, res_d)
			# Make sure it is a unique name for the output
			while lig_name in out_d:
				lig_name = lig_name + "Z"
			# Add it to the dict 
			out_d[lig_name+my_prot] = {}
			for res in new_d:
				# Assign each residue the scores for each molecule
				out_d[lig_name+my_prot][res] = new_d[res]
		# Make the ligand
		ligref = OBMol()
		notatend = sdconv.Read(ligref)
	# Now write the results out
	write_res(out_d, res_d)


def one_to_many():
	"""Function to take multiple confs of ONE ligand and generate their PLIFS against one template protein"""
	# Set up the OpenBaebel conversion modules
	sdconv = OBConversion()
	ligref = OBMol()
	# Define the residues and the proteisn to analyse
	res_d, prot_list = get_dict("myFirstFile.txt")
	# Now read in the ligand
	sdconv.SetInFormat("sdf")
	notatend = sdconv.ReadFile(ligref,"../out.sdf")
	out_d = {}
	counter = 0
	my_prot = "1qmz"
	protref = read_prot(r"C:\www\Protoype\media_coninchi\pdb" + "\\" + my_prot + "al.pdb", res_d)
	# Now read the ligand file
	while notatend:
		lig_name = ligref.GetTitle().strip(",")
		prot_name = lig_name.split("_")[0]
		ligref.AddHydrogens()
		counter +=1
		print counter
		# Get the reference dictionary
		refresdict = pp.getresiduedict(protref, res_d)
		# Update this dict, to only residues in the binding site
		new_d = get_fp(protref,ligref, res_d)
		# Add it to the dict 
		out_d[lig_name+str(counter)] = {}
		for res in new_d:
			# Assign each residue the scores for each molecule
			out_d[lig_name+str(counter)][res] = new_d[res]
		# Make the ligand
		ligref = OBMol()
		notatend = sdconv.Read(ligref)
	# Now write the results out
	write_res(out_d, res_d)	

# Script to read in a series of ligands (as an SD file) and output a CSV file of the bit vectors
# Define the OB objects
if __name__ == "__main__":
	# Make the error file
	out_err = open("out.std.err","w")
	out_err.write("file,errors\n")
	one_to_many()#all_to_all()



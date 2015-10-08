#!/usr/bin/env python

import getopt, sys, os, copy, glob
from openbabel import OBMol, OBConversion
from optparse import OptionParser
from time import time
from tanimoto_coef import *
from ring import *
from interactions import *


if __name__ == "__main__":
	x = time()
  
	#  Default configuration
	config   = "config.txt"

	#  Identifiers initialization
	protein_reference = ""
	ligand_reference = ""
	residue_of_choice = []
	protein_ligand_folder = ""
	output_file = ""

	parser = OptionParser()
	parser.add_option("-c", "--config", dest="config",
					  help="read config from FILE", metavar="FILE")
	parser.add_option("-o", "--output", dest="output_file",
					  help="write result to FILE", metavar="FILE")
	parser.add_option("-f", "--fixed", action="store_true",
					  dest="fixed", default=False, help="Only using protein reference, ignoring the flexible residue(s)")
	parser.add_option("-d", "--debug_mode", action="store_true",
					  dest="debug_mode", default=False, help="Activate debug mode")
	parser.add_option("-a", "--autocharge", action="store_true",
						dest="autocharge", default=False, help="assign formal charge with open babel")
	parser.add_option("-r", "--protref", action = "store_true",
						dest="protein_reference", default=False, help="Using the protein reference instead of binding_site_fixed.mol2 in PLANTS")
#	parser.add_option("-p", "--plants", action="store_true",
#					  dest="plants", default=False, help="Use confxx_entryxx_protein.mol2 from PLANTS output")
	
	(options_parser, args) = parser.parse_args()
	
	#  Use the config file from the command option instead
	#  when it is supplied through command option.
	if options_parser.config:
		config = options_parser.config
	
	try:
		configread  = open(config, 'r')
	except:
		print "The config file: '", config, "' can not be found"
		sys.exit(1)

	configlines = [line for line in configread]
	configread.close()

	for line in configlines:
		options = line.split()
		if not options:
			continue
		if options[0] == "protein_reference":
			protein_reference = options[1]
		if options[0] == "ligand_reference":
			ligand_reference = options[1]
		if options[0] == "residue_of_choice":
			residue_of_choice = options[1:]
		if options[0] == "protein_ligand_folder":
			protein_ligand_folder = options[1]
		if options[0] == "output_file":
			output_file = options[1]
	
	if options_parser.output_file:
		output_file = options_parser.output_file
	
	print protein_ligand_folder
	os.chdir(protein_ligand_folder)
	print "IN THE FOLDER"
	conflist = open('features.csv', 'r')

	firstline = conflist.readline()
	mollisttemp = [line for line in conflist]
	mollist   = []
	scorelist = []
	for mol in mollisttemp:
		mollist.append(mol.split(',')[0])
		scorelist.append(mol.split(',')[1])
	os.chdir('..')

	# opening the molecule files
	pbf = protein_ligand_folder + '/protein_bindingsite_fixed.mol2'
	conv = OBConversion()
	conv.SetInFormat("mol2")

	protfix  = OBMol()
	protref  = OBMol()
	ligref   = OBMol()
	docklig  = OBMol()
	dockprot = OBMol()

	conv.ReadFile(protfix, pbf)
	conv.ReadFile(protref, protein_reference)
	conv.ReadFile(ligref, ligand_reference)
	
	# creating backbone atom list
	protrefhandle = open(protein_reference, 'r')
	protreflines = [line for line in protrefhandle]
	pbfhandle = open(pbf, 'r')
	pbflines = [line for line in pbfhandle]
	backboneatom = ['C', 'N', 'O', 'H', 'CA', 'HA']
	backboneatomlistref = []
	backboneatomlist = []
	for line in protreflines:
		if len(line.split()) > 1:
			if line.split()[1] in backboneatom:
				backboneatomlistref.append(int(line.split()[0]))
	for line in pbflines:
		if len(line.split()) > 1:
			if line.split()[1] in backboneatom:
				backboneatomlist.append(int(line.split()[0]))
    
	refresdict    = getresiduedict(protref, residue_of_choice)
	refringdict   = getringdict(protref)
	fixringdict   = getringdict(protfix)

	ringinteraction(refresdict, refringdict, residue_of_choice, protref, ligref)
	otherinteractions(refresdict, residue_of_choice, protref, ligref, backboneatomlistref) 

	cvsoutdict = {}
	bitarraydict = {}
	filelist = {}
	fixresdict_template = getresiduedict(protfix, residue_of_choice)
	for compound in mollist:
		print compound
		fixresdict = copy.deepcopy(fixresdict_template)
		os.chdir(protein_ligand_folder)
		molname = glob.glob('*' + compound.split('_entry')[1] + '*')
		os.chdir('..')
		molname.sort()
		ligand_file = protein_ligand_folder + molname[0]
		protein_file = protein_ligand_folder + 'protein_bindingsite_fixed.mol2'#molname[1]
		conv.ReadFile(docklig, ligand_file)
		if options_parser.protein_reference:
			conv.ReadFile(dockprot, protein_file)
			ringinteraction(fixresdict, refringdict, residue_of_choice, protref, docklig)
			otherinteractions2(fixresdict, residue_of_choice, protref, docklig, backboneatomlistref)
			hbonddockprot(fixresdict, residue_of_choice, dockprot, docklig)
		elif options_parser.fixed:
			ringinteraction(fixresdict, refringdict, residue_of_choice, protref, docklig)
			otherinteractions(fixresdict, residue_of_choice, protref, docklig, backboneatomlistref)
		else:
			conv.ReadFile(dockprot, protein_file)
			ringinteraction(fixresdict, fixringdict, residue_of_choice, protfix, docklig)
			otherinteractions(fixresdict, residue_of_choice, protfix, docklig, backboneatomlist)
			hbonddockprot(fixresdict, residue_of_choice, dockprot, docklig)
			
		cvsoutdict[compound]   = gettcfromdict(refresdict, fixresdict, residue_of_choice)
		bitarraydict[compound] = collectbit(fixresdict, residue_of_choice)
		filelist[compound]     = ligand_file
    
	refstringbit = collectbit(refresdict, residue_of_choice)
	
	outfile = open(output_file, 'w')
	outfile.write("".ljust(78))
	for res in residue_of_choice:
		outfile.write(res.ljust(7))
	outfile.write('\n')
	outfile.write(ligand_reference.ljust(77))
	outfile.write(" %s" % refstringbit)
	outfile.write('\n')
	for compound, score in zip(mollist, scorelist):
		outfile.write(filelist[compound].ljust(60))
		outfile.write(" %s" % score.ljust(9) )
		outfile.write(" %.3f " % cvsoutdict[compound])
		outfile.write(" %s" % bitarraydict[compound])
		outfile.write('\n')
                    
	outfile.close()
	y = time()
	print 'Total time taken %.3f s.' % (y-x)


#!/usr/bin/env python

import numpy as np
from openbabel import OBMolAtomIter, OBResidueIter, OBResidueAtomIter, OBAtomAtomIter
from ring import *
from atom_property import *
from bitarray import bitarray

def getresiduedict(protein, residuechoice):
	resdict = {}
	for residue in OBResidueIter(protein):
		resdict[residue.GetName()] = bitarray('0000000')
	for residue in residuechoice:
		resdict[residue] = bitarray('0000000')
	return resdict

def ringinteraction(resdict, ringdict, residuechoice, protein, ligand):
	RADTODEG = 180 / np.pi
	for ring in ringdict:
		if ringdict[ring][0] in residuechoice:
			coords = []
			numcoord = 0
			while numcoord < 3:
				coords.append(np.array([protein.GetAtom(ringdict[ring][1][numcoord]).x(),
				protein.GetAtom(ringdict[ring][1][numcoord]).y(), protein.GetAtom(ringdict[ring][1][numcoord]).z()]))
				numcoord += 1
			for ringligand in ligand.GetSSSR():
				if ringligand.IsAromatic():
					ligcoords = []
					numcoord = 0
					while numcoord < 3:
						ligcoords.append(np.array([ligand.GetAtom(ringligand._path[numcoord]).x(),
						ligand.GetAtom(ringligand._path[numcoord]).y(), ligand.GetAtom(ringligand._path[numcoord]).z()]))
						numcoord += 1
	
					ringpath = ringdict[ring][1]
					inrange = ringdistance(ringpath, ringligand, protein, ligand)
					if inrange:
						ringcross = getringcross(coords)
						ligringcross = getringcross(ligcoords)
						dot = np.dot(ringcross, ligringcross)
						cross_modulus = np.sqrt((ringcross*ringcross).sum())
						ligcross_modulus = np.sqrt((ligringcross*ligringcross).sum())
						cos_angle = dot / cross_modulus / ligcross_modulus
						ring_angle    = np.arccos(cos_angle) * RADTODEG
						#  The result of arccos is ranging from 0 to 180 degrees.
						#  Ring angle of 0 deg = 180 deg = parallel to each other
						#  Ring angle of 90 deg = perpendicular to each other

						#  identifying aromatic face to face
						if (30.0 >= ring_angle) or (150.0 <= ring_angle):
							resdict[ringdict[ring][0]] |= bitarray('0100000')
						#  identifying aromatic edge to face
						if (30.0 <= ring_angle <= 150.0):
							resdict[ringdict[ring][0]] |= bitarray('0010000')
    
def otherinteractions(resdict, residuechoice, protein, ligand, backboneatomlist):
	for residue in OBResidueIter(protein):
		residuename = residue.GetName()
		if residuename in residuechoice:	
			for atom in OBResidueAtomIter(residue):
				if not atom.IsHydrogen():
					for atomlig in OBMolAtomIter(ligand):
						if not atomlig.IsHydrogen():
							distance = atom.GetDistance(atomlig)
							if distance <= 4.5:
								if isnonpolar(atomlig) & isnonpolar(atom):
									#print "POLAR"
									resdict[residuename] |= bitarray('1000000')
								if distance <= 4.0:
									setformalcharge(atom, backboneatomlist)
									fakebbatomlist = []
									setformalcharge(atomlig, fakebbatomlist)
									if (atom.GetFormalCharge()>0) & (atomlig.GetFormalCharge()<0):
										#print "NEGATIVE"
										resdict[residuename] |= bitarray('0000010')
									if (atom.GetFormalCharge()<0) & (atomlig.GetFormalCharge()>0):
										#print "POSITIVE"
										resdict[residuename] |= bitarray('0000001')
									if distance <= 3.5:
										if atom.IsHbondDonor() & atomlig.IsHbondAcceptor():
											for neighborDon in OBAtomAtomIter(atom):
												if neighborDon.IsHydrogen():
													angle = atom.GetAngle(neighborDon, atomlig)
													if angle>135.0:
														resdict[residuename] |= bitarray('0001000')
#											donorresidue = atom.GetResidue()
#											for atomres in OBResidueAtomIter(donorresidue):
#												if atomres.IsHydrogen() & atomres.IsConnected(atom):
#													angle = atom.GetAngle(atomres, atomlig)
#													if angle>135.0:
#														resdict[residuename] |= bitarray('0001000')
										if atom.IsHbondAcceptor() & atomlig.IsHbondDonor():
											for neighborDon in OBAtomAtomIter(atomlig):
												if neighborDon.IsHydrogen():
													#print True
													angle = atomlig.GetAngle(neighborDon, atom)
													#print angle
													if angle>135.0:
														#print "ACCEPTOR"
														resdict[residuename] |= bitarray('0000100')
											
#											for atomres in OBMolAtomIter(ligand):
#												if atomres.IsHydrogen() & atomres.IsConnected(atomlig):
#													angle = atomlig.GetAngle(atomres, atom)
#													if angle>135.0:
#														resdict[residuename] |= bitarray('0000100')

def otherinteractions2(resdict, residuechoice, protein, ligand, backboneatomlist):
	RigidHdonors = ['HIS', 'GLN', 'ASN', 'ARG', 'TRP']
	for residue in OBResidueIter(protein):
		residuename = residue.GetName()
		if residuename in residuechoice:
			for atom in OBResidueAtomIter(residue):
				if not atom.IsHydrogen():
					for atomlig in OBMolAtomIter(ligand):
						print atomlig
						if not atomlig.IsHydrogen():
							distance = atom.GetDistance(atomlig)
							if distance <= 4.5:
								if isnonpolar(atomlig) & isnonpolar(atom):
									resdict[residuename] |= bitarray('1000000')
								if distance <= 4.0:
									setformalcharge(atom, backboneatomlist)
									fakebbatomlist = []
									setformalcharge(atomlig, fakebbatomlist)
									if (atom.GetFormalCharge()>0) & (atomlig.GetFormalCharge()<0):
										resdict[residuename] |= bitarray('0000010')
									if (atom.GetFormalCharge()<0) & (atomlig.GetFormalCharge()>0):
										resdict[residuename] |= bitarray('0000001')
									if distance <= 3.5:
										if residuename[:3] in RigidHdonors:
											if atom.IsHbondDonor() & atomlig.IsHbondAcceptor():
												for neighborDon in OBAtomAtomIter(atom):
													if neighborDon.IsHydrogen():
														angle = atom.GetAngle(neighborDon, atomlig)
														if angle>135.0:
															resdict[residuename] |= bitarray('0001000')
										elif atom.GetIdx() in backboneatomlist:
											if atom.IsHbondDonor() & atomlig.IsHbondAcceptor():
												for neighborDon in OBAtomAtomIter(atom):
													if neighborDon.IsHydrogen():
														angle = atom.GetAngle(neighborDon, atomlig)
														if angle>135.0:
															resdict[residuename] |= bitarray('0001000')
										if atom.IsHbondAcceptor() & atomlig.IsHbondDonor():
											for neighborDon in OBAtomAtomIter(atomlig):
												if neighborDon.IsHydrogen():
													angle = atomlig.GetAngle(neighborDon, atom)
													if angle>135.0:
														resdict[residuename] |= bitarray('0000100')

def hbonddockprot(resdict, residuechoice, protein, ligand):
	for residue in OBResidueIter(protein):
		residuename = residue.GetName()
		if residuename in residuechoice:
			for atom in OBResidueAtomIter(residue):
				if not atom.IsHydrogen():
					for atomlig in OBMolAtomIter(ligand):
						if not atomlig.IsHydrogen():
							distance = atom.GetDistance(atomlig)
							if distance <= 3.5:
								if atom.IsHbondDonor() & atomlig.IsHbondAcceptor():
									for neighborDon in OBAtomAtomIter(atom):
										if neighborDon.IsHydrogen():
											angle = atom.GetAngle(neighborDon, atomlig)
											if angle>135.0:
												resdict[residuename] |= bitarray('0001000')
								if atom.IsHbondAcceptor() & atomlig.IsHbondDonor():
									for neighborDon in OBAtomAtomIter(atomlig):
										if neighborDon.IsHydrogen():
											angle = atomlig.GetAngle(neighborDon, atom)
											if angle>135.0:
												resdict[residuename] |= bitarray('0000100')


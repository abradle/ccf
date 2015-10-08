from IOhandle.models import Compound,Molecule
from rdkit import Chem
from rdkit.Chem import Lipinski
import sys
from rdkit.Chem import Descriptors,AllChem
from django.core.exceptions import ValidationError
import os
import sys
import tempfile
import subprocess
from django.db import connections


# Contribution to the RDKit from Hans de Winter
def _InitialiseNeutralisationReactions():
    """Contribution from Hans de Winter"""
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None


def NeutraliseCharges(smiles, reactions=None):
    """Contribution from Hans de Winter"""
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions = _InitialiseNeutralisationReactions()
        reactions = _reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i, (reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return (Chem.MolToSmiles(mol, True), True)
    else:
        return (smiles, False)


# OOMMPPAA functions
def set_pH(in_smi):
    """Function to set the pH of a molecule
    Currently uses OpenBabel to perform protonation (to pH 7.4)
    Takes a smiles string
    Returns the protonated smiles string"""
    # Attempt to use the babel that has been  included in the distro
    try:
        d = sys._MEIPASS
        babel_path = os.path.join(d, "babel")
    # If not use the linux babel or the windows babel in defined locations
    except AttributeError:
        if "linux" in sys.platform:
            babel_path = "/usr/bin/babel"
        # Check for mac babel
        elif "darwin" in sys.platform:
            babel_path = "/usr/local/bin/babel"
            if os.path.isfile(babel_path):
                pass
            else:
                print "PY2APP FIX THIS"
                # Find the babel path
        else:
            sys.stderr.write("USING SYSTEM BABEL")
            babel_path = r"C:\Program Files (x86)\OpenBabel-2.3.2\babel.exe"

    in_f = tempfile.NamedTemporaryFile("w", delete=False)
    in_f.write(in_smi)
    in_f.close()
    out_f = tempfile.NamedTemporaryFile("w", delete=False)
    out_f.close()
    my_list = [babel_path, "-ismi", in_f.name,  "-p", "7.4", "-osmi", out_f.name]
    subprocess.call(my_list, stderr=tempfile.NamedTemporaryFile())
    out_smi = open(out_f.name).read().rstrip()
    return out_smi


# OOMMPPAA functions
def set_pH_sdf(in_sdf):
    """Function to set the pH of a molecule
    Currently uses OpenBabel to perform protonation (to pH 7.4)
    Takes a smiles string
    Returns the protonated smiles string"""
    # Attempt to use the babel that has been  included in the distro
    try:
        d = sys._MEIPASS
        babel_path = os.path.join(d, "babel")
    # If not use the linux babel or the windows babel in defined locations
    except AttributeError:
        if "linux" in sys.platform:
            babel_path = "/usr/bin/babel"
        # Check for mac babel
        elif "darwin" in sys.platform:
            babel_path = "/usr/local/bin/babel"
            if os.path.isfile(babel_path):
                pass
            else:
                print "PY2APP FIX THIS"
                # Find the babel path
        else:
            sys.stderr.write("USING SYSTEM BABEL")
            babel_path = r"C:\Program Files (x86)\OpenBabel-2.3.2\babel.exe"
    in_f = tempfile.NamedTemporaryFile("w", delete=False)
    in_f.write(in_sdf)
    in_f.close()
    out_f = tempfile.NamedTemporaryFile("w", delete=False)
    out_f.close()
    my_list = [babel_path, "-isdf", in_f.name, "-osdf", out_f.name, "-p", "7.4"]
    subprocess.call(my_list, stderr=tempfile.NamedTemporaryFile())
    out_smi = open(out_f.name).read().rstrip()
    return out_smi


def add_h_pdb(in_pdb):
    """Function to set the pH of a molecule
    Currently uses OpenBabel to perform protonation (to pH 7.4)
    Takes a smiles string
    Returns the protonated smiles string"""
    # Attempt to use the babel that has been  included in the distro
    try:
        d = sys._MEIPASS
        babel_path = os.path.join(d, "babel")
    # If not use the linux babel or the windows babel in defined locations
    except AttributeError:
        if "linux" in sys.platform:
            babel_path = "/usr/bin/babel"
        # Check for mac babel
        elif "darwin" in sys.platform:
            babel_path = "/usr/local/bin/babel"
            if os.path.isfile(babel_path):
                pass
            else:
                print "PY2APP FIX THIS"
                # Find the babel path
        else:
            sys.stderr.write("USING SYSTEM BABEL")
            babel_path = r"C:\Program Files (x86)\OpenBabel-2.3.2\babel.exe"

    in_f = tempfile.NamedTemporaryFile("w", delete=False)
    in_f.write(in_pdb)
    in_f.close()
    out_f = tempfile.NamedTemporaryFile("w", delete=False)
    out_f.close()
    my_list = [babel_path, "-ipdb", in_f.name, "-omol2", out_f.name, "-p", "7.4", "-h"]
    subprocess.call(my_list, stderr=tempfile.NamedTemporaryFile())
    out_mol2 = open(out_f.name).read()
    return out_mol2


def desalt_compound(smiles):
    """Function to desalt compound a given smiles string
    Takes a smiles string
    Returns a desalted smiles string."""
    # Chose the biggest fragment, after splitting into fragments
    return sorted([(x, Lipinski.HeavyAtomCount(Chem.MolFromSmiles(x))) for x in smiles.split(".")], key=lambda x: x[1], reverse=True)[0][0]


def add_new_comp(mol, option=None, comp_id=None):
    """Function to add a new compound to the database given an RDKit molecule
    Takes an RDKit molecule. Option of LIG to return original smiles with the Compound object
    Returns a compound object for the RDKit molecule."""
    # Neutralise and desalt compound the compound
    s_store_mol = NeutraliseCharges(desalt_compound(Chem.MolToSmiles(mol, isomericSmiles=True)))[0]
    store_mol = Chem.MolFromSmiles(s_store_mol)
    if store_mol is None:
        sys.stderr.write("NEUTRALISING MADE NONE MOL" + " " + s_store_mol + " " + Chem.MolToSmiles(mol, isomericSmiles=True))
        return None
    # Store the isomeric smiles
    smiles = Chem.MolToSmiles(store_mol, isomericSmiles=True)
    # The inchi string is used for unique identification
    inchi = Chem.MolToInchi(store_mol)
    # Now convert back to inchi to canonicalise
    tmp_mol = Chem.MolFromInchi(inchi)
    if tmp_mol is None:
        # If error in INNCHI READ -> NOT NECCESARILY A KILLER
        sys.stderr.write("INCHI ERROR: " + inchi)
    else:
        inchi = Chem.MolToInchi(tmp_mol)
# Now attribute all this meta-deta to the compound object
    new_comp = Compound()
    new_comp.smiles = smiles
    if len(smiles) > Compound._meta.get_field('smiles').max_length:
        print "SMILES TOO LONG"
        return None
    new_comp.inchi = inchi
    m = store_mol
    try:
        new_comp.mol_log_p = Chem.Crippen.MolLogP(m)
    except:
        sys.stderr.write("NONE MOLECULE PRODUCED\n" + smiles + "\n" + inchi)
        return None
    new_comp.mol_wt = float(Chem.rdMolDescriptors.CalcExactMolWt(m))
    new_comp.heavy_atom_count = Chem.Lipinski.HeavyAtomCount(m)
    new_comp.heavy_atom_mol_wt = float(Descriptors.HeavyAtomMolWt(m))
    new_comp.nhoh_count = Chem.Lipinski.NHOHCount(m)
    new_comp.no_count = Chem.Lipinski.NOCount(m)
    new_comp.num_h_acceptors = Chem.Lipinski.NumHAcceptors(m)
    new_comp.num_h_donors = Chem.Lipinski.NumHDonors(m)
    new_comp.num_het_atoms = Chem.Lipinski.NumHeteroatoms(m)
    new_comp.num_rot_bonds = Chem.Lipinski.NumRotatableBonds(m)
    new_comp.num_val_electrons = Descriptors.NumValenceElectrons(m)
    new_comp.ring_count = Chem.Lipinski.RingCount(m)
    new_comp.tpsa = Chem.rdMolDescriptors.CalcTPSA(m)
    # Validate that the compound is unique
    try:
        new_comp.validate_unique()
        new_comp.save()
        if option is None:
            return new_comp
        elif option == "LIG":
            return (new_comp, smiles)
    except ValidationError:
        if option is None:
            return Compound.objects.get(inchi=inchi)
        elif option == "LIG":
            return (Compound.objects.get(inchi=inchi), smiles)


def centre_of_mass_from_SD_block(sdf_text):
    """Function to return the unweighted Centre of Mass from an SD block.
    Takes an SD block
    Returns three floats, the x,y and z coordinate for the centre of mass"""
    # Convert it to an RDMol
    rdmol = Chem.MolFromMolBlock(str(sdf_text))
    # Gives the atoms
    atoms = rdmol.GetAtoms()
    conf = rdmol.GetConformer()
    x_coord = y_coord = z_coord = 0.0
    numatoms = 0.0
    # Assume all heavy atoms have the same mass
    for atom in atoms:
        if atom.GetAtomicNum() == 1 or atom.GetSmarts() == "[*]":
            continue
        numatoms += 1.0
        coords = conf.GetAtomPosition(atom.GetIdx())
        if coords.x == 0.0 and coords.y == 0.0 and coords.z == 0.0:
            print sdf_text
            sys.exit()
        x_coord += float(coords.x)
        y_coord += float(coords.y)
        z_coord += float(coords.z)
    # Now we have all the coords -> we want to loop through
    if numatoms == 0:
        sys.stderr.write("Zero atoms error...")
        return None, None, None
    return x_coord / numatoms, y_coord / numatoms, z_coord / numatoms


def make_smarts_from_frag(smi):
    """Function to convert a fragmented smiles string to a valid smarts
    Takes a fragmented smiles
    Returns a smarts string."""
    # Make the RDMol
    my_mol = Chem.MolFromSmiles(smi)
    if my_mol is None:
        print "NONE MOL:\n", smi
    # Get the atoms
    atms = Chem.Mol.GetAtoms(my_mol)
    # For each one set the isotope to one -> MARKS IT AS SPECIAL
    for a in atms:
        a.SetIsotope(1)
    # Now create the smiles string -> ISOMERICALLY. FOR WHATEVER REASON THIS MAKES EXPLICIT H's
    new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(my_mol, isomericSmiles=True))
    # REPLACE THE ATACCHEMENT POINTS WITH BLANK SPAEC
    return Chem.MolToSmarts(new_mol).replace("(-[*])", "").replace("[*]-", "").replace("-[*]", "")

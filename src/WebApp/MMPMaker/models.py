from django.db import models
from rdkit import Chem
from rdkit.Chem import AllChem,Lipinski,Descriptors
from IOhandle.models import Molecule,ActivityPoint,Target,Compound
from Pharmacophore.models import PharmaPoint
from django.core.exceptions import ValidationError
from MMPMaker.helpers import *


class MMP(models.Model):
    """Django model to store the "context". The key of the data model"""
    # Canonical smiles for tracking and tings
    context = models.CharField(max_length=500, unique=True, db_index=True)
    # Links back to all the compound's that own it
    compound_id = models.ManyToManyField(Compound)


class MMPFrag(models.Model):
    """Django model to store the fragments. The value of the data model"""
    # The smiles that is stored
    smistore = models.CharField(max_length=500, db_index=True)
    # The 3D information for the fragment if it exists
    sdf_info = models.TextField()
    # The molecule it is associated to
    mol_id = models.ForeignKey(Molecule, null=True, db_index=True)
    # Or if it's an activity point the point it is associated to
    act_id = models.ForeignKey(ActivityPoint, null=True, db_index=True)
    # Link back to the appropriate context
    mmp_link = models.ForeignKey(MMP)
    # To give the core size
    core_size = models.FloatField(db_index=True)
    # To give the core to MMP ratio
    core_ratio = models.FloatField(db_index=True)
    # Define the centre of mass
    x_com = models.FloatField(null=True, db_index=True)
    y_com = models.FloatField(null=True, db_index=True)
    z_com = models.FloatField(null=True, db_index=True)

    class Meta:
        unique_together = ('sdf_info', 'mmp_link', )


class Act3DFrag(models.Model):
    """Django model of 3D fragments found by aligning ActivityPoint
    to Molecule"""
    # Link to original fragment
    mmp_frag_id = models.ForeignKey(MMPFrag)
    # Link to overlapped molecule
    over_mol_id = models.ForeignKey(Molecule, related_name="%(app_label)s_%(class)s_overlapmol")
    # Link to the molecule it is part of
    mol_id = models.ForeignKey(Molecule)
    # Link to the activity point
    act_id = models.ForeignKey(ActivityPoint)
    # Own sdf_info
    sdf_info = models.TextField()
    # Heavy atom count
    heavy_atom_count = models.IntegerField(null=True)
    # Molecular weight
    mol_wt = models.FloatField(null=True)

    class Meta:
        unique_together = ('over_mol_id', 'mmp_frag_id', 'act_id', )


class ActPharmaPoint(models.Model):
    """Django model of pharmacophore differences between molecules.
    Relates them to their activity change and their respective differences."""
    # Direct link to the target for this
    target_id = models.ForeignKey(Target)
    # Define all the pharmacophore related stuff
    pharma_id = models.ForeignKey(PharmaPoint)
    # Distance from nearest identical point
    dist_same = models.FloatField(null=True)
    # The Act3DFrag it is associated to, gives it's activity and the other activity and what not
    act_3d_frag_id = models.ForeignKey(Act3DFrag)
    # Whether it is xtal or ChEMBL activity
    from_xtal = models.BooleanField()
    # type of comparison, xtal -> xtal or xtal -> activity
    type = models.CharField(max_length=50, null=True)
    # Number using 3.5
    num_diff_10 = models.IntegerField(null=True, db_index=True)
    # Is it in this point?
    in_diff_10 = models.BooleanField(default=False, db_index=True)
    # Number using 3.5
    num_diff_15 = models.IntegerField(null=True, db_index=True)
    # Is it in this point?
    in_diff_15 = models.BooleanField(default=False, db_index=True)
    # Number of points using 2.5,
    num_diff_25 = models.IntegerField(null=True, db_index=True)
    # Is it in this point?
    in_diff_25 = models.BooleanField(default=False, db_index=True)
    # Number using 3.0
    num_diff_30 = models.IntegerField(null=True, db_index=True)
    # Is it in this point?
    in_diff_30 = models.BooleanField(default=False, db_index=True)
    # Number using 3.5
    num_diff_35 = models.IntegerField(null=True, db_index=True)
    # Is it in this point?
    in_diff_35 = models.BooleanField(default=False, db_index=True)
    # Log difference in activity
    act_change = models.FloatField(db_index=True)

    class Meta:
        unique_together = ('pharma_id', 'act_3d_frag_id', )


class MMPComparison(models.Model):
    """Django model of comparisons between two molecules in 3D space"""
    sdf_info = models.TextField()
    int_id = models.IntegerField(null=True)
    title = models.TextField()
    # Is it an act or inact or non point??
    type = models.TextField()
    # Two foreign key relations to the molecules it relates to
    xtal_mol = models.ForeignKey(Molecule)
    xtal_act = models.ForeignKey(ActivityPoint)
    chembl_mol = models.ForeignKey(Molecule, related_name="%(app_label)s_%(class)s_chemblmol")
    chembl_act = models.ForeignKey(ActivityPoint, related_name="%(app_label)s_%(class)s_chemblact")
    target_id = models.ForeignKey(Target)
    num_change = models.IntegerField(null=True)
    activity_change = models.FloatField(null=True)
    activity_type = models.TextField()
    # A text field to hold the list of things it is linked to
    map_inds = models.TextField()

    class Meta:
        unique_together = ('target_id', 'activity_type', 'xtal_mol', 'chembl_mol', 'type',)


class MMPDiffMap(models.Model):
    """Django model to hold information about differences for a target.
    e.g. Activity improving pharmacophore points will be held in one map for one target."""
    # Hold the PDB text
    pdb_info = models.TextField(null=True)
    # The type of map (shape,act and inact)
    type = models.TextField()
    # The Target used
    target_id = models.ForeignKey(Target)
    # The number of changes
    num_change = models.IntegerField()
    # The activity change
    activity_change = models.FloatField()
    # The type of activity change
    activity_type = models.TextField()

    class Meta:
        unique_together = ('target_id', 'activity_type', 'type', )


class ActMapPoint(models.Model):
    """Django model to relate points in a map to their respective information.
    Used to relate points to compound information in the viewer."""
    # Foreign keys to relate to toher Information
    map_id = models.ForeignKey(MMPDiffMap)
    mmp_comparsion_id = models.ForeignKey(MMPComparison)
    pharma_id = models.ForeignKey(PharmaPoint)
    # Point to make unique
    point_id = models.IntegerField(db_index=True)

    class Meta:
        unique_together = ('map_id', 'point_id')


def make_new_3d_frag(mmpfrag, act, overmol, rdmol, prot, new_rdmol):
    """Function to make an Act3DFrag object
    Takes an MMPFrag, the activity point and the molecule it relates to for the overlay,
    a protein and the new coordinates of its overlayed molecule.
    Returns None"""
    # So now make this new molecule which is associated to the MMPProtein
    # and the activity point associated to it
    new_mol = Molecule()
    new_mol.prot_id = prot
    new_mol.sdf_info = Chem.MolToMolBlock(new_rdmol)
    new_mol.cmpd_id = act.cmpd_id
    new_mol.rscc = 0.0
    new_mol.lig_id = ""
    new_mol.chain_id = ""
    new_mol.occupancy = 0.0
    try:
        new_mol.validate_unique()
        new_mol.save()
    except ValidationError:
        new_mol = Molecule.objects.get(prot_id=prot, cmpd_id=act.cmpd_id, sdf_info=Chem.MolToMolBlock(new_rdmol))
    new_frag = Act3DFrag()
    new_frag.mol_id = new_mol
    new_frag.act_id = act
    new_frag.mmp_frag_id = mmpfrag
    new_frag.sdf_info = Chem.MolToMolBlock(rdmol)
    new_frag.over_mol_id = overmol
    new_frag.mol_wt = Descriptors.MolWt(rdmol)
    new_frag.heavy_atom_count = Lipinski.HeavyAtomCount(rdmol)
    # Make non fragment
    em = AllChem.EditableMol(rdmol)
    i = 0
    for atom in rdmol.GetAtoms():
        if "*" in atom.GetSymbol():
            em.RemoveAtom(atom.GetIdx() - i)
            i += 1
    out_mol = em.GetMol()
    Chem.SanitizeMol(out_mol)
    try:
        new_frag.validate_unique()
        new_frag.save()
    except ValidationError:
        pass


def make_new_mmp(context, my_id, core, rdmol, attachments, smi, num):
    """Function to make a new MMP and MMPFrag object
    Takes context - a string of the index,
    my_id - an id for the activity point or molecule or compound,
    core - a smiles string of the fragment,
    rdmol - an RDKit molecule of the molecule,
    attachments - an integer of the number of attachments,
    smi - a smiles string of the whole compound,
    num - an int denoting which fragment this one is
    Returns None"""
    # If this is just a compound then return the id
    # Make the new index
    my_type = my_id[0:3]
    my_id = my_id[3:]
    # Now make the MMP frag and associate it with this
    new_frag = MMPFrag()
    new_frag.smistore = core.replace("1*", "*:1")
    # Now associate all of the extra informations
    if  my_type == "mol":
        my_mol = Molecule.objects.get(pk=my_id)
        new_frag.mol_id = my_mol
        new_frag.act_id = ActivityPoint.objects.get(units="DUMMY")
        heavy_weight = my_mol.cmpd_id.heavy_atom_mol_wt
    elif my_type == "act":
        my_mol = ActivityPoint.objects.get(pk=my_id)
        new_frag.act_id = my_mol
        rdmol.SetProp("_Name", str(my_mol.target_id.title) +"_" + str(my_mol.pk) + "_" + str(num))
        new_frag.mol_id = Molecule.objects.get(lig_id="DUMMY")
        heavy_weight = my_mol.cmpd_id.heavy_atom_mol_wt
    elif my_type == "cmp":
        my_mol = Compound.objects.get(pk=my_id)
        heavy_weight = my_mol.heavy_atom_mol_wt
    new_frag.core_size = Chem.Mol.GetNumAtoms(Chem.MolFromSmiles(smi)) - attachments
    new_frag.core_ratio = float(new_frag.core_size) / float(heavy_weight)
    if type == "mol":
        new_frag.act_id = None
    elif type == "act":
        new_frag.mol_id = None
    new_frag.sdf_info = Chem.MolToMolBlock(rdmol)
    return context.replace("1*", "*:1"), new_frag


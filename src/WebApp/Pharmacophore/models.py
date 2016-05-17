
from django.db import models
from django.core.exceptions import ValidationError
from rdkit import Chem
from IOhandle.models import Molecule,Protein,Project,Target


class UniquePharma(models.Model):
    """Django model to hold information about unique Pharmacophore types.
    e.g. H-bond donor"""
    # They type
    smiles = models.CharField(max_length=500,unique=True,db_index=True)
    num_added = models.IntegerField(null=True)
    perc_targ = models.FloatField(null=True)
    perc_tot = models.FloatField(null=True)


class PharmaPoint(models.Model):
    """Django model to hold information about each pharmacophore point"""
    # The type
    smiles = models.CharField(max_length=500, db_index=True)
    init_date = models.DateTimeField(auto_now_add=True)
    # The method is the unique thing here
    method_used = models.CharField(max_length=500)
    git_version = models.CharField(max_length=100)
    # The molecule it relates to
    mol_id = models.ForeignKey(Molecule, null=True)
    uniq_id = models.ForeignKey(UniquePharma)
    pharma_id = models.IntegerField()
    # Value to demonstrate the uncertainty associated to a blob -> 0.0 is default
    # Change this to be correctly spelt at the next go round
    uncertainity = models.FloatField(default=0.0)
    std_x = models.FloatField(default=0.0)
    std_y = models.FloatField(default=0.0)
    std_z = models.FloatField(default=0.0)
    # Value to demonstrate the weight of a feature
    Weight = models.FloatField(default=1.0)
    x_com = models.FloatField()
    y_com = models.FloatField()
    z_com = models.FloatField()

    class Meta:
        unique_together = ('pharma_id', 'mol_id', 'method_used', )


def make_unique_pharma(smiles):
    """Function to make or update a pharmacophore type
    Takes a pharmacophore type as string
    Returns a UniquePharma object"""
    new_uni_frag = UniquePharma()
    new_uni_frag.smiles = smiles
    try:
        new_uni_frag.validate_unique()
        new_uni_frag.num_added = 1
        new_uni_frag.save()
        return new_uni_frag
    except ValidationError:
        newer_frag = UniquePharma.objects.get(smiles=smiles)
        newer_frag.num_added += 1
        newer_frag.save(update_fields=["num_added"])
        return UniquePharma.objects.get(smiles=smiles)


def  find_pharmacophore_points(mol, factory, fdefname):
    """Function to find pharmacophore points for a molecule and store them.
    Takes an RDKit molecule, an RDKit feature factory and an fdefname
    Return None"""
    rdmol = Chem.MolFromMolBlock(str(mol.sdf_info))
    feats = factory.GetFeaturesForMol(rdmol)
    counter = 0
    for feat in feats:
        new_fragment = PharmaPoint()
        new_fragment.method_used = fdefname
        new_fragment.smiles = str(feat.GetType())
        coord = feat.GetPos()
        # Assogm the co-ordinates
        new_fragment.x_com = coord.x
        new_fragment.y_com = coord.y
        new_fragment.z_com = coord.z
        new_fragment.mol_id = mol
        uni_frag = make_unique_pharma(new_fragment.smiles)
        new_fragment.uniq_id = uni_frag
        new_fragment.pharma_id = counter
        new_fragment.git_version = "GIT"
        counter += 1
        try:
            new_fragment.validate_unique()
            new_fragment.save()
        except ValidationError:
            print str(feat.GetType())
            pass

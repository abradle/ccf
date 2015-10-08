from django.db import models
from django.core.exceptions import ValidationError
from IOhandle.models import Molecule, ActivityPoint, Target


class Transformation(models.Model):
    """Django model of comparisons between two molecules in 3D space"""
    # The SMARTS of the transformation
    trans_smirks = models.TextField()
    # The before frag smiles
    before_smiles = models.TextField()
    # The before frag SDF
    before_sdf = models.TextField()
    # The after frag smiles
    after_smiles = models.TextField()
    # The after SDF
    after_sdf = models.TextField()
    # The centre of mass of the source
    x_com = models.FloatField()
    y_com = models.FloatField()
    z_com = models.FloatField()
    # The centre of mass of the dest
    x_com_dest = models.FloatField()
    y_com_dest = models.FloatField()
    z_com_dest = models.FloatField()
    # Two foreign key relations to the molecules it relates to
    mol_one = models.ForeignKey(Molecule, related_name="%(app_label)s_%(class)s_mol_one")
    act_one = models.ForeignKey(ActivityPoint, related_name="%(app_label)s_%(class)s_act_one")
    mol_two = models.ForeignKey(Molecule, related_name="%(app_label)s_%(class)s_mol_two")
    act_two = models.ForeignKey(ActivityPoint, related_name="%(app_label)s_%(class)s_act_two")
    target_id = models.ForeignKey(Target)
    # The ratio
    ha_ratio = models.FloatField(null=True)
    # The difference in Heavy atoms
    ha_diff = models.IntegerField(null=True)
    # The difference in molwt
    mol_wt_diff = models.FloatField(null=True)
    # The difference in Heavy atoms
    ha_change = models.IntegerField(null=True)
    # The difference in molwt
    mol_wt_change = models.FloatField(null=True)
    # The change in activity
    activity_change = models.FloatField(null=True)
    # The type of activity change
    activity_type = models.TextField()
    # The type - "ACTIVE__INACTIVE" etc
    change_type = models.TextField()
    # How the molecules were derived - XTAL_XTAL, ACT_ACT, XTAL_ACT, ACT_XTAL all valid
    mols_derived = models.TextField()
    # The Shape distance between the two -> this can then be correlated with the activtiy distance
    shape_dist = models.FloatField()

    class Meta:
        unique_together = ('target_id', 'activity_type', 'mol_one', 'mol_two',)


class TransGroups(models.Model):
    """TRANSFORMATION GROUPS -> AS IN PARTICULAR DATA ONE WANTS TO COMPARE.... (E.G KINETICS)"""
    # Link to the trans
    trans_id = models.ManyToManyField(Transformation)
    # PDB info to represent
    pdb_info = models.TextField()
    # Type of data used 
    data_type = models.TextField()
    # Target it relates to
    target_id = models.ForeignKey(Target)

    class Meta:
        unique_together = ('data_type', 'target_id')


class TransCluster(models.Model):
    """CONFLICTING SAR AND BIG WINS AND LOSSES"""
    # The CofM
    x_com = models.FloatField()
    y_com = models.FloatField()
    z_com = models.FloatField()
    # The mean, the min and the max -> from which we can find the conflicts
    min = models.FloatField()
    max = models.FloatField()
    diff = models.FloatField()
    mean = models.FloatField()
    # Type of data used 
    data_type = models.TextField()
    # The SMILES or name of pharmacophore
    trans_type = models.TextField()

    class Meta:
        unique_together = ('data_type', 'trans_type', 'x_com', 'y_com', 'z_com')

# NON ADDITIVTY -> a+b > a and a+c > a but a+b+c < a HOW TO CODE THIs????

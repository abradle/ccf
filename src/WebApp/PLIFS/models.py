from django.db import models
from IOhandle.models import Target, Molecule, Protein
from WONKA.models import Residue


class PlifProtein(models.Model):
    """Model to add Mol2 version of a protein"""
    # The protein it relates to
    prot_id = models.ForeignKey(Protein)
    # The data that is held
    mol2_data = models.TextField()


class PlifProbe(models.Model):
    """Model to find molecueles that displace conserved waters"""
    # The target it relates to
    target_id = models.ForeignKey(Target)
    # The molecule and protein foreign keys
    mol_id = models.ForeignKey(Molecule, null=True)
    prot_id = models.ForeignKey(Protein, null=True)
    # Need this to be linked to the residue or molecule
    res_id = models.ForeignKey(Residue, null=True)
    # Coordinates
    x_com = models.FloatField(db_index=True)
    y_com = models.FloatField(db_index=True)
    z_com = models.FloatField(db_index=True)
    # Type of probe
    type = models.CharField(max_length=30, db_index=True)
    # Extra information to derive planes between features
    x_extent = models.FloatField(null=True)
    y_extent = models.FloatField(null=True)
    z_extent = models.FloatField(null=True)
    x_extent_2 = models.FloatField(null=True)
    y_extent_2 = models.FloatField(null=True)
    z_extent_2 = models.FloatField(null=True)
    x_extent_3 = models.FloatField(null=True)
    y_extent_3 = models.FloatField(null=True)
    z_extent_3 = models.FloatField(null=True)
    unique_id = models.CharField(max_length=100)
    # Atom ids - the atom ids for the molecule it relates to
    atom_ids = models.CharField(max_length=100, null=True)

    class Meta:
        unique_together = ('x_com', 'y_com', 'z_com', 'unique_id', 'type',)


class PlifProbeBit(models.Model):
    """Model to store a protein-ligand interaction"""
    # The target
    target_id = models.ForeignKey(Target)
    # The molecule it relates to
    mol_id = models.ForeignKey(Molecule)
    # The protein it relates to
    prot_id = models.ForeignKey(Protein)
    # The source probe
    probe_source_id = models.ForeignKey(PlifProbe, related_name="%(app_label)s_%(class)s_source")
    # The destination probe
    probe_dest_id = models.ForeignKey(PlifProbe, related_name="%(app_label)s_%(class)s_dest")
    # The interaction type
    type = models.CharField(max_length=300)
    # The interaction distance
    dist = models.FloatField()
    # The interaction angle
    angle_1 = models.FloatField(null=True)
    # The interaction angle
    angle_2 = models.FloatField(null=True)
    # The third interaction angle
    angle_3 = models.FloatField(null=True)
    # The definition of the interaction, e.g. edge-to-face
    definition = models.CharField(max_length=300, null=True)

    class Meta:
        unique_together = ('definition', 'probe_source_id', 'probe_dest_id', )


class InteractionScheme(models.Model):
    """Model to store a scheme of interactions"""
    # The name of the scheme
    scheme_name = models.CharField(max_length=300, unique=True)
    # The dict of the scheme
    scheme_dict = models.TextField()


class Interaction(models.Model):
    """Model to store a particular interaction"""
    # The link to the scheme
    scheme_id = models.ForeignKey(InteractionScheme)
    # The link to the interaction
    interaction_id = models.ForeignKey(PlifProbeBit)
    # The target
    target_id = models.ForeignKey(Target)
    # The molecule it relates to
    mol_id = models.ForeignKey(Molecule)
    # The protein it relates to
    prot_id = models.ForeignKey(Protein)

    class Meta:
        unique_together = ('scheme_id', 'interaction_id',)

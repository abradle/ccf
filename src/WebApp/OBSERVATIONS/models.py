from django.db import models
from Pharmacophore.models import PharmaPoint
from MMPMaker.models import MMPFrag, ActPharmaPoint, Act3DFrag
from IOhandle.models import Target, Protein, Molecule, Project, Compound
from django.contrib.auth.models import User
from WONKA.models import Residue, KeyCluster


class UserData(models.Model):
    """Model to link the current state to a user"""
    # The targets this user likes
    cool_targs = models.ManyToManyField(Target, related_name="%(app_label)s_%(class)s_cool")
    # The target this user hasn't  seen
    new_targs = models.ManyToManyField(Target, related_name="%(app_label)s_%(class)s_new")
    # The compounds this user is watching
    cool_cmps = models.ManyToManyField(Compound, related_name="%(app_label)s_%(class)s_cool")
    # The projects this user likes
    project_id = models.ManyToManyField(Project, related_name="%(app_label)s_%(class)s_proj")
    # The user this user is
    user_id = models.ForeignKey(User, unique=True, related_name="%(app_label)s_%(class)s_user")
    # The new mols
    new_mols = models.ManyToManyField(Molecule, related_name="%(app_label)s_%(class)s_new")
    # The cool mols
    cool_mols = models.ManyToManyField(Molecule, related_name="%(app_label)s_%(class)s_cool")
    # Registered


class Observation(models.Model):
    """A model to describe an observation"""
    # Comments on the observation
    comments = models.TextField()
    # The HTML of the page
    html_text = models.TextField()
    # The ICM blob of the page
    icm_blob = models.TextField()
    # The image of the page
    image_page = models.TextField()
    # The author of the work - at the moment just text - could be linked to users
    author = models.TextField()
    # The user it is linked to
    user_id = models.ForeignKey(User, db_index=False, related_name="%(app_label)s_%(class)s_user")
    # The number of likes it has
    likes = models.IntegerField()
    # Its Target
    target_id = models.ForeignKey(Target, related_name="%(app_label)s_%(class)s_targ")
    # Anything else associated to it?????
    uuid = models.TextField()


class ObservationChoices(models.Model):
    """A model to append information regarding to an observation"""
    # First the Observation it relates to
    obs_id = models.ForeignKey(Observation)
    # Next the molecules it relates to
    mol_id = models.ManyToManyField(Molecule)
    # Next the compounds it relates to
    cmpd_id = models.ManyToManyField(Compound)
    # Next the water / ph4 clusters it relates to
    kc_id = models.ManyToManyField(KeyCluster)
    # Next the residue clusters it relates to
    res_id = models.ManyToManyField(Residue)
    # Finally the proteins it relates to
    prot_id = models.ManyToManyField(Protein)


class Comment(models.Model):
    """A model consisting of a comment on an observation"""
    # The comment
    comment = models.TextField()
    # The user
    user_id = models.ForeignKey(User, related_name="%(app_label)s_%(class)s_user")
    # The observation
    obs_id = models.ForeignKey(Observation)

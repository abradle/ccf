from django.db import models
from PLIFS.models import PlifProbeBit, PlifProbe
from IOhandle.models import Target


class PlifVis(models.Model):
    """Model to hold the JSON for a PLIF string"""
    # The target it relates to
    target_id = models.ForeignKey(Target, unique=True)
    # The JSON of the vals
    json_text = models.TextField()


class PlifVisGrid(models.Model):
    """Model to hold the JSON for a PLIF string"""
    # The target it relates to
    target_id = models.ForeignKey(Target)
    # The JSON of the vals
    json_text = models.TextField()
    # The spacings of the grid
    grid_space = models.FloatField()

    class Meta:
        unique_together = ('grid_space', 'target_id', )


class PlifProbeScore(models.Model):
    """Model to hold a score for a PLIF probe"""
    # the score
    score = models.FloatField()
    # The item it links to
    plif_probe = models.ForeignKey(PlifProbe, unique=True)


class PlifProbeGridScoreNew(models.Model):
    """Model to hold a score for a PLIF probe - with different grid spacing"""
    # the score
    score = models.FloatField()
    # The item it links to
    plif_probe = models.ForeignKey(PlifProbe)
    # The grid spacing
    grid_space = models.FloatField()

    class Meta:
        unique_together = ('grid_space', 'plif_probe', )
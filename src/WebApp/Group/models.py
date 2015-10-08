from django.db import models
from django.core.exceptions import ValidationError


class Group(models.Model):
    """Django model to hold information about groups"""
# This is how the group was constructed
    method_used = models.CharField(max_length=50, db_index=True)
# The version of git checkout used
    git_version = models.CharField(max_length=100)
# The number of the group
    group_number = models.IntegerField()
# The special identifier
    group_text = models.CharField(max_length=100, db_index=True)
# A group should be a unique combination of these

    class Meta:
        unique_together = ('method_used', 'group_number', 'group_text', )

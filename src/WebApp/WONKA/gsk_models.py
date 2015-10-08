# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#     * Rearrange models' order
#     * Make sure each model has one field with primary_key=True
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin.py sqlcustom [appname]'
# into your database.
from __future__ import unicode_literals

from django.db import models


class Jcstruct(models.Model):
    cd_id = models.IntegerField(primary_key=True)
    cd_structure = models.TextField()
    cd_smiles = models.TextField(blank=True)
    cd_formula = models.CharField(max_length=100L, blank=True)
    cd_sortable_formula = models.CharField(max_length=255L, blank=True)
    cd_molweight = models.FloatField(null=True, blank=True)
    cd_hash = models.IntegerField()
    cd_flags = models.CharField(max_length=20L, blank=True)
    cd_timestamp = models.DateTimeField()
    cd_pre_calculated = models.IntegerField()
    cd_fp1 = models.IntegerField()
    cd_fp2 = models.IntegerField()
    cd_fp3 = models.IntegerField()
    cd_fp4 = models.IntegerField()
    cd_fp5 = models.IntegerField()
    cd_fp6 = models.IntegerField()
    cd_fp7 = models.IntegerField()
    cd_fp8 = models.IntegerField()
    cd_fp9 = models.IntegerField()
    cd_fp10 = models.IntegerField()
    cd_fp11 = models.IntegerField()
    cd_fp12 = models.IntegerField()
    cd_fp13 = models.IntegerField()
    cd_fp14 = models.IntegerField()
    cd_fp15 = models.IntegerField()
    cd_fp16 = models.IntegerField()

    class Meta:
        db_table = 'JCStruct'


class JcstructUl(models.Model):
    update_id = models.IntegerField(primary_key=True)
    update_info = models.CharField(max_length=120L)
    cache_id = models.CharField(max_length=32L)

    class Meta:
        db_table = 'JCStruct_UL'


class Jchemproperties(models.Model):
    prop_name = models.CharField(max_length=200L, primary_key=True)
    prop_value = models.CharField(max_length=200L, blank=True)
    prop_value_ext = models.TextField(blank=True)

    class Meta:
        db_table = 'JChemProperties'


class JchempropertiesCr(models.Model):
    cache_id = models.CharField(max_length=32L, primary_key=True)
    registration_time = models.CharField(max_length=30L)
    is_protected = models.IntegerField()
 
    class Meta:
        db_table = 'JChemProperties_CR'


class Activity(models.Model):
    structure_id = models.IntegerField(primary_key=True)
    projectname = models.CharField(max_length=30L)
    ligandid = models.CharField(max_length=10L, blank=True)
    assayname = models.CharField(max_length=30L, blank=True)
    assaysource = models.CharField(max_length=10L)
    assaydata = models.CharField(max_length=30L, blank=True)

    class Meta:
        db_table = 'activity'


class Editors(models.Model):
    projectname = models.CharField(max_length=30L)
    editor = models.CharField(primary_key=True,max_length=30L)

    class Meta:
        db_table = 'editors'


class Junk(models.Model):
    junkval = models.CharField(max_length=80L, blank=True)
    id = models.IntegerField(primary_key=True, blank=True)

    class Meta:
        db_table = 'junk'


class Linkers(models.Model):
    projectname = models.CharField(max_length=30L)
    hierarchy = models.CharField(max_length=30L)
    parent = models.CharField(max_length=30L)
    child = models.CharField(max_length=30L)

    class Meta:
        db_table = 'linkers'


class Projects(models.Model):
    projectname = models.CharField(primary_key=True, max_length=30L)
    owner = models.CharField(max_length=30L)
    projectdisplayname = models.CharField(max_length=30L)
    location = models.CharField(max_length=100L, blank=True)
    portfolioid = models.CharField(max_length=30L, blank=True)
    katedataquery = models.CharField(max_length=4000L, blank=True)
    katedataqueryid = models.CharField(max_length=30L, blank=True)
    katedataquerydisplayname = models.CharField(max_length=30L, blank=True)
    cornershopid = models.CharField(max_length=30L, blank=True)
    cornershopdataquery = models.CharField(max_length=4000L, blank=True)
    cornershopdataqueryid = models.CharField(max_length=30L, blank=True)
    cornershopdataquerydisplayname = models.CharField(max_length=400L, blank=True)
    projectcomments = models.CharField(max_length=4000L, blank=True)
    targetpedia = models.CharField(max_length=200L, blank=True)
    goid = models.CharField(max_length=200L, blank=True)
    ec = models.CharField(max_length=200L, blank=True)
    cath = models.CharField(max_length=200L, blank=True)
    pfam = models.CharField(max_length=200L, blank=True)
    scop = models.CharField(max_length=200L, blank=True)
    interpro = models.CharField(max_length=200L, blank=True)
    uniprot = models.CharField(max_length=200L, blank=True)
    project_doi = models.CharField(max_length=100L, blank=True)
    modified_by = models.CharField(max_length=30L, blank=True)
    modified_on = models.DateField(null=True, blank=True)

    class Meta:
        db_table = 'projects'


class Structures(models.Model):
    projectname = models.CharField(max_length=30L)
    structure_id = models.IntegerField(primary_key=True)
    hierarchy = models.CharField(max_length=30L, blank=True)
    parent = models.CharField(max_length=30L, blank=True)
    pdbfilename = models.CharField(max_length=100L, blank=True)
    pdbid = models.CharField(max_length=20L, blank=True)
    displayname = models.CharField(max_length=100L, blank=True)
    densityfile = models.CharField(max_length=100L, blank=True)
    strucsource = models.CharField(max_length=13L, blank=True)
    str_type = models.CharField(max_length=5L, blank=True)
    resol = models.CharField(max_length=10L, blank=True)
    rfactor = models.CharField(max_length=10L, blank=True)
    rfree = models.CharField(max_length=10L, blank=True)
    structurecomments = models.CharField(max_length=4000L, blank=True)
    templateid = models.CharField(max_length=30L, blank=True)
    siteid = models.CharField(max_length=200L, blank=True)
    ligandid = models.CharField(max_length=10L, blank=True)
    ligandresname = models.CharField(max_length=3L, blank=True)
    ligandchainname = models.CharField(max_length=1L, blank=True)
    ligandresnum = models.CharField(max_length=10L, blank=True)
    ligandregno = models.CharField(max_length=20L, blank=True)
    ligandname = models.CharField(max_length=20L, blank=True)
    ligandclass = models.CharField(max_length=20L, blank=True)
    ligandsmiles = models.CharField(max_length=1000L, blank=True)
    ligandsite = models.CharField(max_length=4000L, blank=True)
    ligandcomments = models.CharField(max_length=4000L, blank=True)
    extref_doi = models.CharField(max_length=100L, blank=True)
    extref_pp = models.CharField(max_length=20L, blank=True)
    extref_vol = models.CharField(max_length=6L, blank=True)
    extref_year = models.CharField(max_length=5L, blank=True)
    extref_jrnl = models.CharField(max_length=100L, blank=True)
    extref_linktext = models.CharField(max_length=100L, blank=True)
    added_by = models.CharField(max_length=30L, blank=True)
    added_on = models.DateField(null=True, blank=True)
    modified_by = models.CharField(max_length=30L, blank=True)
    modified_on = models.DateField(null=True, blank=True)
    class Meta:
        db_table = 'structures'


class Surfaces(models.Model):
    projectname = models.CharField(max_length=30L)
    pdbfilename = models.CharField(max_length=100L)
    ligandid = models.CharField(max_length=10L)
    surfacename = models.CharField(primary_key=True,max_length=30L)
    surfacetype = models.CharField(max_length=30L)
    surfacecomment = models.CharField(max_length=4000L, blank=True)
    surfacefile = models.CharField(max_length=4000L, blank=True)
    extsource = models.CharField(max_length=20L, blank=True)
    extlink = models.CharField(max_length=200L, blank=True)
    class Meta:
        db_table = 'surfaces'


class Templates(models.Model):
    templateid = models.CharField(primary_key=True,max_length=30L)
    template_chainid = models.CharField(max_length=1L)
    template_type = models.CharField(max_length=5L)
    projectname = models.CharField(max_length=30L)
    templatepdbfilename = models.CharField(max_length=100L)
    templatecomments = models.CharField(max_length=4000L, blank=True)
    template_doi = models.CharField(max_length=100L, blank=True)
    template_site = models.CharField(max_length=200L, blank=True)
    class Meta:
        db_table = 'templates'


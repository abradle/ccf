# Models for the LLOOMPPA app

from django.db import models
# Imports for foreign key links
from MMPMaker.models import Act3DFrag,MMPFrag,MMPComparison,ActPharmaPoint,MMP
from Pharmacophore.models import PharmaPoint
from IOhandle.models import Molecule,Target,Protein,Compound
from WONKA.models import Residue
# These models are specific to LLOOMMPPAA to be used to speed up data access for the LLOOMMPPAA APP


class Cmpd_Complexity(models.Model):
    """Model to contain complexity information about a given compound"""
    # Compound ID
    cmpd_id = models.ForeignKey(Compound)
    # Synthetic accesibilty of a compound
    sa_score = models.FloatField(null=True)
    # The number of heavy atoms
    size = models.IntegerField(null=True)
    # The number of on bits for MF
    num_morg = models.IntegerField(null=True)
    # The number of on bits for MACCS keys
    num_maccs = models.IntegerField(null=True)
    # The number of chiral centers
    num_chir = models.IntegerField(null=True)
    # The number of rings
    num_rings = models.IntegerField(null=True)
    # The number of rotatable bonds
    num_rot_bonds = models.IntegerField(null=True)
    # The number of HBA
    num_hba = models.IntegerField(null=True)
    # The number of HBD
    num_hbd = models.IntegerField(null=True)


class LLConf(models.Model):
    """Model to contain information about a new molecule created"""
    # The target
    target_id = models.ForeignKey(Target)
    # The protein it refers to
    prot_id = models.ForeignKey(Protein)
    # The closest two units come together
    clash = models.FloatField()
    # The molecule it is based upon
    mol_ref = models.ForeignKey(Molecule, related_name="%(app_label)s_%(class)s_mol_ref")
    # The molecule conformation it is
    mol_id = models.ForeignKey(Molecule, related_name="%(app_label)s_%(class)s_mol_id")
    # The strain of the molecule - used to prioritise compoudns
    strain = models.FloatField(null=True)
    # The min RMSD of this CONF
    rmsd = models.FloatField(null=True)
    # The conformation number
    conf_num = models.IntegerField()
    # The shape distance
    shape_dist = models.FloatField()

    class Meta:
        # Unique constraints
        unique_together = ('mol_id', 'mol_ref', 'target_id', 'prot_id', 'conf_num')


class Reaction(models.Model):
    """A model to store the details of an individual reaction"""
    # The name of the reaction
    name = models.CharField(max_length=100, unique=True)
    # The smarts for the reaction
    react_smarts = models.CharField(max_length=100)
    # The smarts for finding the functional group that can be split at this place
#    funct_smarts = models.CharField(max_length=100)
    retro_smarts = models.CharField(max_length=100)
    cont_smarts = models.CharField(max_length=100)
#    # The indices for where the cut will go
#    coord_one = models.IntegerField()
#    coord_two = models.IntegerField()


class PossReact(models.Model):
    """Model to store a potential reaction for a molecule"""
    # The reaction
    react_id = models.ForeignKey(Reaction)
    # The molecule
    mol_id = models.ForeignKey(Molecule)
    # Frag one - the bit that will be retained
    retained_frag = models.CharField(max_length=100)
    retained_frag_frag = models.CharField(max_length=100)
    retained_frag_context = models.CharField(max_length=100)
    # Frag two - the bit that will be replaced
    replaced_frag = models.CharField(max_length=100)
    replaced_frag_frag = models.CharField(max_length=100)
    replaced_frag_context = models.CharField(max_length=100)
    # Get  the SDF data for one frag
    ret_frag_sdf_info = models.TextField()
    # Get the SDF info of the other one
    rep_frag_sdf_info = models.TextField()

    class Meta:
        # Unique constraints
        unique_together = ('react_id', 'mol_id', 'retained_frag', 'replaced_frag',)


class Reactant(models.Model):
    """A model to add extra information to a compound"""
    # The compound id
    cmpd_id = models.ForeignKey(Compound, unique=True)
    # An optional flag to determine if these compounds are available
    is_avaialable = models.BooleanField(default=True)
    # The type of reaction they can be part of
    react_id = models.ManyToManyField(Reaction)


class ReactantLib(models.Model):
    """A model to store the reactants for a library"""
    # The link to compounds
    reactant_id = models.ManyToManyField(Reactant)
    # The library they come from -> so we don't always use all of the
    library_name = models.CharField(max_length=100, unique=True)
    # The type they come under
    reactant_type = models.CharField(max_length=100)


class Process(models.Model):
    """A model to determine how a reactant was made"""
    # Whether this compound was made or not
    is_made_lloommppaa = models.BooleanField()
    # The reactant id
    reactant_id = models.ForeignKey(Reactant, null=True)
    # The reaction id
    reaction_id = models.ForeignKey(Reaction)
    # The molecule being reacted
    mol_id = models.ForeignKey(Molecule)


class Product(models.Model):
    """A model to store the products of a library - either made here
    or given directly as an SDF file"""
    # The link to the compound
    cmpd_id = models.ForeignKey(Compound, unique=True)
    process_id = models.ManyToManyField(Process)


class ProductLib(models.Model):
    """A model to store a library of finished products"""
    # The link to the compounds
    product_id = models.ManyToManyField(Product)
    # The name of the lib
    lib_name = models.CharField(max_length=100, unique=True)


class ReactionProcess(models.Model):
    """Model to store the reactions being used"""
    # The title of the reaction process
    title = models.CharField(max_length=100)
    # The context
    context = models.CharField(max_length=100)
    # The fragment to react
    react_frag = models.CharField(max_length=100)
    # The molecule being used
    mol_id = models.ForeignKey(Molecule)
    # The proteins used as a conflict
    prot_id = models.ManyToManyField(Protein)
    # The name of the reaction
    react_id = models.ForeignKey(Reaction)
    # The list of products used
    products_id = models.ManyToManyField(ProductLib)
    # The list of products left
    product_queue = models.ManyToManyField(Compound, related_name="%(app_label)s_%(class)s_prod_q")
    # The list of reactants left
    reactant_queue = models.ManyToManyField(Compound, related_name="%(app_label)s_%(class)s_react_q")
    # The list of reactants used
    reactants_id = models.ManyToManyField(ReactantLib)
    # Whether or not this is a reaction
    is_reaction = models.BooleanField()
    # Whether or not this has been completed
    finished = models.BooleanField()
    # The LLCOnfs associated to this reaction
    ll_conf = models.ManyToManyField(LLConf)
    # The stage that the process is in
    proc_stage = models.CharField(max_length=100)
    # The completion of the current stage
    stage_completion = models.IntegerField()


class Histogram(models.Model):
    """Model to store a histogram"""
    # The unique elements for this table
    test_made = models.CharField(max_length=1000)
    # The number of samples
    num_samps = models.CharField(max_length=1000)
    # The title of the hist
    hist_title = models.CharField(max_length=1000)

    class Meta:
        # Unique constraints
        unique_together = ('test_made', 'num_samps', 'hist_title',)


class HistBin(models.Model):
    """Model to store a histogram bin 
    No target id here - since the Histograms can be compared between 
    targets"""
    # The min
    min_char = models.CharField(max_length=30)
    # The max
    max_char = models.CharField(max_length=30)
    # The graph it relates to
    hist_id = models.ForeignKey(Histogram)

    class Meta:
        # Unique constraints - the charred Float (taking the first 6 characters)
        unique_together = ('min_char', 'max_char', 'hist_id', )


class GridPoint(models.Model):
    """Model to store a grid point for a target"""
    # The x_coord
    x_char = models.CharField(max_length=30)
    # The y_coord
    y_char = models.CharField(max_length=30)
    # The z_coord
    z_char = models.CharField(max_length=30)
    # The target it relates to
    target_id = models.ForeignKey(Target)

    class Meta:
        # Unique constraints - the charred Float (taking the first 6 characters)
        unique_together = ('x_char', 'y_char', 'z_char', 'target_id')


class CmpdScore(models.Model):
    """Model containing the picks of list"""
    # The score
    score = models.FloatField()
    # The list of compounds
    cmpd_id = models.ForeignKey(Compound)


class RunAnalysis(models.Model):
    """Model to store the analysis run this time round"""
    # The clash
    clash = models.FloatField()
    # The 
    rmsd = models.FloatField()
    # The shape dist
    shape_dist = models.FloatField()
    # The ReactionProcess id
    react_proc_id = models.ForeignKey(ReactionProcess)
    # The number of trials
    num_trials = models.FloatField()
    # The number of iterations
    num_iters = models.FloatField()

    class Meta:
        # Unique constraints - the charred Float (taking the first 6 characters)
        unique_together = ('clash', 'rmsd', 'shape_dist', 'react_proc_id', 'num_trials', 'num_iters')


class DivList(models.Model):
    """Model containing the picks of list"""
    # The method used
    method_id = models.CharField(max_length=100)
    # The ReactionProcess id
    react_proc_id = models.ForeignKey(ReactionProcess)
    # The number of samples
    num_samps = models.IntegerField()
    # The list of compounds
    score_id = models.ManyToManyField(CmpdScore)
    # The link to the run used in this analysis
    run_id = models.ForeignKey(RunAnalysis)


class GPVal(models.Model):
    """Model to contain the value for this bin of this graph"""
    # The value
    value = models.FloatField(null=True)
    # The histogram it's for e.g. Acceptor
    type_id = models.CharField(max_length=100)
    # The bin it's for
    gp_id = models.ForeignKey(GridPoint)
    # The target it's for
    target_id = models.ForeignKey(Target)
    # The reaction process it relates to
    my_anal_id = models.ForeignKey(RunAnalysis)
    # The method it relates to
    method_id = models.CharField(max_length=100)
    # The number of samples made
    num_samps = models.FloatField()
    # Some text that represents this as PDB data
    pdb_info = models.TextField(null=True)
    # Some text that represents this as SDF data
    sdf_info = models.TextField(null=True)

    class Meta:
        # Unique constraints
        unique_together = ('gp_id', 'target_id', 'type_id', 'my_anal_id',
                           'method_id', 'num_samps')


class HistBinVal(models.Model):
    """Model to contain the value for this bin of this graph"""
    # The value
    value = models.FloatField(null=True)
    # The histogram it's for e.g. Acceptor
    hist_id = models.ForeignKey(Histogram)
    # The bin it's for
    bin_id = models.ForeignKey(HistBin)
    # The target it's for
    target_id = models.ForeignKey(Target)
    # The value
    value = models.FloatField(null=True)
    # The histogram it's for - like for example 
    type_id = models.CharField(max_length=100)
    # The reaction process it relates to
    my_anal_id = models.ForeignKey(RunAnalysis)
    # The method it relates to
    method_id = models.CharField(max_length=100)
    # The number of samples 
    num_samps = models.FloatField()

    class Meta:
        # Unique constraints
        unique_together = ('hist_id', 'bin_id', 'target_id', 'type_id',
                           'my_anal_id', 'method_id', 'num_samps',)

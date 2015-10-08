import os,sys
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "WebApp.settings")
from WebApp import settings#DATABASES  as dbs
from LLOOMMPPAA.functions import make_synth_points, do_diff_analysis, find_clashes
from PLIFS.functions import restore_probs
from PLIFS.interaction_defs import init_schemes
from IOhandle.models import Target, Molecule


def do_tot_analys(target_id=2):
#    make_synth_points(target_id, mmpdbs=["MMPDB"], sdf_flag=True)
### HERE WE SET THE ALLOWED RMSD / CLASH / SHAPE DIST CONSTRAINTS ---->>>> WE NEED SOME RULES FOR EACH  ONE
    rmsd = 0.5## Seems to allow for on average 2.0 conformations per compound - which feels reasonable (both options)
    shape_dist = 0.2#1.0## NORMALLY DO NO FILTERING ON SHAPE DIST - VERY SPECIFIC VALUE WHEN NEEDED e.g. BRD1
    ## The most controversial value could be set to go below - we find closest approach of known ligand set and allow a furhter -0.5 A of leeway
    clash = find_clashes(target_id)
    #Or at the max for a current ligand + 10% or something else equally abitrary......
#    restore_probs(target_id, "LLOOMMPPAA", clash=clash, rmsd=rmsd, shape_dist=shape_dist)
    target = Target.objects.get(pk=target_id)
    mol_ids = Molecule.objects.filter(prot_id__code=target.title + "SYNTH")
    init_schemes(refresh=True, mol_ids=mol_ids)
    do_diff_analysis(target_id, 1000, [30])#, 50, 100, 150])


if os.path.split(settings.DATABASES["default"]["NAME"])[1] != "WONKA_db_LL_NEW":
    print "NOT CONNECTED TO CORRECT DB - CHANGE IT UP - WONKA_db_LL_NEW"
else:
    pass
do_tot_analys()

import os,sys
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "WebApp.settings")
from WebApp import settings#DATABASES  as dbs
from LLOOMMPPAA.reactions import find_procs_to_run, find_dead_procs
find_procs_to_run()
find_dead_procs()
find_procs_to_run()
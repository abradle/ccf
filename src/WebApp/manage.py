#!/usr/bin/env python
# All these imports are required for pyinstaller
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "WebApp.settings")

### ADD THIS
import djutils
# Now import some python modules
import numpy
import subprocess
import ctypes
import sqlite3
import StringIO
import random
import math
from argparser import myargparse
import getopt
import tempfile
import sys
import re
import string
import ast
import WebApp.settings
import Cookie
import uuid
# Import the core django stuff
from django.contrib import messages
import django.core.exceptions
from django.core.exceptions import ValidationError
from django.db.models import Count
import django.test
import HTMLParser
# Import the RDKit stuff
from rdkit import Chem
from rdkit.ML.Cluster import Butina
from rdkit.Chem.rdShapeHelpers import  ShapeProtrudeDist,ComputeConfBox,ComputeUnionBox,EncodeShape
from rdkit.Geometry.rdGeometry import UniformGrid3D
from rdkit.Chem import MCS
from rdkit.Chem import AllChem
# Now import the Django apps stuff
import IOhandle.models
import IOhandle.functions
import Pharmacophore.models
import Pharmacophore.functions
import MMPMaker.models
import MMPMaker.functions
import MMPMaker.helpers
import Viewer
try:
    import gunicorn
except:
    pass
import OOMMPPAA
import jfu
import IOhandle
import MMPMaker
import WebApp.urls
import Group
from loading import load_mols, load_activity_data, do_oommppaa_proc, initialise_dummys

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "WebApp.settings")
    print "HELLO"
    myargparse()

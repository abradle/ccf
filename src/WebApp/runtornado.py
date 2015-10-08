#!/usr/bin/env python

try:
    import sys
    me = sys._MEIPASS
    import os
    import django.core.handlers.wsgi
    from tornado import wsgi,web,httpserver,ioloop
    import time, signal
except AttributeError:
    import os
    print "IMPORTED OS"
    import sys
    print "IMPORTED SYS"
    import re
    print "IMPORTED RE"
    import string
    print "IMPORTED STRING"
    import ast
    print  "IMPORTED AST"
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "WebApp.settings")
    try:
        import django.templatetags.jfutags
        import django.templatetags.disqus_tags
    except ImportError:
        pass
    import WebApp.settings
    print "IMPORTED SETTINGS"
    import Cookie
    import sqlite3
    print "IMPORTED COOKIE"
    from django.contrib import messages
    import django.core.exceptions
    from django.core.exceptions import ValidationError
    from django.contrib.sessions import *
    from django.db.models import Count
    from django.core.management import execute_from_command_line
    print "IMPORTED DJANGO STUFF"
    import django.test
    print "IMPORTED TEST"
    import HTMLParser
    print "IMPORTED PARSER"
    from rdkit import Chem,RDConfig
    from rdkit.ML.Cluster import Butina
    from rdkit.Chem.rdShapeHelpers import ShapeProtrudeDist,ComputeConfBox,ComputeUnionBox,EncodeShape
    from rdkit.Geometry.rdGeometry import UniformGrid3D
    from rdkit.Chem import MCS
    from rdkit.Chem import AllChem
    from rdkit import RDConfig,rdBase
    from rdkit import DataStructs
    from rdkit.DataStructs import cDataStructs
    print "IMPORTED RDKIT"
    import IOhandle.models
    print "IMPORTED IOHANDLE.MODELS"
    import Pharmacophore.models
    print "IMPORTED PHARMACOPHORE.MODELS"
    import Pharmacophore.functions
    print "IMPORTED IOHANDLE.FUNCTIONS"
    import MMPMaker.models
    print "IMPORTED MMPMAKER.MODELS"
    import MMPMaker.functions
    print "IMPORTED MMPMAKER.FUNCTIONS"
    import MMPMaker.helpers
    print "IMPORTED MMPMAKER.HELPERS"
    import numpy
    print "IMPORTED NUMPY"
    import ctypes
    print "IMPORTED CTYPES"
    import StringIO
    print "IMPORTED STRINGIO"
    import random
    print "IMPORTED RANDOM"
    import math
    print "IMPORTED MATH"
    from tornado import wsgi
    print "DONE TORNADO WSGI"
    from tornado import web
    print "DONE TORNADO WEB"
    import django.core.handlers.wsgi
    import django.views.defaults
    print "DONE OPTPARSE STUFF"
    from tornado import httpserver
    print "DONE TORNADO HTTPSERVER"
    from tornado import ioloop
    import time, signal
    print "DONE TORNADO IOLOOP"
    import webbrowser
    print "IMPORTED WEB BROWSER"
    import Viewer
    try:
        import gunicorn
    except:
        pass
    import OOMMPPAA
    import jfu
    import IOhandle
    import MMPMaker
    import WebApp
    import Group

from optparse import OptionParser
parser = OptionParser()
parser.add_option('--port', dest='port', default="9020")
parser.add_option('--address', dest='address', default="127.0.0.1")
options, args = parser.parse_args()

app_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.dirname(app_dir))
wsgi_app = wsgi.WSGIContainer(django.core.handlers.wsgi.WSGIHandler())


def start_tornado():
    application = web.Application([
        (r"/static/(.*)", web.StaticFileHandler, {"path": os.path.join(os.path.split(sys.argv[0])[0],'data/static')}),
        (r".*", web.FallbackHandler, dict(fallback=wsgi_app)),])
    server = httpserver.HTTPServer(application)
    server.listen(options.port, options.address)
    print "Starting Tornado"
    try:
        ioloop.IOLoop.instance().start()
    except KeyboardInterrupt:
        print "Tornado finished"


def stop_tornado():
    from tornado import ioloop
    ioloop = ioloop.IOLoop.instance()
    ioloop.add_callback(lambda x: x.stop(), ioloop)
    print "Asked Tornado to exit"


def runserver():
    os.environ['DJANGO_SETTINGS_MODULE'] = 'WebApp.settings'
    # open the URL
    print "Starting server..."
    print "Serving at " + options.address + " on port " + options.port
    start_tornado()


if __name__ == '__main__':
    runserver()

import os, sys
name = os.path.join(os.path.split(sys.argv[0])[0], 'data/wonka_db')
backend = 'django.db.backends.sqlite3'
user = ''
password = ''
host = ''            # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
port = ''
extra_apps = ['IOhandle',
    'Pharmacophore',
    'MMPMaker',
    'Group',
    'Viewer',
    'gunicorn',
    'OBSERVATIONS',
    'OOMMPPAA',
    'jfu',
    'south',
    "WONKA",
    "LLOOMMPPAA",
    "PLIFS",]

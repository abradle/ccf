import os,sys
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "WebApp.settings")
import WebApp.settings
from django.core.management import execute_from_command_line

if __name__ == "__main__":
    ### Need to write my own parser -> WRITE PARSER
    execute_from_command_line(sys.argv)

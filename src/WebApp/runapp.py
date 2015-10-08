import os
import sys
import time
import subprocess
from time import sleep
from PySide.QtCore import *
from PySide.QtGui import *
from PySide.QtWebKit import *
import urllib2
import signal
import psutil
import pkg_resources


def onClose():
    PROCNAME = "OOMMPPAAserver.exe"
    for proc in psutil.process_iter():
        if proc.name == PROCNAME:
            proc.kill()

if __name__ == "__main__":
    preamble = """
 _______  _______  _______  _______  _______  _______  _______  _______
(  ___  )(  ___  )(       )(       )(  ____ )(  ____ )(  ___  )(  ___  )
| (   ) || (   ) || () () || () () || (    )|| (    )|| (   ) || (   ) |
| |   | || |   | || || || || || || || (____)|| (____)|| (___) || (___) |
| |   | || |   | || |(_)| || |(_)| ||  _____)|  _____)|  ___  ||  ___  |
| |   | || |   | || |   | || |   | || (      | (      | (   ) || (   ) |
| (___) || (___) || )   ( || )   ( || )      | )      | )   ( || )   ( |
(_______)(_______)|/     \||/     \||/       |/       |/     \||/     \|
"""
    print preamble
    print "Loading please wait...."
    os.environ['path'] = ""
    p = subprocess.Popen([r'OOMMPPAAserver.exe', '--port', '8020', '--address', '127.0.0.1'])
    # Now try for 30 seconds to get the server
    i = 0
    while True:
        i += 1
        try:
            req = urllib2.urlopen("http://127.0.0.1:8020/")
            break
        except:
          # Wait for 45 seconds
            if i == 45:
                print "SERVER NOT RESPONDING!!!"
                sys.exit()
            time.sleep(1)
            pass

    app = QApplication(sys.argv)
    app.addLibraryPath("")
    app.aboutToQuit.connect(onClose)
    web = QWebView()
    mSettings = web.settings()
    mSettings.setAttribute(QWebSettings.LocalStorageEnabled, True)
    mSettings.setAttribute(QWebSettings.OfflineStorageDatabaseEnabled, True)
    mSettings.setAttribute(QWebSettings.DeveloperExtrasEnabled, True)
    mSettings.setAttribute(QWebSettings.PluginsEnabled,  True)
    mSettings.setAttribute(QWebSettings.OfflineWebApplicationCacheEnabled, True)
    mSettings.enablePersistentStorage(path=os.getenv("HOME"))
    web.load(QUrl("http://127.0.0.1:8020/OOMMPPAA/"))
    web.show()
    sys.exit(app.exec_())


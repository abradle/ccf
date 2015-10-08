import subprocess
import sys

args = sys.argv

try:
    import sys
    me = sys._MEIPASS
    new_args = []
    args[0] = "winmanage.exe"
except:
    new_args = ["python"]
    args[0]= "run.py"
new_args.extend(args)
out_std = open("out1.std","w")
p = subprocess.Popen(new_args,stderr=out_std)
p.wait()

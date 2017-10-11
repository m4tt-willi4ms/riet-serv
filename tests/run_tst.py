import subprocess, os

# p = subprocess.Popen(["setpaths"], stdout=subprocess.PIPE,shell=True)
# os.environ['PATH'] = ':'.join([os.getenv('PATH'), os.getenv('PYTHONPATH')])
p2 = subprocess.Popen([r"C:\cctbx\build\bin\scitbx.python.bat","tst_RietveldPhases.py"], stdout=subprocess.PIPE,shell=True)
print p2.communicate()[0]
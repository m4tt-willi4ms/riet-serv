from distutils.core import setup
import numpy
import py2exe
import matplotlib
from distutils.core import setup
import py2exe
from distutils.filelist import findall
import os

matplotlibdatadir = matplotlib.get_data_path()
matplotlibdata = findall(matplotlibdatadir)
matplotlibdata_files = []
for f in matplotlibdata:
    dirname = os.path.join('matplotlibdata', f[len(matplotlibdatadir)+1:])
    matplotlibdata_files.append((os.path.split(dirname)[0], [f]))

setup(
    options = {
            "py2exe":{
            "dll_excludes": ["MSVCP90.dll", "HID.DLL", "w9xpopen.exe"],
            # "excludes": ['_gtkagg', '_tkagg'],
            "includes": ["matplotlib.backends.backend_tkagg"],
            "packages": ["FileDialog"],
        }
    },
    data_files =  matplotlibdata_files,
    console = [{'script': 'tst_Rietveld_Cement.py'}]
)
# setup(console=['tst_Rietveld_Cement.py'])
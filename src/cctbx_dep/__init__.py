r"""
A short python script to set the paths needed in order to use the cctbx library.
"""
import sys, os

#Path to the cctbx build directory goes here
os.environ["LIBTBX_BUILD"] = r"C:/rietveld2_cctbx/build"

#These commands update automatically
sys.path.append(os.path.abspath(os.environ["LIBTBX_BUILD"]
   + r"/../modules/cctbx_project"))
sys.path.append(os.path.abspath(os.environ["LIBTBX_BUILD"]
   + r"/../modules/cctbx_project/boost_adaptbx"))
sys.path.append(os.path.abspath(os.environ["LIBTBX_BUILD"]
   + r"/../modules/cctbx_project/libtbx/pythonpath"))
sys.path.append(os.path.abspath(os.environ["LIBTBX_BUILD"]
   + r"/lib"))
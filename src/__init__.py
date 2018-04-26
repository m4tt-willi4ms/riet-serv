r"""
A short python script to set the paths needed in order to use the cctbx codebase.
"""
import sys, os

os.environ["LIBTBX_BUILD"] = r"C:/cctbx/build" #os.path.abspath(r"../../../build")

sys.path.append(os.path.abspath(os.environ["LIBTBX_BUILD"] + r"/../modules/cctbx_project"))
sys.path.append(os.path.abspath(os.environ["LIBTBX_BUILD"] + r"/../modules/cctbx_project/boost_adaptbx"))
sys.path.append(os.path.abspath(os.environ["LIBTBX_BUILD"] + r"/../modules/cctbx_project/libtbx/pythonpath"))
sys.path.append(os.path.abspath(os.environ["LIBTBX_BUILD"] + r"/lib"))
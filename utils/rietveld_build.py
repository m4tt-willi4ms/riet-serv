from __future__ import division, print_function, absolute_import
import os, sys
import shutil

INPUT_DIR = "C:\\cctbx\\build"
OUTPUT_DIR = "C:\\rietveld_cctbx\\build"

DIRS_TO_CREATE = [
    "\\..\\modules",
    "\\lib",
    ]

PACKAGES_TO_CREATE = [
    "\\..\\modules\\cctbx_project\\iotbx",
    "\\..\\modules\\cctbx_project\\iotbx\\pdb",
    "\\..\\modules\\cctbx_project\\libtbx",
    "\\..\\modules\\cctbx_project\\libtbx\\auto_build",
    "\\..\\modules\\cctbx_project\\libtbx\\phil",
    "\\..\\modules\\cctbx_project\\libtbx\\queuing_system_utils",
    "\\..\\modules\\cctbx_project\\libtbx\\test_utils",
    "\\..\\modules\\cctbx_project\\boost_adaptbx",
    "\\..\\modules\\cctbx_project\\cctbx",
    "\\..\\modules\\cctbx_project\\cctbx\\array_family",
    "\\..\\modules\\cctbx_project\\cctbx\\miller",
    "\\..\\modules\\cctbx_project\\cctbx\\sgtbx",
    "\\..\\modules\\cctbx_project\\cctbx\\uctbx",
    "\\..\\modules\\cctbx_project\\cctbx\\xray",
    "\\..\\modules\\cctbx_project\\cctbx\\maptbx",
    "\\..\\modules\\cctbx_project\\cctbx\\crystal",
    "\\..\\modules\\cctbx_project\\cctbx\\covariance",
    "\\..\\modules\\cctbx_project\\cctbx\\geometry",
    "\\..\\modules\\cctbx_project\\cctbx\\eltbx",
    "\\..\\modules\\cctbx_project\\cctbx\\eltbx\\xray_scattering",
    "\\..\\modules\\cctbx_project\\scitbx",
    "\\..\\modules\\cctbx_project\\scitbx\\array_family",
    "\\..\\modules\\cctbx_project\\scitbx\\matrix",
    "\\..\\modules\\cctbx_project\\scitbx\\math",
    "\\..\\modules\\cctbx_project\\scitbx\\lbfgs",
    "\\..\\modules\\cctbx_project\\scitbx\\linalg",
    "\\..\\modules\\cctbx_project\\scitbx\\python_utils",

]

FILES_TO_COPY = [
    "\\libtbx_env",
    "\\..\\modules\\cctbx_project\\libtbx\\forward_compatibility.py",
    "\\..\\modules\\cctbx_project\\libtbx\\cpp_function_name.py",
    "\\..\\modules\\cctbx_project\\libtbx\\easy_run.py",
    "\\..\\modules\\cctbx_project\\libtbx\\load_env.py",
    "\\..\\modules\\cctbx_project\\libtbx\\env_config.py",
    "\\..\\modules\\cctbx_project\\libtbx\\path.py",
    "\\..\\modules\\cctbx_project\\libtbx\\str_utils.py",
    "\\..\\modules\\cctbx_project\\libtbx\\utils.py",
    "\\..\\modules\\cctbx_project\\libtbx\\option_parser.py",
    "\\..\\modules\\cctbx_project\\libtbx\\introspection.py",
    "\\..\\modules\\cctbx_project\\libtbx\\auto_build\\regenerate_module_files.py",
    "\\..\\modules\\cctbx_project\\libtbx\\auto_build\\installer_utils.py",
    "\\..\\modules\\cctbx_project\\libtbx\\queuing_system_utils\\sge_utils.py",
    "\\..\\modules\\cctbx_project\\libtbx\\queuing_system_utils\\pbs_utils.py",
    "\\..\\modules\\cctbx_project\\libtbx\\phil\\tokenizer.py",
    "\\..\\modules\\cctbx_project\\libtbx\\containers.py",
    "\\..\\modules\\cctbx_project\\libtbx\\math_utils.py",
    "\\..\\modules\\cctbx_project\\libtbx\\assert_utils.py",
    "\\..\\modules\\cctbx_project\\libtbx\\smart_open.py",
    "\\..\\modules\\cctbx_project\\libtbx\\complex_math.py",
    "\\..\\modules\\cctbx_project\\libtbx\\table_utils.py",
    "\\..\\modules\\cctbx_project\\iotbx\\pdb\\secondary_structure.py",
    "\\..\\modules\\cctbx_project\\iotbx\\pdb\\records.py",
    "\\..\\modules\\cctbx_project\\iotbx\\pdb\\hierarchy.py",
    "\\..\\modules\\cctbx_project\\iotbx\\pdb\\atom_name_interpretation.py",
    "\\..\\modules\\cctbx_project\\iotbx\\pdb\\hybrid_36.py",
    "\\..\\modules\\cctbx_project\\cctbx\\array_family\\flex.py",
    "\\..\\modules\\cctbx_project\\cctbx\\eltbx\\neutron.py",
    "\\..\\modules\\cctbx_project\\cctbx\\eltbx\\attenuation_coefficient.py",
    "\\..\\modules\\cctbx_project\\cctbx\\eltbx\\wavelengths.py",
    "\\..\\modules\\cctbx_project\\cctbx\\eltbx\\tiny_pse.py",
    "\\..\\modules\\cctbx_project\\cctbx\\crystal\\find_best_cell.py",
    "\\..\\modules\\cctbx_project\\cctbx\\sgtbx\\lattice_symmetry.py",
    "\\..\\modules\\cctbx_project\\cctbx\\sgtbx\\subgroups.py",
    "\\..\\modules\\cctbx_project\\cctbx\\sgtbx\\bravais_types.py",
    "\\..\\modules\\cctbx_project\\cctbx\\xray\\ext.py",
    "\\..\\modules\\cctbx_project\\cctbx\\xray\\observation_types.py",
    "\\..\\modules\\cctbx_project\\cctbx\\xray\\scatterer.py",
    "\\..\\modules\\cctbx_project\\cctbx\\xray\\structure.py",
    "\\..\\modules\\cctbx_project\\cctbx\\xray\\target_functors.py",
    "\\..\\modules\\cctbx_project\\cctbx\\xray\\weighting_schemes.py",
    "\\..\\modules\\cctbx_project\\cctbx\\xray\\minimization.py",
    "\\..\\modules\\cctbx_project\\cctbx\\adptbx.py",
    "\\..\\modules\\cctbx_project\\cctbx\\math_module.py",
    "\\..\\modules\\cctbx_project\\cctbx\\r_free_utils.py",
    "\\..\\modules\\cctbx_project\\scitbx\\array_family\\flex.py",
    "\\..\\modules\\cctbx_project\\scitbx\\array_family\\shared.py",
    "\\..\\modules\\cctbx_project\\scitbx\\python_utils\\dicts.py",
    "\\..\\modules\\cctbx_project\\scitbx\\math\\ext.py",
    "\\..\\modules\\cctbx_project\\scitbx\\math\\gaussian.py",
    "\\..\\modules\\cctbx_project\\scitbx\\linalg\\eigensystem.py",
    "\\..\\modules\\cctbx_project\\scitbx\\linalg\\ext.py",
    "\\..\\modules\\cctbx_project\\scitbx\\linalg\\householder.py",
    "\\..\\modules\\cctbx_project\\scitbx\\python_utils\\misc.py",
    "\\..\\modules\\cctbx_project\\scitbx\\cubicle_neighbors.py",
    "\\..\\modules\\cctbx_project\\scitbx\\fftpack.py",
    ".\\lib\\boost_python.dll",
    ".\\lib\\boost_python_meta_ext.pyd",
    ".\\lib\\iotbx_cif_ext.pyd",
    ".\\lib\\iotbx_pdb_ext.pyd",
    ".\\lib\\iotbx_pdb_hierarchy_ext.pyd",
    ".\\lib\\boost_optional_ext.pyd",
    ".\\lib\\boost_rational_ext.pyd",
    ".\\lib\\std_pair_ext.pyd",
    ".\\lib\\scitbx_array_family_flex_ext.pyd",
    ".\\lib\\scitbx_array_family_shared_ext.pyd",
    ".\\lib\\scitbx_stl_set_ext.pyd",
    ".\\lib\\scitbx_stl_vector_ext.pyd",
    ".\\lib\\scitbx_stl_map_ext.pyd",
    ".\\lib\\scitbx_random_ext.pyd",
    ".\\lib\\scitbx_math_ext.pyd",
    ".\\lib\\scitbx_lbfgs_ext.pyd",
    ".\\lib\\scitbx_fftpack_ext.pyd",
    ".\\lib\\scitbx_linalg_ext.pyd",
    ".\\lib\\scitbx_cubicle_neighbors_ext.pyd",
    ".\\lib\\cctbx_array_family_flex_ext.pyd",
    ".\\lib\\cctbx_uctbx_ext.pyd",
    ".\\lib\\cctbx_sgtbx_ext.pyd",
    ".\\lib\\cctbx_math_ext.pyd",
    ".\\lib\\cctbx_eltbx_neutron_ext.pyd",
    ".\\lib\\cctbx_eltbx_attenuation_coefficient_ext.pyd",
    ".\\lib\\cctbx_eltbx_wavelengths_ext.pyd",
    ".\\lib\\cctbx_eltbx_xray_scattering_ext.pyd",
    ".\\lib\\cctbx_eltbx_tiny_pse_ext.pyd",
    ".\\lib\\cctbx_xray_ext.pyd",
    ".\\lib\\cctbx_covariance_ext.pyd",
    ".\\lib\\cctbx_geometry_ext.pyd",
    ".\\lib\\cctbx_miller_ext.pyd",
    ".\\lib\\cctbx_maptbx_ext.pyd",
    ".\\lib\\cctbx_adptbx_ext.pyd",
    ".\\lib\\cctbx_crystal_ext.pyd",
    ".\\lib\\cctbx_asymmetric_map_ext.pyd",
    ]


SUB_PACKAGES_TO_COPY = [
    "\\..\\modules\\cctbx_project\\iotbx\\cif",
    "\\..\\modules\\cctbx_project\\boost_adaptbx\\boost",
    "\\..\\modules\\cctbx_project\\scitbx\\stl",
    "\\..\\modules\\cctbx_project\\scitbx\\random",
    "\\..\\modules\\cctbx_project\\cctbx\\xray\\structure_factors",

]

def create_dir(rel_path):
    absdir = os.path.abspath(OUTPUT_DIR + rel_path)
    if not os.path.exists(absdir):
        print("Creating directory:", absdir)
        os.makedirs(absdir)
    return absdir

def copy_file(rel_path):
    absinfile = os.path.abspath(INPUT_DIR + rel_path)
    absoutfile = os.path.abspath(OUTPUT_DIR + rel_path)
    if not os.path.exists(absoutfile):
        print("Copying file from: {0}\n\
               to: {1}".format(
            absinfile, absoutfile))
        shutil.copy2(absinfile, absoutfile)

def main():
    if not os.path.exists(OUTPUT_DIR):
        print("Creating output directory:", OUTPUT_DIR)
        os.makedirs(OUTPUT_DIR)

    for dir in DIRS_TO_CREATE:
        create_dir(dir)

    for package in PACKAGES_TO_CREATE:
        abspackage = create_dir(package)
        copy_file(package + '\\__init__.py')

    for file in FILES_TO_COPY:
        copy_file(file)

        # absinit = os.path.join(abspackage, '__init__.py')
        # if not os.path.exists(absinit):
        #     print("Copying into file:", absinit)
        #     shutil.copy2(
        #         os.path.abspath(INPUT_DIR + package) + '\\__init__.py',
        #         os.path.join(OUTPUT_DIR, absinit))

    for subpackage in SUB_PACKAGES_TO_COPY:
        abssubpackage = os.path.abspath(OUTPUT_DIR + subpackage)
        if not os.path.exists(abssubpackage):
            print("Copying into package:", abssubpackage)
            shutil.copytree(
                os.path.abspath(INPUT_DIR + subpackage), abssubpackage)

if __name__ == "__main__":
    main()

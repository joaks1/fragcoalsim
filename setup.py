import os
import re
from setuptools import setup, Extension

#######################################################################
# Hacky way of building mspi in place
#
# import os
# import sys
# import subprocess
#
# BASE_DIR = os.path.abspath(os.path.dirname(__file__))
# MSPI_DIR = os.path.join(BASE_DIR, "fragcoalsim", "libmspi")
# DEV_MODE = False
# INSTALL_MODE = False
# if 'develop' in sys.argv:
#     DEV_MODE = True
# if 'install' in sys.argv:
#     INSTALL_MODE = True

# if INSTALL_MODE:
#     raise Exception("Please install as dev (python setup.py develop)")

# def compile_mspi():
#     current_dir = os.getcwd()
#     os.chdir(MSPI_DIR)
#     cmd = "gcc -fPIC -shared -o libmspi.so mspi.c streec.c rand2t.c -lm"
#     args = cmd.split()
#     try:
#         process = subprocess.Popen(
#                 args,
#                 stdout = subprocess.PIPE,
#                 stderr = subprocess.PIPE,
#                 shell = False)
#         exit_code = process.wait()
#     except:
#         os.chdir(current_dir)
#         sys.stderr.write("ERROR: Problem compiling mspi\n")
#         raise
#     os.chdir(current_dir)

# compile_mspi()
#######################################################################


def _get_version():
    _version_pattern = re.compile(r"^\s*__version__\s*=\s*['\"](?P<version>\S+)['\"]\s*$")
    path = os.path.join(
            os.path.dirname(__file__),
            "fragcoalsim",
            "__init__.py")
    with open(path, "r") as stream:
        for line in stream:
            m = _version_pattern.match(line)
            if m:
                return m.group("version")
    return None


setup(
        name = "fragcoalsim",
        version = _get_version(),
        description = "Package for coalescent simulations of fragmentation",
        author = "Jamie Oaks",
        author_email = "joaks1@gmail.com",
        license = "GPL",
        packages = [
                "fragcoalsim",
                "fragcoalsim.cli",
                "fragcoalsim.test",
                ],
        include_package_data = True,
        zip_safe = False,
        test_suite = "fragcoalsim.test.get_unittest_suite",
        # test_suite = "fragcoalsim.test",
        # test_loader = "unittest:TestLoader",
        install_requires = [
                # 'matplotlib'
        ],
        ext_modules = [
                Extension(name = "fragcoalsim.libmspi.libmspi",
                        sources = [
                                "fragcoalsim/libmspi/mspi.c",
                                "fragcoalsim/libmspi/streec.c",
                                "fragcoalsim/libmspi/rand2t.c",
                                ],
                        libraries = ["m"],
                        language = "c")
                        ],
        entry_points = {
            'console_scripts': [
                'pitracer = fragcoalsim.cli.pitracer:main',
            ],
        },
)

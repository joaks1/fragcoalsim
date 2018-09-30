import os
import sys
import subprocess
from setuptools import setup

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
MSPI_DIR = os.path.join(BASE_DIR, "fragcoalsim", "mspi")
DEV_MODE = False
INSTALL_MODE = False
if 'develop' in sys.argv:
    DEV_MODE = True
if 'install' in sys.argv:
    INSTALL_MODE = True

if INSTALL_MODE:
    raise Exception("Please install as dev (python setup.py develop)")

def compile_mspi():
    current_dir = os.getcwd()
    os.chdir(MSPI_DIR)
    cmd = "gcc -fPIC -shared -o mspi.so mspi.c streec.c rand2t.c -lm"
    args = cmd.split()
    try:
        process = subprocess.Popen(
                args,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                shell = False)
        exit_code = process.wait()
    except:
        os.chdir(current_dir)
        sys.stderr.write("ERROR: Problem compiling mspi\n")
        raise
    os.chdir(current_dir)

compile_mspi()


from fragcoalsim import __version__

setup(
        name = "fragcoalsim",
        version = __version__,
        description = "Package for coalescent simulations of fragmentation",
        author = "Jamie Oaks",
        author_email = "joaks1@gmail.com",
        license = "GPL",
        packages = ["fragcoalsim"],
        include_package_data = True,
        zip_safe = False,
        test_suite = "fragcoalsim.test.get_unittest_suite",
        # test_suite = "fragcoalsim.test",
        # test_loader = "unittest:TestLoader",
        install_requires = [
            # 'matplotlib'
        ],
        entry_points = {
            'console_scripts': [
                'pitracer = fragcoalsim.cli.pitracer:main',
            ],
        },
)

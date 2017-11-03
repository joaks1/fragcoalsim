from setuptools import setup

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

from setuptools import setup, find_packages

def get_version(path):
    """ Parse the version number variable __version__ from a script. """
    import re
    string = open(path).read()
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str


setup(
    name = 'CNVRepeat',
    version = get_version("src/CNVRepeat/__init__.py"),

    packages = find_packages('src'),
    package_dir = {"": "src"},

    entry_points = {
        'console_scripts' : ["CNVRepeat = CNVRepeat.main:main"]
    },

    install_requires = ["h5py", "networkx", "pysam>=0.10.0", "scipy", "traitlets", "ipython-cluster-helper", "pyfaidx", "pandas"],

)

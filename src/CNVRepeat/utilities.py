import collections
import errno
import os
import re
import string 
import sys

class BinaryNotFoundError(Exception):
    pass

def ensure_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise

class cd: 
  def __init__(self, newPath):
    self.newPath = newPath

  def __enter__(self):
    self.savedPath = os.getcwd()
    os.chdir(self.newPath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self.savedPath)
    
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]    


def get_key(options_dict, key, type_=basestring, default="error", error_msg="configuration"):
    if default == "error" and key not in options_dict:
        print "CONFIG ERROR: {} key '{}' is missing".format(error_msg, key)
        sys.exit(1)
    value = options_dict.get(key, default)
    if type_ is not None and not isinstance(value, type_):
        print "CONFIG ERROR: {} key '{}' should be type '{}', not '{}'".format(
            error_msg, key, type_.__name__, type(value).__name__)
        sys.exit(1)
    return value


comp = string.maketrans('ATCGatcg','TAGCtagc')
def revcomp(seq):
    return seq[::-1].translate(comp)

def cpu_count_physical():
    """
    tries to get the number of physical (ie not virtual) cores
    """
    try:
        import psutil
        return psutil.cpu_count(logical=False)
    except:
        import multiprocessing
        return multiprocessing.cpu_count()

def check_memory(logger, min_memory=16):
    try:
        import psutil
        physical_mem_gb = psutil.virtual_memory().total / (1000.**3)
        if physical_mem_gb < min_memory:
            logger.log("WARNING: GROC-SVs typically requires ~16 GB of memory to run; "
                       "you appear to have only {:.1f}GB".format(physical_mem_gb))
    except:
        pass


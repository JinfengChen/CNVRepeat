import collections
import json
import os
import sys

from CNVRepeat import utilities
from CNVRepeat.reference import Reference
from CNVRepeat.utilities import get_key


class ClusterSettings(object):
    def __init__(self):
        self.cluster_type = "local"
        self.processes = utilities.cpu_count_physical()
        self.cluster_options = {}

    @staticmethod
    def deserialize(options_dict):
        settings = ClusterSettings()

        if "processes" in options_dict:
            settings.processes = options_dict["processes"]
        if "cluster_type" in options_dict:
            settings.cluster_type = options_dict["cluster_type"]
        if "cluster_options" in options_dict:
            settings.cluster_options = options_dict["cluster_options"]

        return settings

    def serialize(self):
        return {
            "processes": self.processes,
            "cluster_type": self.cluster_type,
            "cluster_options": self.cluster_options
        }

    
class Options(object):
    def __init__(self, options_path, debug=False):
        self.options_path = options_path

        self.ref_fasta = None
        self.gtf       = None
        self.bed       = None
        self.bam       = None
        self.fastq1    = None
        self.fastq2    = None
        self.repeat    = None
        self.black     = None
        self.binaries  = {}
        self._reference = None
        self._constants = None
        self.output = None

        self.cluster_settings = ClusterSettings()

        self.debug  = debug
        self.method = None 
        self.random_dna_length = None
        self.random_dna_number = None

    def serialize(self, ):

        d = {"ref_fasta": self.ref_fasta,
             "gtf": self.gtf,
             "bed": self.bed,
             "bam": self.bam,
             "fastq1": self.fastq1,
             "fastq2": self.fastq2,
             "repeat": self.repeat,
             "black": self.black,
             "output": self.output,
             "cluster_settings": self.cluster_settings.serialize(),
             "binaries": self.binaries
        }

        return d

    @staticmethod
    def deserialize(options_dict, options_path):
        options = Options(options_path)
        options.ref_fasta = get_key(options_dict, "ref_fasta")
        options.gtf       = get_key(options_dict, "gtf")
        options.bed       = get_key(options_dict, "bed")
        options.bam       = get_key(options_dict, "bam")
        options.fastq1    = get_key(options_dict, "fastq1")
        options.fastq2    = get_key(options_dict, "fastq2")
        options.repeat    = get_key(options_dict, "repeat")
        options.black     = get_key(options_dict, "black") 
        options.binaries  = get_key(options_dict, "binaries", dict, default={})
        options.output    = get_key(options_dict, "output")

        if not os.path.exists(options.output):
            os.mkdir(options.output)

        options.cluster_settings = ClusterSettings.deserialize(
            options_dict.get("cluster_settings", {}))

        return options

    @property
    def output_dir(self):
        return self.output   
 
    @property
    def results_dir(self):
        return os.path.join(self.output_dir, "results")

    @property
    def working_dir(self):
        return os.path.join(self.output_dir, "working")

    @property
    def log_dir(self):
        return os.path.join(self.output_dir, "logs")

    @property
    def reference(self):
        if self._reference is None:
            self._reference = Reference(self.ref_fasta, self.debug)
        return self._reference

    def binary(self, name):
        """
        Checks to see if a path has been specified for an external binary,
        otherwise just return the name of the binary to try running it
        if it's in $PATH
        """
        
        bin_path = self.binaries.get(name, name)
        if utilities.which(bin_path) is None:
            raise utilities.BinaryNotFoundError(
                "Failed to locate binary '{}'; please make sure it is in ".format(name) + 
                "your $PATH or add it to the configuration.json file")
        return bin_path


    @property
    def debug(self):
        return self._debug
    
    @debug.setter
    def debug(self, mode=True):
        self._reference = None
        self._debug = mode

    def __str__(self):
        d = self.serialize()
        d["debug"] = self.debug
        return json.dumps(d, sort_keys=True, indent=4)

    def __getstate__(self):
        """
        allows pickling of Options instances, necessary for ipyparallel
        """
        state = self.__dict__.copy()
        state["_reference"] = None
        state["_constants"] = None

        return state    

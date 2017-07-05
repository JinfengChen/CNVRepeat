import collections
import json
import os
import sys

from CNVRepeat import utilities
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
        self.bam       = None
        self.repeat    = None
        self.binaries  = {}
        self._reference = None
        self._constants = None
        self.output = None

        self.cluster_settings = ClusterSettings()

        self.debug = debug


    def serialize(self, ):

        d = {"ref_fasta": self.ref_fasta,
             "gtf": self.gtf,
             "bam": self.bam,
             "repeat": self.repeat,
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
        options.bam       = get_key(options_dict, "bam")
        options.repeat    = get_key(options_dict, "repeat") 
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


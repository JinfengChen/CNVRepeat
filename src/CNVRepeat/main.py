#the package use modules from grocsvs (https://github.com/grocsvs/grocsvs) as basis for designing the pipeline.
from __future__ import print_function

import argparse
import collections
import json
import logging
import sys

from CNVRepeat import log
from CNVRepeat import options as opts
from CNVRepeat import pipeline
from CNVRepeat import analysis

def load_config(config_path):
    try:
        config = json.load(open(config_path))
    except ValueError as err:
        print("Error parsing configuration file '{}': '{}'\n  Check that this is a properly formatted JSON file!".format(config_path, err))
        sys.exit(1)

    options = opts.Options.deserialize(config, config_path)
    return options

def run(options):
    analysis_steps = prepare_analysis(options)
    runner = pipeline.Runner(options)

    print("Running")
    for analysis_name, analysis_step in analysis_steps.items():
        print ('Running analysis: "{}"'.format(analysis_name))
        runner.run_stage(analysis_step, analysis_name)

def prepare_analysis(options):
    analysis_steps = collections.OrderedDict()
    if not os.path.exists(options.bed):
        analysis_steps["Single Copy Exon"] = analysis.single_copy_exon.SingleCopyExonStep

    return analysis_steps

def main():
    parser = argparse.ArgumentParser(description="CNVRepeat: estimate copy number of repeat sequence in the genome")
    parser.add_argument("--config", help="Path to configuration.json file")
    parser.add_argument("--local", action="store_true", help="run job locally in multiprocess mode")
    parser.add_argument("--scheduler", help="run job using scheduler, SLURM, SGE, PBS/Torque")
    parser.add_argument("--cpu", default=1, help="number of cpu")
    parser.add_argument("--debug", action="store_true", help="run in debug mode")
    args = parser.parse_args()
     
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    options = load_config(args.config)
    options.debug = args.debug

    log.log_command(options, sys.argv)

    run(options)



if __name__ == '__main__':
    main()

from __future__ import print_function

import argparse
import collections
import json
import logging
import sys

from CNVRepeat import log
from CNVRepeat import options as opts


def load_config(config_path):
    try:
        config = json.load(open(config_path))
    except ValueError as err:
        print("Error parsing configuration file '{}': '{}'\n  Check that this is a properly formatted JSON file!".format(config_path, err))
        sys.exit(1)

    options = opts.Options.deserialize(config, config_path)
    return options

def run(options):
    #stages = get_stages()
    #runner = pipeline.Runner(options)

    #for stage_name, stage in stages.items():
    #    print ('Running state: "{}"'.format(stage_name))
    #    runner.run_stage(stage, stage_name)
    


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

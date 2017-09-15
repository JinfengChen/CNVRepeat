#the package use modules from grocsvs (https://github.com/grocsvs/grocsvs) as basis for designing the pipeline.
from __future__ import print_function

import argparse
import collections
import json
import logging
import sys
import os

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
   
    if options.method == 'single_copy_exon':
        if not os.path.exists(options.bed) or not os.path.splitext(options.bed)[1] == '.bed':
            analysis_steps["Single Copy Exon"]      = analysis.single_copy_exon.SingleCopyExonStep
        analysis_steps["Genome Coverage Estimator"] = analysis.estimate_genome_coverage_bed.EstimateGenomeCoverageStep
        analysis_steps["Genome Coverage Merger"]    = analysis.estimate_genome_coverage_bed.CombineGenomeCoverageStep
    elif options.method == 'single_copy_exon':
        if not os.path.exists(options.bed) or not os.path.splitext(options.bed)[1] == '.bed':
            analysis_steps["Random Region"]      = analysis.single_copy_exon.RandomRegionStep
        analysis_steps["Genome Coverage Estimator"] = analysis.estimate_genome_coverage_bed.EstimateGenomeCoverageStep
        analysis_steps["Genome Coverage Merger"]    = analysis.estimate_genome_coverage_bed.CombineGenomeCoverageStep
    elif options.method == 'goleft':
        analysis_steps["Genome Coverage Estimator Goleft"] = analysis.estimate_genome_coverage_bed.EstimateGenomeCoverageGoleftStep
    analysis_steps["Repaet Coverage Estimator"] = analysis.estimate_repeat_coverage.EstimateRepeatCoverageStep
    analysis_steps["Repeat Copy Number"]        = analysis.estimate_repeat_copy_number.EstimateRepeatCopyNumberStep

    return analysis_steps

def main():
    parser = argparse.ArgumentParser(description="CNVRepeat: estimate copy number of repeat sequence in the genome")
    parser.add_argument("--config", help="Path to configuration.json file")
    parser.add_argument("--local", action="store_true", help="run job locally in multiprocess mode")
    parser.add_argument("--scheduler", help="run job using scheduler, SLURM, SGE, PBS/Torque")
    parser.add_argument("--cpu", default=1, help="number of cpu")
    parser.add_argument("--method", default='goleft', help="method for estimation of genome coverage: goleft, single_copy_exon, random_region")
    parser.add_argument("--random_dna_length", default=1000, help="length of DNA for random selection of method random_region")
    parser.add_argument("--random_dna_number", default=100000, help="number of DNA for random selection of method random_region")
    parser.add_argument("--debug", action="store_true", help="run in debug mode")
    args = parser.parse_args()
     
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    options = load_config(args.config)
    options.debug  = args.debug
    options.method = args.method
    options.random_dna_length = args.random_dna_length
    options.random_dna_number = args.random_dna_number

    log.log_command(options, sys.argv)

    run(options)



if __name__ == '__main__':
    main()

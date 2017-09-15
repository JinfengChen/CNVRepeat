import numpy as np
import os
import re
import sys
import subprocess
import pandas
from collections import defaultdict

from CNVRepeat import step 
from CNVRepeat.analysis import single_copy_exon

CHUNKSIZE = 1e7

def chunks_for_chrom(options, chrom):
    return int(np.ceil(options.reference.chrom_lengths[chrom]/CHUNKSIZE))

class EstimateGenomeCoverageGoleftStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        yield EstimateGenomeCoverageGoleftStep(options)

    def __init__(self, options):
        self.options = options

    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])

    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        genomecov_name = "genomecov.goleft.covstats"

        paths = {
            "genomecov": os.path.join(directory, genomecov_name)
        }

        return paths
    

    def run(self):
        outpaths = self.outpaths(final=False)
        self.logger.log("Calculating genome coverage using goleft-covstats ...")
        command  = '{} covstats {} > {}'.format(self.options.binaries['goleft'], self.options.bam, outpaths['genomecov'])
        cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) 
        cmd_error = cmd.wait()
        if cmd_error != 0:
            self.logger.log("goleft genome coverage error code: {}\n{}".format(command, cmd_error)) 

class EstimateGenomeCoverageStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        for chrom in options.reference.chroms:
            for chunk in range(chunks_for_chrom(options, chrom)):
                yield EstimateGenomeCoverageStep(options, chrom, chunk)

    def __init__(self, options, chrom, chunk):
        self.options = options
        self.chrom = chrom
        self.chunk = chunk
    
    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.chrom,
                         str(self.chunk)])

    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        bed_name      = "genomecov.{}.{}.bed".format(self.chrom, self.chunk)

        bedgraph_name = "genomecov.{}.{}.bedgraph".format(self.chrom, self.chunk)

        paths = {
            "exonbed"  : os.path.join(directory, bed_name),
            "genomebedgraph": os.path.join(directory, bedgraph_name)
        }

        return paths

    def run(self):
        outpaths = self.outpaths(final=False)
        self.logger.log("Calculating genome coverage ...")
     
        chrom_length = self.options.reference.chrom_lengths[self.chrom]
        start = int(self.chunk*CHUNKSIZE)
        end = int(min((self.chunk+1)*CHUNKSIZE, chrom_length))    

        self.logger.log("Running chunk: {}:{:,}-{:,}".format(self.chrom, start, end))

        self.run_samtools_bedcov(self.options, self.chrom, start, end, outpaths["exonbed"], outpaths["genomebedgraph"])
        
    def run_samtools_bedcov(self, options, chrom, start, end, bed, bedgraph): 
        if not os.path.exists(options.bed):
            if self.options.method == 'single_copy_exon':
                single_exon_bed_out = single_copy_exon.SingleCopyExonStep(options).outpaths(final=True)["bed_out"]
                if not os.path.exists(single_exon_bed_out):
                    self.logger.log("bed file error: single exon bed file is not found: exiting pipeline \n")
                    sys.exit(1)
                else:
                    options.bed = single_exon_bed_out  
            elif self.options.method == 'random_region':
                random_region_bed_out = single_copy_exon.RandomRegionStep(options).outpaths(final=True)["bed_out"]
                if not os.path.exists(random_region_bed_out):
                    self.logger.log("bed file error: random region bed file is not found: exiting pipeline \n") 
                    sys.exit(1)
                else:
                    options.bed = random_region_bed_out          

        commands = []
        commands.append('awk \'$1~/^{}$/ && $2>={} && $3<={}\' {} > {}'.format(chrom, start, end, options.bed, bed))
        commands.append('{} bedcov {} {} > {}'.format(options.binaries['samtools'], bed, options.bam, bedgraph))
        for command in commands:
            cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
            cmd_error = cmd.wait()
            if cmd_error != 0:
                self.logger.log("genome coverage error code: {}\n{}".format(command, cmd_error))


class CombineGenomeCoverageStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        yield CombineGenomeCoverageStep(options)
    
    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])

    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        genomecov_name = "genomecov.bedgraph"

        paths = {
            "genomecov": os.path.join(directory, genomecov_name)
        }

        return paths

    def run(self):
        outpaths = self.outpaths(final=False)
        self.logger.log("Merging genome coverage ...")
   
        genomecov_bedgraph = []
        ofile = open(outpaths["genomecov"], 'w')
        header = ['chrom','start','end','strand','bedcov','genomecov']
        print >> ofile, '\t'.join(header)
        for i, inpath in enumerate(self.get_input_paths()):
            try:
                genomecov_bedgraph_temp = pandas.read_csv(inpath, sep="\t", header=None, names=header)
                genomecov_bedgraph_temp['genomecov'] = genomecov_bedgraph_temp['bedcov']/(genomecov_bedgraph_temp['end']-genomecov_bedgraph_temp['start']+1)
                genomecov_bedgraph_temp.to_csv(ofile, sep="\t", columns=header, header=False, index=False)
            except pandas.io.common.EmptyDataError:
                self.logger.log("No genome coverage found in {}; skipping".format(inpath))
        ofile.close() 

    def get_input_paths(self):
        paths = []
        for chrom in self.options.reference.chroms:
            for chunk in range(chunks_for_chrom(self.options, chrom)):
                input_step = EstimateGenomeCoverageStep(self.options, chrom, chunk)
                paths.append(input_step.outpaths(final=True)["genomebedgraph"])

        return paths 

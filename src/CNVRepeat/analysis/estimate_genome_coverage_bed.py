import numpy as np
import os
import re
import sys
import subprocess
from collections import defaultdict

from CNVRepeat import step 

CHUNKSIZE = 1e7

def chunks_for_chrom(options, chrom):
    return int(np.ceil(options.reference.chrom_lengths[chrom]/CHUNKSIZE))

class EstimateGenomeCoverageStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        for chrom in options.reference.chroms:
            for chunk in range(chunks_for_chrom(options, chrom)):
                yield EstimateGenomeCoverageStep(options, chrom, chunk)
        #for chrom in ["chr1", "chr2"]:
        #    for chunk in range(chunks_for_chrom(options, chrom)):
            #for chunk in range(1):
        #        yield EstimateGenomeCoverageStep(options, chrom, chunk)

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
        options.bed
        commands = []
        commands.append('awk \'$1~/^{}$/ && $2>={} && $3<={}\' {} > {}'.format(chrom, start, end, options.bed, bed))
        commands.append('{} bedcov {} {} > {}'.format(options.binaries['samtools'], bed, options.bam, bedgraph))
        for command in commands:
            print(command)
            cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
            cmd_error = cmd.wait()
            if cmd_error != 0:
                self.logger.log("genome coverage error code: {}\n{}".format(command, cmd_error))


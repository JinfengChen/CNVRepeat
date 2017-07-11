import os
import re
import sys
import subprocess
import pandas
import numpy as np
from collections import defaultdict

from CNVRepeat import step 
from CNVRepeat.reference import Reference
from CNVRepeat.analysis import estimate_genome_coverage_bed 
from CNVRepeat.analysis import estimate_repeat_coverage

class EstimateRepeatCopyNumberStep(step.StepChunk):
   
    @staticmethod
    def get_steps(options):
        yield EstimateRepeatCopyNumberStep(options)

    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])  
 
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir
        paths = {
            "repeat_cnv" : os.path.join(directory, '{}.repeat_cnv.csv'.format(os.path.splitext(os.path.split(self.options.gtf)[1])[0])),
        }

        return paths

    def run(self):
        self.repeat_cnv = self.outpaths(final=False)["repeat_cnv"]
        repeat_depth    = estimate_repeat_coverage.EstimateRepeatCoverageStep(self.options).outpaths(final=True)["repeat_depth_out"]
        single_copy_exon_bedgraph = estimate_genome_coverage_bed.CombineGenomeCoverageStep(self.options).outpaths(final=True)["genomecov"]

        self.run_estimate_repeat_copy_number(single_copy_exon_bedgraph, repeat_depth, self.options.repeat, self.repeat_cnv)
       
    def run_estimate_repeat_copy_number(self, single_copy_exon_bedgraph, repeat_depth, repeat, repeat_cnv):
        repeat_seq = Reference(repeat)
        single_copy_exon_cov = self.parse_bedgraph(single_copy_exon_bedgraph)
        repeat_cov           = self.parse_depth(repeat_depth, repeat_seq)
        ofile = open(repeat_cnv, 'w')
        print >> ofile, 'Repeat\tRepeat_Coverage\tSingle_Copy_Exon_Genome_Coverage\tRepeat_Copy_Number'
        for repeat in sorted(repeat_cov.keys()):
            '''print >> ofile, '\t'.join(map(str, [repeat, repeat_cov[repeat], single_copy_exon_cov, repeat_cov[repeat]/single_copy_exon_cov]))'''
            print >> ofile, '{}\t{:.2f}\t{:.2f}\t{:.2f}'.format(repeat, repeat_cov[repeat], single_copy_exon_cov, repeat_cov[repeat]/single_copy_exon_cov)
        ofile.close()

    '''
    chrom   start   end     strand  bedcov  genomecov
    Chr1    2903    3268    +       11508   31.4426229508
    Chr1    3354    3616    +       5743    21.8365019011
    '''
    def parse_bedgraph(self, single_copy_exon_bedgraph):
        data_bedgraph  = pandas.read_csv(single_copy_exon_bedgraph, sep="\t", header=0) 
        genomecov_mean = data_bedgraph["genomecov"].mean()
        return genomecov_mean
    
    '''
    mPing   1       2158
    mPing   2       2184
    '''
    def parse_depth(self, repeat_depth, repeat_seq):
        data_depth     = pandas.read_csv(repeat_depth, sep="\t", header=None)
        repeat_dict    = defaultdict(lambda : list())
        repeat_cov     = defaultdict(lambda : float())
   
        for i in range(len(data_depth)):
            repeat_dict[data_depth[0][i]].append(data_depth[2][i])    

        for repeat in repeat_dict.keys():
            repeat_cov[repeat] = np.mean(repeat_dict[repeat])
        return repeat_cov

 

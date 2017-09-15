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
            "repeat_cnv" : os.path.join(directory, '{}.repeat_cnv.csv'.format(os.path.splitext(os.path.split(self.options.ref_fasta)[1])[0])),
        }

        return paths

    def run(self):
        self.repeat_cnv = self.outpaths(final=False)["repeat_cnv"]
        repeat_depth    = estimate_repeat_coverage.EstimateRepeatCoverageStep(self.options).outpaths(final=True)["repeat_depth_out"]
        black_list = defaultdict(lambda : int())
        try:
            self.logger.log("black list {}".format(self.options.black))
            black_list = self.parse_black_list(self.options.black)
            self.logger.log("parse black list: {}".format(self.options.black)) 
        except:
            self.logger.log("no black list found")
            pass

        if self.options.method == 'single_copy_exon' or self.options.method == 'random_region':
            region_bedgraph = estimate_genome_coverage_bed.CombineGenomeCoverageStep(self.options).outpaths(final=True)["genomecov"]
            self.run_estimate_repeat_copy_number(region_bedgraph, repeat_depth, self.options.repeat, self.repeat_cnv, black_list, self.options.method)
        elif self.options.method == 'goleft':
            goleft_covstats = estimate_genome_coverage_bed.EstimateGenomeCoverageGoleftStep(self.options).outpaths(final=True)["genomecov"]      
            self.run_estimate_repeat_copy_number_goleft(goleft_covstats, repeat_depth, self.options.repeat, self.repeat_cnv, black_list) 

    def run_estimate_repeat_copy_number_goleft(self, goleft_covstats, repeat_depth, repeat, repeat_cnv, black_list):
        repeat_seq = Reference(repeat)
        goleft_covstats_cov = self.parse_goleft_covstats(goleft_covstats)
        repeat_cov           = self.parse_depth(repeat_depth, repeat_seq, black_list)
        ofile = open(repeat_cnv, 'w')
        print >> ofile, 'Repeat\tRepeat_Coverage\tGoleft_Covstats_Genome_Coverage\tRepeat_Copy_Number'
        for repeat in sorted(repeat_cov.keys()):
            '''print >> ofile, '\t'.join(map(str, [repeat, repeat_cov[repeat], single_copy_exon_cov, repeat_cov[repeat]/single_copy_exon_cov]))'''
            print >> ofile, '{}\t{:.2f}\t{:.2f}\t{:.2f}'.format(repeat, repeat_cov[repeat], goleft_covstats_cov, repeat_cov[repeat]/goleft_covstats_cov)
        ofile.close() 

    def parse_goleft_covstats(self, goleft_covstats):
        data_bedgraph  = pandas.read_csv(goleft_covstats, sep="\t", header=0) 
        genomecov_mean = data_bedgraph["coverage"][0]
        return float(genomecov_mean)

    def run_estimate_repeat_copy_number(self, region_bedgraph, repeat_depth, repeat, repeat_cnv, black_list, method):
        repeat_seq = Reference(repeat)
        region_cov = self.parse_bedgraph(region_bedgraph)
        repeat_cov           = self.parse_depth(repeat_depth, repeat_seq, black_list)
 
        method_title = '' 
        if method == 'single_copy_exon':
            method_title = 'Single_Copy_Exon_Genome_Coverage'
        elif method == 'random_region':
            method_title = 'Random_Region_Genome_Coverage'

        ofile = open(repeat_cnv, 'w')
        print >> ofile, 'Repeat\tRepeat_Coverage\t{}\tRepeat_Copy_Number'.format(method_title)
        for repeat in sorted(repeat_cov.keys()):
            '''print >> ofile, '\t'.join(map(str, [repeat, repeat_cov[repeat], region_cov, repeat_cov[repeat]/region_cov]))'''
            print >> ofile, '{}\t{:.2f}\t{:.2f}\t{:.2f}'.format(repeat, repeat_cov[repeat], region_cov, repeat_cov[repeat]/region_cov)
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
    def parse_depth(self, repeat_depth, repeat_seq, black_list):
        data_depth     = pandas.read_csv(repeat_depth, sep="\t", header=None)
        repeat_dict    = defaultdict(lambda : list())
        repeat_cov     = defaultdict(lambda : float())
   
        for i in range(len(data_depth)):
            self.logger.log('repeat: {} {}'.format(data_depth[0][i], data_depth[1][i]))
            if black_list.has_key(data_depth[0][i]):
                if not black_list[data_depth[0][i]].has_key(str(data_depth[1][i])):
                    self.logger.log('site not in black list')
                    repeat_dict[data_depth[0][i]].append(float(data_depth[2][i]))
                else:
                    self.logger.log('excluding site from black list: {} {}'.format(data_depth[0][i], data_depth[1][i]))
            else:
                self.logger.log('repeat not in black list')
                repeat_dict[data_depth[0][i]].append(float(data_depth[2][i]))    

        for repeat in repeat_dict.keys():
            repeat_cov[repeat] = np.mean(repeat_dict[repeat])
        return repeat_cov

    '''
    pong	1	260
    pong	3260	5166
    '''
    def parse_black_list(self, black_bed):
        data_dict  = defaultdict(lambda : defaultdict(lambda : int()))
        data_black = np.loadtxt(black_bed, dtype=str)
        for line in data_black:
            self.logger.log('line: {}'.format(line))
            for i in range(int(line[1]), int(line[2])+1):
                data_dict[line[0]][str(i)] = 1
                self.logger.log('site: {}'.format(str(i)))
        return data_dict


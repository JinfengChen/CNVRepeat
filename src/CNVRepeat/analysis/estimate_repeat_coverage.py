import os
import re
import sys
import subprocess
from collections import defaultdict

from CNVRepeat import step 
from CNVRepeat.reference import Reference

class EstimateRepeatCoverageStep(step.StepChunk):
   
    @staticmethod
    def get_steps(options):
        yield EstimateRepeatCoverageStep(options)

    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])  
 
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir
        paths = {
            "bwa_bam" : os.path.join(directory, '{}.bwa_sorted.bam'.format(os.path.splitext(os.path.split(self.options.ref_fasta)[1])[0])),
            "repeat_depth_out" : os.path.join(directory, '{}.bwa_sorted.depth'.format(os.path.splitext(os.path.split(self.options.ref_fasta)[1])[0]))
        }

        return paths

    def run(self):
        self.bwa_bam          = self.outpaths(final=False)["bwa_bam"]
        self.repeat_depth_out = self.outpaths(final=False)["repeat_depth_out"]
        self.run_map2repeat_bwa(self.bwa_bam, self.repeat_depth_out)

    def run_map2repeat_bwa(self, bwa_bam, repeat_depth_out):
        bwa_sam         = '{}.sam'.format(os.path.splitext(bwa_bam)[0])
        bwa_filter_sam  = '{}.filter.sam'.format(os.path.splitext(bwa_bam)[0])
        if not os.path.exists('{}.bwt'.format(self.options.repeat)):
            command = '{} index {}'.format(self.options.binaries['bwa'], self.options.repeat)
            cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
            cmd_error = cmd.wait()
            if cmd_error != 0:
                self.logger.log("bwa index error code: {}\n{}".format(command, cmd_error))
            
        if not os.path.exists(bwa_bam):
            command = '{} mem -t {} -O2,2 -a -Y -k 15 -T 10 {} {} {} | {} view -Sh -F 4 - > {}'.format(self.options.binaries['bwa'], self.options.cluster_settings.processes, self.options.repeat, self.options.fastq1, self.options.fastq2, self.options.binaries['samtools'], bwa_sam)
            cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
            cmd_error = cmd.wait()
            if cmd_error != 0:
                self.logger.log("bwa mapping error code: {}\n{}".format(command, cmd_error))
            
            self.parse_sam(bwa_sam, bwa_filter_sam)
            
            command = '{} view -Shb -F 4 {} | {} sort - -o {}'.format(self.options.binaries['samtools'], bwa_filter_sam, self.options.binaries['samtools'], bwa_bam)
            cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
            cmd_error = cmd.wait()
            if cmd_error != 0:
                self.logger.log("bam sorting error code: {}\n{}".format(command, cmd_error))

        command = '{} depth -d 80000 -Q 20 {} > {}'.format(self.options.binaries['samtools'], bwa_bam, repeat_depth_out)
        cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
        cmd_error = cmd.wait()
        if cmd_error != 0:
            self.logger.log("samtools depth error code: {}\n{}".format(command, cmd_error)) 

    def parse_sam(self, sam, filter_sam):
        repeat_seq = Reference(self.options.repeat) 
        ofile = open(filter_sam, 'w')
        right_clip  = re.compile(r'^(\d+)M\d+S$')
        left_clip   = re.compile(r'^\d+S\d+M$')
        perfect     = re.compile(r'^\d+M$')
        mm          = re.compile(r'NM:i:(\d+)')
        xs          = re.compile(r'XS:i:(\d+)')
        with open (sam, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if line.startswith(r'@'):
                    print >> ofile, line
                else:
                    unit = re.split(r'\t',line)
                    mm_m = mm.search(line)
                    mm_n = mm_m.groups(0)[0] if mm_m else 0
                    xs_m = xs.search(line)
                    xs_n = xs_m.groups(0)[0] if xs_m else 100
                    if int(mm_n) > 2:
                        continue
                    if perfect.search(unit[5]):
                        print >> ofile, line
                    elif left_clip.search(unit[5]) and int(unit[3]) < 10:
                        print >> ofile, line
                    elif right_clip.search(unit[5]):
                        match_len  = right_clip.search(unit[5]).groups(0)[0]
                        repeat_len = len(repeat_seq.fasta[unit[2]])
                        if int(unit[3]) + int(match_len) >= repeat_len - 10:
                            print >> ofile, line
        ofile.close()


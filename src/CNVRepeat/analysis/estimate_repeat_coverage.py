import os
import re
import sys
import subprocess
from collections import defaultdict

from CNVRepeat import step 

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
            "bwa_bam" : os.path.join(directory, '{}.bwa_sorted.bam'.format(os.path.splitext(os.path.split(self.options.gtf)[1])[0])),
            "repeat_depth_out" : os.path.join(directory, '{}.bwa_sorted.depth'.format(os.path.splitext(os.path.split(self.options.gtf)[1])[0]))
        }

        return paths

    def run(self):
        self.bwa_bam          = self.outpaths(final=False)["bwa_bam"]
        self.repeat_depth_out = self.outpaths(final=False)["repeat_depth_out"]
        self.run_map2repeat_bwa(self.bwa_bam, self.repeat_depth_out)

    def run_map2repeat_bwa(self, bwa_bam, repeat_depth_out):
        commands = []
        if not os.path.exists('{}.bwt'.format(self.options.repeat)):
            commands.append('{} index {}'.format(self.options.binaries['bwa'], self.options.repeat))
        if not os.path.exists(bwa_bam):
            commands.append('{} mem -t {} -O2,2 -a -Y -k 15 -T 10 {} {} {} | {} view -Shb -F 4 - | {} sort - -o {}'.format(self.options.binaries['bwa'], self.options.cluster_settings.processes, self.options.repeat, self.options.fastq1, self.options.fastq2, self.options.binaries['samtools'], self.options.binaries['samtools'], bwa_bam))
        commands.append('{} depth -Q 30 {} > {}'.format(self.options.binaries['samtools'], bwa_bam, repeat_depth_out))
        for command in commands:
            print(command)
            cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
            cmd_error = cmd.wait()
            if cmd_error != 0:
                self.logger.log("bwa mapping error code: {}\n{}".format(command, cmd_error)) 


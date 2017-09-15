import os
import re
import sys
import subprocess
from collections import defaultdict

from CNVRepeat import step 

class RandomRegionStep(step.StepChunk):
   
    @staticmethod
    def get_steps(options):
        yield RandomRegionStep(options)

    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])  
 
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir
        paths = {
            "bed_out" : os.path.join(directory, '{}.random_region_l{}_n{}.bed'.format(os.path.splitext(os.path.split(self.options.ref_fasta)[1])[0]), self.options.random_dna_length, self.options.random_dna_number)
        }

        return paths

    def run(self):
        self.bed_out = self.outpaths(final=False)["bed_out"]
        if not os.path.exists(self.bed_out):
            self.run_random_region(self.options.ref_fasta, self.bed_out)

    def run_random_region(self, bed_out)
        commands = []
        fai = '{}.fai'.format(self.options.ref_fasta)
        if not os.path.exists(fai):
            commands.append('{} faidx {}'.format(self.options.binaries['samtools'], self.options.ref_fasta))
        commands.append('{} random -l {} -n {} -seed 123 -g {} | cut -f1,2,3,6 > {}'.format(self.options.binaries['bedtools'], self.options.random_dna_length, self.options.random_dna_number, fai, bed_out))
        
        for command in commands:
            cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
            cmd_error = cmd.wait()
            if cmd_error != 0:
                self.logger.log("random region genome coverage error code: {}\n{}".format(command, cmd_error))

class SingleCopyExonStep(step.StepChunk):
   
    @staticmethod
    def get_steps(options):
        yield SingleCopyExonStep(options)

    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])  
 
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir
        paths = {
            "blast_out" : os.path.join(directory, '{}.blast.m8'.format(os.path.splitext(os.path.split(self.options.gtf)[1])[0])),
            "bed_out" : os.path.join(directory, '{}.single_exon.bed'.format(os.path.splitext(os.path.split(self.options.gtf)[1])[0]))
        }

        return paths

    def run(self):
        self.blast_out = self.outpaths(final=False)["blast_out"]
        self.bed_out = self.outpaths(final=False)["bed_out"]
        if not os.path.exists(self.blast_out):
            self.exon_fa   = self.run_gtf2fasta(self.options.gtf) 
            self.run_formatdb(self.exon_fa)
            self.run_blastn(self.exon_fa)
        if not os.path.exists(self.bed_out):
            self.run_parse_single_exon(self.blast_out, self.options.gtf, self.bed_out)

    def run_gtf2fasta(self, gtf):
        exon_fa = '{}.fa'.format(os.path.splitext(self.options.gtf)[0])
        command = '{} getfasta -fi {} -bed {} -fo {}.fa'.format(self.options.binaries['bedtools'], self.options.ref_fasta, self.options.gtf, os.path.splitext(self.options.gtf)[0])
        cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
        bedtools_error = cmd.wait()
        if bedtools_error != 0:
            self.logger.log("bedtools error code: {}".format(bedtools_error))
        return exon_fa

    def run_formatdb(self, fasta):
        command = '{} -i {} -p F'.format(self.options.binaries['formatdb'], fasta)
        cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
        formatdb_error = cmd.wait()
        if formatdb_error != 0:
            self.logger.log("formatdb error code: {}".format(formatdb_error))

    def run_blastn(self, fasta):
        command = '{} -p blastn -i {} -d {} -o {} -m 8 -e 1-e5 -a {}'.format(self.options.binaries['blastall'], fasta, fasta, self.blast_out, self.options.cluster_settings.processes)
        cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
        blastall_error = cmd.wait()
        if blastall_error != 0:
            self.logger.log("blastall error code: {}".format(blastall_error))

    
    def blastm8_reader(self, infile):
        data = defaultdict(lambda : int())
        with open (infile, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if len(line) > 2: 
                    unit = re.split(r'\t',line)
                    if not unit[0] == unit[1]:
                        if float(unit[10]) <= 1e-5:
                            data[unit[0]] = 1
                            data[unit[1]] = 1
        return data

 

    def run_parse_single_exon(self, blast_out, gtf, bed_out):
        duplicate = self.blastm8_reader(blast_out)
        ofile = open(bed_out, 'w') 
        with open (gtf, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if len(line) > 2:
                    unit = re.split(r'\t',line)
                    exon = '{}:{}-{}'.format(unit[0], str(int(unit[3])-1), unit[4]) 
                    if not duplicate.has_key(exon):
                        print >> ofile, '{}\t{}\t{}\t{}'.format(unit[0], unit[3], unit[4], unit[6])
        ofile.close()

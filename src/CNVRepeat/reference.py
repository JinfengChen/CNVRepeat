import collections
import numpy
import pyfaidx
import re

class Reference(object):
    def __init__(self, path, debug=False):
        self.path = path
        self.fasta = pyfaidx.Fasta(path, as_raw=True)
        self.debug = debug

        self.chroms_to_std_chroms = collections.OrderedDict()
        self.std_chroms_to_chroms = collections.OrderedDict()
        self.chrom_lengths = collections.OrderedDict()

        chrom_re = re.compile(r"(chr|Chr|scaffold_|scaffold)?((\d+)|(X)|(Y))$")

        for chrom in self.fasta.keys():
            std_chrom = self.standardize_chrom(chrom, chrom_re)
            if re.match(chrom_re, std_chrom):
                self.chroms_to_std_chroms[chrom] = std_chrom
                self.std_chroms_to_chroms[std_chrom] = chrom
                self.chrom_lengths[chrom] = len(self.fasta[chrom])

    @property
    def chroms(self):
        if self.debug:
            return ["chr20", "chr21"]

        return self.chrom_lengths.keys()

    def standardize_chrom(self, chrom, chrom_re):
        if "chr" in chrom:
            return chrom
        else:
            chrom_num = 0
            if chrom_re.search(chrom):
                chrom_num = chrom_re.search(chrom).groups(0)[1]
                return "chr{}".format(chrom_num)
            else:
                return "chr{}".format(chrom)

    def compare_chroms(self, chromx, chromy):
        return cmp(self.chroms.index(chromx), self.chroms.index(chromy))

def split_genome(chroms, chroms_to_lengths, chunk_size):
    for chrom in chroms:
        cur_length = chroms_to_lengths[chrom]
        if cur_length <= chunk_size:
            yield (chrom, 0, cur_length)
        else:
            cur_chunk_size = int(cur_length/numpy.ceil(cur_length/float(chunk_size))) + 1

            for i in range(0, chroms_to_lengths[chrom], cur_chunk_size):
                yield (chrom, i, min(chroms_to_lengths[chrom], i+cur_chunk_size))


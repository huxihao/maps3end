"""
Utilities
"""
import os
import string   
from StringIO import StringIO

###############################################################################

def dna_reverse_complement(dna):
    """
    Reverse complement a DNA sequence
    """        
    complements = string.maketrans('acgtACGT', 'tgcaTGCA')
    return dna.translate(complements)[::-1]

###############################################################################

def blocks2bed(block_sizes, block_starts, 
    start=0, strand='+', chrom='chr1', name='b', rgb='0,0,0'):
    """
    Make Bed12 record according to blockSizes, blockStarts
    """
    assert start >= 0, 'start needs to be >=0'
    block_sizes = block_sizes.strip(',')
    try:
        sizes = map(int, block_sizes.split(','))
    except ValueError:
        assert False, "block_sizes should be '100,100'"
    
    block_starts = block_starts.strip(',')
    try:
        starts = map(int, block_starts.split(','))
    except ValueError:
        assert False, "block_starts should be '0,300'"
    assert starts[0] == 0, "block_starts[0] should be 0"
     
    num_blocks = len(sizes)
    assert len(starts) == num_blocks, "block_starts #fields != block_sizes"
    
    for i in range(num_blocks-1):
        assert sizes[i] < starts[i+1], "block %d start in last block" % (i+1)
    
    end = start + starts[-1] + sizes[-1] 
       
    return (chrom, start, end, name, 0, strand, start, start, rgb, 
        num_blocks, block_sizes, block_starts)
    
###############################################################################

def generator_bed_from_blocks(fin, col_bsizes=0, col_bstarts=1, col_start=2, 
    strand='+', chrom='chr1', name_prefix='', rgb='0,0,0'):
    """
    Generator of bed records from file
    Each line with block_sizes, block_starts, start
    """
    i = 0
    for line in fin:
        i += 1
        fields = line.rstrip().split("\t")
        assert len(fields) >= 3, 'error at line %d, not enough fields'
        try:
            start = int(fields[col_start])
        except ValueError:
            assert False, "col_start should be integer"
        
        yield blocks2bed(fields[col_bsizes], fields[col_bstarts], 
            start, strand=strand, chrom=chrom, 
            name="%s%d" % (name_prefix, i), rgb=rgb)
        
###############################################################################

#http://stackoverflow.com/questions/3012421/python-lazy-property-decorator
class lazy_property(object):
    """
    meant to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    """
    def __init__(self, fget):
        self.fget = fget
        self.func_name = fget.__name__

    def __get__(self, obj, cls):
        if obj is None:
            return None
        value = self.fget(obj)
        setattr(obj, self.func_name, value)
        return value

###############################################################################

def bed2bam(f_bed, f_chrominfo, bamprefix, is_bed12=False):
    """Convert bed file to BAM file by BEDTools and SAMtools"""
    f_bam = bamprefix + ".bam"
    bam_tmp = bamprefix + "_tmp.bam"
    bed12 = ''
    if is_bed12:
        bed12 = '-bed12'
    os.system("bedToBam %s -i %s -g %s > %s" % (bed12, f_bed, f_chrominfo, bam_tmp))
    os.system("samtools sort %s %s" % (bam_tmp, bamprefix))
    os.system("samtools index %s" % (f_bam, ))
    os.unlink(bam_tmp)

###############################################################################

def string2file(s, filename):
    """Write given string to file with filename"""
    with open(filename, "w") as oh:
        for line in StringIO(s):
            oh.write(line)
            
            
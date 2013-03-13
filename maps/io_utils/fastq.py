#!/usr/bin/env python
"""
Fastq format
"""

################################################################################

class Fastq(object):
    """
    Fastq record
    """
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
    
    def __str__(self):
        return '\n'.join(['@%s' % str(self.name), self.seq, 
            '+%s' % self.name, self.qual])
    
    def as_simple(self):
        return '\n'.join(['@%s' % str(self.name), self.seq, '+', self.qual])
    
    def get_name(self):
        return self.name
    
    def clean_name(self):
        self.name = self.name.split(" ")[0]
    
    def get_seq(self):
        return self.seq
    
    def get_qual(self):
        return self.qual
    
    def get_ascii(self):
        return map(ord, self.qual)
    
    def get_length(self):
        """Return length"""
        return len(self.seq)

    def get_num_N(self, first=28):
        """Return number of Ns"""
        return self.seq.count('N', 0, first)

    def set_name(self, name):
        """Reset name"""
        self.name = name

    def remove_tail_N(self):
        """Remove N at the end"""
        i = len(self.seq) - 1
        while i >= 0:
            if self.seq[i] == 'N':
                i -= 1
            else:
                break
        self.seq = self.seq[:i+1]
        self.qual = self.qual[:i+1]

################################################################################
    
def reader_fastq(infile):
    """Generator of Fastq object from given file"""
    i = 0
    name = None
    seq = None
    qual = None
    for line in open(infile):
        i += 1
        curr_line = line.strip()
        if i % 4 == 1:
            name = curr_line[1:]
        elif i % 4 == 2:
            seq = curr_line
        elif i % 4 == 0:
            qual = curr_line
            yield Fastq(name, seq, qual)
            
################################################################################
    
def rename_fastq(infile, prefix, outfile):
    """Rename tags with given prefix+id (order in original file)"""
    i = 0
    outhandle = open(outfile, 'w')
    for fq in reader_fastq(infile):
        i += 1
        fq.set_name(prefix+str(i))
        outhandle.write(str(fq)+'\n')
    outhandle.close()
    
################################################################################
    
def filter_fastq(infile, outfile, minlen=6, filterN=False):
    """
    Remove tailing Ns;
    Filter those tags having more than Ns in first given_number nucleotides"""
    outhandle = open(outfile, 'w')
    for fq in reader_fastq(infile):
        keep = True
        
        lenseq = len(fq.get_seq())
        lenqual = len(fq.get_qual())
        if lenseq != lenqual or lenseq < minlen:
            keep = False
        
        if filterN:
            if lenseq <= 10:
                if fq.get_num_N() > 0:
                    keep = False
            elif lenseq <= 20:
                if fq.get_num_N() > 1:
                    keep = False
            else:
                if fq.get_num_N(28) > 2:
                    keep = False
             
        if keep:
            outhandle.write(str(fq)+'\n')
    outhandle.close()
    
################################################################################


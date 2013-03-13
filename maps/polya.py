#!/usr/bin/env python
"""
PolyA Signal and A-stretch module
"""

import re
import os
from maps.utils import dna_reverse_complement
from maps.io_utils.bed import BedReader, Bed
from distutils.version import LooseVersion
import HTSeq
assert LooseVersion(HTSeq.__version__) >= LooseVersion('0.5.3p9') 
from maps.utils import lazy_property

############################################################################### 

class ElemScan(object):
    """
    Sequence elements scanner
    """
    def __init__(self, genomedir, fival=None):
        self.genomedir = genomedir
        self.fival = fival
    
    @lazy_property
    def gas(self): 
        """GenomicArray: intervals marked"""
        if self.fival is None:
            return None
        
        gas = HTSeq.GenomicArray([], stranded=True, typecode='b')
        for b in BedReader(open(self.fival)):
            if not isinstance(b, Bed):
                continue
            
            chrom = b.get_chrom()
            if chrom not in gas.chrom_vectors:
                gas.add_chrom(chrom)
            
            iv = HTSeq.GenomicInterval(chrom, b.get_start(), b.get_end(), 
                                       b.get_strand())
            gas[iv] = True
        
        return gas    
    
    def iterfasta(self):
        """Iterate Fasta files by chromosome"""
        for fname in os.listdir(self.genomedir):
            if fname.endswith('.fa') and fname.startswith('chr'):
                chrom = fname[:-3]
                yield chrom, os.path.join(self.genomedir, fname)
    
    def scan(self):
        """Scan for class-specific elements"""
        for chrom, f_fasta in self.iterfasta():
            fin = open(f_fasta)
            fin.next() # skip first line >name
            
            for (s, e, strand) in self.scan_elem(fin):
                # check whether covered by one interval if given
                if self.gas is None:
                    yield chrom, s, e, strand
                else:
                    if chrom not in self.gas.chrom_vectors:
                        continue
                    iv = HTSeq.GenomicInterval(chrom, s, e, strand)
                    if all([step[1] for step in self.gas[iv].steps()]):
                        yield chrom, s, e, strand
                    
            fin.close()
                
    def scan_elem(self, fin):
        """Scan the first Fasta record from given file handle"""
        start = 0
        curr_line = ''
        curr_end = 0
        for line in fin:
            if line.startswith(">"): 
                continue

            curr_line = curr_line[curr_end:] + line.rstrip().upper()
            for strand in ["+", "-"]:
                for s, e in self.itermatch(curr_line, strand):
                    yield (start+s, start+e, strand)    
            
            curr_end = max(0, len(curr_line) - self.ELEM_LENGTH + 1)
            start += curr_end
        
        curr_line = curr_line[curr_end:]
        for strand in ["+", "-"]:
            for s, e in self.itermatch(curr_line, strand):
                yield (start+s, start+e, strand) 

    def itermatch(self, dna, strand):
        """Abstract method to iterate match"""
        raise Exception("Abstract method, not to be called")
        
############################################################################### 

class PasScan(ElemScan):
    """Scan for PAS"""
    PAS_MOTIFS = (
        'AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'CATAAA',
        'GATAAA', 'AATATA', 'AATACA', 'AATAGA', 'AATGAA', 
        'ACTAAA', 'AACAAA','TTTAAA')
    PAS_FWD = re.compile("|".join(PAS_MOTIFS))
    PAS_REV = re.compile("|".join([dna_reverse_complement(m) for m in PAS_MOTIFS]))
    ELEM_LENGTH = max([len(m) for m in PAS_MOTIFS])

    def itermatch(self, dna, strand):
        """Generate matched (s, e) index on given sequence and strand"""
        sidx = None
        eidx = None
        if strand == '+':
            pat = PasScan.PAS_FWD
        else:
            pat = PasScan.PAS_REV
        
        for m in pat.finditer(dna):
            if sidx is None:
                sidx = m.start()
                eidx = m.end()
            elif eidx >= m.start():
                eidx = m.end()
            else:
                yield (sidx, eidx)
                sidx = m.start()
                eidx = m.end()
        
        if sidx:
            yield (sidx, eidx)
            
############################################################################### 
        
class AstretchScan(ElemScan):
    """Scan for A-stretches"""
    ASTRETCH_WINDOW = 10
    ELEM_LENGTH = ASTRETCH_WINDOW
    
    def itermatch(self, dna, strand):
        """Generate matched (s, e) index on given sequence and strand"""
        sidx = None
        eidx = None
        for j in xrange(len(dna) - AstretchScan.ASTRETCH_WINDOW):
            win = dna[j:(j+AstretchScan.ASTRETCH_WINDOW)]
            pos = None
            if strand == '+':
                if "AAAAAAAA" in win or win.count('A') >= 9:
                    pos = j 
            else:
                if "TTTTTTTT" in win or win.count('T') >= 9:
                    pos = j
            
            if pos is None:
                continue

            if sidx is None:
                sidx = pos
                eidx = pos + AstretchScan.ASTRETCH_WINDOW
            elif eidx >= pos:
                eidx = pos + AstretchScan.ASTRETCH_WINDOW
            else:
                yield (sidx, eidx)
                sidx = pos
                eidx = None
        
        if eidx is not None:
            yield (sidx, eidx)

###############################################################################
    
#!/usr/bin/env python
""" 
Module for decoding multiplex barcodes
"""

import sys
from maps.io_utils.fastq import reader_fastq

################################################################################

def read_barcode(file_barcode):
    """
    Return barcode:name dict
    """
    barcode2name = {}
    for line in open(file_barcode):
        fields = line.rstrip().split('\t')
        barcode2name[fields[0]] = fields[1]
    return barcode2name

################################################################################

def barcode_distance(str_a, str_b):
    """
    Compute distance between two barcodes
    """
    assert len(str_a) == len(str_b)
    return sum([1 for i in xrange(len(str_a)) if str_a[i] != str_b[i]])

################################################################################
    
def decode_fastq(infile, barcode2name, mismatch=0, startpos=37, outprefix=None):
    """
    Decode barcodes in given fastq file
    """
    if not outprefix:
        outprefix = infile
    barcode2outhandle = {
        'failed':open(outprefix+'.failed.m'+str(mismatch), 'w'),
        }
    for bc in barcode2name:
        outfile = outprefix+'.'+barcode2name[bc]+'.m'+str(mismatch)
        barcode2outhandle[bc] = open(outfile, 'w')

    for record in reader_fastq(infile):
        try:
            currseq = record[1][startpos-1:startpos+len(bc)-1]
        except:
            continue
        
        barcode2dist = [(bc, barcode_distance(bc, currseq)) 
            for bc in barcode2name]
        
        #print >> sys.stderr, record[0], record[1], barcode2dist
        min_dist = min([barcode2dist[i][1] for i in xrange(len(barcode2dist))])
        bc_mindist = [barcode2dist[i][0] for i in xrange(len(barcode2dist)) 
            if barcode2dist[i][1] == min_dist]
        
        if min_dist > mismatch or len(bc_mindist) > 1:
            barcode2outhandle['failed'].write('\n'.join(record)+'\n')
        else:
            barcode2outhandle[bc_mindist[0]].write('\n'.join(record)+'\n')

################################################################################

class Decoder(object):
    """
    Decoding sequencing data given barcode file
    Use DecorderExact if no mismatch allowed
    """
    def __init__(self, infile, barcode2name, 
        mismatch=1, startpos=37, outprefix=None):
        """
        barcode2name: a dict for barcode:samplename
        mismatch: maximum allowed mismatch
        startpos: 1-based start position
        """
        self.infile = infile
        self.barcode2name = barcode2name
        self.mismatch = mismatch
        self.startpos = startpos - 1
        assert self.startpos >= 0, "given start position should be 1-based"
        
        self.outprefix = outprefix
        if not outprefix:
            self.outprefix = infile
            
        self.lenbc = self.check_barcode_len()
        assert not self.lenbc is None, "No barcode"
            
        self.barcode2outhandle = {
            'failed':open(self.outprefix+'.failed.m'+str(mismatch), 'w'),
        }
        
        for bc in self.barcode2name:
            outfile = self.outprefix+'.'+self.barcode2name[bc]+'.m'+str(self.mismatch)
            self.barcode2outhandle[bc] = open(outfile, 'w')
        #capture IOError if too many out-handles at the same time (>=1021) 
        
    def check_barcode_len(self):
        """
        Check whether all barcodes have same length and return it
        Return None if no barcode
        """
        lens = [len(k) for k in self.barcode2name]
        if len(lens) == 0:
            return None
        bclen = lens[0]
        for bl in lens[1:]:
            assert bl == bclen
        return bclen

    def decode(self):
        """decode and write results"""
        for record in reader_fastq(self.infile):
            try:
                currseq = record[1][self.startpos:(self.startpos+self.lenbc)]
                if len(currseq) < self.lenbc:
                    sys.stderr.write("%s not have enough length" % currseq)
            except:
                continue
            
            bc = self.decode_one(currseq)
            self.barcode2outhandle[bc].write('\n'.join(record)+'\n')
        
        for oh in self.barcode2outhandle.values():
            oh.close()
    
    def decode_one(self, seq):
        """
        Return the barcode to the given sequence, 'failed' if no one found
        """
        bc_decode = "failed"
        barcode2dist = []
        for bc in self.barcode2name:
            dist = 0
            for i in xrange(self.lenbc):
                if bc[i] != seq[i]:
                    dist += 1
                if dist > self.mismatch:
                    break
            if dist > self.mismatch:
                continue
            barcode2dist.append((bc, dist))
        
        if len(barcode2dist) == 0:
            return bc_decode
        
        min_dist = min([barcode2dist[i][1] for i in xrange(len(barcode2dist))])
        bc_mindist = [barcode2dist[i][0] for i in xrange(len(barcode2dist)) 
            if barcode2dist[i][1] == min_dist]
        
        if len(bc_mindist) == 1:
            bc_decode = bc_mindist[0]
        
        return bc_decode

################################################################################

class DecoderExact(Decoder):
    """
    Decoding without mismatch
    """ 
    def decode_one(self, seq):
        """
        Return the barcode to the given sequence, 'failed' if no one found
        """
        if seq in self.barcode2name:
            return seq
        return "failed"
    
################################################################################

def create_decoder(infile, barcode2name, mismatch=0, startpos=37, outprefix=None):
    """Factor of decoder"""
    if mismatch == 0:
        dclass = DecoderExact
    else:
        dclass = Decoder
        
    return dclass(infile, barcode2name, mismatch=mismatch, startpos=startpos, 
                  outprefix=outprefix)

################################################################################

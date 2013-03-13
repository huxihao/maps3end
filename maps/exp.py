#!/usr/bin/env python
"""
Expression calculation for 3'-end RNA-seq 
"""
import os
import sys
from distutils.version import LooseVersion
import pysam
import HTSeq
assert LooseVersion(HTSeq.__version__) >= LooseVersion('0.5.3p9') 
from maps.utils import lazy_property

################################################################################

class Counter(object):
    """Counter unit:"""
    # see doc/design
    def __init__(self, name):
        self._name = name
        self._unique = 0
        self._shared = 0.0 
        self._id_features = set()

    @property
    def name(self):
        return self._name
        
    @property
    def unique(self):
        return self._unique
    
    @unique.setter
    def unique(self, value):
        self._unique = value
    
    @property
    def shared(self):
        return self._shared
    
    @shared.setter
    def shared(self, value):
        self._shared = value
    
    @property
    def num(self):
        return self._unique + self._shared
        
    def __str__(self):
        return "%s\t%f" % (self.name, self.num)
    
    def add_idfeature(self, id_feature):
        self._id_features.add(id_feature)
    
    def get_idfeatures(self):
        return self._id_features
   
    def rpkm(self, length, numread):
        return 10.0**9 * self.num / (length * numread)
    
################################################################################

class DGE(object):
    """DGE"""
    def __init__(self, tag_file, gtf_file, feature_type, id_feature, id_count):
        assert os.path.exists(tag_file), '%s not exists' % tag_file 
        assert os.path.exists(gtf_file), '%s not exists' % gtf_file
        self.tag_file = tag_file
        self.gtf_file = gtf_file
        self.feature_type = feature_type
        self.id_feature = id_feature
        self.id_count = id_count
        self.numread = 0
        self.nshared = 0
        self.nempty = 0
        self.nunique = 0
    
    def __str__(self):
        return "\n".join(["#%s" % self.tag_file, 
                          "#%s" % self.gtf_file,
                          "#read: %d" % self.numread, 
                          "#unique: %d" % self.nunique,
                          "#shared: %d" % self.nshared, 
                          "#empty: %d" % self.nempty])
        
    def count(self, outhandle=sys.stdout):
        """count and write output"""
        sys.stderr.write("counting unique ...\n")
        self.count_unique()
        
        sys.stderr.write("counting shared ...\n")
        self.count_shared()
        
        sys.stderr.write("writing result ...\n")
        outhandle.write("%s\n" % str(self))
        for cid in sorted(self.cid2counter):
            fields = []
            for v in self.output(cid):
                if isinstance(v, float):
                    v = "%.4f" % v
                fields.append(str(v))
                
            outhandle.write("\t".join(map(str, fields))+"\n")
    
    def output(self, cid):
        """output fields for one Counter object"""
        c = self.cid2counter[cid]
        clen = self.cid2length[cid]
        rpkm = c.rpkm(clen, self.numread)
        return [c.name, c.unique, c.shared, c.num, self.numread, clen, rpkm]
    
    @lazy_property
    def gas(self): 
        """GenomicArrayOfSets: intervals associated id_feature"""
        gas = HTSeq.GenomicArrayOfSets([], stranded=True)
        for f in HTSeq.GFF_Reader(self.gtf_file):
            if f.iv.chrom not in gas.chrom_vectors:
                gas.add_chrom(f.iv.chrom)
            if f.type == self.feature_type:
                feature_id = f.attr[self.id_feature]
                gas[f.iv] += feature_id
        
        return gas

    @lazy_property
    def fid2cid(self):
        """dict of feature_id to count_id"""
        fid2cid = {}
        for f in HTSeq.GFF_Reader(self.gtf_file):
            if f.type == self.feature_type:
                feature_id = f.attr[self.id_feature]
                count_id = f.attr[self.id_count]
                if feature_id not in fid2cid:
                    fid2cid[feature_id] = count_id
        
        return fid2cid
    
    @lazy_property
    def cid2counter(self):
        """dict of count_id to Counter object"""
        cid2counter = {}
        for fid, cid in self.fid2cid.iteritems():
            if cid not in cid2counter:
                cid2counter[cid] = Counter(cid)
            cid2counter[cid].add_idfeature(fid)
                    
        return cid2counter

    @lazy_property
    def cid2length(self):
        """dict of count_id to length of merged interval"""
        cid2length = {}
        for iv, fs in self.gas.steps():
            cids = set([self.fid2cid[f] for f in fs])
            for cid in cids:
                if cid not in cid2length:
                    cid2length[cid] = 0
                cid2length[cid] += iv.length
        
        return cid2length
     
    def iter_reads(self):
        """Generator of reads as tuple of GenomicInterval object"""
        bam = pysam.Samfile(self.tag_file, "rb")
        for hit in bam.fetch():
            if hit.is_unmapped:
                continue
            
            chrom = bam.references[hit.rname]
            strand = '+'
            if hit.is_reverse:
                strand = '-'

            iv_list = []
            start = hit.positions[0]
            for ctype, clen in hit.cigar: # support only match or skipped region
                if ctype == 0:
                    iv_list.append(HTSeq.GenomicInterval(chrom, start, start+clen, strand))
                elif ctype == 3:
                    pass
                else:
                    sys.exit("unknown CIGAR type: %d, now only support[0|3]" % ctype)
                start += clen
                                
            yield tuple(iv_list)
    
    def count_unique(self):
        """Count reads that is uniquely assignable"""
        for read in self.iter_reads():
            self.numread += 1
            fs = self.read2features(read)
            if len(fs) > 1: 
                self.nshared += 1
            elif len(fs) == 0:
                self.nempty += 1
            else: 
                self.nunique += 1
                cid = self.fid2cid[fs[0]]
                cnter = self.cid2counter[cid]
                cnter.unique += 1
    
    def count_shared(self):
        """Assign ambiguous/shared reads to counter unit"""
        for read in self.iter_reads():
            fs = self.read2features(read)
            if len(fs) > 1:
                cnters = [self.cid2counter[self.fid2cid[f]] for f in fs]
                uniqcnts = [ct.unique + 1 for ct in cnters] # with prior 1
                uniqsum = float(sum(uniqcnts))
                for i in range(len(cnters)):
                    cnters[i].shared += uniqcnts[i] / uniqsum
    
    def read2features(self, iv_tuple):
        """
        Find features assignable to given read
        iv_tuple: representation of mapped read as tuple of GenomicInterval
        intersection-strict mode
        """
        fs = None
        for iv in iv_tuple:
            if iv.chrom not in self.gas.chrom_vectors:
                continue
            for step in self.gas[iv].steps():
                fs2 = step[1]
                if fs is None:
                    fs = fs2.copy()
                else:
                    fs = fs.intersection(fs2)
        
        if fs is None: 
            return []
        else:
            return list(fs)
        
################################################################################
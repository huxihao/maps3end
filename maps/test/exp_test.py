#!/usr/bin/env python

import unittest
import HTSeq
from maps.exp import Counter, DGE

###############################################################################

class Counter_Test(unittest.TestCase):
    def setUp(self):
        self.c = Counter('g1')
    
    def test_name(self):
        assert self.c.name == 'g1'
        
    def test_unique(self):
        assert self.c.unique == 0
        self.c.unique += 2
        assert self.c.unique == 2
    
    def test_shared(self):
        assert self.c.shared == 0
        self.c.shared += 0.5
        assert self.c.shared == 0.5
    
    def test_num(self):
        self.c.unique += 2
        self.c.shared += 0.5
        assert self.c.num == 2.5
        
    def test_rpkm(self):
        assert self.c.rpkm(1000, 1000000) == 0
        self.c.unique += 2
        self.c.shared += 1
        assert self.c.rpkm(1000, 1000000) == 3
        
###############################################################################    
    
class DGE_Test(unittest.TestCase):
    def setUp(self):
        self.dg = DGE("data/read.bam", "data/gene.gtf", feature_type="exon", 
                       id_feature="gene_id", id_count="gene_id")
        self.dc = DGE("data/read.bam", "data/gene.gtf", feature_type="exon", 
                       id_feature="gene_id", id_count="cluster_id")
        
    def test_gas(self):
        assert self.dg.gas.stranded 
    
    def test_fid2cid(self):
        assert self.dg.fid2cid['g2'] == 'g2'
        assert self.dc.fid2cid['g2'] == '1'
        
    def test_cid2counter(self):
        assert 'g2' in self.dg.cid2counter
        assert '1' in self.dc.cid2counter
        
    def test_cid2length(self):
        assert self.dg.cid2length['g2'] == 300
        assert self.dc.cid2length['1'] == 400
        
    def test_iter_reads(self):
        assert len([r for r in self.dg.iter_reads()]) == 16
        
    def test_count_unique(self):
        self.dg.count_unique()
        assert self.dg.nunique == 5
        assert self.dg.nshared == 1
        
    def test_count_shared(self):    
        self.dg.count_unique()
        self.dg.count_shared()
        assert self.dg.cid2counter['g2'].shared == 0.5
        self.dc.count_unique()
        self.dc.count_shared()
        assert self.dc.cid2counter['1'].shared == 1 # 
        
    def test_read2features(self):
        iv = HTSeq.GenomicInterval('chr1', 1600, 1620, '+')
        assert sorted(self.dg.read2features((iv,))) == ['g1', 'g2']
        assert sorted(self.dc.read2features((iv,))) == ['g1', 'g2']
    
###############################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
    

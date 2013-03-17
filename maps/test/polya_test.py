#!/usr/bin/env python

import unittest
from maps.polya import ElemScan, PasScan, AstretchScan

###############################################################################

class ElemScan_Test(unittest.TestCase):
    def setUp(self):
        self.esc = ElemScan('data/hgex')
    
    def test_iterfasta(self):
        assert [fname for fname in self.esc.iterfasta()] == [
            ('chr1', 'data/hgex/chr1.fa'), ('chr2', 'data/hgex/chr2.fa')]
        

class PasScan_Test(unittest.TestCase):
    def setUp(self):
        self.esc = PasScan('data/hgex')
    
    def test_itermatch(self):
        assert [m for m in self.esc.itermatch('CCCAAA', '+')] == [] 
        assert [m for m in self.esc.itermatch('AATAAA', '+')] == [(0, 6)] 
        assert [m for m in self.esc.itermatch('ATTAAA', '+')] == [(0, 6)]
        assert [m for m in self.esc.itermatch('AATAAAAGTAAA', '+')] == [(0, 12)] 
        assert [m for m in self.esc.itermatch('AATAAACCAGTAAA', '+')] == [(0, 6), (8, 14)] 
        assert [m for m in self.esc.itermatch('AATAAATAAA', '+')] == [(0, 6)]
        
        
class AstretchScan_Test(unittest.TestCase):
    def setUp(self):
        self.esc = AstretchScan('data/hgex')
    
    def test_itermatch(self):
        assert [m for m in self.esc.itermatch('AAAAAAAA', '+')] == [] #>=10
        assert [m for m in self.esc.itermatch('AAAAAAAAAA', '+')] == [(0, 10)]
        assert [m for m in self.esc.itermatch('AAAACAAAAA', '+')] == [(0, 10)]
        assert [m for m in self.esc.itermatch('AAAACACAAA', '+')] == []
        assert [m for m in self.esc.itermatch('CCAAAAAAAA', '+')] == [(2, 10)]
        assert [m for m in self.esc.itermatch('AAAAAAAACC', '+')] == [(0, 8)]
        assert [m for m in self.esc.itermatch('CCCCCAAAAAAAAAACCCCC', '+')] == [(5, 15)]
        assert [m for m in self.esc.itermatch('AAAAAAAACCAAAAAAAACC', '+')] == [(0, 18)]
        assert [m for m in self.esc.itermatch('AAAAAAAACCCAAAAAAAACC', '+')] == [(0, 8), (11, 19)]
      
        assert [m for m in self.esc.itermatch('TTTTTTTT', '-')] == [] #>=10
        assert [m for m in self.esc.itermatch('TTTTTTTTTT', '-')] == [(0, 10)]
        assert [m for m in self.esc.itermatch('TTTTCTTTTT', '-')] == [(0, 10)]
        assert [m for m in self.esc.itermatch('TTTTCTCTTT', '-')] == []
        assert [m for m in self.esc.itermatch('CCTTTTTTTT', '-')] == [(2, 10)]
        assert [m for m in self.esc.itermatch('TTTTTTTTCC', '-')] == [(0, 8)]
        assert [m for m in self.esc.itermatch('CCCCCTTTTTTTTTTCCCCC', '-')] == [(5, 15)]
        assert [m for m in self.esc.itermatch('TTTTTTTTCCTTTTTTTTCC', '-')] == [(0, 18)]
        assert [m for m in self.esc.itermatch('TTTTTTTTCCCTTTTTTTTCC', '-')] == [(0, 8), (11, 19)]
              
###############################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
    

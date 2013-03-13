#!/usr/bin/env python
"""
Tests of module utils
"""
import unittest
from maps.utils import *
from StringIO import StringIO

############################################################################### 

class Testblocks2bed(unittest.TestCase):
    """Test blocks2bed function"""
    def test_start(self):
        try:
            blocks2bed('10,10', '0,50', start=-100)
            self.fail("start should >=0")
        except AssertionError:
            pass

    def test_block1_start(self):
        try:
            blocks2bed('10,10', '10,50')
            self.fail("block1 start should =0")
        except AssertionError:
            pass

    def test_size_start_length(self):
        try:
            blocks2bed('10,10,10', '0,50')
            self.fail("block sizes and starts not same length")
        except AssertionError:
            pass
        
    def test_size_start_conflict(self):
        try:
            blocks2bed('50,10', '0,30', start=0)
            self.fail("block2 start in block1")
        except AssertionError:
            pass
    
    def test_blocks2bed(self):
        assert blocks2bed('10,10,', '0,50,', start=10) == (
            'chr1', 10, 70, 'b', 0, '+', 10, 10, '0,0,0', 2, '10,10', '0,50')

        assert blocks2bed('10,10', '0,50', start=10, strand='-') == (
            'chr1', 10, 70, 'b', 0, '-', 10, 10, '0,0,0', 2, '10,10', '0,50')

############################################################################### 

blocks = """\
10,10,\t0,50\t10
20,40,\t0,50\t100
"""

class TestGeneratorBed(unittest.TestCase):
    """Test generator_bed_from_blocks function"""
    def setUp(self):
        self.fin = StringIO(blocks)
        
    def test_generator(self):
        beds = [bed for bed in generator_bed_from_blocks(self.fin)]
        assert beds == [('chr1', 10, 70, '1', 0, '+', 10, 10, '0,0,0', 2, '10,10', '0,50'), 
            ('chr1', 100, 190, '2', 0, '+', 100, 100, '0,0,0', 2, '20,40', '0,50')]
    
    def test_names(self):
        names = [bed[3] for bed in generator_bed_from_blocks(self.fin, name_prefix='g')]
        assert names == ['g1', 'g2']

    def test_strands(self):
        strands = [bed[5] for bed in generator_bed_from_blocks(self.fin, strand='-')]
        assert strands == ['-', '-']        
        
############################################################################### 

if __name__ == '__main__':        
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
    
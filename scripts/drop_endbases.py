#!/usr/bin/env python
"""
Drop end bases of given length for fastq format
"""

import os
import sys

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: program infile cut_len[8]'
        sys.exit()
          
    infile = sys.argv[1]
    assert os.path.exists(infile)
    cut_len = 8
    if len(sys.argv) >= 3:
        cut_len = int(sys.argv[2])
    
    i = 0
    for line in open(infile):
        i += 1
        line = line.rstrip()
        if i % 2 == 0:
            line = line[:-cut_len]
        print line
        

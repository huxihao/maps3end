#!/usr/bin/env python
"""Rename bed file, given each record a unique name"""

import os
import sys
from maps.io_utils.bed import BedReader, Bed

################################################################################

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage: program infile prefix outfile'
        sys.exit()
        
    infile = sys.argv[1]
    assert os.path.exists(infile)
    prefix = sys.argv[2]
    outfile = sys.argv[3]
    
    outhandle = open(outfile, 'w')
    i = 0
    for b in BedReader(open(infile)):
        if not isinstance(b, Bed):
            continue
        i += 1
        
        b.fields[3] = '%s%d' % (prefix, i)
        outhandle.write("\t".join(map(str, b.fields))+"\n")
    outhandle.close()
    
    

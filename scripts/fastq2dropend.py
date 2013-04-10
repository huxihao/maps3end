#!/usr/bin/env python
"""
Drop end bases of given length for Fastq format
"""

import os
import sys
import optparse

################################################################################

def process_command_line(argv):
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    usage = "usage: %prog [options] infile.fastq"
    parser = optparse.OptionParser(usage,
        formatter=optparse.TitledHelpFormatter(width=78),
        add_help_option=None)

    # define options here:
    parser.add_option("-o", "--outfile", dest="outfile",
        help="outfile name", metavar="FILE")

    parser.add_option("-d", "--drop_len", dest="drop_len", default=0,
        help='number of nt to drop(default 0)', type="int")
    
    parser.add_option("-e", "--end", dest="end", default=3,
        help='which end to drop: 5 or 3 (default 3)', type="int")    

    parser.add_option('-h', '--help', action='help',
        help='Show this help message and exit.')

    parser.add_option("-q", "--quiet",
        action="store_false", dest="verbose", default=True,
        help="don't print status messages to stdout")

    (options, args) = parser.parse_args(argv)
    if len(args) < 2:
        parser.error(usage)
    if options.end not in (5, 3):
        parser.error('--end only accepts 5 (5prime end) or 3(3prime end)')
        
    return options, args

################################################################################

if __name__ == '__main__':
    options, args = process_command_line(sys.argv)
    infile = args[1]
    assert os.path.exists(infile)
    cut_len = options.drop_len
    
    outhandle = sys.stdout
    if options.outfile:
        outhandle = open(options.outfile, 'w')
    
    num_short = 0
    i = 0
    for line in open(infile):
        i += 1
        line = line.rstrip()
        
        if i % 2 == 0:
            if len(line) < cut_len:
                num_short += 1
                
            if options.end == 3:
                line = line[:-cut_len]
            elif options.end == 5:
                line = line[cut_len:]
                
        outhandle.write(line+"\n")   
        
    sys.stderr.write("#shorter than cut_len(%d): %d\n" % (cut_len, num_short))

    if options.outfile:
        outhandle.close()
           

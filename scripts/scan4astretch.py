#!/usr/bin/env python
"""
Find PolyA stretches [in given intervals] of genome
The genome_dir contains files of sequences by chromosome in Fasta format 
Example: unzipped *.fa files as in http://hgdownload.cse.ucsc.edu/goldenPath/hg18/chromosomes/ 
"""
import os
import sys
import optparse
from maps.polya import AstretchScan

################################################################################

def process_command_line(argv):
    usage = "%s\nUsage: %s [options] genome_dir" % (__doc__, argv[0])
    parser = optparse.OptionParser(usage,
        formatter=optparse.TitledHelpFormatter(width=78),
        add_help_option=None)
         
    parser.add_option("-o", "--outfile", dest="outfile",
        help="output file name", metavar="FILE")

    parser.add_option("--ival", dest="f_ival", 
        help="find only in the intervals[default None]", metavar="FILE")
    
    parser.add_option('-h', '--help', action='help',
        help='Show this help message and exit.')

    (options, args) = parser.parse_args(argv)
    if len(args) < 2:
        parser.error("genome_dir is required")
        
    return options, args

################################################################################

if __name__ == '__main__':
    options, args = process_command_line(sys.argv)    
    genome_dir = args[1]
    assert os.path.exists(genome_dir)
    
    outhandle = sys.stdout
    if options.outfile:
        outhandle = open(options.outfile, "w")
    
    scanner = AstretchScan(genome_dir, options.f_ival)
    i = 0
    for chrom, start, end, strand in scanner.scan():
        i += 1
        outhandle.write("\t".join(map(str, [
            chrom, start, end, i, 0, strand]))+"\n")
    
    if options.outfile:
        outhandle.close()
    
    

#!/usr/bin/env python
"""
Compute expression in the defined landing intervals associated with counting units
"""

import sys
import optparse
from maps.exp import DGE

################################################################################

def process_command_line(argv):
    if argv is None:
        argv = sys.argv[1:]
           
    usage = "%s\nusage: prog [options] f_read f_gtf" % __doc__
    parser = optparse.OptionParser(usage, 
        formatter=optparse.TitledHelpFormatter(width=178),
        add_help_option=True)
         
    parser.add_option("-t", "--type", type="string", dest="featuretype",
        default = "exon", help = "feature type (3rd column in GTF file)[exon]")
         
    parser.add_option("-u", "--unit", type="string", dest="unit",
        default = "transcript_id", help = "GTF attribute as counting unit[transcript_id]")

    parser.add_option("-o", "--outfile", type="string", dest="outfile",
        help = "out file name")

    (options, args) = parser.parse_args()
   
    if len(args) != 2:
        parser.error('No required parameters')
          
    return options, args

################################################################################

def main():
    options, args = process_command_line(None)
    f_read = args[0]
    f_gtf = args[1]
    
    dge = DGE(f_read, f_gtf, feature_type=options.featuretype, 
        id_feature = "gene_id", id_count=options.unit)

    outhandle = sys.stdout
    if options.outfile:
        outhandle = open(options.outfile, "w")
        
    dge.count(outhandle)
    
    if options.outfile:
        outhandle.close()
        
        
if __name__ == "__main__":
    main()


#!/usr/bin/env python
"""
Decode barcodes in one lane
"""

import os
import sys
import optparse
from maps.barcode import read_barcode, create_decoder

################################################################################

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    'argv' is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    usage = "usage: %s [options] infile file_barcode" % sys.argv[0]
    parser = optparse.OptionParser(usage,
        formatter=optparse.TitledHelpFormatter(width=78),
        add_help_option=None)

    # define options here:
    parser.add_option("-s", "--startpos", dest="startpos", default=37,
        help="1-based start position of barcode[37 default]", type="int")
    
    parser.add_option("-m", "--mismatch", dest="mismatch", default=0,
        help="allowed mismatch[0 default]", type="int")

    parser.add_option("--outprefix", dest="outprefix", 
        help="out prefix", type="str")

    parser.add_option('-h', '--help', action='help',
        help='Show this help message and exit.')

    parser.add_option("-q", "--quiet",
        action="store_false", dest="verbose", default=True,
        help="don't print status messages to stdout")

    (options, args) = parser.parse_args(argv)
    if argv:
        if len(sys.argv) < 3:
            parser.error(usage)

    if not options.outprefix:
        parser.error("--outprefix required")

    # further process settings & args if necessary
    return options, args

################################################################################

if __name__ == '__main__':
    options, args = process_command_line(sys.argv)

    infile = args[1]
    assert os.path.exists(infile)
    file_barcode = args[2]
    assert os.path.exists(file_barcode)
    allowed_mismatch = options.mismatch
    
    barcode2name = read_barcode(file_barcode)
    dc = create_decoder(infile, barcode2name, 
        mismatch=allowed_mismatch, startpos=options.startpos,
        outprefix=options.outprefix)
    
    dc.decode()
    
################################################################################

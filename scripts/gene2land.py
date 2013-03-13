#!/usr/bin/env python
"""Convert gene/clusters to 3'-end exons
Input genes should have unique id, and are clustered by clusterGenes(kent), and input clusters should also have unique id.
3'-ends in a cluster are grouped into unique 3'-ends, which are assigned unique key_id.
Output is in GTF format, feature:exon, attributes: gene_ids, transcript_id set as key_id, cluster_id.
"""

import os
import sys
import optparse
from maps.land import LandMaker

################################################################################

def process_command_line(argv):
    if argv is None:
        argv = sys.argv[1:]
        
    usage = "%s\nusage: prog [options] f_cluster f_gene" % __doc__
    parser = optparse.OptionParser(usage,
        formatter=optparse.TitledHelpFormatter(width=178),
        add_help_option=True)

    parser.add_option('-o', '--outfile', dest='outfile',  
        help='outfile name', type='str')

    parser.add_option('-s', '--source', dest='source', default='refGene', 
        help='gene source[refGene default]', type='str')
    
    parser.add_option('-L', '--lmrna', dest='lmrna', default=300, 
        help='length at mature RNA level to 3-end[300 default]', type='int')    
    (options, args) = parser.parse_args(argv)

    if len(args) < 2:
        parser.error('No required parameters')
    
    return options, args

################################################################################
    
if __name__ == '__main__':
    options, args = process_command_line(None)
    
    f_cluster = args[0]
    f_gene = args[1]
    source = options.source
    assert os.path.exists(f_cluster)
    assert os.path.exists(f_gene)
    outhandle = sys.stdout
    if options.outfile:
        outhandle = open(options.outfile, "w")

    lmk = LandMaker(f_cluster, f_gene, options.lmrna)
    for (chrom, s, e, strand, kid, gids, cid) in lmk.emit():
        outhandle.write("\t".join(map(str, [
            chrom, source, "exon", s + 1, e, 0, strand, '.', "; ".join([
                'gene_id "%s"' % "|".join(gids), 
                'transcript_id "%s"' % kid,
                'cluster_id "%s"' % cid,
                ])
            ]))+"\n")
            
    if options.outfile:
        outhandle.close()
    

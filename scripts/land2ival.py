#!/usr/bin/env python
"""
Merge landing intervals associated with transcript_id/cluster_id
add gene names if a refFlat.bed is given
"""

import os
import sys
import optparse
import HTSeq

################################################################################

def process_command_line(argv):
    if argv is None:
        argv = sys.argv[1:]
        
    usage = "%s\nusage: prog [options] f_gtf refGene.bed" % __doc__
    parser = optparse.OptionParser(usage,
        formatter=optparse.TitledHelpFormatter(width=178),
        add_help_option=True)

    parser.add_option('-o', '--outfile', dest='outfile',  
        help='outfile name', type='str')

    parser.add_option("-t", "--type", type="string", dest="featuretype",
        default = "exon", help = "feature type (3rd column in GTF file)[exon]")
    
    parser.add_option("-u", "--unit", type="string", dest="unit",
        default = "transcript_id", help = "GTF attribute as counting unit[transcript_id]")

    parser.add_option('--refFlat', dest='refFlat',  
        help='refFlat in Bed format', type='str')
    
    (options, args) = parser.parse_args(argv)

    if len(args) < 2:
        parser.error('No required parameters')
    
    return options, args

################################################################################

def gene2name(f_bed, f_flat):
    """
    Match genes by chrom,start,end,strand,cdsStart,cdsEnd,exonStarts,exonEnds
    return dict of f_bed:name to f_flat:name
    require both files are in BED12 format
    """
    key2name = {}
    keycols = [0, 1, 2, 5, 6, 7, 10, 11]
    for line in open(f_flat):
        fields = line.rstrip().split("\t")
        key = tuple([fields[i] for i in range(len(fields)) if i in keycols])
        key2name[key] = fields[3]
    
    gid2name = {}
    for  line in open(f_bed):
        fields = line.rstrip().split("\t")
        key = tuple([fields[i] for i in range(len(fields)) if i in keycols])
        if key in key2name:
            gid2name[fields[3]] = key2name[key]
        else:
            gid2name[fields[3]] = ''
    
    return gid2name    

################################################################################
        
if __name__ == '__main__':
    options, args = process_command_line(None)
    
    f_gtf = args[0]
    f_bed = args[1]
    assert os.path.exists(f_gtf)
    assert os.path.exists(f_bed)
    outhandle = sys.stdout
    if options.outfile:
        outhandle = open(options.outfile, "w")

    gid2name = None
    if options.refFlat:
        assert os.path.exists(options.refFlat)
        gid2name = gene2name(f_bed, options.refFlat)
    
    uid2ivals = {}
    uid2names = {}
    for f in HTSeq.GFF_Reader(f_gtf):
        if f.type == options.featuretype:
            gids = f.attr['gene_id'].split("|")
            uid = f.attr[options.unit]
            if uid not in uid2ivals:
                uid2ivals[uid] = []
            uid2ivals[uid].append(f.iv)
            
            if uid not in uid2names:
                uid2names[uid] = set()
            for gid in gids:
                name = gid2name[gid]
                if name:
                    uid2names[uid].add(name)
    
    for uid in sorted(uid2names):
        ivs = uid2ivals[uid]
        chrom = ivs[0].chrom
        strand = ivs[0].strand
        start = min([iv.start for iv in ivs])
        end = max([iv.end for iv in ivs])
        outhandle.write("\t".join(map(str, [
            chrom, start, end, uid, "|".join(uid2names[uid]), strand]))+"\n")
        
    if options.outfile:
        outhandle.close()
    

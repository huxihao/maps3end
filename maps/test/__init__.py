#!/usr/bin/env python
import os
from StringIO import StringIO
from maps.utils import generator_bed_from_blocks, bed2bam, string2file

test_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(test_dir, "data")
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

 
f_chrominfo = os.path.join(data_dir, "chromInfo")
fn_track = os.path.join(data_dir, "generead.track")
fn_gene_bed = os.path.join(data_dir, "gene.bed")
fn_read_bed = os.path.join(data_dir, "read.bed")
fn_bam = os.path.join(data_dir, "read.bam")

chrominfo = """\
chr1\t100000
chr2\t100000"""

genes = """\
100,100,100,100,100,100,\t0,200,400,600,800,1000\t1000
100,100,100,100,100,100,\t0,200,400,600,800,1300\t1000
"""

reads = """\
20\t0\t2000
20\t0\t2300
10,10\t0,110\t1890
10,10\t0,410\t1890
20\t0\t1400
20\t0\t1600
20\t0\t1750
10,10\t0,80\t1690
"""

if not os.path.exists(f_chrominfo):
    string2file(chrominfo, f_chrominfo)

fwd_gene = ["\t".join(map(str, g)) for g in 
    generator_bed_from_blocks(StringIO(genes), name_prefix='g')]

fwd_read = ["\t".join(map(str, r)) for r in 
    generator_bed_from_blocks(StringIO(reads), name_prefix='r')]

neg_gene = ["\t".join(map(str, g)) for g in 
    generator_bed_from_blocks(StringIO(genes), name_prefix='gn', strand='-')]

neg_read = ["\t".join(map(str, r)) for r in 
    generator_bed_from_blocks(StringIO(reads), name_prefix='rn', strand='-')]

if not os.path.exists(fn_gene_bed):    
    with open(fn_gene_bed, "w") as oh_gene: 
        oh_gene.write("\n".join(fwd_gene + neg_gene)+"\n")

if not os.path.exists(fn_read_bed):
    with open(fn_read_bed, "w") as oh_read: 
        oh_read.write("\n".join(fwd_read + neg_read)+"\n")

if not os.path.exists(fn_bam):
    bed2bam(fn_read_bed, f_chrominfo, os.path.join(data_dir, "read"), 
            is_bed12=True)
    
if not os.path.exists(fn_track):    
    with open(fn_track, "w") as oh_track: 
        oh_track.write("""\
browser position chr1:890-2500
browser hide all
track name="GenesFwd" description="Genes Forward" visibility=pack color=255,0,0
%s
track name="ReadsFwd" description="Reads Forward" visibility=pack color=255,0,0
%s
track name="GenesRev" description="Genes Reverse" visibility=pack color=0,0,255
%s
track name="ReadsRev" description="Reads Reverse" visibility=pack color=0,0,255
%s
""" % ("\n".join(fwd_gene), "\n".join(fwd_read), "\n".join(neg_gene), "\n".join(neg_read), )) 


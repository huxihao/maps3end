#!/bin/bash

f_gene=mm9refGene_uid.bed
../scripts/bed_rename.py mm9refGene.bed '' ${f_gene}
f_cluster=refGene.cluster
bedToGenePred ${f_gene} refGene.genePred
clusterGenes -conflicted ${f_cluster} mm9 refGene.genePred
f_land=refGene.gtf
../scripts/gene2land.py --lmrna 300 -o ${f_land} ${f_cluster} ${f_gene}
#../scripts/land2exp.py infile.bam ${f_land} -u cluster_id -o output_count


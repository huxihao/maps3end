#!/bin/bash
rm -f data/*
./__init__.py
bedToGenePred data/gene.bed data/gene.genePred 
clusterGenes -conflicted data/gene.cluster no data/gene.genePred
../../scripts/gene2land.py --lmrna 300 -o data/gene.gtf data/gene.cluster data/gene.bed


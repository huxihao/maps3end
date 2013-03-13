#!/usr/bin/env python
"""
Cluster of genes like from clusterGenes(kent)
"""

################################################################################

class GeneCluster(object):
    """
    Cluster of genes
    """
    def __init__(self, chrom, start, end, strand, name, genes=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        if genes is None:
            self.genes = list()
        else:
            self.genes = genes
        self.intervals = list()

    def __str__(self):
        return "\t".join(map(str, [self.chrom, self.start, self.end,
            self.name, 0, self.strand]))

    def add_gene(self, gene):
        self.genes.append(gene)
        
    def add_interval(self, interval):
        self.intervals.append(interval)

    def set_start(self, start):
        self.start = min(self.start, start)

    def set_end(self, end):
        self.end = max(self.end, end)

    def get_name(self):
        return self.name
    
    def get_genes(self):
        return self.genes
    
    def get_chrom(self):
        return self.chrom
    
    def get_start(self):
        return self.start

    def get_end(self):
        return self.end
    
    def get_strand(self):
        return self.strand
    
    def overlaps(self, chrom, start, end, strand=None):
        """whether overlaps with other interval"""
        if (self.chrom != chrom 
            or min(self.end, end) - max(self.start, start) <= 0 
            or (strand is not None and self.strand != strand)): 
            return False
        return True
    
    def get_uniqtss_set(self):
        uniqtss = set()
        if self.strand == '+':
            uniqtss.add(self.start)
        else:
            uniqtss.add(self.end) # 1-based
        for (s, e) in self.intervals:
            if self.strand == '+':
                uniqtss.add(s)
            else:
                uniqtss.add(e)
        return uniqtss
        
################################################################################

def generator_cluster(infile):
    """
    Determine the clusters from program clusterGenes in kent.
    Taking the min(starts), max(ends) as the clustered interval.
    """
    genecls = None
    for line in open(infile):
        if line[0] == '#': continue
        
        cluster, table, gene, chrom, start, end, strand = \
            line.rstrip().split('\t')[:7]
        start = int(start)
        end = int(end)
        
        if genecls is None:
            genecls = GeneCluster(chrom, start, end, strand, cluster)

        if genecls.get_name() != cluster:
            yield genecls
            genecls = GeneCluster(chrom, start, end, strand, cluster)    
        else:
            genecls.set_start(start)
            genecls.set_end(end)
        
        genecls.add_gene(gene)
        genecls.add_interval((start, end))

    # the last cluster
    if genecls:
        yield genecls
    
################################################################################

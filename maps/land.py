#!/usr/bin/env python
"""
Extract potential landing intervals of reads from 3'-end RNA-seq from gene annotation
"""

from maps.io_utils.cluster import generator_cluster
from maps.io_utils.gene import reader_gene


class ClusterToLand(object):
    """
    Extract 3'-end exons within a given distance (LM/lmrna) to 3'-end of gene at mature mRNA level.
    3'-ends in a cluster are grouped into unique 3'-ends, which are assigned unique 3end_id.
    Generate exon intervals associated with unique id
    """
    def __init__(self, gcls, genes, lmrna=None):
        """lmrna: length on mRNA level"""
        self.gcls = gcls
        self.genes = genes
        self.lmrna = lmrna
    
    @property
    def strand(self):
        return self.gcls.get_strand()

    def emit(self):
        """Generate exon, 3end_id, concatenated gene names"""
        positions = self.uniq_keys()
        sortedpos = self.sort(positions)
        key2id = self.positions2id(sortedpos)
        for i in range(len(sortedpos)):
            pos = sortedpos[i]
            kid = key2id[pos]
            gids = self.key2gids(pos)
            for (s, e) in self.key2exons(pos):
                yield (s, e, kid, gids)
    
    def get_key(self, gene):
        """Return key position for a given gene (3'-end)"""
        if self.strand == "+":
            return gene.get_txEnd()
        else:
            return gene.get_txStart()

    def uniq_keys(self):
        """Get unique 3'-end positions"""
        return set([self.get_key(g) for g in self.genes])

    def positions2id(self, positions):
        """Group key positions and assign unique id"""
        key2id = {}
        for i in range(len(positions)):
            key2id[positions[i]] = "%s_%s" % (self.gcls.get_name(), i+1)
        
        return key2id

    def sort(self, positions):
        """sort positions in the sense direction"""
        positions = sorted(positions)
        if self.strand == '-':
            positions.reverse()

        return positions

    def key2gids(self, pos):
        """Return gene names/ids tuple associated with given key position"""
        return tuple([g.get_name() for g in self.genes if self.get_key(g) == pos])

    def key2exons(self, pos):
        """Extract exons associated with given key position"""
        uexons = set()
        for g in self.genes:
            if self.get_key(g) != pos:
                continue
            
            e_starts = map(int, g.get_exonStarts().rstrip(",").split(","))
            e_ends = map(int, g.get_exonEnds().rstrip(",").split(","))                
            exons = [(e_starts[i], e_ends[i]) for i in range(len(e_starts))]

            clen = 0
            if self.strand == '+': # from 3'-end to 5'-end
                exons.reverse()
            for i in range(len(exons)):
                s, e = exons[i]
                if self.lmrna is None:
                    uexons.add((s, e))
                    continue
                else:
                    if clen >= self.lmrna:
                        break
                    if self.strand == '+':
                        s = max(s, e - (self.lmrna - clen))
                    else:
                        e = min(e, s + (self.lmrna - clen))
                    
                    uexons.add((s, e))
                    clen += (e - s)

        return uexons


class LandMaker(object):
    """
    For 3'-end RNA-seq, make potential landing intervals of reads
    from clustered genes by kent program clusterGenes
    """
    def __init__(self, f_cluster, f_gene, lmrna=500):
        """lmrna: for reads spaning two exons within given distance to 3'-end"""
        self.f_cluster = f_cluster
        self.f_gene = f_gene
        self.lmrna = lmrna
        self.id2gene = LandMaker.read_gene(open(self.f_gene)) 
    
    @staticmethod
    def read_gene(fh):
        """Return dict of gene_id:gene object, require name fields are unique"""
        id2gene = {}
        for gene in reader_gene(fh):
            id2gene[gene.get_name()] = gene
        
        return id2gene
    
    def emit(self):
        """Emit intervals"""
        for gcls in generator_cluster(self.f_cluster):
            gids = gcls.get_genes()
            genes = [self.id2gene[gid] for gid in gids]
            c2l = ClusterToLand(gcls, genes, self.lmrna)
            for (s, e, kid, gids) in c2l.emit():
                yield (gcls.get_chrom(), s, e, gcls.get_strand(), kid, gids, gcls.get_name())




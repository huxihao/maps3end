#!/usr/bin/env python
"""
Support for reading Gene format
"""
from maps.io_utils.tabular import TableRow, TableReader, Header

################################################################################

class GeneRecord(TableRow):
    """
    A Gene Record having following fields (from UCSC db table)
    'name', 'chrom', 'strand', 'txStart','txEnd',
    'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds'
    """
    header = ['name', 'chrom', 'strand', 'txStart','txEnd',
        'cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
    no = 0
    
    def __init__(self, tbl_row):
        TableRow.__init__(self, tbl_row.reader, tbl_row.fields)
        GeneRecord.no += 1
        
    def get_name(self):
        return self.__getitem__('name')
    
    def get_chrom(self):
        return self['chrom']
    
    def get_txStart(self):
        return int(self['txStart'])

    def get_txEnd(self):
        return int(self['txEnd'])

    def get_strand(self):
        return self['strand']

    def get_cdsStart(self):
        return int(self['cdsStart'])

    def get_cdsEnd(self):
        return int(self['cdsEnd'])

    def get_exonCount(self):
        return int(self['exonCount'])

    def get_exonStarts(self):
        return self['exonStarts']

    def get_exonEnds(self):
        return self['exonEnds']
    
    def get_score(self):
        return GeneRecord.no

    def as_bed12(self):
        exon_starts = map(int, self.get_exonStarts().strip(',').split(','))
        exon_ends = map(int, self.get_exonEnds().strip(',').split(','))
        block_sizes = [exon_ends[i] - exon_starts[i] 
            for i in xrange(self.get_exonCount())]
        block_starts = [exon_starts[i] - self.get_txStart() 
            for i in xrange(self.get_exonCount())]
        
        return (self.get_chrom(), self.get_txStart(), self.get_txEnd(),
            self.get_name(), self.get_score(), self.get_strand(), self.get_cdsStart(),
            self.get_cdsEnd(), '0,0,0', self.get_exonCount(), 
            ",".join(map(str, block_sizes)),
            ",".join(map(str, block_starts)),
            )
        
    def get_exons(self):
        exon_starts = map(int, self.get_exonStarts().strip(',').split(','))
        exon_ends = map(int, self.get_exonEnds().strip(',').split(','))
        return [(exon_starts[i], exon_ends[i]) for i in range(self.get_exonCount())]

    def get_tss(self):
        if self.get_strand() == '+':
            return self.get_txStart()
        else:
            return self.get_txEnd() - 1
        
################################################################################

class BedGeneRecord(GeneRecord):
    """
    Gene record in Bed format
    """
    header = ['chrom', 'txStart','txEnd', 'name', 'score', 'strand',
        'cdsStart','cdsEnd', 'itemRgb', 
        'exonCount', 'blockSizes', 'blockStarts']
    
    def get_exonStarts(self): 
        exon_starts = [int(i) + self.get_txStart() 
            for i in self['blockStarts'].strip(',').split(',')]
        
        return ",".join(map(str, exon_starts))+','
    
    def get_exonEnds(self):
        exon_starts = [int(i) + self.get_txStart() 
            for i in self['blockStarts'].strip(',').split(',')]
        blockSizes = [int(i) 
            for i in self['blockSizes'].strip(',').split(',')]
        exon_ends = [base + offset 
            for base,offset in zip(exon_starts, blockSizes)]
        
        return ",".join(map(str, exon_ends))+','

    def get_score(self):
        return self['score']

################################################################################

class NCGeneRecord(GeneRecord):
    """
    Non-coding gene record
    """
    header = ['name', 'chrom', 'strand', 'txStart','txEnd',
         'exonCount','exonStarts','exonEnds']
    
    def get_cdsStart(self):
        return int(self['txStart'])

    def get_cdsEnd(self):
        return int(self['txStart'])
    
################################################################################

def reader_gene(fhd, format='bed'):
    """Read gene records"""
    if format == 'bed': # change record class according to format
        rowclass = BedGeneRecord
    elif format == 'ncgene':
        rowclass = NCGeneRecord
    else:
        rowclass = GeneRecord
        
    for tblrow in TableReader(fhd, 
         force_header = Header(rowclass.header),
         comment_lines_startswith = ["#", "track ", "browser"],
         return_comments=False, 
        ):
        if isinstance(tblrow, TableRow):
            yield rowclass(tblrow)            

################################################################################


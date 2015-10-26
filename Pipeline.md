Here is a simple list of steps to do MAPS data analysis.

## Decoding sequencing data ##
We generally multiplex 12 or more samples in one lane of HiSeq2000 Ilumina sequencing run. Each sample has a unique barcode sequence. After getting the Fastq file (f\_fastq), we can decode and split the data into multiple Fastq files per sample with script _decode.py_.
  1. Create a barcode configuration file as f\_barcode (see an example below)
  1. Get the start position of barcode in sequencing read (37 by default)
  1. Choose maximum number of mismatches (0 by default)
  1. Run script _decode.py_ to decode samples (see an example below)
  1. Run script _fastq2dropend.py_ to keep only the target sequence by providing the length to drop for each read

An example of f\_barcode: 4 samples with barcode sequence in 1st column, and sample in 2nd column (separated by \t, tab).
```
ACTG	WT1
GGAT	KO1
AGGA	WT2
CGTC	KO2
```
Decoding usage example:
```
decode.py [--startpos=37 –mismatch=0 –outprefix lane1] f_fastq f_barcode
```
This will generate Fastq files named as lane1\_WT1.fastq, etc.
```
fastq2dropend.py -d 7 -o WT1clean.fastq WT1.fastq
```
This will drop last 7 nt of reads in WT1.fastq and write reads to file WT1clean.fastq.


## Mapping sequencing reads ##

The user could choose different mappers to do mapping from reads to the locations where the reads are from. We generally use Bowtie to map reads to the genome, for example human genome hg18.  Here are the basic steps.
  1. Build genome index with program _bowtie-build_ from sample genome sequence
  1. Remove 3'-adaptor and/or polyA sequences at the end of read with program _[cutadapt](https://code.google.com/p/cutadapt/)_
  1. Keep reads longer than 18 nt with simple script
  1. Map reads to the reference and keep unique hits (see example below)

```
bowtie EBWT --trim5 4 -l25 -n2 -k1 --best --strata
```
EBWT is the reference index. We skip the first 4 nt of reads in mapping, as the beginning nucleotides from random primers may not match template, and allow two mismatches in the remaining first 25 nt seed sequence.

## Create reference genes ##

We generally use annotated genes from [UCSC genome browser](http://genome.ucsc.edu/) to compute their expression. As MAPS method theoretically generate exonic reads close to the gene's 3-end, we normally only count reads close to the 3'-end. And also, one gene may have multiple annotated 3'-ends, we count all of them for computing the gene's expression level. Users could compute expression on different 3'-ends separately to investigate APA switch problem. Below are basic steps to create the "landing exons" of reads for reference genes.

  1. Download refGene genomic coordinates in BED format with UCSC table browser
  1. Rename the gene IDs (column 4) to unique names (for example row id)
  1. Cluster isofoms/transcripts of genes with program _clusterGenes_
  1. Create the counting intervals for each gene with script _gene2land.py_

```
gene2land.py --lmrna 300 -o f_land f_cluster f_gene
```

f\_gene is the BED file from step 2 and f\_cluster is from step 3. The intervals ("landing exons" of reads) is set by option -lmrna. The output file f\_land is in GTF format.

## Filtering reads resulting from internal priming ##

We use a heuristic filter based on previous reported motifs to remove potential internal priming events. The downstream sequence (up to 300 nt) of each read is checked on the presence of polyA stretches (consecutive 8 As or 9 As in a 10 nt window). The polyA stretches can be scanned by the script _scan4astretch.py_ (see usage below) across the sequences of all annotated genes or the genome once and saved for later use. The output BED file can be compared with the mapped reads by the _intersectBed_ program in [BEDTools](https://github.com/arq5x/bedtools2) to filter out those potential internal priming events.

```
scan4astretch.py -o polyA.bed --ival gene.bed genome_dir
```
The directory genome\_dir contains the Fasta files of the sample genome sequences. The search can be limited in annotated genes with option --ival by provided the BED file of gene intervals.

## Compute gene expression ##

Finally, we use the script _land2exp.py_ to count reads and calculate gene expression in RPKM. Below is a usage example.

```
land2exp.py f_tagbam f_land -u cluster_id -o f_count
```

f\_tagbam is the file of mapped reads in BAM format. The option -u/--unit sets the computation mode for counting reads for a clustered gene (cluster\_id)  or a transcript (transcript\_id).
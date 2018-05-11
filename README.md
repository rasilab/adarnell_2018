# Code for analysis : Darnell, Subramaniam, and O'Shea 2018

This github repository includes jupyter notebooks (and associated static HTML files) for reproducing all plots in the paper, and R scripts for sequencing data analysis, in [scripts/](scripts)

## Table of Contents

- [Annotation files](#annotation-files)
- [Generating bowtie indices and annotation files](#generating-bowtie-indices-and-annotation-files)
- [Pipeline for processing ribosome profiling sequencing data](#pipeline-for-processing-ribosome-profiling-sequencing-data)

## Annotation files

The following bowtie indices / annotation files should be generated in new folders called `bowtie_indices` and `sequence_annotation_files` before running the script below. These files are not provided in the repository due to their large size, but instructions for generating them are included below.
 
* bowtie index for subtractive alignment of rRNA contaminants: `hg38.rrna`
* transcriptome rsem reference bowtie index for alignment: `genome_rsem/hg38.gencode.v24.rsem/bowtie`
* annotations in fai format for preparing transcriptome rsem reference bowtie index: `hg38.fa.masked`
* canonical CCDS susbet annotations for assigning reads after alignment: `gencode.v24.canonical_ccds_transcripts.20170315.gff3`

## Generating bowtie indices and annotation files

### Required software
- bowtie/1.1.1
- rsem/1.3.0
- samtools/1.4.0
- gzip/1.3.12
- R/3.3.3
- tophat/2.0.14
- python/2.7.13

### Download, unzip, and index hg38 genome annotations (hg38.fa.masked)
```bash
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.masked.gz .
```
```bash
gunzip hg38.fa.masked.gz
```
```bash
samtools faidx hg38.fa.masked
```

### Download Gencode annotations and extract unique CCDS annotations to a gff3 file (gencode.v24.canonical_ccds_transcripts.20170315.gff3)
```bash
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gff3.gz
```
```bash
gunzip gencode.v24.annotation.gff3.gz
```
```bash
Rscript make_unique_ccds_transcripts.R
```

### Download rRNA sequences from NCBI, write to fasta, and create human rRNA bowtie index
```python
from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'rasi1983@gmail.com'

refseqids = {'28S' : 'NR_003287.2',
             '18S' : 'NR_003286.2',
             '5.8S' : 'NR_003285.2',
             '5S' : 'NR_023363.1',
            }
```
```python
outputlist = []
for (name,ID) in refseqids.items():
    # identify ids by searching
    rec = Entrez.read(Entrez.esearch(db="nucleotide", term=ID))
    # retrieve full record for the found id
    handle = Entrez.efetch(db="nucleotide", id=rec["IdList"][0], rettype="fasta")
    output = SeqIO.read(handle, 'fasta')
    # clean up record to be minimalistic
    output.id = ID
    output.description = output.name = ''
    outputlist.append(output)
    
SeqIO.write(outputlist, 'sequence_annotation_files/rrna.fa', 'fasta')
```
```bash
/app/bowtie/1.1.1/bowtie-build sequence_annotation_files/rrna.fa hg38.rrna
```

### Prepare gencode transcriptome index (hg38.gencode.v24.rsem) with rsem for alignment
```bash
rsem-prepare-reference --bowtie-path /app/bowtie/1.1.1 --bowtie --gff3 sequence_annotation_files/gencode.v24.annotation.gff3 sequence_annotation_files/hg38.fa.masked bowtie
```


## Pipeline for processing ribosome profiling sequencing data

Starting with `.fastq` files.

### Required software
- bowtie/1.1.1
- cutadapt/1.4.1
- R/3.3.3
- rsem/1.3.0

### Make a directory for each sample

```bash
for SAMPLE in `ls rawdata/*.fq.gz | cut -f1 -d. | cut -f2 -d/`;
do
    mkdir processeddata/$SAMPLE; 
    echo $SAMPLE; 
done
```

### Trim reads


```bash
for SAMPLE in `ls rawdata/*.fq.gz | cut -f1 -d. | cut -f2 -d/`; 
do 
    cutadapt --adapter=AAAAAAAAAA --minimum-length=13 --discard-untrimmed --output processeddata/$SAMPLE/$SAMPLE.trim.fq rawdata/$SAMPLE.R1.fq.gz; 
done
```

### Subtractive alignment to rRNA


```bash
for SAMPLE in `ls processeddata/*/*.trim.fq | cut -f1 -d. | cut -f2 -d/`; 
do 
    bowtie --seedlen=23 --threads=8 --un processeddata/$SAMPLE/norrna.fq --sam bowtie_indices/hg38.rrna processeddata/$SAMPLE/$SAMPLE.trim.fq; 
done
```

### Align to transcriptome


```bash

for SAMPLE in `ls processeddata/*/*.norrna.fq | cut -f1 -d. | cut -f2 -d/`; 
do
    /n/osheafs1/LAB/alicia/bioinformatics/RSEM-1.3.0/rsem-calculate-expression --num-threads 8 --output-genome-bam --sort-bam-by-coordinate  processeddata/$SAMPLE/$SAMPLE.norrna.fq /n/osheafs1/LAB/alicia/bioinformatics/bowtie_indices/genome_rsem/hg38.gencode.v24.rsem/bowtie processeddata/$SAMPLE/gencode; 
done
```
### Calculate codon-specific ribosome density


```bash

for SAMPLE in `ls processeddata/*/*norrna* | cut -f1 -d. | cut -f2 -d/`; 
do 
    Rscript scripts/calculate_codon_ribo_density_150overhang_20170807.R $SAMPLE; echo $SAMPLE; 
done
```

### Run R script to calculate gene fold changes


```bash
for SAMPLE in `ls processeddata/*/*norrna* | cut -f1 -d. | cut -f2 -d/`; 
do 
    Rscript scripts/calculate_gene_fold_change_20170615.R $SAMPLE; echo $SAMPLE; 
done
```

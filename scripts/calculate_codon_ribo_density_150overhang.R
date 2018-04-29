###Imports
.libPaths(c("/n/home11/adarnell/R/x86_64-unknown-linux-gnu-library/3.3/", .libPaths()))
library(BSgenome.Hsapiens.UCSC.hg38) # genome sequence
library(GenomicFeatures)  # for working with annotations
library(GenomicAlignments)  # for working with bam files
library(rtracklayer)  # read/write gff3 files
source('scripts/functions.R')
library(tidyverse)  # for tab-data manipulations
library(magrittr) #pipe operations

args <- commandArgs(trailingOnly=TRUE)
samplename <- args[1]

### Import sequence, annotations, alignments
genome <- BSgenome.Hsapiens.UCSC.hg38
canonical <- import.gff3('/n/osheafs1/LAB/alicia/bioinformatics/sequence_annotation_files/gencode.v24.canonical_ccds_transcripts.20170315.gff3')
canonicaldf <- canonical %>% tbl_df()
bamfile <- BamFile(paste0('processeddata/',
                         samplename,
                         '/gencode.genome.sorted.bam'))

### Parameter choices
left.trim <- 12  # nt trimmed from 5' side of alignment
right.trim <- 12  # nt trimmed from 3' side of alignment
min.read.density <- 0.33  # reads / nt, only cds above this used for codon density
overhang <- 150 # nt around each codon for codon density calculation


### get length of each cds
cds.exons <- canonical[canonical$type == 'CDS']
# number for individual GRanges within each CDS
names(cds.exons) <- mcols(cds.exons)$exon_number
mcols(cds.exons) <- mcols(cds.exons)['transcript_id']
cds <- split(cds.exons, cds.exons$transcript_id)
cds.mcols <- (
  mcols(cds)
  %>% tbl_df()
  %>% mutate(
    transcript_id=names(cds),
    length=sum(width(cds))
  )
  %>% left_join(
        canonicaldf
        %>% filter(type == 'transcript')
        %>% select(transcript_id, gene_name),
        by='transcript_id'
      )
 %>% as.data.frame()
)
row.names(cds.mcols) <-  cds.mcols$transcript_id

### Calculate trimmed and weighted read coverage over genome
# retrieve alignments along with the RSEM posterior probability
alns <- readGAlignments(bamfile, param=ScanBamParam(tag='ZW'))
# trim reads by 12nt from left and right, see functions.R
# use trim5Prime for trimming to a single nt from 5' side
psites <- trimBothEnds(left.trim, right.trim, alns)
## expand out psites to GRangesList and then to a Dataframe to assign
## correct weights for spliced reads
psites.grl <- grglist(psites, use.mcols=TRUE)
mcols(psites.grl)$grlwidth <- sum(width(psites.grl))
psites.gr <- GRanges(as.data.frame(psites.grl, use.outer.mcols=TRUE))
psites.gr$weight <- psites.gr$ZW / psites.gr$grlwidth
# go to cds coords
psites.cds <- mapToTranscripts(psites.gr, cds)
seqlengths(psites.cds) <- cds.mcols[seqlevels(psites.cds), 'length']
mcols(psites.cds) = mcols(psites.gr)[mcols(psites.cds)$xHits,]
# calculate read coverage with correct weights and full chr lengths
cds.cvg <- coverage(
  psites.cds,
  width=seqlengths(psites.cds),
  weight=psites.cds$weight
)
# write cds counts to a tsv file
counts <- sum(cds.cvg) 
cds.mcols[names(counts), 'counts'] = counts
(cds.mcols
  %>% tbl_df()
  %>% arrange(desc(counts))
  %>% mutate(counts=round(counts, digits=1))
  %>% write_tsv(
          paste0(
            'processeddata/',
            samplename,
            '/cds.counts_overhang150.tsv'
          )
          )
)
# normalize by the mean read count for each CDS
mean.cds.cvg <- mean(cds.cvg)
cds.cvg <- cds.cvg / mean.cds.cvg
# subset to cds that pass min.read.density threshold
cds.cvg <- cds.cvg[names(mean.cds.cvg[mean.cds.cvg > min.read.density])]

# get sequences
cdsseq <- extractTranscriptSeqs(genome, cds)

codons <- names(GENETIC_CODE[GENETIC_CODE != '*'])
calculateCodonDensity <- function (codon) {
  codonlocs <- (
    vmatchPattern(codon, cdsseq)  # find codon location in all frames
    %>% unlist()  # convert from list to iranges
    %>% tbl_df()  # convert to tibble
    %>% rename(seqnames=names)  # rename column to prep for granges
    %>% filter(start %% 3 == 1)  # select only codons in frame 1
    %>% mutate(transcript_id=seqnames)  # add transcript_id
    %>% left_join(cds.mcols, by='transcript_id')
    %>% GRanges()  # convert to GRanges
    %>% resize(., width=2*overhang, fix='center')  # resize to get 150nt each side
  )

  # subset codons
  codonlocs.subset <- codonlocs[
    codonlocs$transcript_id %in% names(cds.cvg)  # only for tx in cvg vector
    & start(codonlocs) > 0  # codons > overhang nt from start
    & end(codonlocs) <= codonlocs$length  # codons > overhang nt from end
  ]
  
  codon.cvg <- cds.cvg[codonlocs.subset]  # extract coverage for subset codons
  codon.cvg.df <- codon.cvg %>% as.matrix() %>% tbl_df() # convert to dplyr
  colnames(codon.cvg.df) = c(seq(-overhang, -1), seq(1, overhang)) # name columns 
  codon.cvg.df %>% summarize_all(mean)  # calculate mean across all codon instances
}

codon.ribo.density <- lapply(codons, calculateCodonDensity)
names(codon.ribo.density) <- codons
codon.ribo.density <- codon.ribo.density %>%
  bind_rows(.id = 'codon') %>%
  mutate_if(is.numeric, . %>% round(4)) %>%
  write_tsv(paste0(
    'processeddata/', samplename, '/codon.ribo.density_overhang150.tsv'))

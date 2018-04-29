library(tidyverse)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(stringr)

sessionInfo()

# Gtf file with protein-coding gene models
gffFile <- '../sequence_annotation_files/gencode.v24.annotation.gff3'

# read in all gencode annotations
annotations <- import.gff3(gffFile)
annotationsdf <- tbl_df(annotations)


# subset transcripts that have ccds and are identified by APPRIS
# for each gene pick lowest number ccds
# for each ccds pick lowest number transcript
ccdstx <- (
  annotationsdf
  %>% filter(!is.na(ccdsid))
  %>% filter(str_detect(tag, "appris_principal"))
  %>% filter(gene_type == 'protein_coding')
  %>% filter(type == 'transcript')
  %>% mutate(ccdsnum = as.numeric(substring(ccdsid, 5)))
  %>% mutate(txnum = as.numeric(substring(ID, 5)))
  %>% group_by(gene_id)
  %>% filter(ccdsnum == min(ccdsnum))
  %>% filter(txnum == min(txnum))
)

### write canonical ccds tx annotations including genes to a new file
### used for rsem
ccds <- annotations[
(annotations$transcript_id %in% ccdstx$transcript_id) |
(annotations$gene_id  %in% ccdstx$gene_id & annotations$type == 'gene')]
export.gff3(ccds, 'gencode.v24.canonical_ccds_transcripts.20170315.gff3')

#####################################################################################################
# extract CDS GRanges for ccdstx
ccds <- annotations[annotations$transcript_id %in% ccdstx$transcript_id
                   & annotations$type == 'exon']
# group exons for each transcript
ccds <- split(ccds, ccds$transcript_id)

# get sequences for the ccds
ccdsseqs <- extractTranscriptSeqs(Hsapiens, ccds)
# give them transcript id as names
names(ccdsseqs) <- names(ccds)
# write to fasta file
writeXStringSet(ccdsseqs, file = 'gencode.v24.canonical_ccds_transcripts.20170315.fa', format='fasta')


#####################################################################################################
# extract CDS GRanges for ccdstx
ccds <- annotations[annotations$transcript_id %in% ccdstx$transcript_id
                   & annotations$type == 'CDS']
# group exons for each transcript
ccds <- split(ccds, ccds$transcript_id)

# get sequences for the ccds
ccdsseqs <- extractTranscriptSeqs(Hsapiens, ccds)
# give them transcript id as names
names(ccdsseqs) <- names(ccds)
# write to fasta file
writeXStringSet(ccdsseqs, file = 'gencode.v24.canonical_ccds.20170315.fa', format='fasta')

######################################################################################################
# extract remaining non-canonical annotations and write them to a file
nonccds <- annotations[!(annotations$transcript_id %in% ccdstx$transcript_id)
                      & annotations$type == 'exon' & annotations$gene_type == 'protein_coding' &
                        !str_detect(annotations$tag, "PAR")
                      ]
# group exons for each transcript
nonccds <- split(nonccds, nonccds$transcript_id)
# get sequences for the nonccds
nonccdsseqs <- extractTranscriptSeqs(Hsapiens, nonccds)
# give them transcript id as names
names(nonccdsseqs) <- names(nonccds)
# write to fasta file
writeXStringSet(nonccdsseqs, file = 'gencode.v24.noncanonical_transcripts.20170315.fa', format = 'fasta')

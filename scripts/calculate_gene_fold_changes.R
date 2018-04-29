suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(DESeq2))

######################################################################

# genes that do have these counts in at least 1 sample are discarded
readcutoff <- 100

# get count files
cds.counts.files <- list.files('processeddata',
                              pattern = 'cds.counts.tsv$',
                              recursive = TRUE,
                              full.names = TRUE)

# read counts
cds.counts <- cds.counts.files %>%
  lapply(. %>%
         read_tsv(col_types = c(gene_name = col_character())))

# assign sample names to each list element
names(cds.counts) <- cds.counts.files %>%
  lapply(. %>%
         str_match('/([^/]+)/') %>%
         extract2(., 2)) %>%
  make.names()

# join all sample data
cds.counts <- bind_rows(cds.counts, .id = 'sample') %>%
  spread(key = 'sample', value = 'counts') %T>%
  write_tsv('tables/cds.counts.20170615.tsv')

######################################################################
## DESeq2 run

samplepairfiles <- list.files('tables/', pattern = 'samplepairs_for_deseq2.+tsv',
                            full.names = TRUE)

samplepairlist <- samplepairfiles %>%
  lapply(. %>% read_tsv(col_types = c(col_character()))
         %>% mutate_all(make.names))

# gets the count matrix for each table of samplepairs
get_count_data <- function(samplepairs) {
  samples <- samplepairs %>%
    gather(key, sample) %>%
    select(sample) %>%
    distinct() %>%
    extract2(1) 

  subset <- cds.counts %>%
    # select only tx that pass the readcutoff in at least 1 sample
    filter_at(vars(-transcript_id, -gene_name, -length),
              any_vars(. > readcutoff)) %>%
    na.omit()

  countdata <- subset %>%
    select(one_of(samples)) %>%
    mutate_if(is.numeric, as.integer) %>% 
    as.data.frame()

  rownames(countdata) <- subset %>% select(transcript_id) %>% extract2(1)
  countdata
}

# gets the sample matrix for each count matrix
get_col_data <- function(countdata) {
  data.frame(row.names = colnames(countdata),
             sample = colnames(countdata))
}

# run DESeq2
run_deseq2 <- function(countdata, coldata) {
  ddsObject <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = coldata,
                                     design = ~sample)
  dds <- DESeq(ddsObject, betaPrior=TRUE)
}

# calculate fold-change between a single sample pair
calculate_fold_change <- function(deseq2object, sample1, sample2) {
  res <- deseq2object %>% 
    results(contrast = c("sample", sample1, sample2), 
            addMLE = TRUE) %>% 
    as_tibble() %>%
    rownames_to_column('transcript_id') %>% 
    select(transcript_id, lfcMLE, baseMean) %>%  
    mutate(lfcMLE = round(lfcMLE, 3),
           baseMean = round(baseMean, 0),
           sample1 = sample1, sample2 = sample2)
}

# repeat above for all samplepairs in a given table
calculate_all_fold_changes <- function(deseq2object, samplepairs) {
  samplelist1 <- samplepairs %>% select(sample1) %>% extract2(.,1)
  samplelist2 <- samplepairs %>% select(sample2) %>% extract2(.,1)

  lfc <- mapply(
    function(x, y) calculate_fold_change(deseq2object, x, y),
   samplelist1, samplelist2, SIMPLIFY = FALSE) %>%  
    bind_rows
}
  
# write fold change tables to output
write_fold_change <- function(n) {
  foldchangelist[[n]] %>%
    mutate(sample = paste(sample1 , "vs", sample2, sep = '.')) %>%
    select(-sample1, -sample2) %>%
    spread(sample, lfcMLE) %>% 
    # add gene name for easy comparison
    left_join((cds.counts %>%
               select(transcript_id, gene_name)), on = transcript_id) %>% 
    select(transcript_id, gene_name, everything()) %>% 
    write_tsv(paste0('tables/foldchange_samplepairs_', n, '.tsv'))
}

countdatalist <- lapply(samplepairlist, get_count_data)
coldatalist <- lapply(countdatalist, get_col_data)
deseq2objectlist <- mapply(run_deseq2, countdatalist, coldatalist,
                          SIMPLIFY = FALSE)
foldchangelist <- mapply(calculate_all_fold_changes,
                        deseq2objectlist, samplepairlist, SIMPLIFY = FALSE)
lapply(seq_along(foldchangelist), write_fold_change)

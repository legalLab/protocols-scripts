#!/usr/bin/env Rscript

# collection of eDNA functions
# Tomas Hrbek starting July 2019

##################################
# process NCBI local blast+ output
# function for extracting taxonomic results from BLAST output
# with sensible defaults
##################################

analyze_edna <- function(infile, path, lookup_table, taxa_to_filter = NA, upid = 1, lpid = .97, len = 99, n_reads = 4) {
  # generate a list of taxa to be removed
  # standard list plus any additional taxa passed to the function
  always_remove <- c("Homo;", "Pan;", "Bos;", "Pseudonovibos;", "Sus;", "Mus;", "Rattus;", "Felis;", "Canis;", "Gallus;", "Meleagris;")
  taxa_to_remove <- if(!is.na(taxa_to_filter[1])) {
    c(always_remove, taxa_to_filter)
  } else {
    always_remove
  }
  
  # check if file is not empty first - empty gz files are <=75 bytes
  # make a fake file with 1 sample "bob"
  if(file.size(paste0(path, infile)) <= 100L) {
    blast_table <- tibble(id = "bob", n = 1, freq = 0, taxonomy = "Bob;Bob;Bob;Bob;Bob;Bob;Bob bob")
    return(blast_table)
    break
  }
  
  blast <- read.table(gzfile(paste0(path, infile)), header = FALSE, sep = "\t")
  colnames(blast) <- c('qacc', 'sacc', 'evalue', 'pident', 'mismatch', 'length', 'qstart', 'qend', 'sstart', 'send')
  
  blast <- as_tibble(blast) %>%
    distinct(qacc, .keep_all = TRUE) %>%
    # making 100% matches 99.99% so that 100% will not be dropped (pident < upid)
    mutate(pident = if_else(pident == 100, 99.99, pident)) %>%
    mutate(pident = pident/100) %>%
    # filter based on a range of values; e.g. 97%-95% representing genus
    filter(pident < upid & pident >= lpid & length > len)
  
  #blast$sacc <- str_extract(blast$sacc, regex("[0-9_]+$"))
  # swarm files end with size=666; and not size=666 like dereped only files
  blast_table <- blast %>%
    mutate(n = as.numeric(str_extract(.$qacc, pattern = "(?<=;size=)[0-9]+"))) %>%
    group_by(sacc) %>%
    summarise(n = sum(n)) %>%
    filter(n > n_reads) %>%
    arrange(desc(n)) %>%
    rename(id = sacc) %>%
    left_join(lookup_table, by = 'id') %>%
    filter(!reduce(map(taxa_to_remove, ~str_detect(taxonomy, .x)), `|`)) %>%
    mutate(freq = n/sum(n)) %>%
    relocate(freq, .after = n) %>%
    add_row(id = "bob", n = 1, freq = 0, taxonomy = "Bob;Bob;Bob;Bob;Bob;Bob;Bob bob")
  
  return(blast_table)
}


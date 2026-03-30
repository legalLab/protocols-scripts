#!/usr/bin/env Rscript

# collection of eDNA functions
# Tomas Hrbek starting July 2019

##################################
# process NCBI local blast+ output
# function for extracting taxonomic results from BLAST output
# with sensible defaults
##################################

##################################
# correct taxonomy
##################################
correct_taxonomy <- function(blast_tbl, dedup_tbl) {
  # parse named dups into long format
  dedup_df <- purrr::map_dfr(seq_along(dedup_tbl), function(i) {
    entries <- str_split(dedup_tbl[i], ",")[[1]] %>%
      str_trim()
    parsed <- str_split(entries, ";", n = 2)
    tibble(
      cluster_id = i,
      id = map_chr(parsed, 1),
      taxonomy  = map_chr(parsed, 2)
    )
  })
  
  normalize_taxonomy_7 <- function(tax_vec) {
    sapply(tax_vec, function(x) {
      if (is.na(x) || x == "") return(NA_character_)
      parts <- strsplit(x, ";")[[1]]
      # ensure length = 7
      if (length(parts) < 7) {
        parts <- c(parts, rep("", 7 - length(parts)))
      } else if (length(parts) > 7) {
        parts <- parts[1:7]
      }
      paste(parts, collapse = ";")
    })
  }
  
  # compute corrected taxonomy per cluster
  dedup_tax <- dedup_df %>%
    group_by(cluster_id) %>%
    summarise(new_taxonomy = collapse_taxa(taxonomy), .groups = "drop")
  
  # map accession → corrected taxonomy
  tax_map <- dedup_df %>%
    left_join(dedup_tax, by = "cluster_id") %>%
    select(id, new_taxonomy)
  
  # replace taxonomy if accession exists in clusters
  blast_tbl_updated <- blast_tbl %>%
    left_join(tax_map, by = "id") %>%
    #mutate(taxonomy = ifelse(!is.na(new_taxonomy), new_taxonomy, taxonomy)) %>%
    mutate(new_taxonomy = normalize_taxonomy_7(new_taxonomy),
           taxonomy = dplyr::coalesce(new_taxonomy, taxonomy)) %>%
    select(-new_taxonomy)
  
  return(blast_tbl_updated)
}


##############################
# extract taxa in a given similarity range
# sensible suggestions
# >98% species
# 95-98% genus
# 90-95% family
##############################
analyze_edna <- function(infile, path, lookup_table, dedup_table = NA, taxa_to_filter = NA, upid = 1, lpid = .98, len = 99, n_reads = 4) {
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
  
  blast_table <- as_tibble(blast) %>%
    distinct(qacc, .keep_all = TRUE) %>%
    # making 100% matches 99.99% so that 100% will not be dropped (pident < upid)
    mutate(pident = if_else(pident == 100, 99.99, pident)) %>%
    mutate(pident = round(pident/100, 4)) %>%
    # filter based on a range of values and length; e.g. 97%-95% representing genus
    filter(pident < upid & pident >= lpid & length > len) %>%
    # swarm files end with size=666;
    # derep files end with size=666
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
  
  # correct taxonomy if dedup table exists
  if(!is.na(dedup_table)[1]) {
    blast_table <- correct_taxonomy(blast_table, dedup_table)
  }
  
  return(blast_table)
}

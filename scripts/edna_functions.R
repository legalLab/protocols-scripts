#!/usr/bin/env Rscript

# collection of eDNA functions
# Tomas Hrbek starting July 2019

##################################
# process NCBI local blast+ output
# function for extracting taxonomic results from BLAST output
# with sensible defaults
# will extract any similarity range
##################################

##################################
# read fasta into a dataframe
##################################
read_fasta_df <- function(file) {
  lines <- readLines(file)
  header_idx <- grep("^>", lines)
  ids <- sub("^>", "", lines[header_idx])
  
  seqs <- mapply(function(start, end) {
    paste(lines[(start + 1):(end - 1)], collapse = "")
  },
  start = header_idx,
  end = c(header_idx[-1], length(lines) + 1), SIMPLIFY = TRUE)
  
  seq_tbl <- tibble::tibble(id = ids, sequence = seqs)
  return(seq_tbl)
}

##################################
# write dataframe as fasta
##################################
write_fasta_df <- function(df, outfile) {
  lines <- c(rbind(paste0(">", df$id), df$sequence))
  writeLines(lines, outfile)
}

##################################
# convert dataframe to DNAbin
##################################
tibble_to_DNAbin <- function(df) {
  seq_list <- lapply(df$sequence, function(x) strsplit(x, "")[[1]])
  names(seq_list) <- df$id
  ape::as.DNAbin(seq_list)
}

##################################
# calculate diversity indices
##################################
get_div_indices <- function(df) {
  all_groups <- list()
  samples <- select(df, -1) %>%
    colnames()
  for (sample in samples) {
    x <- df[[sample]]
    N <- sum(x, na.rm = TRUE)
    S <- sum(x > 0, na.rm = TRUE)
    p <- x / N
    shannon <- -sum(p * log(p), na.rm = TRUE)
    simpson <- 1 - sum(p^2, na.rm = TRUE)
    inv_simpson <- 1 / sum(p^2, na.rm = TRUE)
    pielou <- ifelse(S > 0, shannon / log(S), NA)
    berger <- max(p, na.rm = TRUE)
    all_groups[[sample]] <- list(S = S, Shannon = shannon, Simpson = simpson, InvSimpson = inv_simpson, Pielou = pielou, Berger = berger)
  }
  return(all_groups)
}


##################################
# compute MRCA + collapsed descendant
##################################
collapse_taxa <- function(taxa_vec) {
  split_taxa <- strsplit(taxa_vec, ";")
  
  # find minimum length
  min_len <- min(lengths(split_taxa))
  
  # find longest common prefix index
  prefix_len <- 0
  for (i in seq_len(min_len)) {
    vals <- sapply(split_taxa, `[`, i)
    if (length(unique(vals)) == 1) {
      prefix_len <- i
    } else {
      break
    }
  }
  
  # if fully identical → return original
  if (prefix_len == max(lengths(split_taxa)) && length(unique(taxa_vec)) == 1) {
    return(taxa_vec[1])
  }
  
  # extract prefix
  prefix <- if (prefix_len > 0) {
    paste(split_taxa[[1]][1:prefix_len], collapse = ";")
  } else {
    ""
  }
  
  # extract differing immediate descendants
  suffixes <- sapply(split_taxa, function(x) {
    if (length(x) > prefix_len) x[prefix_len + 1] else NA
  }) %>%
    unique() %>%
    na.omit()
  
  collapsed_suffix <- paste(suffixes, collapse = "/")
  
  if (prefix == "") {
    return(collapsed_suffix)
  } else {
    return(paste(prefix, collapsed_suffix, sep = ";"))
  }
}

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

##############################
# extract all taxonomic levels
# extract all associated FASTAs
# extract taxa (same as analyze_edna)
##############################
extract_edna <- function(infile, path, lookup_table, dedup_table = NA, taxa_to_filter = NA, path_raw = NA, upid = 1, lpid = .98, len = 99, n_reads = 4) {
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
    blast_table_big <- tibble(id = "bob", n = 1, freq = 0, taxonomy = "Bob;Bob;Bob;Bob;Bob;Bob;Bob bob",
                              qacc = "bob", pident = 1)
    seqs_subset <- matrix(c("a", "c", "t", "g"), nrow = 1)
    rownames(seqs_subset) <- "bob"
    seqs_subset <- ape::as.DNAbin(seqs_subset)
    named_seqs <- tibble(id = "bob", taxon = "Bob bob", pident = 1)
    return(blast_table, blast_table_big, seqs_subset, named_seqs)
    break
  }
  
  blast <- read.table(gzfile(paste0(path, infile)), header = FALSE, sep = "\t")
  colnames(blast) <- c('qacc', 'sacc', 'evalue', 'pident', 'mismatch', 'length', 'qstart', 'qend', 'sstart', 'send')
  
  # if not specified, reconstruct raw data folder
  if(is.na(path_raw)) {
    path_raw <- basename(path) %>%
      str_remove("_[^_]+$") %>%
      paste0(dirname(path), '/', ., '/')
  }
  
  # extract base infile
  infile_base <- infile %>%
    str_remove("\\..*")
  
  # big table with all matched qacc
  blast_table_big <- as_tibble(blast) %>%
    distinct(qacc, .keep_all = TRUE) %>%
    # making 100% matches 99.99% so that 100% will not be dropped (pident < upid)
    mutate(pident = if_else(pident == 100, 99.99, pident)) %>%
    mutate(pident = round(pident/100, 4)) %>%
    # filter on minimal length
    filter(length > len) %>%
    select(c(qacc, sacc, pident)) %>%
    relocate(id = sacc, qacc, pident) %>%
    # swarm files end with size=666;
    # derep files end with size=666
    mutate(n = as.numeric(str_extract(.$qacc, pattern = "(?<=;size=)[0-9]+"))) %>%
    arrange(desc(n)) %>%
    left_join(lookup_table, by = 'id') %>%
    filter(!reduce(map(taxa_to_remove, ~str_detect(taxonomy, .x)), `|`)) %>%
    mutate(freq = n/sum(n)) %>%
    relocate(n, freq, taxonomy, .after = id) %>%
    add_row(id = "bob", n = 1, freq = 0, taxonomy = "Bob;Bob;Bob;Bob;Bob;Bob;Bob bob",
            qacc = "bob", pident = 1)
  
  # correct taxonomy if dedup table exists
  if(!is.na(dedup_table)[1]) {
    blast_table_big <- correct_taxonomy(blast_table_big, dedup_table)
  }
  
  # standard table just with all matched taxa based on pid
  blast_table <- blast_table_big %>%
    filter(!id == "bob") %>%
    # filter based on a range of values; e.g. 97%-95% representing genus
    filter(pident < upid & pident >= lpid) %>%
    group_by(id) %>%
    summarise(n = sum(n)) %>%
    filter(n > n_reads) %>%
    arrange(desc(n)) %>%
    left_join(lookup_table, by = 'id') %>%
    mutate(freq = n/sum(n)) %>%
    relocate(freq, .after = n) %>%
    add_row(id = "bob", n = 1, freq = 0, taxonomy = "Bob;Bob;Bob;Bob;Bob;Bob;Bob bob")
  
  # correct taxonomy if dedup table exists
  if(!is.na(dedup_table)[1]) {
    blast_table <- correct_taxonomy(blast_table, dedup_table)
  }
  
  # filter by minimal read number
  blast_table_big <- blast_table_big %>%
    filter(n > n_reads) %>%
    mutate(pident = if_else(pident == 0.9999, 1, pident))
  
  # get FASTA of retained qacc
  # vector of sequence IDs to keep
  ids <- blast_table_big %>%
    filter(!id == "bob") %>%
    pull(qacc)
  
  # read fasta
  seqs <- read_fasta_df(gzfile(paste0(path_raw, infile_base, ".fn.gz")))
  
  # subset sequences
  seqs_subset <- seqs[seqs$id %in% ids, ] %>%
    mutate(id = str_extract(.$id, "^[^;]+"))
  #seqs_subset <- seqs[rownames(seqs) %in% ids, ]
  
  # get named qaccs
  named_seqs <- blast_table_big %>%
    # user input or default to 0.98?
    #filter(pident >= lpid) %>%
    # get the deepest named taxon
    mutate(
      taxon = map_chr(strsplit(taxonomy, ";"), ~ {
        x <- .x[.x != ""]
        if (length(x) == 0) NA_character_ else tail(x, 1)
      })
    ) %>%
    select(c(qacc, taxon, pident)) %>%
    mutate(id = str_extract(.$qacc, "^[^;]+")) %>%
    relocate(id) %>%
    select(-qacc)
  
  return(list(blast_table, blast_table_big, seqs_subset, named_seqs))
}

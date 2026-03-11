#!/usr/bin/env Rscript

# flu_gdb_genbank_pipeline.R ---------------------------------------------------
# v1.1 – 2026-01-27 
# Pipeline to generate:
#   1) discarded.rds
#   2) transposed.rds
#   3) metadata_clusters_segments.rds
#   4) cluster_glue.rds (not required)
#   5) cluster_reference.rds
#
# Required inputs (relative to ROOT_DIR):
#   1) current_version/IAV_DB_matrix.tsv
#   2) current_version/clusters_trees/proteins_AA/sgt_[1-8]_*_AA_cluster_rep.fasta
#   3) current_version/clusters_trees/cluster_members/seg*_clusters.tsv
#   4) flu_gdb_app_data/discarded copy.rds
#   5) flu_gdb_app_data/taxonomic_names.csv
#   6) flu_gdb_app_data/reference_set_list.tsv
#   7) flu_gdb_app_data/strains_info.tsv

# Setup ------------------------------------------------------------------------
options(warn = 1, stringsAsFactors = FALSE)

## Libraries -------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(tools)
library(ape)
library(seqinr)
library(janitor)
library(taxize)
library(phytools)

## Paths ----------------------------------------------------------------------

ROOT_DIR <- "/home/laura/IAV_DB/Flu-Mutation-Explorer/"

DATA_DIR            <- file.path(ROOT_DIR, "flu_gdb_app_data") 

CURRENT_VERSION_DIR <- file.path(ROOT_DIR, "current_version")
MATRIX_FILE         <- file.path(CURRENT_VERSION_DIR, "IAV_DB_matrix_filtered.tsv")
CLUSTERS_DIR        <- file.path(CURRENT_VERSION_DIR, "clusters_trees")
PROTEINS_AA_DIR     <- file.path(CLUSTERS_DIR, "proteins_AA")


## write outputs directly to ROOT_DIR
OUTPUT_DIR <- ROOT_DIR

# OUTPUT_DIR <- file.path(ROOT_DIR, "flu_gdb_app_data")
# if (!dir.exists(OUTPUT_DIR)) {
#   dir.create(OUTPUT_DIR, recursive = TRUE)
# }

message("ROOT_DIR: ", ROOT_DIR)
message("OUTPUT_DIR: ", OUTPUT_DIR)

if (!dir.exists(CLUSTERS_DIR)) {
  stop("CLUSTERS_DIR does not exist: ", CLUSTERS_DIR)
}
if (!file.exists(MATRIX_FILE)) {
  stop("MATRIX_FILE does not exist: ", MATRIX_FILE)
}


## Define helper functions -----------------------------------------------------

# Function to rename sequence lists from file paths to segment names
rename_sequences <- function(seqs, to_replace = "_AA_cluster_rep") {
seqs %<>% 
  names %>% 
  basename %>% 
  stringr::str_extract("\\d+") %>%   
  setNames(seqs, .) 
}

# Function to get amino acids at the consensus position from a gapless position.
consensus <- function(segment, sequence_id, gapless_position) {
  message("Segment: ", segment)
  message("Sequence ID: ", sequence_id)
  message("Gapless position: ", gapless_position)
  
  # get consensus position from gapless sequence position
  consensus_position <- 
    discarded %>% 
    pluck(segment, sequence_id) %>% # get sequence
    .[gapless_position] %>% # get amino acid at gapless position
    names %>% # get consensus position as a string
    as.integer # convert to index
  
  message("Consensus position: ", consensus_position)
  
  # get amino acid at consensus position for all sequences in segment
  transposed %>% 
    pluck(segment, consensus_position) %>% 
    enframe(name = "Sequence ID", value = "Amino Acid")
}


# STEP 1 – Load trees and FASTA alignments -------------------------------------
tree_files <- list.files(
  path       = CLUSTERS_DIR,
  pattern    = "\\.nwk$",
  full.names = TRUE
)

# Midpoint-root + ladderize, then overwrite the original .nwk files
for (i in seq_along(tree_files)) {
  tr <- ape::read.tree(tree_files[i])
  tr_rooted <- phytools::midpoint.root(tr)
  tr_rooted <- ape::ladderize(tr_rooted, right = FALSE)
  tmp <- paste0(tree_files[i], ".tmp")
  ape::write.tree(tr_rooted, file = tmp)
  file.rename(tmp, tree_files[i])
}


fa_files <- list.files(
  path       = PROTEINS_AA_DIR,
  pattern    = "^sgt_[1-8]_AA_cluster_rep\\.fasta$",
  full.names = TRUE
)

# STEP 2 – Filter sequences and build gapless/transposed objects ---------------
##   - Produce:
##       - discarded.rds
##       - transposed.rds

sequences <- 
  fa_files %>% 
  purrr::set_names() %>% 
  purrr::map(function(x) {
    fa <- read.fasta(file = x, seqtype = "AA") 
    
    seqs <- purrr::map(fa, function(y){
      y %>% 
        getSequence %>% 
        setNames(seq_along(.))
    }) 
  })

sequences %<>% rename_sequences

# Create sequence list with gaps discarded retaining the original positions as names
discarded <- sequences %>%
  purrr::map(function(seg) {
    purrr::map(seg, function(seq_vec) {
      seq_vec %>% purrr::discard(grepl("-", .))
    })
  })


transposed <- purrr::map(discarded, function(segment_data) {
  all_positions <- sort(unique(as.integer(unlist(lapply(segment_data, names)))))
  all_positions_chr <- as.character(all_positions)

  result <- lapply(all_positions_chr, function(pos_str) {
    aa_at_pos <- character(0)
    for (seq_name in names(segment_data)) {
      seq_data <- segment_data[[seq_name]]
      if (pos_str %in% names(seq_data)) {
        aa_at_pos[seq_name] <- unname(seq_data[[pos_str]])
      }
    }
    aa_at_pos
  })

  names(result) <- all_positions_chr
  result
})


# Save required outputs
readr::write_rds(discarded,  file.path(OUTPUT_DIR, "discarded.rds"))
readr::write_rds(transposed, file.path(OUTPUT_DIR, "transposed.rds"))

# STEP 3 - Read GenBank matrix, filter to kept sequences, add host taxonomy ----
##   - Produce:
##       - metadata_clusters_segments.rds

genbank_metadata <- readr::read_tsv(MATRIX_FILE, show_col_types = FALSE) %>%
  janitor::clean_names()

unique_ids <- discarded %>%
  purrr::map(names) %>%
  unlist() %>%
  unique()


genbank_filtered <- genbank_metadata %>%
  dplyr::filter(primary_accession %in% unique_ids) %>%
  dplyr::mutate(
    serotype_validated = stringr::str_to_upper(serotype_validated),
    serotype_validated = dplyr::na_if(serotype_validated, "MIXED"),
    serotype_validated = dplyr::na_if(serotype_validated, "UNKNOWN")
  ) %>%
  mutate(
    h_subtype = stringr::str_extract(serotype_validated, "^H[0-9]+"),
    n_subtype = stringr::str_extract(serotype_validated, "N[0-9]+")
  ) %>%
  dplyr::select(
    primary_accession,
    accession_version,
    h_subtype,
    n_subtype,
    dplyr::starts_with("host"),
    parsed_strain,
    segment_validated
  )


## generate taxonomic_names.csv for version updates only ----------------------
# TO-DO: generate separated script
# (run using interactive mode for ncbi taxonomy assignment)

# hosts <- genbank_metadata %>%
#   distinct(host_validated) %>% 
#   arrange(host_validated)
# 
# taxonomic_names <- tax_name(sci = hosts$host_validated, get = c("class", "order"), db = "ncbi")
# # correct mislabel
# taxonomic_names <- taxonomic_names %>%
#   mutate(
#     order = case_when(
#       query == "yellow-legged gull"        ~ "Charadriiformes",
#       query == "common gull"               ~ "Charadriiformes",
#       query == "peacock"                   ~ "Galliformes",
#       query == "unidentified influenza virus" ~ NA_character_,
#       TRUE                                 ~ order
#     ),
#     class = case_when(
#       query %in% c("yellow-legged gull","common gull","peacock") ~ "Aves",
#       query == "unidentified influenza virus"                    ~ NA_character_,
#       TRUE                                                       ~ class
#     )
#   )
#
# # asign host_group
# taxonomic_names <- taxonomic_names %>%
#   mutate(
#     host_group = case_when(
#       query == "Homo sapiens" ~ "Human",
#       order %in% c("Primates","Carnivora","Eulipotyphla","Chiroptera","Artiodactyla",
#                    "Perissodactyla","Rodentia","Pilosa","Lagomorpha") ~ "Other Mammals",
#       class == "Aves" ~ "Birds",
#       query == "environmental samples" ~ "Environment",
#       TRUE ~ "Unknown"
#     )
#   )
# 
# # save 
# taxonomic_names %>% write_csv(file.path(DATA_DIR, "taxonomic_names.csv"))


# Load precomputed host taxonomy
taxonomic_names <- readr::read_csv(file.path(DATA_DIR,"taxonomic_names.csv"),
  show_col_types = FALSE)


genbank_filtered <- genbank_filtered %>%
  left_join(
    taxonomic_names %>% dplyr::select(-db),
    by = c("host_validated" = "query")
  )

# rename column
genbank_filtered <- genbank_filtered %>%
  dplyr::rename(host_order = order)

genbank_filtered <- genbank_filtered %>%
  dplyr::rename(segment = segment_validated)


# Save required output
readr::write_rds(genbank_filtered,  file.path(OUTPUT_DIR, "metadata_clusters_segments.rds"))


# STEP 4 – Build clusters & segments metadata -----------------------------------------------------------
##   - Produce:
##       - cluster_glue.rds


# list files with cluster representatives (e.g. seg1_clusters.tsv)
CLUSTER_MEMBERS_DIR <- file.path(CLUSTERS_DIR, "cluster_members")

cluster_files <- list.files(
  path       = CLUSTER_MEMBERS_DIR,
  pattern    = "^seg\\d+_clusters\\.tsv$",
  full.names = TRUE
)


# df renamed: clusters -> cluster_rep
cluster_rep <-
  purrr::map(cluster_files, function(x){
    # get segment number from file name
    segment <- 
      x %>% 
      basename %>% 
      stringr::str_extract("\\d+") %>% 
      as.integer
    
    # read cluster file and add  segment number as first column
    readr::read_tsv(x, col_names = FALSE, show_col_types = FALSE, progress = FALSE) %>% 
      setNames(c("representative", "accession")) %>% 
      mutate(segment = segment, .before = 1)
  }) %>% list_rbind


# Load list of references sequences (for numbering scheme)

ref_set <- readr::read_delim(file.path(DATA_DIR,"reference_set_list.tsv"),
  show_col_types = FALSE) 


# Filtering join of representative sequences to reference cluster representatives + metadata
ref_set_cluster_rep <- ref_set %>%
  left_join(cluster_rep,
    by = c("accession", "segment")) %>%
  select(-accession) %>%
  distinct()

# add metadata
ref_set_cluster_rep <- ref_set_cluster_rep %>%
  inner_join(
    genbank_filtered %>% select(-segment),
    by = c("representative" = "primary_accession")
  ) %>% # create strain display
  mutate(
    h_num         = as.integer(stringr::str_extract(h_subtype, "\\d+")),
    n_num         = as.integer(stringr::str_extract(n_subtype, "\\d+")),
    strain        = parsed_strain,
    strain_display = stringr::str_c(strain, " (", h_subtype, n_subtype, ")")
  ) %>% # sort by segment, H and N subtypes
  arrange(segment, h_num, n_num) %>%
  # select only relevant columns
  select(
    segment,
    strain,
    strain_display,
    representative,
    h_subtype,
    h_num,
    n_subtype,
    n_num
  ) %>% distinct()

readr::write_rds(ref_set_cluster_rep, file.path(OUTPUT_DIR, "cluster_glue.rds"))


# STEP 5 – Build cluster_reference for display in the app (Subtype dropdown mennu) -------------------------------
##   - Produce:
##       - cluster_reference.rds

# Read reference strains list
reference_strains <- readr::read_delim(
  file.path(DATA_DIR, "strains_info.tsv"),
  show_col_types = FALSE,
  col_select = c(primary_accession, segment_validated, parsed_strain)
) 

# Define reference strains to force-include
reference <- c(
  "A/New_York/392/2004"
  # "A/Puerto_Rico/8/1934", "A/goose/Guangdong/1/1996"
)

# Filter the cluster representatives for strains on the reference list
cluster_reference <- dplyr::bind_rows(
  cluster_rep %>%
    dplyr::semi_join(
      reference_strains,
      by = c("accession" = "primary_accession",
             "segment"  = "segment_validated")
    ),
  reference_strains %>%
    dplyr::filter(parsed_strain %in% reference) %>%
    dplyr::transmute(
      segment        = segment_validated,
      representative = primary_accession,
      accession      = primary_accession
    )
) %>%
  dplyr::distinct(segment, representative, .keep_all = TRUE)

# Add metadata
cluster_reference <- cluster_reference %>%
  dplyr::inner_join(
    genbank_filtered %>% dplyr::select(-segment),
    by = c("representative" = "primary_accession")
  ) %>%
  dplyr::mutate(
    h_num          = as.integer(stringr::str_extract(h_subtype, "\\d+")),
    n_num          = as.integer(stringr::str_extract(n_subtype, "\\d+")),
    strain         = parsed_strain,
    strain_display = dplyr::if_else(is.na(h_subtype) | is.na(n_subtype),
    strain, stringr::str_c(strain, " (", h_subtype, n_subtype, ")"))
  ) %>%
  dplyr::arrange(segment, h_num, n_num) %>%
  dplyr::select(
    segment,
    strain,
    strain_display,
    representative,
    h_subtype,
    h_num,
    n_subtype,
    n_num
  )

readr::write_rds(cluster_reference, file.path(OUTPUT_DIR, "cluster_reference.rds"))


outputs <- c(
  "discarded.rds",
  "transposed.rds",
  "metadata_clusters_segments.rds",
  "cluster_glue.rds",
  "cluster_reference.rds"
)

for (f in outputs) {
  p <- file.path(OUTPUT_DIR, f)
  if (file.exists(p)) {
    sz <- file.info(p)$size / 1024^2
    message(sprintf("[OUTPUT] %s | %.2f MB", f, sz))
  } else {
    message(sprintf("[OUTPUT] MISSING: %s", f))
  }
}

message(sprintf("[COMPLETE] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))




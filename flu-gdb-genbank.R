#!/usr/bin/env Rscript

# flu_gdb_genbank.R ------------------------------------------------------------
# v1.2 – 2026-10-08
#
# Pipeline to generate (into APP_DIR):
#   1) sequence_metadata_full.rds    metadata for every sequence in the database (not used by the app)
#   2) product_positions.rds         residues per sequence (tree tips), all 12 products
#   3) product_columns.rds           the same, one entry per column
#   4) sequence_metadata.rds         metadata for the tree tips
#   5) trees/seg[1-8]_cluster_rep.nwk  midpoint-rooted trees, tips renamed to strains
#   6) subtype_references.rds        Subtype dropdown for the Tree tab
#   7) reference_numbering.rds       position <-> column for every reference and product
#   8) blastdb/sgt_[1-8].*           BLAST protein DBs of the tree tips (primary products only)
#   9) master_blastdb/reference_proteins_db.*  BLAST DB of the master reference proteins
#  10) reference_proteins_db.fasta   master reference proteins (copied from INPUT_DIR)
#  11) alignment_column_counts.rds   residue counts per column and host category, all sequences
#  12) cluster_column_residues.rds   residue counts per column and cluster
#  13) cluster_host_composition.rds  host composition per cluster, preformatted for display
#  14) IAV_DB_summary.log            database summary for the Home tab
#
# Required inputs (relative to ROOT_DIR):
#   1) <CURRENT_VERSION_DIR>/IAV_DB_matrix.tsv
#   2) <CURRENT_VERSION_DIR>/software_info.tsv
#   3) <CURRENT_VERSION_DIR>/alignments/proteins_AA/sgt_*_AA.fasta
#   4) <CURRENT_VERSION_DIR>/clusters_trees/seg[1-8]_cluster_rep.nwk  (tips = accessions)
#   5) <CURRENT_VERSION_DIR>/clusters_trees/cluster_members/seg[1-8]_clusters.tsv
#   6) <INPUT_DIR>/taxonomic_names.csv           (from flu-gdb-taxonomic-names.R)
#   7) <INPUT_DIR>/ref_set_numbering.txt
#   8) <INPUT_DIR>/reference_proteins_db.fasta   (master reference proteins)
#
# Requires makeblastdb (BLAST+) in the PATH.

# Setup ------------------------------------------------------------------------
options(warn = 1, stringsAsFactors = FALSE)

## Libraries -------------------------------------------------------------------
library(tidyverse)
library(tools)
library(ape)
library(janitor)
library(phytools)

## Paths ----------------------------------------------------------------------

ROOT_DIR <- getwd()

CURRENT_VERSION_DIR <- file.path(ROOT_DIR, "current_version")
INPUT_DIR           <- file.path(ROOT_DIR, "flu_gdb_inputs")
APP_DIR             <- file.path(ROOT_DIR, "flu_gdb_app_data")

MATRIX_FILE         <- file.path(CURRENT_VERSION_DIR, "IAV_DB_matrix.tsv")
SOFTWARE_INFO_FILE  <- file.path(CURRENT_VERSION_DIR, "software_info.tsv")

CLUSTERS_DIR        <- file.path(CURRENT_VERSION_DIR, "clusters_trees")
CLUSTER_MEMBERS_DIR <- file.path(CLUSTERS_DIR, "cluster_members")
ALIGNMENT_DIR       <- file.path(CURRENT_VERSION_DIR, "alignments", "proteins_AA")

TAXONOMIC_NAMES_FILE <- file.path(INPUT_DIR, "taxonomic_names.csv")
REF_SET_FILE         <- file.path(INPUT_DIR, "ref_set_numbering.txt")
MASTER_FASTA         <- file.path(INPUT_DIR, "reference_proteins_db.fasta")

# outputs go to a staging folder, published as APP_DIR at the end of the run
OUTPUT_DIR <- stringr::str_c(APP_DIR, "_new")
PREV_DIR   <- stringr::str_c(APP_DIR, "_prev")

TREES_DIR          <- file.path(OUTPUT_DIR, "trees")
BLASTDB_DIR        <- file.path(OUTPUT_DIR, "blastdb")
MASTER_BLASTDB_DIR <- file.path(OUTPUT_DIR, "master_blastdb")

message("ROOT_DIR: ", ROOT_DIR)
message("CURRENT_VERSION_DIR: ", CURRENT_VERSION_DIR)
message("INPUT_DIR: ", INPUT_DIR)
message("APP_DIR: ", APP_DIR)

for (f in c(MATRIX_FILE, SOFTWARE_INFO_FILE, TAXONOMIC_NAMES_FILE, REF_SET_FILE, MASTER_FASTA)) {
  if (!file.exists(f)) stop("Input file does not exist: ", f)
}
if (!dir.exists(CLUSTERS_DIR)) {
  stop("CLUSTERS_DIR does not exist: ", CLUSTERS_DIR)
}

MAKEBLASTDB <- Sys.which("makeblastdb")
if (MAKEBLASTDB == "") {
  stop("makeblastdb (BLAST+) not found in the PATH - install it in this environment ",
       "(e.g. conda install -c bioconda blast)")
}

# fresh staging folder
if (dir.exists(OUTPUT_DIR)) unlink(OUTPUT_DIR, recursive = TRUE)
for (d in c(OUTPUT_DIR, TREES_DIR, BLASTDB_DIR, MASTER_BLASTDB_DIR)) {
  dir.create(d, recursive = TRUE)
}

## Products (proteins) ---------------------------------------------------------

PRODUCT_FILES <- c(seg1 = "sgt_1_PB2_AA.fasta",
                   seg2 = "sgt_2_PB1_AA.fasta", seg2_PB1_F2 = "sgt_2_PB1_F2_AA.fasta",
                   seg3 = "sgt_3_PA_AA.fasta",  seg3_PA_X   = "sgt_3_PA_X_AA.fasta",
                   seg4 = "sgt_4_HA_AA.fasta",
                   seg5 = "sgt_5_NP_AA.fasta",  seg6        = "sgt_6_NA_AA.fasta",
                   seg7 = "sgt_7_M1_AA.fasta",  seg7_M2     = "sgt_7_M2_AA.fasta",
                   seg8 = "sgt_8_NS1_AA.fasta", seg8_NEP    = "sgt_8_NEP_AA.fasta")

PRODUCT_SEGMENT <- c(seg1 = 1,
                     seg2 = 2, seg2_PB1_F2 = 2,
                     seg3 = 3, seg3_PA_X = 3,
                     seg4 = 4,
                     seg5 = 5,
                     seg6 = 6,
                     seg7 = 7, seg7_M2 = 7,
                     seg8 = 8, seg8_NEP = 8)

PRIMARY_PRODUCTS <- stringr::str_c("seg", 1:8)


## Define helper functions -----------------------------------------------------

# Host categories to be excluded form the barplots;
EXCLUDED_HOSTS <- c("Unknown", "Environment", "Unclassified")

# Residues allowed in the counts (20 amino acids + ambiguity/stop)
VALID_RESIDUES <- c(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], "X", "?", "*", "B", "Z", "J")

# Build the position lookups for one product from its alignment:
#   discarded  - each sequence's non-gap residues, named by their alignment column
#   transposed - the same turned inside out, one entry per column
build_product <- function(aligned, accessions = names(aligned)) {
  aligned <- aligned[accessions]

  gapless <- purrr::map(as.character(aligned), function(sequence) {
    residues <- stringr::str_split_1(sequence, "")
    names(residues) <- seq_along(residues)
    residues[residues != "-"]
  })
  names(gapless) <- names(aligned)

  # trailing columns gapped in every sequence produce no entry and are not addressable
  width <- max(as.integer(unlist(purrr::map(gapless, names))))

  columns <- purrr::map(seq_len(width), function(column) {
    residues <- purrr::map_chr(gapless, function(sequence) {
      residue <- sequence[as.character(column)]
      if (is.na(residue)) NA_character_ else unname(residue)
    })
    residues[!is.na(residues)]
  })
  names(columns) <- seq_len(width)

  list(discarded = gapless, transposed = columns)
}

# Read an alignment, dropping repeated accessions
read_alignment <- function(path) {
  aligned <- Biostrings::readAAStringSet(path)
  names(aligned) <- stringr::str_remove(names(aligned), "\\s.*$")
  duplicates <- sum(duplicated(names(aligned)))
  if (duplicates > 0) {
    message(basename(path), ": ", duplicates, " duplicate accession(s) - keeping the first")
    aligned <- aligned[!duplicated(names(aligned))]
  }
  aligned
}

# Host percentages per cluster: <1% shown as "<1%", the rest as integers that add up
# to their rounded total (largest remainder)
host_percent <- function(n) {
  percent <- 100 * n / sum(n)
  shown   <- percent >= 1
  whole   <- floor(percent)
  missing <- round(sum(percent[shown])) - sum(whole[shown])
  top     <- order(ifelse(shown, percent - whole, -1), decreasing = TRUE)[seq_len(missing)]
  whole[top] <- whole[top] + 1
  dplyr::if_else(shown, stringr::str_c(whole, "%"), "<1%")
}

# STEP 1 - Read GenBank matrix, add host taxonomy ------------------------------
##   - Produce:
##       - sequence_metadata_full.rds

# Load precomputed host taxonomy
taxonomic_names <- readr::read_csv(TAXONOMIC_NAMES_FILE, show_col_types = FALSE)

sequence_metadata_full <- readr::read_tsv(MATRIX_FILE, show_col_types = FALSE) %>%
  janitor::clean_names()

# kept for the summary log and the tree tip checks
total_genbank  <- nrow(sequence_metadata_full)
total_curated  <- sum(is.na(sequence_metadata_full$exclusion))
exclusion_flags <- sequence_metadata_full %>%
  dplyr::filter(!is.na(exclusion)) %>%
  dplyr::select(primary_accession, segment = segment_validated, exclusion)

sequence_metadata_full <- sequence_metadata_full %>%
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

sequence_metadata_full <- sequence_metadata_full %>%
  left_join(
    taxonomic_names %>% dplyr::select(-db),
    by = c("host_validated" = "query")
  )

# rename columns
sequence_metadata_full <- sequence_metadata_full %>%
  dplyr::rename(host_order = order)

sequence_metadata_full <- sequence_metadata_full %>%
  dplyr::rename(segment = segment_validated)

readr::write_rds(sequence_metadata_full,  file.path(OUTPUT_DIR, "sequence_metadata_full.rds"))


# STEP 2 – Full protein alignments ---------------------------------------------

alignments <- purrr::map(PRODUCT_FILES, ~ read_alignment(file.path(ALIGNMENT_DIR, .x)))

# STEP 3 – Trees and tree tips -------------------------------------------------
# Raw IQ-TREE trees (tips = accessions): midpoint-root + ladderize. 
# Renamed trees are written in STEP 6.
tree_files <- list.files(
  path       = CLUSTERS_DIR,
  pattern    = "^seg[1-8]_cluster_rep\\.nwk$",
  full.names = TRUE
)
names(tree_files) <- stringr::str_extract(basename(tree_files), "\\d+")

if (!setequal(names(tree_files), as.character(1:8))) {
  stop("Expected seg1-8 trees in ", CLUSTERS_DIR, ", found: ",
       paste(basename(tree_files), collapse = ", "))
}

trees <- purrr::map(tree_files[as.character(1:8)], function(path) {
  tr <- ape::read.tree(path)
  tr <- phytools::midpoint.root(tr)
  ape::ladderize(tr, right = FALSE)
})

# Tree tips (representatives + references) per segment, in tree order. Every product
# of a segment uses this same order.
tree_tips <- purrr::map(trees, "tip.label")

tips_table <- purrr::imap(tree_tips, ~ tibble::tibble(primary_accession = .x,
                                                      segment = as.numeric(.y))) %>%
  purrr::list_rbind()

duplicated_tips <- purrr::keep(tree_tips, ~ any(duplicated(.x)))
if (length(duplicated_tips) > 0) {
  stop("Duplicated tip(s) in tree(s) of segment(s): ", paste(names(duplicated_tips), collapse = ", "))
}

# tips must be accessions of that segment in the matrix (raw, not renamed, trees)
missing_metadata <- tips_table %>%
  dplyr::anti_join(sequence_metadata_full, by = c("primary_accession", "segment"))

if (nrow(missing_metadata) > 0) {
  stop(nrow(missing_metadata), " tree tip(s) absent from the matrix (or under another segment) - ",
       "are the trees the raw IQ-TREE output with accession tips? e.g. ",
       paste(head(missing_metadata$primary_accession, 5), collapse = ", "))
}

for (product in names(PRODUCT_FILES)) {
  tips    <- tree_tips[[as.character(PRODUCT_SEGMENT[[product]])]]
  missing <- setdiff(tips, names(alignments[[product]]))
  if (length(missing) > 0) {
    stop(length(missing), " tree tip(s) of segment ", PRODUCT_SEGMENT[[product]],
         " absent from ", PRODUCT_FILES[[product]], ": ",
         paste(head(missing, 5), collapse = ", "))
  }
}

# tree tips flagged as excluded upstream
excluded_tips <- exclusion_flags %>%
  dplyr::semi_join(tips_table, by = c("primary_accession", "segment"))

if (nrow(excluded_tips) > 0) {
  warning(nrow(excluded_tips), " tree tip(s) carry an exclusion flag: ",
          paste0(excluded_tips$primary_accession, " (", excluded_tips$exclusion, ")", collapse = "; "),
          call. = FALSE)
}

# STEP 4 – Position lookups for every product ----------------------------------
##   - Produce:
##       - product_positions.rds, product_columns.rds (keyed seg1 ... seg8_NEP)
positions_all <- purrr::imap(PRODUCT_FILES, function(file, product) {
  product_data <- build_product(alignments[[product]],
                                tree_tips[[as.character(PRODUCT_SEGMENT[[product]])]])
  message(product, ": ", length(product_data$discarded), " tree tip(s), ",
          length(product_data$transposed), " addressable column(s)")
  product_data
})

product_positions <- purrr::map(positions_all, "discarded")
product_columns   <- purrr::map(positions_all, "transposed")

readr::write_rds(product_positions, file.path(OUTPUT_DIR, "product_positions.rds"))
readr::write_rds(product_columns,   file.path(OUTPUT_DIR, "product_columns.rds"))


# STEP 5 – Metadata for the tree tips ------------------------------------------
##   - Produce:
##       - sequence_metadata.rds
sequence_metadata <- sequence_metadata_full %>%
  dplyr::semi_join(tips_table, by = c("primary_accession", "segment"))

readr::write_rds(sequence_metadata,  file.path(OUTPUT_DIR, "sequence_metadata.rds"))


# STEP 6 – Trees with strain names ---------------------------------------------
##   - Produce:
##       - trees/seg[1-8]_cluster_rep.nwk
for (segment in names(trees)) {
  strains <- sequence_metadata %>%
    dplyr::filter(segment == as.numeric(!!segment)) %>%
    dplyr::select(primary_accession, parsed_strain) %>%
    tibble::deframe()

  tr <- trees[[segment]]
  tr$tip.label <- unname(strains[tr$tip.label])

  if (any(is.na(tr$tip.label)) || any(duplicated(tr$tip.label))) {
    stop("seg", segment, ": missing or duplicated strain names among the tree tips")
  }
  if (any(stringr::str_detect(tr$tip.label, "[(),:;\\s]"))) {
    stop("seg", segment, ": strain names with characters not allowed in newick")
  }

  ape::write.tree(tr, file = file.path(TREES_DIR, stringr::str_c("seg", segment, "_cluster_rep.nwk")))
}


# STEP 7 – Subtype references (Subtype dropdown menu) --------------------------
##   - Produce:
##       - subtype_references.rds

ref_set_numbering <- readr::read_tsv(REF_SET_FILE, show_col_types = FALSE) %>%
  dplyr::mutate(segment = as.numeric(segment))

available <- tips_table %>%
  dplyr::rename(accession = primary_accession)

missing_refs <- ref_set_numbering %>%
  dplyr::anti_join(available, by = c("segment", "accession"))

if (nrow(missing_refs) > 0) {
  warning("References missing from the tree tips: ",
          paste0(missing_refs$accession, " (segment ", missing_refs$segment, ")",
                 collapse = "; "),
          call. = FALSE)
}

cluster_reference <- ref_set_numbering %>%
  dplyr::semi_join(available, by = c("segment", "accession")) %>%
  dplyr::mutate(
    representative = accession, # colname to match expected, real reference accession
    h_subtype      = stringr::str_extract(serotype, "^H\\d+"),
    n_subtype      = stringr::str_extract(serotype, "N\\d+"),
    h_num          = as.integer(stringr::str_extract(h_subtype, "\\d+")),
    n_num          = as.integer(stringr::str_extract(n_subtype, "\\d+")),
    strain_display = dplyr::if_else(is.na(h_subtype) | is.na(n_subtype),
                                    strain,
                                    stringr::str_c(strain, " (", h_subtype, n_subtype, ")"))
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

readr::write_rds(cluster_reference, file.path(OUTPUT_DIR, "subtype_references.rds"))


# STEP 8 – Reference numbering -------------------------------------------------
##   - Produce:
##       - reference_numbering.rds

reference_numbering <- purrr::map(names(PRODUCT_FILES), function(product) {
  refs <- ref_set_numbering %>%
    dplyr::filter(segment == PRODUCT_SEGMENT[[product]],
                  accession %in% names(alignments[[product]]))

  if (nrow(refs) == 0) return(NULL)

  gapless <- build_product(alignments[[product]], refs$accession)$discarded

  purrr::pmap(refs, function(serotype, strain, segment, accession) {
    residues <- gapless[[accession]]
    tibble::tibble(product   = product,
                   segment   = as.integer(segment),
                   accession = accession,
                   strain    = strain,
                   serotype  = serotype,
                   position  = seq_along(residues),
                   column    = as.integer(names(residues)),
                   residue   = unname(residues))
  }) %>%
    purrr::list_rbind()
}) %>%
  purrr::list_rbind()

missing_numbering <- ref_set_numbering %>%
  dplyr::anti_join(reference_numbering %>% dplyr::distinct(segment, accession),
                   by = c("segment", "accession"))

if (nrow(missing_numbering) > 0) {
  warning("References absent from the full alignments: ",
          paste0(missing_numbering$accession, " (segment ", missing_numbering$segment, ")",
                 collapse = "; "),
          call. = FALSE)
}

attr(reference_numbering, "sources") <-
  tools::md5sum(c(REF_SET_FILE, file.path(ALIGNMENT_DIR, PRODUCT_FILES)))
attr(reference_numbering, "created") <- Sys.time()

readr::write_rds(reference_numbering,
                 file.path(OUTPUT_DIR, "reference_numbering.rds"), compress = "gz")


# STEP 9 – BLAST protein DBs of the tree tips ----------------------------------
##   - Produce:
##       - blastdb/sgt_[1-8]_ungapped.fasta + sgt_[1-8] BLAST DB files
for (product in PRIMARY_PRODUCTS) {
  segment <- PRODUCT_SEGMENT[[product]]
  db_name <- stringr::str_c("sgt_", segment)

  tips     <- tree_tips[[as.character(segment)]]
  ungapped <- Biostrings::AAStringSet(stringr::str_remove_all(as.character(alignments[[product]][tips]), "-"))
  names(ungapped) <- tips

  fasta <- file.path(BLASTDB_DIR, stringr::str_c(db_name, "_ungapped.fasta"))
  Biostrings::writeXStringSet(ungapped, fasta)

  status <- system2("makeblastdb",
                    c("-in", shQuote(fasta), "-dbtype", "prot", "-title", db_name,
                      "-out", shQuote(file.path(BLASTDB_DIR, db_name)), "-parse_seqids"),
                    stdout = FALSE)
  if (status != 0) {
    stop("makeblastdb failed for ", db_name)
  }
  message(db_name, ": BLAST DB of ", length(ungapped), " tree tip(s)")
}


# STEP 10 – BLAST DB of the master reference proteins --------------------------
##   - Produce:
##       - master_blastdb/reference_proteins_db.*
##       - reference_proteins_db.fasta
file.copy(MASTER_FASTA, file.path(OUTPUT_DIR, "reference_proteins_db.fasta"))

status <- system2("makeblastdb",
                  c("-in", shQuote(MASTER_FASTA), "-dbtype", "prot",
                    "-out", shQuote(file.path(MASTER_BLASTDB_DIR, "reference_proteins_db")),
                    "-title", "flu_reference_proteins"),
                  stdout = FALSE)
if (status != 0) {
  stop("makeblastdb failed for ", MASTER_FASTA)
}
message("master_blastdb: BLAST DB of ",
        length(Biostrings::readAAStringSet(MASTER_FASTA)), " reference protein(s)")


# STEP 11 – Amino acid counts per alignment column and host category -----------
##   - Produce:
##       - alignment_column_counts.rds

addressable <- product_columns

products <- names(PRODUCT_FILES) %>% purrr::keep(~ .x %in% names(addressable))

hosts <- sequence_metadata_full %>%
  dplyr::transmute(accession = primary_accession, host_group, host_order) %>%
  dplyr::distinct(accession, .keep_all = TRUE)


column_counts <- function(sequences, level, category, product, max_column) {
  counts <- Biostrings::consensusMatrix(sequences)
  counts <- counts[rownames(counts) != "-", , drop = FALSE]
  counts <- counts[, seq_len(min(ncol(counts), max_column)), drop = FALSE]

  nonzero <- which(counts > 0, arr.ind = TRUE)

  tibble::tibble(product       = product,
                 column        = as.integer(nonzero[, "col"]),
                 level         = level,
                 host_category = category,
                 amino_acid    = rownames(counts)[nonzero[, "row"]],
                 n             = as.integer(counts[nonzero]))
}

tally_product <- function(product) {
  sequences  <- alignments[[product]]
  max_column <- length(addressable[[product]])

  annotated <-
    tibble::tibble(accession = names(sequences), index = seq_along(sequences)) %>%
    dplyr::left_join(hosts, by = "accession")

  unplaced <- sum(is.na(annotated$host_group))
  if (unplaced > 0) {
    message(product, ": ", unplaced, " sequence(s) absent from the host table - dropped")
  }

  levels <-
    dplyr::bind_rows(
      annotated %>% dplyr::transmute(index, level = "host_group", host_category = host_group),
      annotated %>% dplyr::transmute(index, level = "host_order", host_category = host_order)
    ) %>%
    dplyr::filter(!is.na(host_category), !host_category %in% EXCLUDED_HOSTS)

  result <-
    levels %>%
    dplyr::group_split(level, host_category) %>%
    purrr::map(function(group) {
      column_counts(sequences[group$index],
                    level      = group$level[1],
                    category   = group$host_category[1],
                    product    = product,
                    max_column = max_column)
    }) %>%
    purrr::list_rbind()

  message(product, ": ", length(sequences), " sequence(s), ", max_column,
          " column(s), ", nrow(result), " count row(s)")

  result
}

alignment_column_counts <- purrr::map(products, tally_product) %>% purrr::list_rbind()

unknown_residues <- setdiff(unique(alignment_column_counts$amino_acid), VALID_RESIDUES)
if (length(unknown_residues) > 0) {
  stop("Unexpected residue(s) in the alignments: ", paste(sort(unknown_residues), collapse = ", "))
}

readr::write_rds(alignment_column_counts,
                 file.path(OUTPUT_DIR, "alignment_column_counts.rds"), compress = "gz")


# STEP 12 – Amino acid counts per alignment column and cluster -----------------
##   - Produce:
##       - cluster_column_residues.rds

members_by_segment <-
  sort(unique(PRODUCT_SEGMENT[products])) %>%
  purrr::set_names() %>%
  purrr::map(function(segment) {
    readr::read_tsv(file.path(CLUSTER_MEMBERS_DIR, stringr::str_c("seg", segment, "_clusters.tsv")),
                    col_names = c("representative", "member"),
                    col_types = "cc", progress = FALSE)
  })

tally_clusters <- function(product) {
  sequences  <- alignments[[product]]
  max_column <- length(addressable[[product]])
  members    <- members_by_segment[[as.character(PRODUCT_SEGMENT[[product]])]]

  indices <- split(match(members$member, names(sequences)), members$representative)
  indices <- purrr::map(indices, ~ .x[!is.na(.x)])

  missing_members <- nrow(members) - sum(purrr::map_int(indices, length))
  if (missing_members > 0) {
    message(product, ": ", missing_members, " cluster member(s) absent from the alignment - not counted")
  }

  indices <- purrr::keep(indices, ~ length(.x) > 0)

  counts <- purrr::imap(indices, function(rows, representative) {
    matrix <- Biostrings::consensusMatrix(sequences[rows])
    matrix <- matrix[rownames(matrix) != "-", , drop = FALSE]
    matrix <- matrix[, seq_len(min(ncol(matrix), max_column)), drop = FALSE]

    nonzero <- which(matrix > 0, arr.ind = TRUE)

    tibble::tibble(representative = representative,
                   column         = as.integer(nonzero[, "col"]),
                   amino_acid     = rownames(matrix)[nonzero[, "row"]],
                   n              = as.integer(matrix[nonzero]))
  })

  tallied <- dplyr::bind_rows(counts) %>% dplyr::mutate(product = product, .before = 1)

  message(product, ": ", nrow(tallied), " row(s) over ", length(indices), " cluster(s)")

  tallied
}

cluster_column_residues <- purrr::map(products, tally_clusters) %>% dplyr::bind_rows()

readr::write_rds(cluster_column_residues,
                 file.path(OUTPUT_DIR, "cluster_column_residues.rds"), compress = "gz")


# STEP 13 – Host composition per cluster ---------------------------------------
##   - Produce:
##       - cluster_host_composition.rds

cluster_hosts <-
  members_by_segment %>%
  purrr::imap(~ dplyr::mutate(.x, segment = as.integer(.y))) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(hosts, by = c("member" = "accession"))

unplaced <- sum(is.na(cluster_hosts$host_group))
if (unplaced > 0) {
  message("cluster hosts: ", unplaced, " member(s) absent from the host table - counted as Unknown")
  cluster_hosts$host_group[is.na(cluster_hosts$host_group)] <- "Unknown"
}

cluster_host_composition <-
  cluster_hosts %>%
  dplyr::count(segment, representative, host_group) %>%
  dplyr::arrange(segment, representative, dplyr::desc(n), host_group) %>%
  dplyr::summarise(
    cluster_size  = sum(n),
    cluster_hosts = stringr::str_flatten(stringr::str_c(host_group, " ", host_percent(n)), collapse = "; "),
    .by = c(segment, representative)
  ) %>%
  dplyr::mutate(cluster_members = stringr::str_c(cluster_size,
                                                 dplyr::if_else(cluster_size == 1, " sequence", " sequences"),
                                                 " in cluster:"),
                .after = cluster_size)

message("cluster host composition: ", nrow(cluster_host_composition), " cluster(s)")

readr::write_rds(cluster_host_composition,
                 file.path(OUTPUT_DIR, "cluster_host_composition.rds"), compress = "gz")


# STEP 14 – Database summary log -----------------------------------------------
##   - Produce:
##       - IAV_DB_summary.log
software_info <- readr::read_tsv(SOFTWARE_INFO_FILE, col_names = c("item", "value"),
                                 skip = 1, show_col_types = FALSE)

last_update <- software_info %>%
  dplyr::filter(stringr::str_trim(item) == "Time of creation") %>%
  dplyr::pull(value) %>%
  stringr::str_trim()

if (length(last_update) != 1) {
  stop("No 'Time of creation' row in ", SOFTWARE_INFO_FILE)
}

segment_lines <- purrr::map_chr(PRIMARY_PRODUCTS, function(product) {
  segment <- PRODUCT_SEGMENT[[product]]
  stringr::str_c("Segment_", segment, "\t", length(alignments[[product]]), "\t",
                 length(tree_tips[[as.character(segment)]]))
})

writeLines(c(stringr::str_c("last_update\t[", last_update, "]"),
             stringr::str_c("total_GenBank_sequences\t", total_genbank),
             stringr::str_c("total_curated_sequences\t", total_curated),
             "segment\ttotal\tclustered",
             segment_lines),
           file.path(OUTPUT_DIR, "IAV_DB_summary.log"))


# STEP 15 – Check outputs and publish APP_DIR ----------------------------------
outputs <- c(
  "sequence_metadata_full.rds",
  "product_positions.rds",
  "product_columns.rds",
  "sequence_metadata.rds",
  stringr::str_c("trees/seg", 1:8, "_cluster_rep.nwk"),
  "subtype_references.rds",
  "reference_numbering.rds",
  stringr::str_c("blastdb/sgt_", 1:8, ".pin"),
  "master_blastdb/reference_proteins_db.pin",
  "reference_proteins_db.fasta",
  "alignment_column_counts.rds",
  "cluster_column_residues.rds",
  "cluster_host_composition.rds",
  "IAV_DB_summary.log"
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

missing_outputs <- outputs[!file.exists(file.path(OUTPUT_DIR, outputs))]
if (length(missing_outputs) > 0) {
  stop("Missing output(s), ", APP_DIR, " left unchanged: ", paste(missing_outputs, collapse = ", "))
}

# current APP_DIR -> <APP_DIR>_prev, staging -> APP_DIR
if (dir.exists(PREV_DIR)) unlink(PREV_DIR, recursive = TRUE)
if (dir.exists(APP_DIR) && !file.rename(APP_DIR, PREV_DIR)) {
  stop("Could not move ", APP_DIR, " to ", PREV_DIR)
}
if (!file.rename(OUTPUT_DIR, APP_DIR)) {
  stop("Could not move ", OUTPUT_DIR, " to ", APP_DIR)
}
message("Outputs saved to: ", APP_DIR, if (dir.exists(PREV_DIR)) str_c(" | backup of previous version: ", PREV_DIR, ")") else "")

message(sprintf("[COMPLETE] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

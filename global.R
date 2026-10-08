library(shiny)
library(shiny.react)
library(shinyjs)
library(jsonlite)
library(tidyverse) # TODO individual packages
library(magrittr)
# library(gt)
library(rBLAST)
library(DT)
library(janitor)
library(ape)
library(Biostrings)
library(pwalign)
library(tools)
library(readxl)
library(ggiraph)
library(bslib)
library(gdtools)
library(config)
library(logger)

# Initialize logger
log_threshold(INFO)
log_layout(layout_glue_generator(format = "{time} [{level}] {msg}"))

# Load configuration
conf <- config::get()

# The tabs address a protein product, not a segment: segments 2, 3, 7 and 8 each encode a
# second one, and the catalogue numbers its M2 and NEP entries against those products
# rather than against M1 and NS1 - position 11 of M2 is a different residue from position
# 11 of M1.
#
# Clustering and the phylogenies are still per segment, so a secondary product shares its
# parent's tree; the Tree tab maps back to the parent for the tree and the metadata and
# uses the product itself for the position search.
#
# PB1-F2 and PA-X have no catalogue entries yet, so their Adaptation Mutations table is
# empty. They are listed all the same, because their position data exists and a catalogue
# row added later is then numbered against the right reading frame.
#
# Each entry names the parent segment, the label used in captions, and the full alignment
# the per-sequence counts are tallied from.
products <- list(
  seg1        = list(segment = 1, label = "PB2",    alignment = "sgt_1_PB2_AA.fasta"),
  seg2        = list(segment = 2, label = "PB1",    alignment = "sgt_2_PB1_AA.fasta"),
  seg2_PB1_F2 = list(segment = 2, label = "PB1-F2", alignment = "sgt_2_PB1_F2_AA.fasta"),
  seg3        = list(segment = 3, label = "PA",     alignment = "sgt_3_PA_AA.fasta"),
  seg3_PA_X   = list(segment = 3, label = "PA-X",   alignment = "sgt_3_PA_X_AA.fasta"),
  seg4        = list(segment = 4, label = "HA",     alignment = "sgt_4_HA_AA.fasta"),
  seg5        = list(segment = 5, label = "NP",     alignment = "sgt_5_NP_AA.fasta"),
  seg6        = list(segment = 6, label = "NA",     alignment = "sgt_6_NA_AA.fasta"),
  seg7        = list(segment = 7, label = "M1",     alignment = "sgt_7_M1_AA.fasta"),
  seg7_M2     = list(segment = 7, label = "M2",     alignment = "sgt_7_M2_AA.fasta"),
  seg8        = list(segment = 8, label = "NS1",    alignment = "sgt_8_NS1_AA.fasta"),
  seg8_NEP    = list(segment = 8, label = "NEP",    alignment = "sgt_8_NEP_AA.fasta")
)

# Segment number and display label for a product key, for captions and the HA checks
product_segment <- function(key) products[[key]]$segment
product_label   <- function(key) products[[key]]$label

# Most hits listed in the BLAST hits table. Defaulted here so a config predating
# the setting still starts rather than failing at the first search.
max_blast_hits <<- conf$blast$max_hits %||% 100L
log_info("BLAST hits table capped at {max_blast_hits} hits")

register_gfont("Open Sans")

# Everything below reads from R/, so it is sourced before any data is loaded
sapply(list.files("R", full.names = TRUE), source)

# The data the app reads, loaded by the functions in R/data_loading.R and copied into the
# global environment under the names the modules and the tests look for: tree_names,
# discarded, transposed, ref_seqs, numbering_refs and so on - see read_app_data() for the
# list. The full list is also kept whole as app_data.
app_data <- read_app_data(conf, products)
list2env(app_data, envir = globalenv())

# ref_set.tsv is no longer read at startup: blast_query_segment() searches the product
# references, where the subject id is the product and no ref_id lookup applies. It stays in
# BLAST_segment_recognizer/ as the labelled panel that replacement was measured against -
# see MIN_SEGMENT_BITSCORE in R/query_segment.R.

critical <- tryCatch(
  validate_data(),
  error = function(e) {
    log_error("Validation failed: {e$message}")
    e$message
  })

# A missing data file or an unnumberable segment produces silently wrong or empty
# results, so refuse to start rather than serve them
if(length(critical) > 0) {
  stop("Data validation failed:\n", str_c("  - ", critical, collapse = "\n"))
}

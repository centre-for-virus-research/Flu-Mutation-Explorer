#!/usr/bin/env Rscript

# flu_gdb_taxonomic_names.R ----------------------------------------------------
# v1.2 – 2026-10-02
# Standalone script to generate taxonomic_names.csv from GenBank matrix.
# Hosts are classified from the matrix tax_id; hosts left without
# class/order are logged to flu_gdb_inputs/unclassified_hosts.txt.
#
# Required inputs (relative to ROOT_DIR):
#   1) current_version/IAV_DB_matrix.tsv
#
# Output:
#   flu_gdb_inputs/taxonomic_names.csv
#   flu_gdb_inputs/unclassified_hosts.txt hosts left without class/order

# Setup ------------------------------------------------------------------------
options(warn = 1, stringsAsFactors = FALSE)

## Libraries -------------------------------------------------------------------
library(tidyverse)
library(janitor)
library(taxize)

## Paths -----------------------------------------------------------------------
ROOT_DIR <- getwd()

INPUT_DIR   <- file.path(ROOT_DIR, "flu_gdb_inputs") # input folder of flu-gdb-genbank.R
MATRIX_FILE <- file.path(ROOT_DIR, "current_version", "IAV_DB_matrix.tsv")

stopifnot(file.exists(MATRIX_FILE))

# Load GenBank matrix ----------------------------------------------------------
genbank_metadata <- readr::read_tsv(MATRIX_FILE, show_col_types = FALSE) %>%
  janitor::clean_names()

# Extract unique host names and their tax_id -----------------------------------
hosts <- genbank_metadata %>%
  distinct(host_validated, tax_id) %>%
  mutate(tax_id = as.character(as.integer(tax_id))) %>% # avoid "1e+05"
  arrange(host_validated)

multiple_ids <- hosts %>% dplyr::count(host_validated) %>% dplyr::filter(n > 1)
if (nrow(multiple_ids) > 0) {
  warning("Host(s) with more than one tax_id (first kept): ",
          paste(multiple_ids$host_validated, collapse = ", "), call. = FALSE)
  hosts <- hosts %>% distinct(host_validated, .keep_all = TRUE)
}

taxids <- hosts %>% dplyr::filter(!is.na(tax_id)) %>% distinct(tax_id) %>% pull(tax_id)

message("Querying NCBI classification for ", length(taxids), " tax_id(s) (",
        nrow(hosts), " unique host names)...")

# Get classification from tax_id -----------------------------------------------
options(taxize_api_sleep = 0.4)
cls <- taxize::classification(taxids, db = "ncbi", messages = FALSE)

# Extract class and order ------------------------------------------------------
ranks <- purrr::map2_dfr(cls, names(cls), function(x, id) {
  if (is.null(x) || inherits(x, "logical") || nrow(x) == 0) {
    return(tibble(tax_id = id, class = NA_character_, order = NA_character_))
  }
  tibble(
    tax_id = id,
    class  = x$name[x$rank == "class"][1],
    order  = x$name[x$rank == "order"][1]
  )
})

taxonomic_names <- hosts %>%
  left_join(ranks, by = "tax_id") %>%
  transmute(db = "ncbi", query = host_validated, class, order)

# Correct known mislabels ------------------------------------------------------
# matrix tax_id points to insects for "common gull" and "peacock", and to the virus
# for "unidentified influenza virus"
taxonomic_names <- taxonomic_names %>%
  mutate(
    order = case_when(
      query == "common gull"                  ~ "Charadriiformes",
      query == "peacock"                      ~ "Galliformes",
      query == "unidentified influenza virus"  ~ NA_character_,
      query == "environmental samples"         ~ "Environment",
      TRUE                                    ~ order
    ),
    class = case_when(
      query %in% c("common gull", "peacock") ~ "Aves",
      query == "unidentified influenza virus" ~ NA_character_,
      query == "environmental samples"        ~ "Environment",
      TRUE                                    ~ class
    )
  )

# Log hosts left without classification ----------------------------------------
unclassified <- taxonomic_names %>%
  dplyr::filter(is.na(class) & is.na(order)) %>%
  left_join(hosts, by = c("query" = "host_validated"))

if (nrow(unclassified) > 0) {
  log_file <- file.path(INPUT_DIR, "unclassified_hosts.txt")
  writeLines(
    c(
      paste("# Unclassified hosts –", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "# host_validated <tab> tax_id",
      "",
      paste(unclassified$query, unclassified$tax_id, sep = "\t")
    ),
    con = log_file
  )
  message("Unclassified hosts: ", nrow(unclassified), " - logged to: ", log_file)
}

# Missing class/order as "Unknown", like host_group (no NA state in the app)
taxonomic_names <- taxonomic_names %>%
  mutate(
    class = replace_na(class, "Unknown"),
    order = replace_na(order, "Unknown")
  )

# Assign host_group ------------------------------------------------------------
taxonomic_names <- taxonomic_names %>%
  mutate(
    host_group = case_when(
      query == "Homo sapiens"      ~ "Human",
      order %in% c(
        "Primates", "Carnivora", "Eulipotyphla", "Chiroptera",
        "Artiodactyla", "Perissodactyla", "Rodentia", "Pilosa", "Lagomorpha"
      )                            ~ "Other Mammals",
      class == "Aves"              ~ "Birds",
      query == "environmental samples" ~ "Environment",
      TRUE                         ~ "Unknown"
    )
  )


# Save -------------------------------------------------------------------------
out_file <- file.path(INPUT_DIR, "taxonomic_names.csv")
taxonomic_names %>% readr::write_csv(out_file)
message("Saved: ", out_file)

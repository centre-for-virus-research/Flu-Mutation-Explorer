#!/usr/bin/env Rscript

# flu_gdb_taxonomic_names.R ----------------------------------------------------
# v1.1 – 2026-03-11
# Standalone script to generate taxonomic_names.csv from GenBank matrix.
# Runs non-interactively — ambiguous NCBI matches are resolved by selecting
# the first result automatically and logged to flu_gdb_app_data/ambiguous_hosts.log.
#
# Required inputs (relative to ROOT_DIR):
#   1) current_version/IAV_DB_matrix_filtered.tsv
#
# Output:
#   flu_gdb_app_data/taxonomic_names.csv

# Setup ------------------------------------------------------------------------
options(warn = 1, stringsAsFactors = FALSE)

## Libraries -------------------------------------------------------------------
library(tidyverse)
library(janitor)
library(taxize)

## Paths -----------------------------------------------------------------------
ROOT_DIR    <- "/home/laura/IAV_DB/Flu-Mutation-Explorer/"
DATA_DIR    <- file.path(ROOT_DIR, "flu_gdb_app_data")
MATRIX_FILE <- file.path(ROOT_DIR, "current_version", "IAV_DB_matrix_filtered.tsv")

stopifnot(file.exists(MATRIX_FILE))

# Load GenBank matrix ----------------------------------------------------------
genbank_metadata <- readr::read_tsv(MATRIX_FILE, show_col_types = FALSE) %>%
  janitor::clean_names()

# Extract unique host names ----------------------------------------------------
hosts <- genbank_metadata %>%
  distinct(host_validated) %>%
  arrange(host_validated)

message("Querying NCBI for ", nrow(hosts), " unique host names...")

# Resolve UIDs non-interactively -----------------------------------------------
options(taxize_api_sleep = 0.4)
uids <- taxize::get_uid(hosts$host_validated, ask = FALSE, messages = FALSE)

# Log ambiguous cases (first result selected automatically) --------------------
ambiguous_idx <- attr(uids, "match") == "NA due to ask=FALSE & > 1 result"
if (any(ambiguous_idx)) {
  log_file <- file.path(DATA_DIR, "ambiguous_hosts.txt")
  writeLines(
    c(
      paste("# Ambiguous hosts –", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("# First NCBI result selected automatically for each"),
      "",
      hosts$host_validated[ambiguous_idx]
    ),
    con = log_file
  )
  message("Ambiguous hosts (first result selected): ",
    paste(hosts$host_validated[ambiguous_idx], collapse = ", "))
  message("Full list logged to: ", log_file)
}

# Manual UID overrides for known problematic queries ---------------------------
# (check ambiguous_hosts.txt) 

uid_overrides <- c(
  "Snow goose"              = "8849",    # Anser caerulescens — species over subspecies
  "environmental samples"   = "61964",
  "greater flamingo"        = "435638",
  "lesser black-backed gull"= "8915",    # Larus fuscus — species over subspecies
  "mink"                    = "452646",
  "pigs"                    = "9821"
)
override_idx <- hosts$host_validated %in% names(uid_overrides)
uids[override_idx] <- uid_overrides[hosts$host_validated[override_idx]]

# Get classification from UIDs -------------------------------------------------
cls <- taxize::classification(uids, db = "ncbi", messages = FALSE)

# Extract class and order ------------------------------------------------------
taxonomic_names <- purrr::map2_dfr(cls, hosts$host_validated, function(x, host) {
  if (is.null(x) || inherits(x, "logical") || nrow(x) == 0) {
    return(tibble(db = "ncbi", query = host, class = NA_character_, order = NA_character_))
  }
  tibble(
    db    = "ncbi",
    query = host,
    class = x$name[x$rank == "class"][1],
    order = x$name[x$rank == "order"][1]
  )
})


# Correct known mislabels ------------------------------------------------------
taxonomic_names <- taxonomic_names %>%
  mutate(
    order = case_when(
      query == "yellow-legged gull"           ~ "Charadriiformes",
      query == "common gull"                  ~ "Charadriiformes",
      query == "peacock"                      ~ "Galliformes",
      query == "unidentified influenza virus"  ~ NA_character_,
      TRUE                                    ~ order
    ),
    class = case_when(
      query %in% c("yellow-legged gull", "common gull", "peacock") ~ "Aves",
      query == "unidentified influenza virus"                       ~ NA_character_,
      TRUE                                                          ~ class
    )
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
out_file <- file.path(DATA_DIR, "taxonomic_names.csv")
taxonomic_names %>% readr::write_csv(out_file)
message("Saved: ", out_file)

#!/usr/bin/env Rscript
# Test script for BLAST segment recognizer
# Tests that each test sequence is correctly identified to its segment

library(Biostrings)
library(rBLAST)
library(dplyr)
library(stringr)
library(logger)

log_threshold(INFO)

# Expected segment mappings from filename
test_files <- list(
  "Blast_test_files/PB2_test.fa" = 1,
  "Blast_test_files/PB1_test.fa" = 2,
  "Blast_test_files/NP_test.fa" = 5,
  "Blast_test_files/NA_test.fa" = 6,
  "Blast_test_files/M1_test.fa" = 7,
  "Blast_test_files/NS1_test.fa" = 8
)

# Load reference set for segment mapping
ref_set <- read.table("BLAST_segment_recognizer/ref_set.tsv", 
                      col.names = c("accession_version", "segment"))

# Create BLAST database
blast_db_ref <- blast(db = "BLAST_segment_recognizer/ref_set_AA", type = "blastp")

log_info("Starting segment recognizer tests...")

all_passed <- TRUE
results <- list()

for (test_file in names(test_files)) {
  expected_segment <- test_files[[test_file]]
  
  if (!file.exists(test_file)) {
    log_warn("Test file not found: {test_file}")
    next
  }
  
  # Read sequence
  seq <- readAAStringSet(test_file)
  
  # Get sequence name and expected segment name
  seq_name <- names(seq)[1]
  segment_name <- str_extract(basename(test_file), "^[^_]+")
  
  log_info("Testing {segment_name} ({test_file})...")
  
  # Perform BLAST search
  hits <- predict(blast_db_ref, seq, BLAST_args = "-max_target_seqs 1") %>%
    arrange(evalue)
  
  if (nrow(hits) == 0) {
    log_error("  ✗ FAILED: No BLAST hits found for {segment_name}")
    all_passed <- FALSE
    results[[segment_name]] <- list(
      expected = expected_segment,
      detected = NA,
      status = "FAILED",
      reason = "No hits"
    )
    next
  }
  
  # Get top hit
  top_hit <- hits %>% slice(1)
  
  # Look up segment from reference set
  detected_segment <- ref_set %>%
    filter(accession_version == top_hit$sseqid) %>%
    pull(segment)
  
  if (length(detected_segment) == 0) {
    log_error("  ✗ FAILED: Top hit {top_hit$sseqid} not found in reference set")
    all_passed <- FALSE
    results[[segment_name]] <- list(
      expected = expected_segment,
      detected = NA,
      top_hit = top_hit$sseqid,
      status = "FAILED",
      reason = "Hit not in ref_set"
    )
    next
  }
  
  # Check if detected segment matches expected
  if (detected_segment == expected_segment) {
    log_success("  ✓ PASSED: {segment_name} correctly identified as segment {detected_segment}")
    log_info("    Top hit: {top_hit$sseqid} (E-value: {top_hit$evalue}, Identity: {top_hit$pident}%)")
    results[[segment_name]] <- list(
      expected = expected_segment,
      detected = detected_segment,
      top_hit = top_hit$sseqid,
      evalue = top_hit$evalue,
      identity = top_hit$pident,
      status = "PASSED"
    )
  } else {
    log_error("  ✗ FAILED: {segment_name} incorrectly identified as segment {detected_segment} (expected {expected_segment})")
    log_error("    Top hit: {top_hit$sseqid} (E-value: {top_hit$evalue})")
    all_passed <- FALSE
    results[[segment_name]] <- list(
      expected = expected_segment,
      detected = detected_segment,
      top_hit = top_hit$sseqid,
      evalue = top_hit$evalue,
      status = "FAILED",
      reason = "Wrong segment"
    )
  }
}

# Summary
log_info("=" %>% rep(60) %>% paste(collapse = ""))
log_info("Test Summary:")
log_info("=" %>% rep(60) %>% paste(collapse = ""))

passed <- sum(sapply(results, function(x) x$status == "PASSED"))
total <- length(results)

for (name in names(results)) {
  result <- results[[name]]
  status_icon <- if (result$status == "PASSED") "✓" else "✗"
  log_info("{status_icon} {name}: Expected seg{result$expected}, Got seg{result$detected %||% 'NA'}")
}

log_info("=" %>% rep(60) %>% paste(collapse = ""))
log_info("Results: {passed}/{total} tests passed")

if (all_passed) {
  log_success("ALL TESTS PASSED!")
  quit(status = 0)
} else {
  log_error("SOME TESTS FAILED")
  quit(status = 1)
}

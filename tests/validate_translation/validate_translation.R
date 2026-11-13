#!/usr/bin/env Rscript
# Validation of the alignment + translation steps of the Flu-ME pipeline.
# Compares the proteins in the protein alignments (proteins_AA) with the CDS translations
# annotated in GenBank (NCBI Entrez efetch, rettype = fasta_cds_aa).
#
#   Task 1: reference strains of ref_set_numbering, all 12 products
#           (skipped if the ref set file is not available).
#   Task 2: random sample of curated DB sequences (IAV_DB_matrix.tsv rows with empty "exclusion"):
#           1. sample accessions and download their GenBank CDS translations
#           2. compare them with the proteins in proteins_AA
#           3. check the origin of each X / mismatch in the codons of the CDS_nt alignment
#
# Settings (paths, Entrez, sampling mode and size) are in tests/config_validation.yml
# (another file can be given with the environment variable VALIDATION_CONFIG).
# Run from the project root:
#   Rscript tests/validate_translation.R

suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(httr)
  library(logger)
  library(config)
})

log_threshold(INFO)

# ---- Configuration ----
config_file <- Sys.getenv("VALIDATION_CONFIG", "tests/config_validation.yml")
conf <- config::get(file = config_file)

matrix_file <- conf$paths$matrix
aln_dir <- conf$paths$proteins_aa
nt_dir <- conf$paths$cds_nt
ref_set_file <- conf$paths$ref_set_numbering
cache_dir <- conf$paths$genbank_cache
product_map_file <- conf$paths$product_map

entrez_email <- conf$entrez$email
entrez_key <- conf$entrez$api_key
entrez_tool <- conf$entrez$tool
entrez_batch_size <- conf$entrez$batch_size

sampling_mode <- match.arg(conf$sampling$mode, c("total", "per_segment"))
sample_n <- if (sampling_mode == "total") conf$sampling$n_total else conf$sampling$n_per_segment
sample_seed <- conf$sampling$seed

max_signal_peptide <- conf$comparison$max_signal_peptide
min_identity_same_orf <- conf$comparison$min_identity_same_orf

# Results of each sampling setting go to their own folder
run_tag <- sprintf("%s_%d_seed%d", sampling_mode, sample_n, sample_seed)
out_dir <- file.path(conf$paths$results, run_tag)
task1_dir <- file.path(conf$paths$results, "task1_ref_set")

for (d in c(out_dir, task1_dir, cache_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
log_info("Sampling: {sampling_mode}, n = {sample_n}, seed = {sample_seed} -> {out_dir}")

products <- data.frame(
  product = c("PB2", "PB1", "PB1-F2", "PA", "PA-X", "HA", "NP", "NA", "M1", "M2", "NS1", "NEP"),
  segment = c(1, 2, 2, 3, 3, 4, 5, 6, 7, 7, 8, 8),
  aln_file = c(
    "sgt_1_PB2_AA.fasta", "sgt_2_PB1_AA.fasta", "sgt_2_PB1_F2_AA.fasta",
    "sgt_3_PA_AA.fasta", "sgt_3_PA_X_AA.fasta", "sgt_4_HA_AA.fasta",
    "sgt_5_NP_AA.fasta", "sgt_6_NA_AA.fasta", "sgt_7_M1_AA.fasta",
    "sgt_7_M2_AA.fasta", "sgt_8_NS1_AA.fasta", "sgt_8_NEP_AA.fasta"
  ),
  stringsAsFactors = FALSE
)
products$nt_file <- sub("_AA\\.fasta$", "_CDS.fasta", products$aln_file)

# ---- Helpers: GenBank ----

# Assign each GenBank CDS to a product with the mapping file (paths$product_map), built from the
# /gene + /product names of all CDS in the DB (validation/collect_cds_product_names.R).
# Lookup key: segment + gene + product + spliced (join in the location); names compared as they
# are, only trimming spaces. Returns, per CDS:
#   product    one of the 12 products when the map gives one, otherwise NA
#   map_status "mapped"              classified_as is one of the 12 products
#              "minor_product"       identified, but not compared (PB1-N40, M42, NS3, ...)
#              "unnamed_cds"         classified_as = unnamed_cds in the map (see classification_note)
#              "no_gene_no_product"  CDS without /gene and /product: discarded
#              "name_not_in_mapping" name combination absent from the map
read_product_map <- function(file) {
  read.delim(file, colClasses = "character", na.strings = character(), quote = "") %>%
    filter(segment != "unknown") %>%
    mutate(gene = str_squish(gene), product = str_squish(product)) %>%
    distinct(segment, gene, product, spliced, .keep_all = TRUE) %>%
    select(segment, gene, product, spliced, classified_as, classification_note)
}

classify_product <- function(gene, protein, location, segment, product_map) {
  key <- data.frame(
    segment = as.character(segment),
    gene = str_squish(ifelse(is.na(gene), "", gene)),
    product = str_squish(ifelse(is.na(protein), "", protein)),
    spliced = ifelse(!is.na(location) & grepl("join", location), "TRUE", "FALSE"),
    stringsAsFactors = FALSE
  ) %>%
    left_join(product_map, by = c("segment", "gene", "product", "spliced"))
  map_status <- case_when(
    key$gene == "" & key$product == "" ~ "no_gene_no_product",
    is.na(key$classified_as) ~ "name_not_in_mapping",
    key$classified_as == "unnamed_cds" ~ "unnamed_cds",
    key$classified_as %in% products$product ~ "mapped",
    TRUE ~ "minor_product"
  )
  data.frame(
    product = ifelse(map_status == "mapped", key$classified_as, NA_character_),
    map_status = map_status,
    classification_note = key$classification_note,
    stringsAsFactors = FALSE
  )
}

# Parse an NCBI fasta_cds_aa file into a data frame (one row per annotated CDS)
parse_cds_fasta <- function(file) {
  empty <- data.frame(
    accession_version = character(), gene = character(), protein_name = character(),
    protein_id = character(), location = character(), translation = character(),
    accession = character(), stringsAsFactors = FALSE
  )
  if (!file.exists(file) || file.size(file) < 2) return(empty)
  aa <- tryCatch(readAAStringSet(file), error = function(e) NULL)
  if (is.null(aa) || length(aa) == 0) return(empty)
  h <- names(aa)
  data.frame(
    accession_version = str_match(h, "^lcl\\|(.+?)_(prot|cds)_")[, 2],
    gene = str_match(h, "\\[gene=([^\\]]*)\\]")[, 2],
    protein_name = str_match(h, "\\[protein=([^\\]]*)\\]")[, 2],
    protein_id = str_match(h, "\\[protein_id=([^\\]]*)\\]")[, 2],
    location = str_match(h, "\\[location=([^\\]]*)\\]")[, 2],
    translation = as.character(aa),
    stringsAsFactors = FALSE
  ) %>%
    mutate(accession = sub("\\.\\d+$", "", accession_version))
}

# GenBank cache: downloaded batches (*.fasta) and the list of ids already requested (fetched_ids.txt)
fetched_ids_file <- file.path(cache_dir, "fetched_ids.txt")

read_fetched_ids <- function() {
  if (file.exists(fetched_ids_file)) unique(readLines(fetched_ids_file)) else character()
}

# Download CDS translations for the ids not fetched before; return the cached CDS of all ids
fetch_cds_aa <- function(ids) {
  ids <- unique(ids)
  done <- read_fetched_ids()
  to_fetch <- setdiff(ids, done)
  log_info("GenBank: {length(ids)} ids needed, {length(ids) - length(to_fetch)} already in cache, {length(to_fetch)} to download")

  if (length(to_fetch) > 0) {
    if (is.null(entrez_email) || entrez_email == "") {
      stop("Entrez email is not set (entrez$email in tests/config_validation.yml / ENTREZ_EMAIL)")
    }
    batches <- split(to_fetch, ceiling(seq_along(to_fetch) / entrez_batch_size))
    n_existing <- length(list.files(cache_dir, pattern = "^cds_aa_batch_.*\\.fasta$"))
    for (i in seq_along(batches)) {
      body <- list(
        db = "nuccore", id = paste(batches[[i]], collapse = ","),
        rettype = "fasta_cds_aa", retmode = "text",
        tool = entrez_tool, email = entrez_email
      )
      if (!is.null(entrez_key) && entrez_key != "") body$api_key <- entrez_key
      log_info("  Entrez efetch batch {i}/{length(batches)} ({length(batches[[i]])} ids)")
      resp <- NULL
      for (attempt in 1:4) {
        resp <- tryCatch(
          POST("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
               body = body, encode = "form", timeout(180)),
          error = function(e) NULL
        )
        if (!is.null(resp) && status_code(resp) == 200) break
        log_warn("    attempt {attempt} failed, retrying")
        Sys.sleep(2^attempt)
      }
      if (is.null(resp) || status_code(resp) != 200) stop("Entrez efetch failed for batch ", i)
      out_file <- file.path(cache_dir, sprintf("cds_aa_batch_%05d.fasta", n_existing + i))
      tmp <- paste0(out_file, ".tmp")
      writeLines(content(resp, as = "text", encoding = "UTF-8"), tmp)
      file.rename(tmp, out_file)
      write(batches[[i]], fetched_ids_file, append = TRUE) # only after the batch is safely written
      Sys.sleep(if (!is.null(entrez_key) && entrez_key != "") 0.15 else 0.4)
    }
  }

  # Keep the requested version (ids without version: whatever version was returned)
  has_version <- grepl("\\.\\d+$", ids)
  ids_noversion <- sub("\\.\\d+$", "", ids)
  cache_files <- list.files(cache_dir, pattern = "\\.fasta$", full.names = TRUE)
  gb <- bind_rows(lapply(cache_files, parse_cds_fasta)) %>%
    filter(accession_version %in% ids[has_version] | accession %in% ids[!has_version]) %>%
    distinct(accession_version, protein_id, location, translation, .keep_all = TRUE)
  missing <- setdiff(ids_noversion, gb$accession)
  if (length(missing) > 0) {
    log_warn("{length(missing)} of {length(ids)} records have no annotated CDS in GenBank")
  }
  gb
}

# ---- Helpers: comparison ----

# Compare one pipeline protein (gapped, from the alignment) with one GenBank translation.
# ref_gapped: gapped protein of the reference the sequence was aligned against (closest_reference),
#             used to check whether differences come from that reference.
# nt_gapped:  gapped CDS of the sequence (CDS_nt alignment; codon i = column i of the protein
#             alignment), used to tell why the pipeline has an X.
compare_protein <- function(pipe_gapped, gb, product, ref_gapped = NA_character_, nt_gapped = NA_character_) {
  p <- gsub("-", "", pipe_gapped, fixed = TRUE)
  stop_pos <- regexpr("*", p, fixed = TRUE)[1]
  has_stop <- stop_pos > 0
  p_core <- if (has_stop) substr(p, 1, stop_pos - 1) else p
  gb <- sub("\\*$", "", gb)
  n_term_allowed <- if (product == "HA") max_signal_peptide else 0

  # Exact identity (HA: our protein = GenBank translation without the signal peptide)
  identical_seq <- if (product == "HA") {
    nchar(p_core) > 0 && endsWith(gb, p_core) && (nchar(gb) - nchar(p_core)) <= n_term_allowed
  } else {
    p_core == gb
  }

  res <- data.frame(
    identical = identical_seq, len_pipeline = nchar(p_core), len_genbank = nchar(gb),
    pipeline_stop = has_stop, n_aa_compared = nchar(p_core), n_aa_identical = nchar(p_core),
    identity = 1, n_term_missing = 0L, n_term_extra = 0L,
    c_term_missing = 0L, c_term_extra = 0L,
    n_x_ambiguous_nt = 0L, n_x_split_codon = 0L, n_x_terminal = 0L, n_x_other = 0L, n_gb_ambiguous = 0L,
    n_subst = 0L, n_subst_eq_ref = 0L, n_ins = 0L, n_del = 0L, n_del_ref_insertion = 0L,
    flag_x = FALSE, flag_subst = FALSE, flag_indel = FALSE, flag_premature_stop = FALSE,
    flag_c_trunc = FALSE, flag_n_term = FALSE, flag_extension = FALSE, flag_orf_annotation = FALSE,
    details = "", stringsAsFactors = FALSE
  )
  if (identical_seq) {
    res$n_term_missing <- nchar(gb) - nchar(p_core)
    return(res)
  }
  if (nchar(p_core) == 0 || nchar(gb) == 0) {
    res$c_term_missing <- nchar(gb)
    res$n_aa_compared <- 0L
    res$n_aa_identical <- 0L
    res$flag_premature_stop <- has_stop
    res$flag_c_trunc <- !has_stop
    res$details <- "empty pipeline protein"
    return(res)
  }

  # Residues outside BLOSUM62 (J, U, O) are masked as X for the alignment only
  mask <- function(s) gsub("[^ACDEFGHIKLMNPQRSTVWYBZX]", "X", s)
  aln <- pairwiseAlignment(
    AAString(mask(p_core)), AAString(mask(gb)),
    type = "overlap", substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 1
  )
  ap <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  ag <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  p0 <- start(pattern(aln)); p1 <- end(pattern(aln))
  g0 <- start(subject(aln)); g1 <- end(subject(aln))
  p_pos <- p0 - 1 + cumsum(ap != "-")
  g_pos <- g0 - 1 + cumsum(ag != "-")
  pc <- strsplit(p_core, "")[[1]]
  gc <- strsplit(gb, "")[[1]]
  ap[ap != "-"] <- pc[p_pos[ap != "-"]] # restore original residues
  ag[ag != "-"] <- gc[g_pos[ag != "-"]]
  both <- ap != "-" & ag != "-"
  res$n_aa_compared <- sum(both)
  res$n_aa_identical <- sum(both & ap == ag)
  res$identity <- sum(both & ap == ag) / max(1, sum(both))

  res$n_term_missing <- g0 - 1
  res$n_term_extra <- p0 - 1
  res$c_term_missing <- nchar(gb) - g1
  res$c_term_extra <- nchar(p_core) - p1

  # GenBank annotates a different ORF (e.g. PB1-F2 from a downstream ATG): not comparable residue by residue.
  # Detected by low identity or by staggered ends (each protein overhangs the other at opposite ends).
  staggered <- (res$n_term_missing > n_term_allowed + 5 && res$c_term_extra > 5) ||
    (res$n_term_extra > 5 && res$c_term_missing > 5)
  if (res$identity < min_identity_same_orf || staggered) {
    res$flag_orf_annotation <- TRUE
    res$details <- sprintf("different_ORF: pipeline %d aa%s, GenBank %d aa, identity %.0f%%",
                           nchar(p_core), if (has_stop) " (stop)" else "", nchar(gb), 100 * res$identity)
    return(res)
  }

  mm <- both & ap != ag
  ins <- ag == "-" # residue only in pipeline
  del <- ap == "-" # residue only in GenBank
  x_pipe <- (mm | ins) & ap == "X" # X produced by the pipeline
  gb_amb <- mm & ap != "X" & ag %in% c("X", "B", "Z", "J") # ambiguous in GenBank, resolved in pipeline
  subst <- mm & !x_pipe & !gb_amb
  ins_real <- ins & ap != "X"
  extra_pos <- c(seq_len(p0 - 1), if (p1 < nchar(p_core)) (p1 + 1):nchar(p_core))
  extra_x <- extra_pos[pc[extra_pos] == "X"]
  extra_other <- extra_pos[pc[extra_pos] != "X"]

  # Why does the pipeline have an X? Look at the codon in the CDS alignment
  cols <- which(strsplit(pipe_gapped, "")[[1]] != "-")
  x_idx <- p_pos[x_pipe]
  if (length(x_idx) > 0) {
    if (!is.na(nt_gapped) && nchar(nt_gapped) == 3 * nchar(pipe_gapped)) {
      cod <- substring(nt_gapped, 3 * cols[x_idx] - 2, 3 * cols[x_idx])
      xc <- ifelse(grepl("-", cod), "split_codon", ifelse(grepl("[^ACGTacgt]", cod), "ambiguous_nt", "other"))
    } else {
      # CDS alignment not codon-matched to the protein alignment (PB1-F2, PA-X, M1, NS1):
      # infer from the protein pattern. Ambiguity code in GenBank (B/Z/J/X) = ambiguous nucleotide;
      # pipeline X next to an extra pipeline X = codon split by a gap.
      ins_x_pos <- p_pos[ins & ap == "X"]
      gb_at_x <- ag[x_pipe]
      split_like <- x_idx %in% ins_x_pos | (x_idx + 1) %in% ins_x_pos | (x_idx - 1) %in% ins_x_pos
      xc <- ifelse(gb_at_x %in% c("B", "Z", "J", "X"), "ambiguous_nt",
                   ifelse(split_like, "split_codon", "other"))
    }
    res$n_x_ambiguous_nt <- sum(xc == "ambiguous_nt")
    res$n_x_split_codon <- sum(xc == "split_codon")
    res$n_x_other <- sum(xc == "other")
  }
  res$n_x_terminal <- length(extra_x)
  res$n_gb_ambiguous <- sum(gb_amb)
  res$n_subst <- sum(subst)
  res$n_ins <- sum(ins_real)
  res$n_del <- sum(del)

  run_bounds <- function(flag) {
    r <- rle(flag)
    ends <- cumsum(r$lengths)
    data.frame(start = ends - r$lengths + 1, end = ends, len = r$lengths)[r$values, , drop = FALSE]
  }

  # Differences that come from the reference used for the alignment
  if (!is.na(ref_gapped)) {
    ref_chars <- strsplit(ref_gapped, "")[[1]]
    if (any(subst)) {
      rc <- ref_chars[cols[p_pos[subst]]]
      res$n_subst_eq_ref <- sum(!is.na(rc) & rc == ap[subst])
    }
    # Residues missing in the pipeline where the reference has no residue either
    # (insertion relative to the reference, without columns in the alignment)
    if (any(del)) {
      rb <- run_bounds(del)
      for (k in seq_len(nrow(rb))) {
        pb <- p_pos[rb$start[k]] # last pipeline residue before the missing run
        if (pb >= 1 && pb < nchar(p_core)) {
          between <- seq_len(cols[pb + 1] - cols[pb] - 1) + cols[pb]
          if (all(ref_chars[between] == "-")) res$n_del_ref_insertion <- res$n_del_ref_insertion + rb$len[k]
        }
      }
    }
  }

  res$flag_x <- (res$n_x_ambiguous_nt + res$n_x_split_codon + res$n_x_other + res$n_x_terminal + res$n_gb_ambiguous) > 0
  res$flag_subst <- res$n_subst > 0
  res$flag_indel <- (res$n_ins + res$n_del) > 0
  res$flag_premature_stop <- has_stop && res$c_term_missing > 0
  res$flag_c_trunc <- !has_stop && res$c_term_missing > 0
  res$flag_n_term <- res$n_term_missing > n_term_allowed
  res$flag_extension <- length(extra_other) > 0

  # Human readable details: GenBank residue + GenBank position + pipeline residue (pipeline position)
  det <- character()
  if (any(mm)) det <- c(det, paste0(ag[mm], g_pos[mm], ap[mm], "(p", p_pos[mm], ")"))
  rb <- run_bounds(del)
  if (nrow(rb) > 0) det <- c(det, sprintf("missing_in_pipeline_gb:%d-%d", g_pos[rb$start], g_pos[rb$end]))
  rb <- run_bounds(ins)
  if (nrow(rb) > 0) det <- c(det, sprintf("extra_in_pipeline_p:%d-%d", p_pos[rb$start], p_pos[rb$end]))
  if (res$flag_n_term) det <- c(det, sprintf("Nterm_missing:%d", res$n_term_missing))
  if (res$n_term_extra > 0) det <- c(det, sprintf("Nterm_extra:%s", substr(p_core, 1, p0 - 1)))
  if (res$flag_premature_stop) {
    det <- c(det, sprintf("premature_stop:p%d(gb_len=%d)", nchar(p_core) + 1, nchar(gb)))
  }
  if (res$flag_c_trunc) det <- c(det, sprintf("Cterm_missing:%d", res$c_term_missing))
  if (res$c_term_extra > 0) det <- c(det, sprintf("Cterm_extra:%s", substr(p_core, p1 + 1, nchar(p_core))))
  res$details <- paste(det, collapse = ";")
  res
}

# Compare against all GenBank candidates for one product and keep the best one
compare_best <- function(pipe_gapped, gb_candidates, product, ref_gapped = NA_character_, nt_gapped = NA_character_) {
  cmp <- lapply(gb_candidates$translation, function(g) compare_protein(pipe_gapped, g, product, ref_gapped, nt_gapped))
  score <- sapply(cmp, function(r) {
    (!r$identical) * 1e6 + r$flag_orf_annotation * 1e5 + r$n_x_ambiguous_nt + r$n_x_split_codon +
      r$n_x_terminal + r$n_x_other + r$n_gb_ambiguous + r$n_subst + r$n_ins + r$n_del +
      r$c_term_missing + r$c_term_extra + r$n_term_extra + r$flag_n_term * r$n_term_missing
  })
  best <- which.min(score)
  cbind(
    cmp[[best]],
    protein_id = gb_candidates$protein_id[best],
    genbank_accession_version = gb_candidates$accession_version[best],
    n_genbank_candidates = nrow(gb_candidates),
    translation = gb_candidates$translation[best],
    stringsAsFactors = FALSE
  )
}

read_alignment_subset <- function(file, accs, reader = readAAStringSet) {
  accs <- unique(accs[!is.na(accs)])
  if (length(accs) == 0) return(character())
  aln <- reader(file)
  nm <- sub("\\s.*$", "", names(aln))
  idx <- match(accs, nm)
  keep <- !is.na(idx)
  setNames(as.character(aln[idx[keep]]), accs[keep])
}

# Width of the first sequence of a FASTA alignment (reads only the first record)
first_record_width <- function(file) {
  con <- file(file, "r")
  on.exit(close(con))
  width <- 0
  seen_header <- FALSE
  repeat {
    l <- readLines(con, n = 500)
    if (length(l) == 0) break
    for (x in l) {
      if (startsWith(x, ">")) {
        if (seen_header) return(width)
        seen_header <- TRUE
      } else {
        width <- width + nchar(x)
      }
    }
  }
  width
}

# Clopper-Pearson exact binomial proportion
prop_ci <- function(k, n) {
  if (n == 0) return(c(est = NA, lo = NA, hi = NA))
  ci <- binom.test(k, n)$conf.int
  c(est = 100 * k / n, lo = 100 * ci[1], hi = 100 * ci[2])
}

write_tsv_out <- function(x, file) write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Step 0: inputs ----
log_info("Reading {matrix_file}")
mat <- fread(
  matrix_file, sep = "\t", quote = "",
  select = c("primary_accession", "accession_version", "segment_validated", "closest_reference", "exclusion"),
  colClasses = "character", na.strings = NULL
) %>%
  as.data.frame() %>%
  mutate(segment = as.integer(as.numeric(segment_validated)))

curated <- mat %>% filter((is.na(exclusion) | exclusion == "") & !is.na(segment))
log_info("  {nrow(mat)} rows; {nrow(curated)} curated (no exclusion, segment assigned)")
pop_sizes <- curated %>% count(segment, name = "N_population")

run_task1 <- file.exists(ref_set_file)
if (run_task1) {
  log_info("Reading reference set: {ref_set_file}")
  ref_set <- read.delim(ref_set_file, colClasses = "character", na.strings = character()) %>%
    mutate(segment = as.integer(segment)) %>%
    left_join(mat %>% select(primary_accession, accession_version) %>% distinct(primary_accession, .keep_all = TRUE),
              by = c("accession" = "primary_accession")) %>%
    # fetch the same version used in the DB when available
    mutate(fetch_id = ifelse(is.na(accession_version) | accession_version == "", accession, accession_version))
  log_info("  {n_distinct(ref_set$strain)} reference strains, {nrow(ref_set)} segment records")
} else {
  log_warn("Reference set not found ({ref_set_file}): Task 1 is skipped")
  ref_set <- data.frame(serotype = character(), strain = character(), segment = integer(),
                        accession = character(), fetch_id = character(), stringsAsFactors = FALSE)
}

# ---- Step 1: random sample of accessions ----
sample_file <- file.path(out_dir, "sample.tsv")
if (file.exists(sample_file)) {
  log_info("Re-using existing sample: {sample_file}")
  smp <- read.delim(sample_file, colClasses = "character", na.strings = character()) %>%
    mutate(segment = as.integer(segment))
} else {
  set.seed(sample_seed)
  smp <- if (sampling_mode == "total") {
    curated %>% slice_sample(n = min(sample_n, nrow(curated)))
  } else {
    curated %>% group_by(segment) %>% slice_sample(n = sample_n) %>% ungroup()
  }
  smp <- smp %>%
    select(primary_accession, accession_version, segment, closest_reference) %>%
    arrange(segment, primary_accession)
  write_tsv_out(smp, sample_file)
  log_info("Sample drawn: {nrow(smp)} sequences -> {sample_file}")
}
log_info("Sample per segment: {paste(names(table(smp$segment)), table(smp$segment), sep = '=', collapse = ', ')}")

# ---- Step 1: GenBank CDS translations ----
gb <- fetch_cds_aa(c(ref_set$fetch_id, smp$accession_version))
seg_by_acc <- bind_rows(
  ref_set %>% select(accession, segment),
  smp %>% select(accession = primary_accession, segment)
) %>% distinct(accession, .keep_all = TRUE)
log_info("Reading product map: {product_map_file}")
product_map <- read_product_map(product_map_file)
gb <- gb %>% left_join(seg_by_acc, by = "accession")
gb <- bind_cols(gb, classify_product(gb$gene, gb$protein_name, gb$location, gb$segment, product_map))

not_mapped <- gb %>%
  filter(map_status != "mapped") %>%
  select(accession_version, segment, gene, protein_name, protein_id, location, map_status, classification_note)
write_tsv_out(not_mapped, file.path(out_dir, "genbank_cds_not_mapped.tsv"))
for (s in setdiff(unique(gb$map_status), "mapped")) {
  log_info("GenBank CDS {s} (not compared): {sum(gb$map_status == s)}")
}

# Candidates indexed by accession + product (fast lookup in the loops)
gb_mapped <- gb %>% filter(map_status == "mapped")
gb_index <- split(gb_mapped, paste(gb_mapped$accession, gb_mapped$product))

# Why a record has no CDS for a product: no CDS at all, or CDS not mapped / not annotated
status_without_cds <- function(acc) {
  st <- gb$map_status[gb$accession == acc]
  if (length(st) == 0) return("no_genbank_cds")
  for (s in c("name_not_in_mapping", "unnamed_cds", "no_gene_no_product")) if (s %in% st) return(s)
  "product_not_annotated"
}

# ---- Steps 2-3: compare with proteins_AA, check the codons in CDS_nt ----
task1_rows <- list()
task2_rows <- list()

for (i in seq_len(nrow(products))) {
  prod <- products$product[i]
  seg <- products$segment[i]
  aln_path <- file.path(aln_dir, products$aln_file[i])
  nt_path <- file.path(nt_dir, products$nt_file[i])

  refs_seg <- ref_set %>% filter(segment == seg)
  smp_seg <- smp %>% filter(segment == seg)
  if (nrow(refs_seg) + nrow(smp_seg) == 0) next

  log_info("Product {prod}: reading {aln_path}")
  aln <- read_alignment_subset(aln_path, c(refs_seg$accession, smp_seg$primary_accession, smp_seg$closest_reference))

  # The CDS alignment is only useful when codon i = column i of the protein alignment
  nt <- character()
  if (nrow(smp_seg) > 0 && file.exists(nt_path)) {
    if (first_record_width(nt_path) == 3 * first_record_width(aln_path)) {
      log_info("  reading {nt_path}")
      nt <- read_alignment_subset(nt_path, smp_seg$primary_accession, reader = readDNAStringSet)
    } else {
      log_info("  {basename(nt_path)} is not codon-matched to the protein alignment: X origin inferred from the protein")
    }
  }

  # Task 1
  for (j in seq_len(nrow(refs_seg))) {
    acc <- refs_seg$accession[j]
    cand <- gb_index[[paste(acc, prod)]]
    base <- data.frame(serotype = refs_seg$serotype[j], strain = refs_seg$strain[j], segment = seg,
                       accession = acc, product = prod, stringsAsFactors = FALSE)
    task1_rows[[length(task1_rows) + 1]] <- if (!acc %in% names(aln)) {
      cbind(base, status = "not_in_alignment")
    } else if (is.null(cand)) {
      cbind(base, status = status_without_cds(acc))
    } else {
      cbind(base, status = "compared", compare_best(aln[[acc]], cand, prod))
    }
  }

  # Task 2
  for (j in seq_len(nrow(smp_seg))) {
    acc <- smp_seg$primary_accession[j]
    ref_acc <- smp_seg$closest_reference[j]
    cand <- gb_index[[paste(acc, prod)]]
    base <- data.frame(accession = acc, accession_version = smp_seg$accession_version[j], segment = seg,
                       product = prod, closest_reference = ref_acc, stringsAsFactors = FALSE)
    task2_rows[[length(task2_rows) + 1]] <- if (!acc %in% names(aln)) {
      cbind(base, status = "not_in_alignment")
    } else if (is.null(cand)) {
      cbind(base, status = status_without_cds(acc))
    } else {
      ref_gapped <- if (!is.na(ref_acc) && ref_acc %in% names(aln)) aln[[ref_acc]] else NA_character_
      nt_gapped <- if (acc %in% names(nt)) nt[[acc]] else NA_character_
      cbind(base, status = "compared", compare_best(aln[[acc]], cand, prod, ref_gapped, nt_gapped))
    }
  }
  rm(aln, nt)
  gc(verbose = FALSE)
}

task1 <- bind_rows(task1_rows)
task2 <- bind_rows(task2_rows)

# ---- Task 1 outputs ----
if (run_task1) {
  write_tsv_out(task1 %>% select(-any_of("translation")), file.path(task1_dir, "comparison.tsv"))

  task1_summary <- task1 %>%
    group_by(product) %>%
    summarise(
      n_refs = n(),
      n_compared = sum(status == "compared"),
      n_identical = sum(status == "compared" & identical %in% TRUE),
      n_differ = sum(status == "compared" & identical %in% FALSE),
      n_no_genbank_cds = sum(status == "no_genbank_cds"),
      n_product_not_annotated = sum(status == "product_not_annotated"),
      n_cds_not_mapped = sum(status %in% c("name_not_in_mapping", "unnamed_cds", "no_gene_no_product")),
      n_not_in_alignment = sum(status == "not_in_alignment"),
      n_aa_genbank = sum(len_genbank[status == "compared"]),
      n_aa_pipeline = sum(len_pipeline[status == "compared"]),
      n_aa_compared = sum(n_aa_compared[status == "compared"]),
      n_aa_identical = sum(n_aa_identical[status == "compared"]),
      pct_aa_identical = 100 * n_aa_identical / n_aa_compared,
      differences = paste(
        paste0(strain, " (", accession, "): ", details)[status == "compared" & identical %in% FALSE],
        collapse = " | "
      ),
      .groups = "drop"
    ) %>%
    arrange(match(product, products$product))
  write_tsv_out(task1_summary, file.path(task1_dir, "summary_by_product.tsv"))
}

# ---- Task 2 outputs ----
write_tsv_out(task2 %>% select(-any_of("translation")), file.path(out_dir, "comparison_by_product.tsv"))
write_tsv_out(task2 %>% count(segment, product, status), file.path(out_dir, "status_counts.tsv"))

task2_seq <- task2 %>%
  filter(status == "compared") %>%
  group_by(accession, accession_version, segment, closest_reference) %>%
  summarise(
    n_products = n(),
    n_aa_compared = sum(n_aa_compared),
    n_aa_identical = sum(n_aa_identical),
    any_diff = any(!identical),
    any_diff_excl_orf_annotation = any(!identical & !flag_orf_annotation),
    cause_residue = any(flag_x | flag_subst | flag_indel),
    cause_x = any(flag_x),
    cause_x_ambiguous_nt = any(n_x_ambiguous_nt + n_gb_ambiguous > 0),
    cause_x_split_codon = any(n_x_split_codon > 0),
    cause_x_terminal_codon = any(n_x_terminal > 0),
    cause_x_other = any(n_x_other > 0),
    cause_subst = any(flag_subst),
    cause_indel = any(flag_indel),
    cause_premature_stop = any(flag_premature_stop),
    cause_orf_annotation = any(flag_orf_annotation),
    cause_other = any(flag_c_trunc | flag_n_term | flag_extension),
    n_subst = sum(n_subst),
    n_subst_eq_ref = sum(n_subst_eq_ref),
    n_del = sum(n_del),
    n_del_ref_insertion = sum(n_del_ref_insertion),
    products_differing = paste(product[!identical], collapse = ","),
    details = paste(paste0(product[!identical], ":", details[!identical]), collapse = " | "),
    .groups = "drop"
  ) %>%
  # All substitutions / missing residues reproduce the alignment reference
  mutate(explained_by_reference = (n_subst + n_del) > 0 & (n_subst_eq_ref + n_del_ref_insertion) == (n_subst + n_del))
write_tsv_out(task2_seq, file.path(out_dir, "comparison_by_sequence.tsv"))

metrics <- c("any_diff", "any_diff_excl_orf_annotation", "cause_residue", "cause_x",
             "cause_x_ambiguous_nt", "cause_x_split_codon", "cause_x_terminal_codon", "cause_x_other",
             "cause_subst", "cause_indel", "cause_premature_stop", "cause_orf_annotation",
             "cause_other", "explained_by_reference")
summ_rows <- list()
for (m in metrics) {
  for (s in sort(unique(task2_seq$segment))) {
    x <- task2_seq[[m]][task2_seq$segment == s]
    ci <- prop_ci(sum(x), length(x))
    summ_rows[[length(summ_rows) + 1]] <- data.frame(group = paste0("segment_", s), metric = m,
      n = length(x), k = sum(x), pct = ci["est"], ci95_low = ci["lo"], ci95_high = ci["hi"])
  }
  x <- task2_seq[[m]]
  ci <- prop_ci(sum(x), length(x))
  summ_rows[[length(summ_rows) + 1]] <- data.frame(group = "all_pooled", metric = m,
    n = length(x), k = sum(x), pct = ci["est"], ci95_low = ci["lo"], ci95_high = ci["hi"])
  if (sampling_mode == "per_segment") {
    # Stratified sample: estimate for the whole DB weighting each segment by its size (normal approx. CI)
    st <- task2_seq %>% group_by(segment) %>% summarise(n_h = n(), p_h = mean(.data[[m]]), .groups = "drop") %>%
      inner_join(pop_sizes, by = "segment") %>% mutate(w = N_population / sum(N_population))
    p_w <- sum(st$w * st$p_h)
    se_w <- sqrt(sum(st$w^2 * st$p_h * (1 - st$p_h) / st$n_h))
    summ_rows[[length(summ_rows) + 1]] <- data.frame(group = "all_weighted_by_segment_size", metric = m,
      n = length(x), k = NA, pct = 100 * p_w, ci95_low = 100 * max(0, p_w - 1.96 * se_w),
      ci95_high = 100 * min(1, p_w + 1.96 * se_w))
  }
}
task2_summary <- bind_rows(summ_rows)
rownames(task2_summary) <- NULL
write_tsv_out(task2_summary, file.path(out_dir, "summary.tsv"))

task2_by_ref <- task2_seq %>%
  group_by(segment, closest_reference) %>%
  summarise(
    n = n(), n_diff = sum(any_diff), n_residue = sum(cause_residue),
    n_premature_stop = sum(cause_premature_stop), n_explained_by_reference = sum(explained_by_reference),
    n_subst = sum(n_subst), n_subst_eq_ref = sum(n_subst_eq_ref),
    n_del = sum(n_del), n_del_ref_insertion = sum(n_del_ref_insertion),
    .groups = "drop"
  ) %>%
  arrange(desc(n_diff))
write_tsv_out(task2_by_ref, file.path(out_dir, "by_reference.tsv"))

writeLines(c(paste("Run:", format(Sys.time())), paste("Sampling:", run_tag),
             paste("Config:", config_file), capture.output(sessionInfo())),
           file.path(out_dir, "run_info.txt"))

# ---- Summary ----
log_info("=" %>% rep(60) %>% paste(collapse = ""))
if (run_task1) {
  log_info("Task 1 - reference set")
  for (k in seq_len(nrow(task1_summary))) {
    r <- task1_summary[k, ]
    log_info("  {r$product}: {r$n_identical}/{r$n_compared} identical, {r$n_aa_identical}/{r$n_aa_compared} aa (no GenBank CDS: {r$n_no_genbank_cds}; product not annotated: {r$n_product_not_annotated}; CDS not mapped: {r$n_cds_not_mapped}; not in alignment: {r$n_not_in_alignment})")
  }
}
log_info("Task 2 - random sample ({run_tag})")
log_info("  residues: {sum(task2_seq$n_aa_identical)}/{sum(task2_seq$n_aa_compared)} identical")
for (m in metrics) {
  r <- task2_summary %>% filter(group == "all_pooled", metric == m)
  log_info("  {m}: {r$k}/{r$n} = {sprintf('%.2f', r$pct)}% (95% CI {sprintf('%.2f', r$ci95_low)}-{sprintf('%.2f', r$ci95_high)})")
}
log_info("Results in {out_dir}")

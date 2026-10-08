#' What each file in R/ is expected to define
#'
#' global.R sources R/ by directory listing, so a file missing from a deployment is not an
#' error - it is a function that never gets defined. The app starts cleanly and dies on the
#' first click that reaches it, which the user sees as "Disconnected from server" with
#' nothing on screen to say why. It has happened: R/cluster_residues.R was once left out,
#' and the Position search killed the session as soon as it called
#' cluster_residue_summary().
#'
#' Keyed by file rather than listed flat, so the report can name the file to copy.
#' validate_data.R is absent from the list on purpose: without it there is nothing here to
#' run.
#' @export
REQUIRED_MODULE_OBJECTS <- list(
  "R/aa_palette.R"          = c("get_aa_palette_js", "get_aa_levels"),
  "R/adaptation.R"          = c("adaptation_product_hits", "adaptation_in_play",
                                "adaptation_query_residues", "adaptation_column_frequencies",
                                "adaptation_coverage", "adaptation_coverage_caption",
                                "adaptation_pooled_hosts", "adaptation_plot_data",
                                "adaptation_pooled_caption", "adaptation_plot"),
  "R/alignment_positions.R" = "build_product",
  "R/cluster_residues.R"    = c("cluster_residue_summary", "residue_share"),
  "R/components.R"          = "TaxoniumComponent",
  "R/data_loading.R"        = c("read_app_data", "load_trees", "load_full_metadata", "load_positions",
                                "load_alignment_counts", "load_cluster_residues",
                                "load_adaptation_mutations", "build_metadata_csv",
                                "load_subtype_references", "load_reference_numbering",
                                "load_reference_sequences", "load_status"),
  "R/host_categories.R"     = "EXCLUDED_HOSTS", # a character vector, not a function
  "R/mod_adaptation.R"      = c("adaptation_sidebar_ui", "adaptation_tab_ui", "adaptation_server"),
  "R/mod_screening.R"       = c("screening_sidebar_ui", "screening_tab_ui", "screening_server"),
  "R/mod_tree.R"            = c("tree_sidebar_ui", "tree_tab_ui", "tree_server"),
  "R/screening.R"           = c("screen_sequences", "screen_blast", "screen_in_play",
                                "screen_walk", "screen_catalogue", "screen_table",
                                "screen_product_label", "hsp_query_residues", "screen_download"),
  "R/numbering_reference.R" = c("canonical_protein", "reference_accession", "product_reference",
                                "resolve_numbering_reference", "resolve_numbering_references"),
  "R/query_segment.R"       = c("blast_query_segment", "resolve_query_segment",
                                "parse_sequence_records", "is_nucleotide",
                                "PROTEIN_ONLY_LETTERS", "MIN_SEGMENT_BITSCORE"),
  "R/utils.R"               = c("consensus", "gapless_to_consensus", "html_attribute",
                                "subtype_references", "papers_dois", "clean_query_sequence",
                                "aligned_query_residues", "record_product", "record_scheme",
                                "convert_reference_position")
)

#' Validate Data Integrity on Startup
#'
#' Checks that every file in R/ was sourced, that the required data files exist, and
#' that every segment has a usable reference for numbering adaptation mutation
#' positions.
#'
#' Returns the critical problems rather than stopping, so the caller can report
#' all of them at once. An empty result means the app is usable.
#'
#' All three inputs are arguments defaulting to their globals, so a test can validate a
#' deliberately broken resolution without having to break the real data first.
#'
#' @param config the configuration naming the required data files
#' @param resolution the numbering references as resolved by resolve_numbering_references()
#' @param required file-to-objects map the sourcing of R/ is expected to have produced
#' @return character vector of critical problems (empty if none)
#' @export
validate_data <- function(config = conf, resolution = numbering_resolution,
                          required = REQUIRED_MODULE_OBJECTS) {
  log_info("Starting data validation...")
  critical <- character()

  # 1. Every file in R/ was sourced
  # inherits = FALSE is what makes this a real check: source() puts these in the global
  # environment, and a search that inherited would let an attached package supply the name
  # and pass - ape exports consensus(), so a missing R/utils.R would report as fine.
  defined <- function(x) exists(x, envir = globalenv(), inherits = FALSE)

  absent <- map(required, ~ .x[!map_lgl(.x, defined)])
  absent <- absent[lengths(absent) > 0]

  if(length(absent) > 0) {
    critical <- c(critical, imap_chr(absent, function(objects, file) {
      if(length(objects) == length(required[[file]])) {
        str_c(file, " was not sourced - the file is missing from this deployment ",
              "(defines ", str_flatten_comma(objects), ")")
      } else {
        # the file is there but does not define what it used to: a rename that got
        # halfway, which breaks the callers just as thoroughly
        str_c(file, " no longer defines ", str_flatten_comma(objects))
      }
    }))
  }

  # 2. File Existence Checks
  #
  # The product reference database is named in config by its prefix, as BLAST wants it, so
  # the index is what gets tested. All three tabs identify a query against it, so a missing
  # one is three broken tabs and is worth catching here rather than at the first submission.
  #
  # NULL when a caller hands in a partial config, and c() drops that, as it does for every
  # other key here - str_c() would turn it into the bare ".pin" instead.
  reference_db <- config$paths$data$reference_protein_db
  reference_index <- if(!is.null(reference_db)) str_c(reference_db, ".pin")

  required_files <- c(
    config$paths$data$discarded,
    config$paths$data$transposed,
    config$paths$data$metadata,
    config$paths$data$subtype_references,
    config$paths$data$cluster_composition,
    config$paths$data$status,
    reference_index
  )

  missing_files <- required_files[!file.exists(required_files)]
  if(length(missing_files) > 0) {
    # The reference database is the one of these with a one-line fix, so name it
    hint <- if(!is.null(reference_index) && reference_index %in% missing_files) {
      " - the reference protein database is expected at the reference_protein_db path in config.yml"
    } else ""

    critical <- c(critical, str_c("Missing core data files: ",
                                  str_flatten_comma(missing_files), hint))
  }

  # 3. Numbering reference integrity
  # numbering_refs is resolved in global.R against the canonical protein for each
  # segment; NA means no representative could number that segment, which would leave
  # the Adaptation Mutations module returning nothing for it.
  unresolved <- names(resolution)[map_lgl(resolution, ~
    is.null(.x$accession) || is.na(.x$accession) || !is.null(.x$error))]

  if(length(unresolved) > 0) {
    critical <- c(critical, str_c("No usable numbering reference for: ", str_flatten_comma(unresolved)))
  }

  # the accessions, taken from the resolution passed in rather than from the global
  # numbering_refs, which is that same map() and would ignore the argument
  accessions <- map(resolution, "accession")
  resolved <- accessions[setdiff(names(accessions), unresolved)]
  log_info("Numbering references resolved for {length(resolved)} of {length(accessions)} segments/schemes.")

  # 4. Report which reference each segment ended up using, so the log records the
  # provenance of every position reported by the app for this database version
  iwalk(resolved, ~ log_info("  {.y}: {.x}"))

  if(length(critical) == 0) {
    log_info("Data validation complete.")
  } else {
    walk(critical, ~ log_error("CRITICAL: {.x}"))
  }

  critical
}

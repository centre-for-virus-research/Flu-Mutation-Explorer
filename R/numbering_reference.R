#' Resolve the reference sequences used to number adaptation mutation positions
#'
#' Catalogue positions are indices into the canonical full-length protein, so the
#' reference used to convert a position into an alignment column must number
#' identically to that protein. A reference carrying an indel shifts every position
#' after it, which would give confidently wrong answers rather than an error.
#'
#' One strain numbers the whole workflow: A/New_York/392/2004, the H3N2 RefSeq set,
#' present for all eight segments and all four secondary reading frames, with no residue
#' dropped from its HA. A poultry-adapted strain will not do - the stalk deletion makes
#' A/turkey/England/50-92/1991's NA 446 aa against a canonical 470, while the NA catalogue
#' runs to position 468. H5 numbering is the one exception: positions reported in mature
#' H5 numbering are numbered against A/goose/Guangdong/1/1996, which is what that means.
#'
#' The accessions themselves are not written down here or in config.yml. They come from
#' subtype_references.rds, which is the pipeline's own record of which accession each
#' reference strain is in this database version, so the two cannot drift apart.

#' The strains every catalogued position is numbered against.
#'
#' Named rather than accessioned on purpose: an accession is a fact about one database
#' version, a strain is the scientific choice, and subtype_references.rds turns one into
#' the other.
NUMBERING_STRAIN    <- "A/New_York/392/2004"      # H3N2; every product, and H3 numbering
NUMBERING_STRAIN_H5 <- "A/goose/Guangdong/1/1996" # H5N1; mature H5 numbering only

#' The reference a product's positions are numbered against
#'
#' The resolved accession for a product, out of what resolve_numbering_references() settled
#' at startup. HA is the one product with two, so its scheme has to be named; everything
#' else has one numbering and ignores the argument.
#'
#' @param product a product key
#' @param scheme "H5" for mature H5 numbering, "H3" or NA for everything else
#' @param references the resolved accessions, defaulting to the global numbering_refs
#' @return the accession, or NA_character_ where the product has no resolved reference
#' @export
product_reference <- function(product, scheme = NA_character_, references = numbering_refs) {
  if(is.null(product) || length(product) != 1 || is.na(product)) return(NA_character_)

  if(product != "seg4") {
    reference <- references[[product]]
    return(if(is.null(reference)) NA_character_ else reference)
  }

  if(identical(scheme, "H5")) references$ha_h5 else references$ha_h3
}

#' Gapless protein string for a stored sequence, without the terminal stop codon
canonical_protein <- function(sequence) {
  str_remove(str_c(sequence, collapse = ""), "\\*$")
}

#' Accession of a reference strain for one segment
#'
#' @param references the subtype reference table
#' @param strain the strain name, as it appears in that table's `strain` column
#' @param segment the segment number
#' @return the accession, or NA_character_ when the strain has no entry for that segment
reference_accession <- function(references, strain, segment) {
  # Base subsetting rather than filter(): `strain` and `segment` are both an argument
  # and a column here, and the data mask would resolve them to the column either way.
  found <- references$accession[references$strain == strain & references$segment == segment]

  if (length(found) == 0) NA_character_ else found[[1]]
}

#' Resolve one numbering reference
#'
#' The strain's accession for this product's segment is used when it is present in the
#' current database and carries exactly the canonical protein. Identity rather than
#' length, since a different sequence can match the canonical length - the turkey
#' references did for six of the eight segments. Where the two disagree the reference
#' FASTA is stale and the fix is to regenerate it from the same references as the database, not to
#' substitute another accession.
#'
#' @return list(accession, error) - accession is NA_character_ when the product cannot
#'   be numbered, and error then says why
resolve_numbering_reference <- function(seg_key, canonical, strain,
                                        references = subtype_reference_table,
                                        label = seg_key,
                                        segment = product_segment(seg_key),
                                        sequence_positions = discarded) {
  sequences <- sequence_positions[[seg_key]]

  if (is.null(sequences)) {
    log_error("No sequence data for '{seg_key}'; cannot number {label}.")
    return(list(accession = NA_character_, error = "sequence data missing"))
  }

  accession <- reference_accession(references, strain, segment)

  if (is.na(accession)) {
    error <- str_c("Numbering strain '", strain, "' has no entry for segment ", segment,
                   " in the subtype references; ", label, " cannot be numbered.")
    log_error("CRITICAL: {error}")
    return(list(accession = NA_character_, error = error))
  }

  if (!accession %in% names(sequences)) {
    error <- str_c("Numbering reference ", accession, " (", strain, ") for ", label,
                   " is not in the current database.")
    log_error("CRITICAL: {error}")
    return(list(accession = NA_character_, error = error))
  }

  stored <- canonical_protein(sequences[[accession]])

  if (!identical(stored, canonical)) {
    error <- str_c("Numbering reference ", accession, " (", strain, ") for ", label,
                   " is ", nchar(stored), " aa in the database but the reference protein ",
                   "on disk is ", nchar(canonical), " aa",
                   if (nchar(stored) == nchar(canonical)) " and a different sequence" else "",
                   " - regenerate the reference proteins so they match the database")
    log_error("CRITICAL: {error}")
    return(list(accession = NA_character_, error = error))
  }

  log_success("Numbering reference for {label}: {accession} ({strain}).")
  list(accession = accession, error = NULL)
}

#' Resolve the numbering reference for every product
#'
#' Segment 4 is keyed by numbering scheme rather than by segment, because HA
#' positions are reported in mature H3 or mature H5 numbering.
#'
#' Every input is an argument defaulting to its global, so a test can resolve against a
#' handful of sequences and an inline reference table - including the failure paths, which
#' the real data cannot reach without corrupting an artefact first.
#'
#' @param sequence_positions per-sequence gapless positions, keyed by product
#' @param references canonical proteins, keyed by product (AAStringSet each)
#' @param product_keys the products to resolve, HA excluded and handled by scheme
#' @param subtype_refs the subtype reference table, strain to accession per segment
#' @param segments segment number for each product key
#' @return named list: one entry per product bar seg4, plus ha_h3 and ha_h5
resolve_numbering_references <- function(sequence_positions = discarded,
                                         references         = ref_seqs,
                                         product_keys       = available_products,
                                         subtype_refs       = subtype_reference_table,
                                         segments           = product_segment) {
  # Every product but HA, which is keyed by numbering scheme below. The secondary reading
  # frames resolve the same way as any other: their own canonical protein and catalogue
  # rows, their parent segment's row of the reference table.
  keys <- setdiff(product_keys, "seg4")

  resolved <-
    keys %>%
    purrr::set_names() %>%
    map(~ resolve_numbering_reference(.x,
                                      canonical = as.character(references[[.x]][[1]]),
                                      strain = NUMBERING_STRAIN,
                                      references = subtype_refs,
                                      segment = segments(.x),
                                      sequence_positions = sequence_positions))

  if ("seg4" %in% product_keys) {
    resolved$ha_h3 <- resolve_numbering_reference("seg4",
                                                  canonical = as.character(references$seg4[[1]]),
                                                  strain = NUMBERING_STRAIN,
                                                  references = subtype_refs,
                                                  label = "HA (mature H3 numbering)",
                                                  segment = segments("seg4"),
                                                  sequence_positions = sequence_positions)
    resolved$ha_h5 <- resolve_numbering_reference("seg4",
                                                  canonical = as.character(references$seg4[[2]]),
                                                  strain = NUMBERING_STRAIN_H5,
                                                  references = subtype_refs,
                                                  label = "HA (mature H5 numbering)",
                                                  segment = segments("seg4"),
                                                  sequence_positions = sequence_positions)
  }

  resolved
}

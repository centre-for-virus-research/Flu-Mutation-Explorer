#' Build the position lookups for one product from its alignment
#'
#' Reproduces what `gapless_genbank.Rmd` produces for the primary products: `discarded`
#' keeps each sequence's non-gap residues named by their alignment column, and `transposed`
#' is that turned inside out, one entry per column holding only the sequences with a
#' residue there. Rebuilding all eight existing segments this way returns objects identical
#' to the stored ones.
#'
#' Lives here rather than inside secondary_orf_positions.R (now in stale/) so the tests can build a small
#' alignment's position data the same way the pipeline builds the real one; a fixture made
#' any other way would be testing the fixture rather than the app.
#'
#' @param aligned AAStringSet of aligned sequences, named by accession
#' @param accessions which of them to keep, in the order they should be stored
#' @return list(discarded, transposed) for one product
#' @export
build_product <- function(aligned, accessions = names(aligned)) {
  aligned <- aligned[accessions]

  gapless <- map(as.character(aligned), function(sequence) {
    residues <- str_split_1(sequence, "")
    names(residues) <- seq_along(residues)
    residues[residues != "-"] # a gap is not a residue, exactly as discard(grepl("-", .))
  })
  names(gapless) <- names(aligned)

  # Trailing columns gapped in every representative produce no entry and so are not
  # addressable, which is the same reason transposed's NS1 is 240 wide against a
  # 252 column alignment.
  width <- max(as.integer(unlist(map(gapless, names))))

  columns <- map(seq_len(width), function(column) {
    residues <- map_chr(gapless, function(sequence) {
      residue <- sequence[as.character(column)]
      if(is.na(residue)) NA_character_ else unname(residue)
    })
    residues[!is.na(residues)]
  })
  names(columns) <- seq_len(width)

  list(discarded = gapless, transposed = columns)
}

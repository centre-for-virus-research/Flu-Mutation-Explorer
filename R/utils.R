# Amino acid carried by every cluster representative of a product at the alignment
# column that `gapless_position` occupies in `sequence_id`.
#
# `product` is a key into `alignment_columns` - "seg1", "seg7_M2" - not a segment
# number. See gapless_to_consensus() below for why that distinction is worth enforcing.
#
# The two data sets are arguments defaulting to the globals, so a test can pass a small
# fixture instead of loading the real ones. R/alignment_positions.R builds such a fixture
# the same way the pipeline builds the real thing.
consensus <- function(product, sequence_id, gapless_position,
                      sequence_positions = discarded,
                      alignment_columns  = transposed) {
  message("Product: ", product)
  message("Sequence ID: ", sequence_id)
  message("Gapless position: ", gapless_position)

  consensus_position <- gapless_to_consensus(product, sequence_id, gapless_position,
                                             sequence_positions = sequence_positions)

  message("Consensus position: ", consensus_position)
  if(is.na(consensus_position)) {
    return(NULL)
  }

  # every representative's residue at that column
  alignment_columns %>%
    pluck(product, consensus_position) %>%
    unlist %>%
    enframe(name = "Node", value = "amino_acid")
}

# Alignment column occupied by a position in one sequence's own, gapless numbering.
#
# `product` is a key - "seg1", "seg7_M2" - and never a segment number, which is enforced
# below. A number indexes `discarded` by position, and the secondary reading frames are
# appended to that list, so discarded[[7]] stopped agreeing with discarded[["seg7"]]. The
# failure is silent: a valid column of the wrong protein, residues plotted from the wrong
# reading frame.
gapless_to_consensus <- function(product, sequence_id, gapless_position,
                                 sequence_positions = discarded){
  if (!is.character(product)) {
    stop("gapless_to_consensus() takes a product key such as \"seg7\" or \"seg7_M2\", ",
         "not a segment number (got ", class(product)[1], ": ", product, ")")
  }

  if (is.null(sequence_positions[[product]]) ||
      !sequence_id %in% names(sequence_positions[[product]])) {
    return(NA_integer_)
  }

  seq_data <- sequence_positions[[product]][[sequence_id]]

  if (gapless_position > length(seq_data) || gapless_position < 1) {
    return(NA_integer_)
  }

  consensus_position <-
    seq_data[gapless_position] %>%
    names %>% # the alignment column, carried as the element name
    as.integer
  
  if (length(consensus_position) == 0) {
    return(NA_integer_)
  }
  
  consensus_position
}

#' The product a reference record belongs to
#'
#' The reference protein database holds HA twice - seg4 and seg4_h5, the mature H3 and H5
#' numbering references - and they are one product in two numbering schemes rather than two
#' products. Every tab has to collapse them the same way: a position resolved against the
#' wrong one of the pair is a residue reported from the wrong sequence, silently, so the
#' rule is written once here rather than three times.
#'
#' @param record a subject id from the reference database, or an ordinary product key
#' @return the product key
#' @export
record_product <- function(record) if_else(record == "seg4_h5", "seg4", record)

#' The numbering scheme a reference record expresses positions in
#'
#' @param record a subject id from the reference database, or an ordinary product key
#' @return "H3" or "H5" for HA, NA_character_ for everything else, which has one numbering
#' @export
record_scheme <- function(record) {
  dplyr::case_when(record == "seg4"    ~ "H3",
                   record == "seg4_h5" ~ "H5",
                   .default = NA_character_)
}

#' A single-quoted JavaScript string literal, escaped
#'
#' The tables build their own anchors, so values interpolated into an onclick handler are
#' quoted by hand. A backslash or an apostrophe in one would end the literal early and
#' leave the rest of the handler as syntax. Backslash first, or the escapes added after it
#' would themselves be escaped.
#'
#' Sequence names come from a FASTA header the user supplied, so they are the case this
#' exists for; a product key or an accession could not contain either character.
#'
#' @param x length-one character vector
#' @return the same text, quoted and escaped, ready to interpolate
#' @export
js_string <- function(x) {
  x %>%
    str_replace_all("\\\\", "\\\\\\\\") %>%
    str_replace_all("'", "\\\\'") %>%
    str_c("'", ., "'")
}

#' Escape text for use inside a single-quoted HTML attribute
#'
#' The adaptation table builds its own anchors, so tooltip text goes into a title=''
#' attribute by hand. An apostrophe in it would close the attribute early and the rest
#' of the sentence would be parsed as markup.
#'
#' @export
html_attribute <- function(x) {
  x %>%
    str_replace_all("&", "&amp;") %>%   # first, or it would double-escape the rest
    str_replace_all("'", "&#39;") %>%
    str_replace_all('"', "&quot;") %>%
    str_replace_all("<", "&lt;") %>%
    str_replace_all(">", "&gt;")
}

#' Subtype references the Tree tab's Position search can be pointed at
#'
#' The contents of that tab's Subtype dropdown for a product's segment: named by the
#' strain label shown to the user, valued by the accession positions are numbered against.
#' Per segment, since that is how cluster_glue is keyed - both products of segment 7 offer
#' the same list, and both can use it, because M2 carries the parent's representatives.
#'
#' Shared with the Adaptation Mutations tab, which checks that a position's numbering
#' reference is in here before offering to show it on the tree. A missing one would not
#' fail: updateSelectInput() falls back silently to the first choice, and the tree would
#' be coloured by that position in a different sequence.
#'
#' @param product a product key, "seg1" or "seg7_M2"
#' @param references the subtype reference table, defaulting to the global cluster_glue
#' @return named character vector, suitable for updateSelectInput(choices = )
#' @export
subtype_references <- function(product, references = cluster_glue) {
  references %>%
    filter(segment == product_segment(product)) %>%
    select(strain_display, accession) %>% # "representative" in the file; see global.R
    deframe
}

#' Carry a position from one reference strain to another through the alignment column
#'
#' The Tree tab colours by alignment column, so a Position search keeps its colouring when
#' the reference changes: the column is the same, only the number the user reads off differs.
#' `reference_numbering` has one row per residue of every reference in every product, so the
#' old position gives the column and the column gives the new position.
#'
#' Either end can fail without it being an error, and the caller words each differently:
#' `column` is NA when the old position is not in the old reference, and `position` is NA
#' when the new reference has a gap at the column - a shorter HA, an NS1 deletion.
#'
#' @param product a product key, "seg6" or "seg7_M2"
#' @param from,to reference accessions
#' @param position the position in `from`'s own numbering. Ignored when `column` is given.
#' @param column an alignment column, to convert one already known rather than a position
#' @param numbering the reference numbering table, defaulting to the global
#' @return list(column, position), integers, NA where there is none
#' @export
convert_reference_position <- function(product, from, to, position = NA_integer_,
                                       column = NA_integer_,
                                       numbering = reference_numbering) {
  rows <- numbering[numbering$product == product, ]

  if(is.na(column)) {
    column <- rows$column[rows$accession == from & rows$position == position][1]
  }

  new_position <- if(is.na(column)) NA_integer_ else rows$position[rows$accession == to & rows$column == column][1]

  list(column = as.integer(column), position = as.integer(new_position))
}

# Semicolon-separated papers and DOIs as HTML links, one per line
papers_dois <- function(paper, doi) {
  # A product with no catalogue rows - PB1-F2 and PA-X today - reaches here through
  # rowwise() + mutate(), which evaluates the expression once on zero-length columns. The
  # `&` chain below is then logical(0), which if() rejects, failing the whole View mode
  # observer rather than rendering an empty table.
  if (length(paper) == 0 || length(doi) == 0) {
    return("")
  }

  if (!is.na(paper) &
      !is.na(doi) &
      is.character(paper) &
      is.character(doi)) {
    paper_list <- str_split_1(paper, ";") %>% str_trim()
    doi_list <- str_split_1(doi, ";") %>% str_trim()
    
    map2_chr(paper_list, doi_list, function(paper_s, doi){
      pubmed_id <- str_match(
        doi,
        regex("(?:PMID|PubMed(?: ID)?)\\s*[-:]?\\s*(\\d+)", ignore_case = TRUE)
      )[, 2]

      if(!is.na(pubmed_id)) {
        str_c("<a href='https://pubmed.ncbi.nlm.nih.gov/", pubmed_id,
              "/' target='_blank'>", paper_s, "</a>")
      } else {
        str_c("<a href='https://doi.org/", doi, "' target='_blank'>", paper_s, "</a>")
      }
    }) %>% str_flatten(collapse = "<br>")
  } else {
    ""
  }
}

#' Clean a pasted amino acid sequence for BLAST or alignment
#'
#' Accepts what users actually paste: a bare sequence, a sequence wrapped over several
#' lines, or FASTA. A header left in the residues makes AAStringSet fail with "key 62
#' (char '>') not in lookup table", which reaches the user as "Invalid query sequence".
#'
#' Only the first record of a multi-record FASTA is used; these inputs align a single
#' query against a single reference.
#'
#' @param x raw text from the query textarea
#' @return list(sequence, records) - cleaned residues, and how many FASTA
#'   records were present (0 when the input carried no header)
#' @export
clean_query_sequence <- function(x) {
  lines <- x %>% str_split_1("\r?\n") # split before stripping, so headers are whole lines

  headers <- which(str_detect(lines, "^\\s*>"))

  if(length(headers) > 0) {
    # keep the lines after the first header, up to the next one
    first <- headers[1]
    nxt   <- headers[headers > first]
    last  <- if(length(nxt) > 0) nxt[1] - 1 else length(lines)
    lines <- lines[seq.int(first + 1, length.out = max(0, last - first))]
  }

  sequence <-
    lines %>%
    str_c(collapse = "") %>%
    str_remove_all("[\\s\\-\\*]") %>% # gaps, stop codons and any internal whitespace
    str_to_upper

  list(sequence = sequence, records = length(headers))
}

#' Residue the query carries at each reference position
#'
#' Reports, for every position of the reference, which amino acid of the query is aligned
#' to it. That lets the adaptation module ask "does the query carry this catalogued
#' residue" at every catalogued position, rather than only where query and reference
#' differ: a residue the reference also carries never appears in a mismatch table, but is
#' present in the query all the same.
#'
#' Columns where the reference has a gap are insertions in the query and carry no
#' reference position, so they are dropped - full-length HA against the mature HA
#' reference opens with 16 such columns. Columns where the query has a gap keep their
#' position but report NA: not covered is not the same as carrying a different residue.
#'
#' @param pairwise a global PairwiseAlignmentsSingleSubject, query as pattern
#' @return tibble(position, query_amino_acid), one row per reference residue
#' @export
aligned_query_residues <- function(pairwise) {
  query_column     <- pairwise %>% alignedPattern %>% as.character %>% str_split_1("")
  reference_column <- pairwise %>% alignedSubject %>% as.character %>% str_split_1("")

  tibble(reference = reference_column, query = query_column) %>%
    filter(reference != "-") %>% # insertions in the query have no reference position
    mutate(position         = row_number(), # reference positions, gaps now removed
           query_amino_acid = na_if(query, "-")) %>% # a gap is "not covered", not a residue
    select(position, query_amino_acid)
}

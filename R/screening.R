#' Screen many sequences against the adaptation mutation catalogue
#'
#' The Batch Screening tool's engine. Pure - no Shiny, no notifications - so it can be
#' exercised from a script against tests/fixtures/batch_screening_100.fasta.
#'
#' One BLAST call answers three questions at once, which is why there is no translation
#' code here and no pairwise alignment: the subject id of an HSP is the product, the
#' subject coordinates say which catalogued positions the query reaches, and qseq/sseq are
#' the alignment, so the residue at a position is read off rather than computed.
#'
#' Nucleotide input goes through blastx, whose six-frame translation is what reaches the
#' secondary reading frames: PB1-F2, PA-X, M2 and NEP come back as their own HSPs in their
#' own frames, so a genomic segment resolves them with no splice junctions supplied. All
#' 100 of the measured whole segments were recovered.

#' Most sequences accepted in one submission
#'
#' One BLAST call serves the whole submission whatever its size, so the cost is linear in
#' the sequences with no step to fall off. Timed on whole nucleotide segments through
#' blastx, the slower of the two programs, against the 13-protein reference database:
#'
#'   screen_sequences()  100 seqs 1.3 s | 500 5.1 s | 1,000 10.0 s | 2,000 19.9 s
#'   the tab, end to end  100 seqs 3.6 s |                         | 2,000 63.9 s
#'
#' The database is 5,344 residues, which is why the search itself is quick - it is bounded
#' by the queries, not by what they are searched against. The rest is the tab building a
#' View link for each of the 37,000 hit rows 2,000 whole segments produce.
#'
#' The limit is here to keep one submission from occupying the server for an unbounded
#' time, not because anything breaks above it.
#' @export
MAX_SCREEN_SEQUENCES <- 2000L

#' E-value a hit must clear to identify a product
#' @export
SCREEN_STRICT_EVALUE <- 1e-5

#' E-value the search itself runs at, so weaker HSPs exist for the second pass
#' @export
SCREEN_SEARCH_EVALUE <- 10

#' Catalogue key behind a product key
#'
#' seg4 and seg4_h5 are one catalogue read in two numbering schemes, not two catalogues -
#' see record_product() in R/utils.R, which every tab collapses the pair with.
#' @noRd
screen_catalogue_key <- function(product) record_product(product)

#' Which numbering reference a product key resolves to
#'
#' NULL rather than NA where there is none, which is what this tab's callers test for.
#' @noRd
screen_reference <- function(product) {
  reference <- product_reference(record_product(product), record_scheme(product))
  if(is.na(reference)) NULL else reference
}

#' Catalogue positions for a product key, in that key's own numbering
#'
#' @param product product key, including seg4_h5
#' @param catalogues the per-product catalogues, defaulting to the global
#' @return tibble(position, wt, mutant, mutation, paper, doi), empty where none
#' @export
screen_catalogue <- function(product, catalogues = adaptation_mutations_all) {
  rows <- catalogues[[screen_catalogue_key(product)]]

  # Same columns as the full return below, since screen_sequences() uses this as the
  # shape of a screen in which nothing hit
  if(is.null(rows) || nrow(rows) == 0) {
    return(tibble(position = integer(), plot_position = integer(), scheme = character(),
                  wt = character(), mutant = character(),
                  paper = character(), doi = character()))
  }

  # HA carries both numberings on every row; the rest carry one. position is the one the
  # query is aligned and matched in; plot_position is the one the Adaptation Mutations tab
  # wants, which for HA is always H3 - in its view mode numbering_reference() resolves to
  # ha_h3 whatever else is set. The two are the same physical site and resolve to the same
  # alignment column, so the plot is identical either way; only the label differs.
  position <- if(product == "seg4_h5") rows$position_h5 else rows$position

  tibble(position      = as.integer(position),
         plot_position = as.integer(if(product == "seg4_h5") rows$position else position),
         scheme        = record_scheme(product),
         wt       = rows$wt,
         mutant   = rows$new_amino_acid,
         paper    = if("paper_s" %in% names(rows)) rows$paper_s else NA_character_,
         doi      = if("doi" %in% names(rows)) rows$doi else NA_character_) %>%
    dplyr::filter(!is.na(position))
}

#' Run the BLAST search behind a screen
#'
#' Split out so a caller can supply HSPs from elsewhere - a fixture, a cached run - and
#' exercise everything downstream without BLAST installed.
#'
#' @param records tibble(name, sequence) from parse_sequence_records()
#' @param nucleotide TRUE to run blastx, FALSE for blastp
#' @param db the product-keyed reference protein database
#' @return tibble of HSPs, one row each
#' @export
screen_blast <- function(records, nucleotide, db = conf$paths$data$reference_protein_db) {
  query <- tempfile(fileext = ".fasta")
  out   <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(query, out)), add = TRUE)

  writeLines(str_c(">", records$name, "\n", records$sequence), query)

  fields <- c("qseqid", "sseqid", "sstart", "send", "evalue", "bitscore", "qseq", "sseq")

  status <- system2(Sys.which(if(nucleotide) "blastx" else "blastp"),
                    c("-db", shQuote(db), "-query", shQuote(query),
                      "-outfmt", shQuote(str_c("6 ", str_flatten(fields, " "))),
                      "-evalue", SCREEN_SEARCH_EVALUE,
                      "-max_target_seqs", 50,
                      "-out", shQuote(out)),
                    stdout = FALSE, stderr = FALSE)

  if(status != 0 || !file.exists(out)) {
    stop("BLAST failed (exit ", status, ")")
  }

  if(file.size(out) == 0) {
    return(tibble(qseqid = character(), sseqid = character(), sstart = integer(),
                  send = integer(), evalue = numeric(), bitscore = numeric(),
                  qseq = character(), sseq = character()))
  }

  readr::read_tsv(out, col_names = fields, show_col_types = FALSE,
                  col_types = readr::cols(qseqid = "c", sseqid = "c", sstart = "i",
                                          send = "i", evalue = "d", bitscore = "d",
                                          qseq = "c", sseq = "c"))
}

#' Products in play for each query, after the two passes and the HA collapse
#'
#' Pass one identifies products on a strict E-value. Pass two then accepts every other
#' HSP against a subject already identified, however weak, because the question has
#' changed: it is no longer which product this is but which residues it covers. That is
#' what recovers the first exon of a spliced product - NEP's exon 1 returns at E = 0.11,
#' far below any threshold that keeps junk out, and NEP position 7 lives in it.
#'
#' The collapse is the other half. An HA query matches seg4 and seg4_h5 both, and they are
#' one product in two numbering schemes, so keeping both matches the catalogue twice -
#' measured at 810 hit rows where 427 were real.
#'
#' @param hsps as returned by screen_blast()
#' @return the HSPs that survive
#' @export
screen_in_play <- function(hsps) {
  if(nrow(hsps) == 0) return(hsps)

  identified <- hsps %>%
    dplyr::filter(evalue <= SCREEN_STRICT_EVALUE) %>%
    dplyr::distinct(qseqid, sseqid)

  kept <- hsps %>% dplyr::semi_join(identified, by = c("qseqid", "sseqid"))

  # one HA reference per query, by best bitscore
  ha <- kept %>%
    dplyr::filter(str_starts(sseqid, "seg4")) %>%
    dplyr::slice_max(bitscore, n = 1, by = qseqid, with_ties = FALSE) %>%
    dplyr::select(qseqid, chosen = sseqid)

  kept %>%
    dplyr::left_join(ha, by = "qseqid") %>%
    dplyr::filter(!str_starts(sseqid, "seg4") | sseqid == chosen) %>%
    dplyr::select(-chosen)
}

#' Walk an HSP's alignment into one row per subject position
#'
#' cumsum over the ungapped subject columns turns the alignment into subject coordinates,
#' so a catalogued position is a join rather than a search. Columns where the subject is
#' gapped carry no subject position and drop out.
#'
#' @param hsps as returned by screen_in_play()
#' @return tibble(qseqid, product, position, residue)
#' @export
screen_walk <- function(hsps) {
  if(nrow(hsps) == 0) {
    return(tibble(qseqid = character(), product = character(),
                  position = integer(), residue = character()))
  }

  purrr::pmap(list(hsps$qseqid, hsps$sseqid, hsps$sstart, hsps$qseq, hsps$sseq),
    function(query, product, start, qseq, sseq) {
      subject <- str_split_1(sseq, "")
      residue <- str_split_1(qseq, "")
      ungapped <- subject != "-"

      tibble(qseqid = query, product = product,
             position = start + cumsum(ungapped) - 1L,
             residue = residue)[ungapped, ]
    }) %>%
    purrr::list_rbind() %>%
    dplyr::distinct(qseqid, product, position, .keep_all = TRUE)
}

#' Residue a query carries at each position of a reference, read off its HSPs
#'
#' The sibling of aligned_query_residues() in R/utils.R, which walks a pairwise alignment.
#' Same output, same meaning - one row per reference position, NA where the query does not
#' reach - taken from a BLAST alignment instead, so it needs no protein query. That is what
#' lets the Adaptation Mutations tab accept nucleotide: blastx has already translated and
#' aligned it, and aligning a nucleotide query against a protein reference is not something
#' pairwiseAlignment() can be asked to do.
#'
#' It is also the more honest reading for a genomic segment. A spliced product's exons are
#' separate HSPs in their own frames, and the walk takes each on its own; translating the
#' segment in one frame and aligning that would run off the end of the splice junction.
#'
#' @param hsps HSPs as blast_query_segment() returns them, best first, any subject
#' @param sseqid the one subject to read, since hsps carries every product the query hit
#' @param reference_length positions the reference has
#' @return tibble(position, query_amino_acid) with one row per reference position
#' @export
hsp_query_residues <- function(hsps, sseqid, reference_length) {
  # The subject is an argument rather than left to the caller to filter on: hsps holds
  # every product the query hit, and walking them together would read M2's residues onto
  # M1's positions
  mine <- hsps %>% dplyr::filter(sseqid == !!sseqid)

  # screen_walk() keeps the first row it sees for a position and the HSPs arrive sorted by
  # bitscore, so where HSPs overlap the residue comes from the better alignment
  walked <- screen_walk(mine) %>%
    dplyr::select(position, query_amino_acid = residue) %>%
    dplyr::filter(dplyr::between(position, 1L, reference_length))

  tibble(position = seq_len(reference_length)) %>%
    dplyr::left_join(walked, by = "position")
}

#' Screen a set of sequences
#'
#' @param records tibble(name, sequence)
#' @param nucleotide whether the set is nucleotide
#' @return list(hits, summary) - the long table and one row per query
#' @export
screen_sequences <- function(records, nucleotide) {
  hsps   <- screen_blast(records, nucleotide)
  inplay <- screen_in_play(hsps)
  walked <- screen_walk(inplay)

  products <- inplay %>%
    dplyr::slice_max(bitscore, n = 1, by = c(qseqid, sseqid), with_ties = FALSE) %>%
    dplyr::select(qseqid, product = sseqid, evalue)

  # ptype so a screen in which nothing hit still has the columns everything below reads.
  # unnest() over zero products drops them, and the summarise that follows then fails
  # on a missing `position`.
  catalogues <- unique(products$product) %>%
    purrr::map(~ dplyr::mutate(screen_catalogue(.x), product = .x, .before = 1)) %>%
    purrr::list_rbind(ptype = dplyr::mutate(screen_catalogue("", catalogues = list()),
                                            product = character(), .before = 1))

  # assessed-of-total, per query and product. The total counts positions resolvable
  # against the reference this query aligned to, not catalogue rows: a position outside
  # that reference can never be assessed and does not belong in the denominator.
  totals <- catalogues %>%
    dplyr::summarise(catalogued = dplyr::n_distinct(position), .by = product)

  assessed <- walked %>%
    dplyr::inner_join(catalogues %>% dplyr::distinct(product, position),
                      by = c("product", "position")) %>%
    dplyr::summarise(assessed = dplyr::n_distinct(position), .by = c(qseqid, product))

  hits <- walked %>%
    dplyr::inner_join(catalogues, by = c("product", "position"),
                      relationship = "many-to-many") %>%
    dplyr::filter(residue == mutant)

  # Where each query stood in the submitted FASTA. Everything the tab shows and downloads
  # is ordered by this rather than by name, so the results read down beside the file they
  # came from instead of putting sample_10 between sample_1 and sample_2.
  #
  # By first appearance, so a FASTA naming the same sequence twice gets one position rather
  # than fanning the join out. BLAST sees one query id there in any case.
  submission_order <- records %>%
    dplyr::distinct(qseqid = name) %>%
    dplyr::mutate(submitted = dplyr::row_number())

  hits <- hits %>%
    dplyr::left_join(assessed, by = c("qseqid", "product")) %>%
    dplyr::left_join(totals,   by = "product") %>%
    dplyr::left_join(submission_order, by = "qseqid") %>%
    dplyr::arrange(submitted, product, position)

  # Every submitted sequence gets a row, including those that resolved to nothing.
  summary <- records %>%
    dplyr::select(qseqid = name) %>%
    dplyr::left_join(submission_order, by = "qseqid") %>%
    dplyr::left_join(products %>% dplyr::slice_min(evalue, n = 1, by = qseqid, with_ties = FALSE),
                     by = "qseqid") %>%
    dplyr::left_join(assessed, by = c("qseqid", "product")) %>%
    dplyr::left_join(totals,   by = "product") %>%
    dplyr::left_join(hits %>% dplyr::summarise(hits = dplyr::n(), .by = qseqid),
                     by = "qseqid") %>%
    dplyr::mutate(
      hits     = tidyr::replace_na(hits, 0L),
      assessed = tidyr::replace_na(assessed, 0L),
      catalogued = tidyr::replace_na(catalogued, 0L),
      # Four things produce no hit and only one means the sequence came back clean, so the
      # outcome is named rather than left as an empty cell. Partial is the dangerous one:
      # "none found" over a fragment is not the claim "none found" over a whole protein is.
      outcome = dplyr::case_when(
        is.na(product)                  ~ "Not screened",
        catalogued == 0                 ~ "Nothing catalogued",
        hits > 0 & assessed < catalogued ~ "Screened (partial)",
        hits > 0                        ~ "Screened",
        assessed < catalogued           ~ "None found (partial)",
        .default                        = "None found"),
      segment = dplyr::if_else(is.na(product), NA_integer_,
                               as.integer(str_extract(product, "(?<=^seg)\\d"))),
      product_label = dplyr::case_when(
        is.na(product)       ~ NA_character_,
        product == "seg4_h5" ~ "HA (H5 numbering)",
        product == "seg4"    ~ "HA (H3 numbering)",
        .default = purrr::map_chr(product, ~ products_label_safe(.x)))) %>%
    dplyr::arrange(submitted)

  list(hits = hits, summary = summary)
}

#' The one table the tab displays
#'
#' One row per catalogued position found, plus one row for every query that found none.
#' A screening tool that drops a sequence from its own results reports a failure as an
#' all-clear, so a query is never absent - if it has nothing to report, the Outcome column
#' says which of the four reasons applies, and only one of them means it came back clean.
#'
#' @param screened as returned by screen_sequences()
#' @return tibble of display columns, in the order the sequences were submitted, then by
#'   product and position within each
#' @export
screen_table <- function(screened) {
  found <- screened$hits %>%
    dplyr::transmute(qseqid, submitted, product,
                     # WT, position and the residue found, as the catalogue writes them:
                     # T160A rather than three columns the reader has to reassemble
                     mutation = str_c(wt, position, residue),
                     position, plot_position, scheme,
                     outcome = "Found", assessed, catalogued,
                     reference = paper, doi)

  # Query, Product and Assessed are repeated on the empty rows rather than blanked: the
  # table is sorted, filtered and downloaded, so a row has to say what it is on its own.
  empty <- screened$summary %>%
    dplyr::filter(hits == 0) %>%
    dplyr::transmute(qseqid, submitted, product, mutation = NA_character_,
                     position = NA_integer_, plot_position = NA_integer_,
                     scheme = NA_character_,
                     outcome, assessed, catalogued,
                     reference = NA_character_, doi = NA_character_)

  # Submission order first, so a query that found nothing sits among the ones that did,
  # where it was submitted, rather than after them. Then product and position within a
  # query, which is the order the catalogue reads in.
  dplyr::bind_rows(found, empty) %>%
    dplyr::arrange(submitted, product, position)
}

#' The screen as a flat table for download
#'
#' The same rows the tab shows, under names that survive leaving the app: the View column
#' is dropped, being two links, and Assessed is split back into the two numbers the table
#' renders as "12 of 32". Everything the screen concluded is here.
#'
#' Its own function rather than a transform inside the download handler, so the column
#' names it depends on are exercised by the tests rather than only on a user's click.
#'
#' @param screened as returned by screen_sequences(), or NULL before any screen
#' @return tibble of download columns, empty but correctly shaped when screened is NULL
#' @export
screen_download <- function(screened = NULL) {
  rows <- if(is.null(screened)) {
    tibble(qseqid = character(), product = character(), mutation = character(),
           position = integer(), plot_position = integer(), scheme = character(),
           outcome = character(), assessed = integer(), catalogued = integer(),
           reference = character(), doi = character())
  } else {
    screen_table(screened)
  }

  rows %>%
    dplyr::transmute(
      query = qseqid,
      # read off the product key, which is why this is taken before product is relabelled
      segment    = as.integer(str_extract(product, "(?<=^seg)\\d")),
      product    = screen_product_label(product),
      # NA for everything but HA, which is catalogued in two numbering schemes and whose
      # position means nothing without the one it is expressed in
      numbering  = scheme,
      mutation, position, outcome, assessed, catalogued,
      paper = reference, doi)
}

#' Display label for a product key, with HA's two references reading as one product
#'
#' seg4 and seg4_h5 are one protein in two numbering schemes, so both read "HA". Which
#' scheme a row is numbered in is said on the position, where it belongs.
#' @export
screen_product_label <- function(product) {
  dplyr::case_when(is.na(product) ~ NA_character_,
                   str_starts(product, "seg4") ~ "HA",
                   .default = purrr::map_chr(product, products_label_safe))
}

#' A product's display label, tolerating a key the products list does not hold
#'
#' NA in, NA out: a query that resolved to nothing has no product, and if_else() evaluates
#' both its branches, so this is called on those rows whatever the guard around it says.
#' @noRd
products_label_safe <- function(key) {
  if(length(key) != 1 || is.na(key)) return(NA_character_)
  entry <- products[[key]]
  if(is.null(entry)) key else entry$label
}

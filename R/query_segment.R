#' Segment detection for a submitted sequence
#'
#' The Tree and Adaptation Mutations tabs open the same way: take what the user pasted,
#' decide whether it is protein or nucleotide, BLAST it against the product references, and
#' read the segment off the top hit.
#'
#' Split in two, so the BLAST and the lookup are pure and can be tested without a session,
#' and a thin wrapper adds the notifications the two tabs would otherwise duplicate.
#'
#' is_nucleotide() decides which of blastp and blastx runs, not which queries to refuse:
#' all three tabs search the same reference proteins, and blastx does the six-frame
#' translation.

#' Letters that cannot be a nucleotide under any IUPAC code
#'
#' Presence of any of these proves a sequence is protein, whatever its composition.
#' U is deliberately absent: it is uracil, so vetoing it would read every RNA sequence
#' as protein. BLAST treats U as T by itself, so RNA needs no other handling.
#'
#' These nine are 33.6% of residues in the reference proteins, so a peptide avoiding all
#' of them runs at 1.67% for a 10-mer and 0.028% for a 20-mer.
#' @export
PROTEIN_ONLY_LETTERS <- c("E", "F", "I", "J", "L", "P", "Q", "X", "Z")

#' Smallest fraction of nucleotide letters for a sequence to be read as nucleotide
#'
#' Measured separation on real data is enormous - 100 sampled segments run 0.9988 to
#' 1.000, the reference proteins top out at 0.3055 - so anything from 0.4 to 0.9 works
#' and the exact value is not delicate.
#' @export
NUCLEOTIDE_FRACTION <- 0.9

#' Is a sequence nucleotide?
#'
#' Two tests, because neither alone is enough. The veto is exact but only fires when one
#' of nine letters happens to be present; the composition test always answers but is a
#' threshold. A closed "only ACGTN" alphabet is not usable: real GenBank sequences carry
#' ambiguity codes - 5 of 100 sampled segments hold Y, K or N - and every ambiguity code
#' is also a valid amino acid letter, so those five would be read as protein and BLASTed
#' as nonsense rather than rejected.
#'
#' @param x character vector of sequences
#' @return logical vector
#' @export
is_nucleotide <- function(x) {
  x <- str_to_upper(x)
  vetoed <- str_detect(x, str_c("[", str_flatten(PROTEIN_ONLY_LETTERS), "]"))
  fraction <- str_count(x, "[ACGTUN]") / pmax(nchar(x), 1L)

  !vetoed & fraction >= NUCLEOTIDE_FRACTION
}

#' Parse pasted text or an uploaded file into named records
#'
#' The multi-record sibling of clean_query_sequence(), which keeps only the first record
#' because its callers align one query against one reference. This keeps all of them.
#'
#' Sequence text with no header at all is one unnamed record, so a user who pastes a bare
#' sequence is not told their input is empty.
#'
#' @param text raw FASTA or bare sequence
#' @return tibble(name, sequence), residues upper-cased with gaps and whitespace removed
#' @export
parse_sequence_records <- function(text) {
  lines <- text %>% str_split_1("\r?\n") %>% str_trim()
  lines <- lines[nzchar(lines)]

  if(length(lines) == 0) {
    return(tibble(name = character(), sequence = character()))
  }

  headers <- str_starts(lines, ">")

  if(!any(headers)) {
    return(tibble(name = "query",
                  sequence = str_to_upper(str_remove_all(str_flatten(lines), "[^A-Za-z*]"))))
  }

  # Lines before the first header are not part of any record and are dropped
  record <- cumsum(headers)
  keep   <- record > 0

  tibble(line = lines[keep], record = record[keep], header = headers[keep]) %>%
    summarise(
      # the first token of the header, as every FASTA reader takes the id
      name = str_remove(word(line[header][1], 1), "^>"),
      sequence = str_to_upper(str_remove_all(str_flatten(line[!header]), "[^A-Za-z*]")),
      .by = record) %>%
    mutate(name = if_else(is.na(name) | !nzchar(name), str_c("record_", record), name)) %>%
    select(name, sequence)
}

#' Smallest bitscore a hit must reach to identify a query
#'
#' A bitscore rather than an E-value floor, because E-value scales with the size of the
#' database: the 1e-3 used against the 748-protein set's 319,693 residues would be some
#' sixty times more permissive against the 13 product references' 5,344. A bitscore is a
#' property of the alignment alone, so this number means the same whatever the database is
#' next rebuilt to hold.
#'
#' Measured with the 748-protein Nextalign set as signal and three kinds of junk as noise
#' - nucleotide pasted into the protein box, homopolymers, and shuffled flu proteins.
#' Noise tops out at 28 bits; no real query of 100 aa or more scores below 52. Real
#' queries rejected, by query length, against junk admitted out of 300:
#'
#'   floor   full    200 aa   100 aa   50 aa    junk
#'      30   0/748    0/608    0/660    4/748   0/300
#'      40   0/748    0/608    0/660   23/748   0/300
#'      50   0/748    0/608    0/660   65/748   0/300
#'
#' 40 keeps a 1.4x margin over the worst junk and rejects nothing of 100 aa or more. The
#' 23 it turns away are all 50 aa fragments, where declining to answer is the honest
#' result rather than a loss.
#' @export
MIN_SEGMENT_BITSCORE <- 40

#' Segment and product of a query sequence, by BLAST against the product references
#'
#' Pure - no Shiny, no notifications. Failures come back as a result rather than as a
#' condition, since both are ordinary outcomes of a user-supplied sequence: nothing
#' matched well enough, or the best match is not a product the app knows.
#'
#' The database is the one the Batch Screening tab searches, holding one record per
#' product keyed by the product key itself, so the subject id *is* the answer and no
#' reference-set table is needed to look the hit up in. Against the 748-protein recognizer
#' set as a labelled panel, the 13 references name the product correctly for 748/748
#' full-length proteins, 608/608 200 aa fragments and 652/660 at 100 aa, and never name
#' the wrong segment at any length tested.
#'
#' HA is two records, seg4 and seg4_h5, being the mature H3 and H5 numbering references.
#' They are one product in two numbering schemes, so the product collapses to seg4 and
#' which record won is reported as `scheme` - the fact the Adaptation Mutations tab
#' otherwise needs a second BLAST against Ref_H3_H5/ for.
#'
#' Nucleotide input goes through blastx, which does the six-frame translation itself. The
#' returned HSPs carry the alignment, so a caller wanting the query's residue at a position
#' of the reference reads it off rather than aligning again - see hsp_query_residues().
#'
#' @param query AAStringSet or DNAStringSet holding the one query
#' @param nucleotide TRUE to run blastx, FALSE for blastp
#' @param db BLAST database of the product references
#' @param known product key -> segment/label list, defaulting to the global products
#' @param min_bitscore floor a hit must reach to identify the query
#' @return list(ok, segment, product, scheme, sseqid, identified, hsps) on success,
#'   list(ok = FALSE, message) otherwise. scheme is "H3" or "H5" for HA and NA otherwise;
#'   identified names every subject that cleared the floor, best first, which for a whole
#'   genomic segment is both of its products; hsps are every HSP of the search, best first,
#'   for a caller that wants the alignment as well as the answer.
#' @export
blast_query_segment <- function(query, nucleotide = FALSE,
                                db = conf$paths$data$reference_protein_db,
                                known = products,
                                min_bitscore = MIN_SEGMENT_BITSCORE) {
  # The fields screen_walk() reads, so its walk works on these HSPs unchanged
  fields <- "qseqid sseqid evalue bitscore sstart send qseq sseq"

  hits <-
    predict(blast(db = db, type = if(nucleotide) "blastx" else "blastp"), query,
            custom_format = fields, verbose = TRUE) %>%
    # Bitscore first, E-value only to break its ties. The secondary products overlap their
    # primary - PA-X is the first 191 aa of PA plus 40-60 of its own - so a PA-X query hits
    # both references at the same E-value, and sorting on E-value alone left the winner to
    # whatever order BLAST happened to return.
    arrange(desc(bitscore), evalue)

  # Without a floor the top hit is taken however bad it is, and junk is reported as a
  # confident answer: 52 of 100 whole nucleotide segments read as protein find a hit here,
  # the best of them 27 bits. See MIN_SEGMENT_BITSCORE for why the floor is a bitscore.
  identified <- hits %>% dplyr::filter(bitscore >= min_bitscore)

  if(nrow(identified) == 0) {
    return(list(ok = FALSE, message = "Unable to identify segment"))
  }

  sseqid <- identified %>% dplyr::slice(1) %>% pull(sseqid) # top hit

  # seg4 and seg4_h5 are HA in two numbering schemes, not two products
  product <- record_product(sseqid)
  entry   <- known[[product]]

  # A subject the products list does not hold would leave segment empty and fail as
  # "argument is of length zero" further down. Reachable if the database is rebuilt from a
  # FASTA carrying keys the app does not know.
  if(is.null(entry)) {
    return(list(ok = FALSE,
                message = str_c("Top BLAST hit ", sseqid,
                                " is not a known product - unable to identify segment")))
  }

  # Every HSP is handed back, not only those clearing the floor. The floor answers "which
  # product is this"; once that is settled the question is "which residues does the query
  # cover", and a weak HSP still answers it - NEP's exon 1 comes back far below any
  # threshold that keeps junk out, and NEP position 7 lives in it.
  #
  # identified is every subject that cleared the floor, not only the winner: a whole
  # genomic segment of 2, 3, 7 or 8 is two products, and the Adaptation Mutations tab
  # reports both.
  list(ok         = TRUE,
       segment    = entry$segment,
       product    = product,
       scheme     = record_scheme(sseqid),
       sseqid     = sseqid,
       identified = identified %>% dplyr::distinct(sseqid) %>% pull(sseqid),
       hsps       = hits)
}

#' Read a submitted query and resolve the segment it came from, reporting to the user
#'
#' Everything the Tree and Adaptation Mutations tools do in common with a submitted
#' sequence, notifications included. Returns NULL when the query cannot be used; the
#' reason has already been shown, so the caller only has to stop.
#'
#' Call from an observer - it notifies, so it needs a reactive context.
#'
#' @param text raw contents of the query textarea
#' @return list(seqs, nucleotide, segment, product, scheme, sseqid, identified, hsps),
#'   or NULL
#' @export
resolve_query_segment <- function(text) {
  query <- clean_query_sequence(text) # accepts bare sequence or FASTA

  if(query$records > 1) {
    showNotification("Multiple FASTA records found - using the first", type = "message")
  }

  if(!nzchar(query$sequence)) {
    showNotification("No sequence found in the query", type = "error")
    return(NULL)
  }

  # Nucleotide is not refused here: these tabs search the same reference proteins Batch
  # Screening does, so it goes to blastx and the six-frame translation is BLAST's to do.
  nucleotide <- is_nucleotide(query$sequence)

  seqs <- tryCatch(
    if(nucleotide) {
      # U is uracil and DNAStringSet has no letter for it, though BLAST reads it as T.
      # An RNA sequence is otherwise an ordinary query and is not worth refusing.
      DNAStringSet(str_replace_all(query$sequence, "U", "T"))
    } else {
      AAStringSet(query$sequence)
    },
    error = function(e) {
      log_warn("Query sequence rejected: {e$message}")
      showNotification("Invalid query sequence", type = "error")
      NULL
    })

  if(is.null(seqs)) {
    return(NULL)
  }

  log_info("Query read as {if(nucleotide) 'nucleotide' else 'protein'} ",
           "({nchar(query$sequence)} {if(nucleotide) 'nt' else 'aa'})")

  detected <- blast_query_segment(seqs, nucleotide = nucleotide)

  if(!detected$ok) {
    showNotification(detected$message, type = "error")
    return(NULL)
  }

  log_info("Query matched {detected$sseqid} - segment {detected$segment}, ",
           "product {detected$product}")
  showNotification(str_c("Detected segment ", detected$segment), type = "message")

  list(seqs = seqs, nucleotide = nucleotide, segment = detected$segment,
       product = detected$product, scheme = detected$scheme,
       sseqid = detected$sseqid, identified = detected$identified,
       hsps = detected$hsps)
}

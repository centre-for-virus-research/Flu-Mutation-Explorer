#' Loading the app's data
#'
#' Everything global.R used to build at top level, as functions of the config and of the
#' objects each step depends on, so a test can load a small fixture through the same code
#' the app runs rather than sourcing the real data. read_app_data() runs them in order and
#' returns one named list; global.R copies its entries into the global environment, which is
#' where the modules and the older tests still look for them.
#'
#' Each loader says what it needs as an argument and returns what it built. None of them
#' reads a global, so the order they run in is the order in read_app_data() and nothing else.

#' The trees, as the app holds them
#'
#' @param trees_dir directory of `seg<n>_cluster_rep.nwk` files
#' @return list of tree_names (`seg1` ... `seg8`), tree_data (a Newick string per name) and
#'   trees_multi (the same trees as an ape MultiPhylo, for the tip labels)
load_trees <- function(trees_dir) {
  tree_files <- list.files(path = trees_dir,
                           pattern = "\\.nwk$",
                           full.names = TRUE)

  tree_names <-
    tree_files %>%
    basename %>%
    map_vec(tools::file_path_sans_ext) %>%
    map_vec(tools::file_path_sans_ext) %>%
    map_vec(tools::file_path_sans_ext) %>%
    str_replace("_cluster_rep", "")

  log_info("Tree names: {paste(tree_names, collapse=', ')}")

  # Newick string per tree name
  tree_data <-
    tree_files %>%
    map(read_file) %>%
    setNames(tree_names)

  list(tree_names = tree_names,
       tree_data  = tree_data,
       trees_multi = tree_data %>%
         paste(collapse = "") %>%
         read.tree(text = ., tree.names = tree_data %>% names))
}

#' One row per cluster representative, with the host composition of its cluster
#'
#' host_group arrives curated - Birds, Human, Other Mammals, Environment, Unknown - so
#' there is nothing to recode here; the taxonomy stage of the pipeline does that work.
#'
#' The composition is precomputed by cluster_host_composition.R. A tip is a cluster
#' representative standing for up to tens of thousands of members, so its own host says
#' little about the cluster.
load_full_metadata <- function(conf) {
  full_metadata <-
    read_rds(conf$paths$data$metadata) %>%
    select(parsed_strain, primary_accession, h_subtype, n_subtype, host_group, host_order, segment) %>%
    dplyr::rename(strain = parsed_strain)

  if(!file.exists(conf$paths$data$cluster_composition)) {
    stop("Missing ", conf$paths$data$cluster_composition,
         " - generate it with: Rscript cluster_host_composition.R")
  }

  cluster_composition <- read_rds(conf$paths$data$cluster_composition)

  full_metadata <-
    full_metadata %>%
    left_join(cluster_composition %>% select(-cluster_size),
              by = c("primary_accession" = "representative", "segment"))

  missing_composition <- sum(is.na(full_metadata$cluster_hosts))
  if(missing_composition > 0) {
    # Rather than a blank pop-up field, say the composition is unavailable
    log_warn("{missing_composition} representative(s) have no cluster composition - rerun cluster_host_composition.R")
    full_metadata <-
      full_metadata %>%
      mutate(cluster_members = replace_na(cluster_members, "cluster size unavailable"),
             cluster_hosts   = replace_na(cluster_hosts, "host composition unavailable"))
  }

  full_metadata
}

#' The per-sequence and per-column position data for every product
#'
#' product_positions.rds and product_columns.rds hold all twelve products - the eight
#' primary ones and the secondary reading frames PB1-F2, PA-X, M2 and NEP - keyed by product
#' ("seg1", "seg7_M2"), which is how everything in the app is addressed. Checked against the
#' tree names rather than trusted: a tree without position data would offer a segment that
#' cannot be searched.
#'
#' @param products the product table from global.R, to find each product's parent segment
#' @return list of discarded (residues, named by alignment column) and transposed (the same
#'   data, per alignment column)
load_positions <- function(conf, tree_names, products) {
  check_product_keys <- function(data, what) {
    missing <- setdiff(tree_names, names(data))
    if(length(missing) > 0) {
      stop(what, " has no entry for tree(s) ", str_flatten_comma(missing),
           " - it holds ", str_flatten_comma(names(data) %||% "nothing named"))
    }

    data
  }

  discarded <- read_rds(conf$paths$data$discarded) %>%
    check_product_keys(conf$paths$data$discarded)

  transposed <- read_rds(conf$paths$data$transposed) %>%
    check_product_keys(conf$paths$data$transposed)

  # The secondary frames share their parent segment's representatives, so one whose
  # representatives had drifted would report a different population under the same tree.
  # Checked on discarded, which is keyed by sequence; transposed is keyed by column.
  iwalk(discarded, function(sequences, key) {
    if(is.null(products[[key]])) return(invisible()) # not a product this app offers
    parent <- str_c("seg", products[[key]]$segment) # by key, never by list position
    if(!setequal(names(sequences), names(discarded[[parent]]))) {
      stop(conf$paths$data$discarded, ": ", key, " representatives do not match those of ", parent)
    }
  })

  # Both are addressed by the same product key, and gapless_to_consensus() resolves a
  # column in one that consensus() then reads out of the other. Keyed differently they
  # would return a column of the wrong protein rather than fail.
  if(!identical(names(discarded), names(transposed))) {
    stop("discarded and transposed are keyed differently - discarded: ",
         str_flatten_comma(names(discarded)), "; transposed: ",
         str_flatten_comma(names(transposed)),
         ". Both are named from the tree files, so a tree without position data (or the ",
         "reverse) will do this - rebuild them together.")
  }

  list(discarded = discarded, transposed = transposed)
}

#' Amino acid counts per alignment column over every sequence in the database
#'
#' Built by alignment_column_counts.R. transposed holds only the cluster representatives -
#' one vote per cluster however many sequences it stands for - so this table lets the
#' adaptation bar plot count them all. It shares transposed's column coordinates, so
#' gapless_to_consensus() addresses both.
#'
#' Optional on purpose: a deployment that has not run that pipeline stage still starts, and
#' offers the cluster representative plot on its own.
#'
#' @return list of alignment_counts, alignment_counts_source, alignment_counts_partial and
#'   alignment_counts_totals; the first and last NULL when the file is absent
load_alignment_counts <- function(conf, available_products, transposed) {
  alignment_counts_path <- conf$paths$data$alignment_column_counts

  if(is.null(alignment_counts_path) || !file.exists(alignment_counts_path)) {
    log_warn("No alignment column counts found at {alignment_counts_path %||% '<unset>'} - ",
             "the Adaptation Mutations plot will offer cluster representatives only. ",
             "Build them with: Rscript alignment_column_counts.R")

    return(list(alignment_counts         = NULL,
                alignment_counts_source  = NA_character_,
                alignment_counts_partial = FALSE,
                alignment_counts_totals  = NULL))
  }

  alignment_counts <- read_rds(alignment_counts_path)

  # A table built against a different alignment would silently plot the wrong column
  overrun <-
    alignment_counts %>%
    dplyr::summarise(max_column = max(column), .by = product) %>%
    dplyr::filter(max_column > map_int(product, ~ length(transposed[[.x]]) %||% 0L))

  if(nrow(overrun) > 0) {
    stop(conf$paths$data$alignment_column_counts, " holds columns beyond the alignment ",
         "for product(s) ", str_flatten_comma(overrun$product),
         " - rebuild it with alignment_column_counts.R")
  }

  # A product with no counts silently falls back to the cluster plot in that mode
  uncounted <- setdiff(available_products, unique(alignment_counts$product))
  if(length(uncounted) > 0) {
    log_warn("No All sequences counts for {str_flatten_comma(uncounted)} - ",
             "rebuild with: Rscript alignment_column_counts.R")
  }

  # The counts are only as complete as the host taxonomy they were built from. On the
  # fallback a few percent of sequences have a host that could not be placed and are
  # absent from the bars, so the plot says so rather than presenting it as everything.
  alignment_counts_source  <- attr(alignment_counts, "host_source") %||% "unknown"
  alignment_counts_partial <- !str_detect(alignment_counts_source, "host_groups\\.rds$")

  # Sequences that entered the tally, per product and level - the denominator the plot
  # caption reports a column's coverage against.
  #
  # A table built before that attribute existed still has to give the caption a number, so
  # fall back to the best covered column. That is a lower bound, since every sequence is
  # gapped somewhere but not all at the same column, so it understates how much of the
  # alignment is missing at the position on screen.
  alignment_counts_totals <- attr(alignment_counts, "product_totals")

  if(is.null(alignment_counts_totals)) {
    alignment_counts_totals <-
      alignment_counts %>%
      dplyr::summarise(counted = sum(n), .by = c(product, level, column)) %>%
      dplyr::summarise(sequences = max(counted), .by = c(product, level))

    log_warn("{alignment_counts_path} predates the product_totals attribute - the plot ",
             "caption will report coverage against the best covered column, which ",
             "understates it. Rebuild with: Rscript alignment_column_counts.R")
  }

  log_info("Alignment column counts: {nrow(alignment_counts)} rows for ",
           "{n_distinct(alignment_counts$product)} product(s), host taxonomy from {alignment_counts_source}")

  if(alignment_counts_partial) {
    log_warn("Alignment column counts were built without flu_gdb_app_data/host_groups.rds - ",
             "sequences whose host could not be placed are excluded from the All sequences plot")
  }

  list(alignment_counts         = alignment_counts,
       alignment_counts_source  = alignment_counts_source,
       alignment_counts_partial = alignment_counts_partial,
       alignment_counts_totals  = alignment_counts_totals)
}

#' What each cluster carries at each alignment column
#'
#' Built by cluster_column_residues.R. A tree tip stands for a whole cluster, so after a
#' Position search the pop-up reports the cluster's residues rather than only the
#' representative's - the same reason cluster_host_composition.rds exists for hosts.
#'
#' Optional, like the other two: without it the tree still colours by residue and the
#' pop-up carries no breakdown.
#'
#' @return the table, or NULL when the file is absent
load_cluster_residues <- function(conf, transposed) {
  cluster_residues_path <- conf$paths$data$cluster_residues

  if(is.null(cluster_residues_path) || !file.exists(cluster_residues_path)) {
    log_warn("No cluster residue breakdown at {cluster_residues_path %||% '<unset>'} - ",
             "tree pop-ups will show the representative's own residue only. ",
             "Build it with: Rscript cluster_column_residues.R")
    return(NULL)
  }

  cluster_residues <- read_rds(cluster_residues_path)

  # Built against a different alignment it would report another column's residues
  residue_overrun <-
    cluster_residues %>%
    dplyr::summarise(max_column = max(column), .by = product) %>%
    dplyr::filter(max_column > map_int(product, ~ length(transposed[[.x]]) %||% 0L))

  if(nrow(residue_overrun) > 0) {
    stop(cluster_residues_path, " holds columns beyond the alignment for product(s) ",
         str_flatten_comma(residue_overrun$product),
         " - rebuild it with: Rscript cluster_column_residues.R")
  }

  log_info("Cluster residue breakdown: {nrow(cluster_residues)} rows for ",
           "{n_distinct(cluster_residues$product)} product(s), ",
           "{n_distinct(cluster_residues$representative)} cluster(s)")

  cluster_residues
}

#' One product's rows out of a catalogue covering more than one, keyed on the label in the
#' segment column. Trailing blank rows carry no label and drop out with everything else
#' that is not this product.
#'
#' The files differ in which trailing columns they carry - extra_notes on PA, M and NS,
#' further_notes plus unnamed spreadsheet columns on PB1 - so the drop list is the union
#' and any_of() ignores the names a given file does not have.
orf_rows <- function(path, orf) {
  rows <-
    read_csv(path, show_col_types = FALSE) %>%
    select(-1) %>%
    clean_names() %>%
    dplyr::filter(segment == orf) %>%
    dplyr::select(-segment, -experimentally_verified,
                  -any_of(c("extra_notes", "further_notes")),
                  -matches("^x(_\\d+)?$")) %>% # unnamed spreadsheet columns
    dplyr::rename(new_amino_acid = mutant)

  if(nrow(rows) == 0) {
    log_warn("No {orf} rows in {path} - that product will show an empty table")
  }

  rows
}

#' The curated adaptation mutation catalogues, one table per product
#'
#' @param available_products products with position data; the rest are dropped, since a
#'   table with no positions behind it cannot plot
#' @return named list keyed by product, each with the Paper column as the tables show it
load_adaptation_mutations <- function(conf, available_products) {
  paths <- conf$paths$adaptation_mutations

  adaptation_mutations_all <- list(
    seg1 = read_csv(paths$seg1, show_col_types = FALSE) %>%
      select(-1) %>% clean_names() %>%
      dplyr::select(-segment, -further_notes, -experimentally_verified) %>%
      dplyr::rename(doi = x, new_amino_acid = mutant),

    # The PB1 and PA files can hold their segment's second reading frame as well, exactly
    # as the M and NS files do. Those rows are taken out here and read back in as their own
    # products below; rows with no label at all are kept, since they are PB1 and PA entries
    # whose segment column was simply left blank rather than another product's.
    seg2 = read_csv(paths$seg2, show_col_types = FALSE) %>%
      select(-1) %>% clean_names() %>%
      dplyr::filter(is.na(segment) | segment != "PB1-F2") %>%
      dplyr::select(-segment, -further_notes, -experimentally_verified, -starts_with("x")) %>%
      dplyr::rename(new_amino_acid = mutant),

    seg3 = read_csv(paths$seg3, show_col_types = FALSE) %>%
      select(-1) %>% clean_names() %>%
      dplyr::filter(is.na(segment) | segment != "PA-X") %>%
      dplyr::select(-segment, -extra_notes, -experimentally_verified) %>%
      dplyr::rename(new_amino_acid = mutant),

    seg4 = read_csv(paths$seg4, show_col_types = FALSE) %>%
      select(-1) %>% clean_names() %>%
      mutate(position_h5 = as.integer(str_extract(mutation_h5_numbering, "-?\\d+")), .after = position) %>%
      dplyr::select(-segment, -experimentally_verified) %>%
      dplyr::rename(new_amino_acid = mutation, mutation = mutation_h3_numbering),

    seg5 = read_csv(paths$seg5, show_col_types = FALSE) %>%
      select(-1) %>% clean_names() %>%
      dplyr::select(-segment, -extra_notes, -experimentally_verified) %>%
      dplyr::rename(new_amino_acid = mutant),

    seg6 = read_csv(paths$seg6, show_col_types = FALSE) %>%
      select(-1) %>% clean_names() %>%
      dplyr::select(-segment, -x, -experimentally_verified) %>%
      dplyr::rename(new_amino_acid = mutation_1),

    # Segments 7 and 8 hold two products each in one file, told apart by the segment column
    # the other catalogues drop. It has to be read here: the M2 and NEP positions all fall
    # inside M1's 252 and NS1's 230, so rows left unsplit resolve silently to the wrong
    # residue rather than failing.
    seg7     = orf_rows(paths$seg7, "M1"),
    seg7_M2  = orf_rows(paths$seg7, "M2"),
    seg8     = orf_rows(paths$seg8, "NS1"),
    seg8_NEP = orf_rows(paths$seg8, "NEP"),

    # PB1-F2 and PA-X are read the same way, out of the PB1 and PA files. Neither has any
    # entry yet, so both tables are empty and orf_rows() says so in the log. Reading them now
    # rather than when the first row appears keeps that row from being numbered against PB1
    # or PA, which it would be silently: PB1-F2's 90 positions all fall inside PB1's 757.
    seg2_PB1_F2 = orf_rows(paths$seg2, "PB1-F2"),
    seg3_PA_X   = orf_rows(paths$seg3, "PA-X")
  )

  # Trailing blank rows out of the spreadsheets - one each at the end of the PB2, PB1, PA
  # and NP files. They carry no position, so they can never be screened or plotted, and they
  # survive the filters above because seg2 and seg3 deliberately keep rows with a blank
  # segment column (unlabelled PB1 and PA entries). Left in, they sit in every denominator:
  # PB2 would report "158 of 159" for a sequence covering the whole protein.
  adaptation_mutations_all <- map(adaptation_mutations_all, ~ remove_empty(.x, "rows"))

  # Products with no position data loaded would offer a table that cannot plot
  adaptation_mutations_all <- adaptation_mutations_all[available_products]

  # The Paper column as the tables show it: a row's papers as links, built once here rather
  # than per render. papers_dois() takes one row at a time, so the display code ran it under
  # rowwise() over a product's whole catalogue every time one was selected.
  map(adaptation_mutations_all,
      ~ .x %>% mutate(paper_html = map2_chr(paper_s, doi, papers_dois)))
}

#' The metadata CSV for each tree, precalculated for faster session startup
#'
#' One string per segment rather than one for all eight: the CSV is keyed on strain to
#' match the tip labels, and hundreds of strains are a cluster representative in more than
#' one segment. Subtype and host are the same row either way, but the cluster a strain
#' represents is not, so a combined table would let Taxonium match the wrong segment's
#' composition. Per-segment also cuts the payload from ~300 KB to ~65 KB per tree.
build_metadata_csv <- function(full_metadata, tree_names) {
  initial_metadata_csv <-
    full_metadata %>%
    select(-primary_accession) %>%
    group_split(segment) %>%
    map(~ .x %>% select(-segment) %>% format_csv()) %>%
    setNames(full_metadata %>% distinct(segment) %>% arrange(segment) %>% pull(segment) %>% str_c("seg", .))

  # Taxonium parses this CSV with a bare split(",") that strips quote characters and
  # ignores quoting (processNewick.ts), so a comma in any value silently shifts every
  # following column - the pop-up would show the wrong field values with no error.
  metadata_commas <-
    full_metadata %>%
    select(-segment, -primary_accession) %>%
    summarise(across(everything(), ~ sum(str_detect(replace_na(.x, ""), ","), na.rm = TRUE))) %>%
    unlist %>%
    keep(~ .x > 0)

  if(length(metadata_commas) > 0) {
    stop("Metadata column(s) contain commas, which Taxonium's CSV parser cannot handle: ",
         str_c(names(metadata_commas), " (", metadata_commas, " value(s))", collapse = ", "))
  }

  # A tree with no metadata entry would hand Taxonium a NULL and render bare tip labels
  missing_trees <- setdiff(tree_names, names(initial_metadata_csv))
  if(length(missing_trees) > 0) {
    stop("No metadata for tree(s): ", str_c(missing_trees, collapse = ", "))
  }

  log_info("Metadata prepared for trees: {paste(names(initial_metadata_csv), collapse=', ')}")

  initial_metadata_csv
}

#' Curated reference strains for each segment and subtype, pre-sorted by segment and by H
#' and N number
#'
#' The accession column is still called "representative" in the file, for readers that
#' predate the change, but it holds the reference's own accession rather than a cluster
#' representative. Renamed on the way in, so nothing downstream has to remember that.
#'
#' Only segment, strain_display and accession build the Subtype dropdown, so a reference
#' with an unassigned H or N is kept and shown as H? or N? rather than hidden. One is
#' excluded only when it has no sequence data, which would give an empty Position search.
#'
#' Two readers, one table: the Tree tab's Subtype dropdown (through subtype_references()
#' in R/utils.R) and the numbering reference resolution, which looks up the accession the
#' numbering strain has for each segment.
load_subtype_references <- function(conf, discarded) {
  subtype_references_all <-
    read_rds(conf$paths$data$subtype_references) %>%
    dplyr::rename(any_of(c(accession = "representative"))) %>%
    mutate(usable =
             !is.na(strain_display) &
             !is.na(accession) &
             map2_lgl(segment, accession, # reference must have amino acid data for its segment
                      ~ .y %in% names(discarded[[str_c("seg", .x)]])))

  if(any(!subtype_references_all$usable)) {
    log_warn("Excluded {sum(!subtype_references_all$usable)} reference(s) without sequence data: {paste(subtype_references_all$accession[!subtype_references_all$usable], collapse = ', ')}")
  }

  subtype_reference_table <-
    subtype_references_all %>%
    filter(usable) %>%
    select(-usable) %>%
    mutate(strain_display = if_else( # label references whose subtype is only partly assigned
      (is.na(h_subtype) | is.na(n_subtype)) & !str_detect(strain_display, "\\)$"),
      str_c(strain_display, " (", coalesce(h_subtype, "H?"), coalesce(n_subtype, "N?"), ")"),
      strain_display))

  log_info("Subtype references available: {nrow(subtype_reference_table)} of {nrow(subtype_references_all)}")

  subtype_reference_table
}

#' Every residue of every reference strain in every product, with the alignment column it
#' occupies
#'
#' It lets the Tree tab carry a Position search from one reference to another: the column
#' is the same, so the tree colouring is too, and only the number differs. See
#' convert_reference_position() in R/utils.R.
#'
#' Tree tab only. Adaptation Mutations and Batch Screening keep the H3/H5 numbering.
#' Optional: without it a change of reference leaves the Position box as it was.
#'
#' @return the table, or NULL when the file is absent
load_reference_numbering <- function(conf, transposed, subtype_reference_table, products) {
  reference_numbering_path <- conf$paths$data$reference_numbering

  if(is.null(reference_numbering_path) || !file.exists(reference_numbering_path)) {
    log_warn("No reference numbering at {reference_numbering_path %||% '<unset>'} - ",
             "the Tree tab's Position box will keep its number when the reference changes")
    return(NULL)
  }

  reference_numbering <-
    read_rds(reference_numbering_path) %>%
    select(product, accession, position, column, residue)

  # Built for different alignments it would convert to the wrong column, silently
  numbering_overrun <-
    reference_numbering %>%
    dplyr::summarise(max_column = max(column), .by = product) %>%
    dplyr::filter(max_column > map_int(product, ~ length(transposed[[.x]]) %||% 0L))

  if(nrow(numbering_overrun) > 0) {
    stop(reference_numbering_path, " holds columns beyond the alignment for product(s) ",
         str_flatten_comma(numbering_overrun$product), " - it was built against a different alignment")
  }

  # Every reference in the Subtype dropdown has to be convertible, or switching to it
  # would report "no equivalent position" for a strain that has one.
  unnumbered <-
    subtype_reference_table %>%
    distinct(segment, accession) %>%
    tidyr::crossing(product = names(products)) %>%
    dplyr::filter(segment == map_dbl(product, ~ products[[.x]]$segment)) %>%
    dplyr::anti_join(reference_numbering, by = c("product", "accession"))

  if(nrow(unnumbered) > 0) {
    log_warn("{nrow(unnumbered)} reference/product pair(s) are absent from reference_numbering.rds - ",
             "the Position box will not follow a switch to them")
  }

  log_info("Reference numbering: {nrow(reference_numbering)} residues over ",
           "{n_distinct(reference_numbering$accession)} reference(s)")

  reference_numbering
}

#' The canonical protein for each product
#'
#' The sequence a submitted query is aligned against, and the sequence catalogue positions
#' are numbered in. One record per product key, so the primary products and the secondary
#' reading frames are read exactly the same way. HA's second numbering scheme is the
#' seg4_h5 record; it is filed under seg4 after the H3 one, so seg4 holds the two schemes
#' in that order.
#'
#' @return named list of AAStringSet, keyed by product
load_reference_sequences <- function(conf, discarded) {
  reference_records <- readAAStringSet(conf$paths$data$reference_proteins)
  ref_seqs <- split(reference_records, str_remove(names(reference_records), "_h5$"))

  missing_references <- setdiff(names(discarded), names(ref_seqs))
  if(length(missing_references) > 0) {
    log_warn("No reference protein for {str_flatten_comma(missing_references)} in ",
             "{conf$paths$data$reference_proteins} - those products cannot be numbered.")
  }

  ref_seqs
}

#' The Home tab's database summary
#'
#' Three header lines ahead of the per-segment table: last_update, total_GenBank_sequences
#' and total_curated_sequences.
#'
#' @return list of status (the header lines) and status_segments (the per-segment table)
load_status <- function(path) {
  status <- read_tsv(path, col_names = FALSE, n_max = 3)

  status_segments <-
    read_tsv(path, col_names = TRUE, skip = 3,
             col_types = list(segment = col_character(),
                              total = col_integer(),
                              clustered = col_integer())) %>%
    mutate(segment = case_when(
      segment == "Segment_1" ~ "PB2",
      segment == "Segment_2" ~ "PB1",
      segment == "Segment_3" ~ "PA",
      segment == "Segment_4" ~ "HA",
      segment == "Segment_5" ~ "NP",
      segment == "Segment_6" ~ "NA",
      segment == "Segment_7" ~ "M",
      segment == "Segment_8" ~ "NS",
      TRUE ~ segment
    )) %>%
    rename_with(~ str_to_title(.), everything())

  list(status = status, status_segments = status_segments)
}

#' Load everything the app reads
#'
#' @param conf the config
#' @param products the product table defined in global.R: parent segment, label and
#'   alignment file for each product key
#' @return named list; global.R copies each entry into the global environment under the
#'   same name
read_app_data <- function(conf, products) {
  trees <- load_trees(conf$paths$data$trees)

  full_metadata <- load_full_metadata(conf)

  positions <- load_positions(conf, trees$tree_names, products)
  discarded  <- positions$discarded
  transposed <- positions$transposed

  # Only the products whose position data actually loaded can be selected
  available_products <- names(products) %>% keep(~ .x %in% names(discarded))
  product_choices <-
    set_names(names(products), map_chr(products, ~ str_c(.x$segment, "(", .x$label, ")"))) %>%
    keep(~ .x %in% available_products)

  # The products that are a second reading frame of their segment. A primary product is
  # keyed exactly "seg<n>", so anything else is a secondary frame and no list has to be kept
  # in step by hand.
  #
  # They share the parent's tree, representatives and clusters, differing only in the frame
  # the residues are read in. The Tree tab says so on screen, or its phylogeny caption reads
  # as a claim that this is a tree of M2.
  secondary_products <-
    available_products %>% keep(~ .x != str_c("seg", products[[.x]]$segment))

  alignment_counts <- load_alignment_counts(conf, available_products, transposed)

  subtype_reference_table <- load_subtype_references(conf, discarded)

  ref_seqs <- load_reference_sequences(conf, discarded)

  # The reference each segment's adaptation mutation positions are numbered against, each
  # verified against the canonical protein rather than trusted to still be in the clustered
  # set
  numbering_resolution <- resolve_numbering_references(
    sequence_positions = discarded,
    references         = ref_seqs,
    product_keys       = available_products,
    subtype_refs       = subtype_reference_table,
    segments           = function(key) products[[key]]$segment)

  c(trees,
    list(full_metadata            = full_metadata,
         discarded                = discarded,
         transposed               = transposed,
         available_products       = available_products,
         product_choices          = product_choices,
         secondary_products       = secondary_products),
    alignment_counts,
    list(cluster_residues         = load_cluster_residues(conf, transposed),
         adaptation_mutations_all = load_adaptation_mutations(conf, available_products),
         initial_metadata_csv     = build_metadata_csv(full_metadata, trees$tree_names),
         subtype_reference_table  = subtype_reference_table,
         cluster_glue             = subtype_reference_table, # the name the Tree tab and the tests still use
         reference_numbering      = load_reference_numbering(conf, transposed, subtype_reference_table, products),
         ref_seqs                 = ref_seqs,
         numbering_resolution     = numbering_resolution,
         numbering_refs           = map(numbering_resolution, "accession")),
    load_status(conf$paths$data$status))
}

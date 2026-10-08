#' The Adaptation Mutations tool's engine
#'
#' Pure - no Shiny, no notifications - so the arithmetic behind the table and the plot can
#' be exercised without a session, the way R/screening.R is for Batch Screening. The module
#' in R/mod_adaptation.R holds the reactives and the rendering and calls these.
#'
#' Everything a function needs arrives as an argument, defaulting to the global it would
#' otherwise read, so a test can pass a fixture in its place.

#' Catalogued mutations a query carries in one product
#'
#' The per-product half of a query submission, split out because a query is not one
#' product: a whole genomic segment covers both products of segments 2, 3, 7 and 8, M1 and
#' M2 being one piece of RNA, and reporting only the better-scoring one would answer half
#' the question while looking like the whole of it.
#'
#' @param product product key
#' @param scheme "H3" or "H5" for HA, NA for everything else
#' @param residues tibble(position, query_amino_acid) over that product's reference
#' @param catalogues per-product catalogues, defaulting to the global
#' @return list(found, catalogued, assessable, unassessed) - found carries the display
#'   columns plus .product and .scheme, which say which product each row belongs to once
#'   the products are bound into one table
#' @export
adaptation_product_hits <- function(product, scheme, residues,
                                    catalogues = adaptation_mutations_all) {
  catalogue <- catalogues[[product]]

  if(is.null(catalogue)) {
    catalogue <- tibble(position = integer(), mutation = character(), wt = character(),
                        new_amino_acid = character(), which_virus = character(),
                        paper_s = character(), doi = character(), paper_html = character(),
                        how_was_the_mutation_discovered = character(), notes = character())
  }

  # HA is catalogued in both H3 and H5 numbering and the recognizer's hit fixes which one
  # applies. Chosen for every query, not only for one carrying a mismatch, so the numbering
  # cannot be left over from an earlier submission.
  if(product == "seg4") {
    catalogue <-
      if(identical(scheme, "H5")) {
        catalogue %>%
          filter(position_h5 > 0) %>% # the signal peptide positions are negative
          mutate(position = position_h5, mutation = mutation_h5_numbering, .keep = "unused")
      } else {
        catalogue %>% select(-position_h5, -mutation_h5_numbering)
      }
  }

  # Catalogue entries whose position is not in the reference drop out here: entries with no
  # position recorded, and the HA signal peptide positions, which are negative and have no
  # column in the mature HA reference.
  assessable <- catalogue %>% inner_join(residues, by = "position")

  # Present means the query carries the catalogued residue at the catalogued position,
  # whatever the reference has there
  found <-
    assessable %>%
    filter(!is.na(query_amino_acid), query_amino_acid == new_amino_acid) %>%
    select(-c(query_amino_acid, new_amino_acid, wt,
              how_was_the_mutation_discovered, notes)) %>%
    dplyr::rename(Position = position, Mutation = mutation, Subtype = which_virus) %>%
    mutate(.product = product, .scheme = scheme)

  list(found      = found,
       catalogued = nrow(catalogue),
       assessable = nrow(assessable),
       # A position the query does not reach is not an absent mutation - it is one that
       # could not be assessed, and saying so avoids a false negative
       unassessed = sum(is.na(assessable$query_amino_acid)))
}

#' Products of one segment a query reaches, best-scoring first
#'
#' A protein query is almost always one product; a whole genomic segment of 2, 3, 7 or 8 is
#' two, M1 and M2 being one piece of RNA read in two frames.
#'
#' @param identified subject ids that cleared the floor, best first, from
#'   blast_query_segment()
#' @param segment the segment number those products must belong to
#' @param available the products whose position data loaded
#' @return tibble(sseqid, product, scheme), one row per product
#' @export
adaptation_in_play <- function(identified, segment, available = available_products) {
  tibble(sseqid = identified) %>%
    mutate(product = record_product(sseqid),
           scheme  = record_scheme(sseqid)) %>%
    filter(product %in% available,
           map_lgl(product, ~ product_segment(.x) == segment)) %>%
    # identified arrives best first, so this keeps HA's better-scoring numbering scheme
    dplyr::slice_head(n = 1, by = product)
}

#' Residue the query carries at each position of one product's reference
#'
#' So that every catalogued position can be checked - a mismatch table would report an
#' adaptation residue the reference also carries as absent.
#'
#' Two routes to it. The product a protein query is mainly of is aligned against its
#' reference here. Everything else reads its residues off the BLAST alignment the
#' recognizer already has:
#'
#'   a nucleotide query, because pairwiseAlignment() has no way to score a codon against an
#'   amino acid. BLAST did the six-frame translation, and a spliced product's exons arrive
#'   as separate HSPs in their own frames rather than as one reading frame run past its
#'   splice junction;
#'
#'   any further product, because that product is a region of the query rather than the
#'   whole of it. Aligning a 716 aa PA query end to end against the 252 aa PA-X reference is
#'   not the comparison being asked for; the HSP is.
#'
#' The two are not interchangeable, which is why the primary protein path is left alone:
#' over the 748-protein reference set they place gaps differently on 185 of them, mostly in
#' the NA stalk, by up to 64 positions.
#'
#' @param detected as resolve_query_segment() returns it
#' @param product the product to read the query against
#' @param scheme "H5" to number against the mature H5 reference, "H3" or NA otherwise
#' @param sseqid the recognizer record that product was hit as
#' @param primary the product the query is mainly of, which is the one aligned here
#' @param references the canonical proteins, keyed by product
#' @return list(residues, error) - residues is NULL when the alignment failed, and error
#'   then carries the message for the caller to report
#' @export
adaptation_query_residues <- function(detected, product, scheme, sseqid, primary,
                                      references = ref_seqs) {
  index     <- if(identical(scheme, "H5")) 2 else 1 # AAStringSet index into references
  reference <- references[[product]][[index]]

  if(isTRUE(detected$nucleotide) || !identical(product, primary)) {
    return(list(residues = hsp_query_residues(detected$hsps, sseqid, length(reference)),
                error = NULL))
  }

  # TODO display alignment TODO refine alignment parameters gap penalty?
  pairwise <- tryCatch(pairwiseAlignment(detected$seqs, reference,
                                         substitutionMatrix = "BLOSUM80"),
                       error = function(e) e)

  if(inherits(pairwise, "error")) {
    return(list(residues = NULL, error = conditionMessage(pairwise)))
  }

  list(residues = aligned_query_residues(pairwise), error = NULL)
}

#' An empty frequency table, in the shape the plot expects
#' @export
ADAPTATION_NO_FREQUENCIES <- data.frame(amino_acid = character(),
                                        host_group = character(),
                                        n = integer())

#' Amino acid counts by host at one alignment column
#'
#' @param product product key
#' @param column alignment column, as gapless_to_consensus() resolves it
#' @param reference the accession positions are numbered against, dropped from the cluster
#'   tally so the plot reports the sequences being numbered rather than the yardstick
#' @param grouping "host_group" or "host_order"
#' @param counting_sequences TRUE to count every sequence, FALSE for cluster representatives
#' @param counts the per-column tally, defaulting to the global
#' @param columns the per-column residues of the representatives, defaulting to transposed
#' @param metadata the sequence metadata, for the representatives' hosts
#' @return tibble(amino_acid, host_category, n), empty where there is nothing to count
#' @export
adaptation_column_frequencies <- function(product, column, reference, grouping,
                                          counting_sequences = FALSE,
                                          counts   = alignment_counts,
                                          columns  = transposed,
                                          metadata = full_metadata) {
  # Every sequence in the database, rather than one vote per cluster. The counts were
  # tallied per alignment column by alignment_column_counts.R, which reads the same
  # coordinate frame transposed uses, so the column addresses it unchanged. The residues are
  # the sequences' own - nothing is inferred from a representative.
  if(counting_sequences) {
    return(
      counts %>%
        filter(product == !!product, column == !!column, level == grouping) %>%
        # already excluded when the table was built; repeated so the two modes cannot
        # drift apart if that ever changes
        filter(!is.na(host_category),
               !host_category %in% EXCLUDED_HOSTS) %>%
        mutate(amino_acid = factor(amino_acid, levels = get_aa_levels())) %>%
        select(amino_acid, host_category, n)
    )
  }

  amino_acids <- columns %>% pluck(product, column) %>% unlist

  if(is.null(amino_acids) || length(amino_acids) <= 1) {
    return(ADAPTATION_NO_FREQUENCIES)
  }

  # Drop the sequence whose numbering defines the position, so the plot reports the
  # sequences being numbered rather than counting the yardstick among them. By accession,
  # not by position: the reference is not first in the alignment for every product. For HA
  # the accession changes with the numbering scheme, which is why the caller resolves it
  # rather than passing a fixed per-segment value.
  amino_acids[names(amino_acids) != reference] %>%
    enframe(name = "primary_accession", value = "amino_acid") %>%
    inner_join(metadata, by = join_by(primary_accession)) %>%
    mutate(host_category = .data[[grouping]]) %>% # host group or host order
    filter(!is.na(host_category)) %>% # remove sequences with no host at this level
    filter(!host_category %in% EXCLUDED_HOSTS) %>% # the same list the other mode uses
    select(primary_accession, amino_acid, host_category) %>%
    mutate(amino_acid = factor(amino_acid, levels = get_aa_levels())) %>%
    dplyr::count(amino_acid, host_category)
}

#' How much of the data the bars actually rest on, in the unit being counted
#'
#' Most positions are covered by about 90% of the cluster representatives, so the shortfall
#' is unremarkable and the same at every position. A handful are not: the mature HA
#' N-terminus and the NA stalk are gapped in most of the alignment. The plot fills each bar
#' to full height whatever it rests on, so without this a position covered by a fifth of the
#' data reads exactly like one covered by all of it.
#'
#' @inheritParams adaptation_column_frequencies
#' @param totals sequences that entered the tally per product and level
#' @param sequences the per-sequence positions, for the count of representatives
#' @return list(covered, total, unit, reference_counted), or NULL where there is nothing to
#'   report against
#' @export
adaptation_coverage <- function(product, column, reference, grouping,
                                counting_sequences = FALSE,
                                counts    = alignment_counts,
                                totals    = alignment_counts_totals,
                                columns   = transposed,
                                sequences = discarded) {
  if(is.null(product) || is.na(column)) {
    return(NULL)
  }

  if(counting_sequences) {
    if(is.null(counts) || is.null(totals)) {
      return(NULL)
    }

    covered <-
      counts %>%
      filter(product == !!product, column == !!column, level == grouping) %>%
      pull(n) %>%
      sum()

    total <-
      totals %>%
      filter(product == !!product, level == grouping) %>%
      pull(sequences)

    if(length(total) != 1 || total == 0) {
      return(NULL)
    }

    # Both figures come from the same table under the same host filters, so the bars
    # already sum to the numerator and there is nothing to reconcile.
    return(list(covered = covered, total = total, unit = "sequences",
                reference_counted = FALSE))
  }

  # The representatives carrying a residue here, against every representative of the
  # product. Counted before the host filters, which is what makes it a statement about the
  # alignment rather than about the taxonomy.
  residues <- columns %>% pluck(product, column) %>% unlist()
  total    <- length(sequences[[product]])

  if(is.null(total) || total == 0) {
    return(NULL)
  }

  # The frequencies drop the numbering reference, but only where it has a residue here at
  # all. Whether it does decides the arithmetic in the caption, so read it rather than
  # assuming it is always one of the covered.
  list(covered = length(residues), total = total, unit = "cluster representatives",
       reference_counted = reference %in% names(residues))
}

#' The coverage line of the plot caption, given the total the bars actually draw
#'
#' @param counted as adaptation_coverage() returns it, or NULL
#' @param plotted what the bars sum to
#' @return the caption line, or NULL when there is no coverage to report
#' @export
adaptation_coverage_caption <- function(counted, plotted) {
  if(is.null(counted)) {
    return(NULL)
  }

  # In the cluster scope the numerator counts the alignment while the bars count the
  # taxonomy, so the two differ by the numbering reference and the representatives with no
  # usable host. Reconcile them rather than leaving the reader to discover that the bars do
  # not sum to the number above them. In the sequence scope both figures come from
  # alignment_column_counts.rds under the same filters, so they agree and nothing is added.
  reconciliation <- if(plotted < counted$covered) {
    hostless <- counted$covered - plotted - as.integer(counted$reference_counted)

    str_c(" ", scales::comma(plotted), " are plotted, setting aside ",
          str_flatten_comma(c(
            if(counted$reference_counted) "the numbering reference",
            if(hostless > 0) str_c(scales::comma(hostless),
                                   " with no host at this level")),
            last = " and "),
          ".")
  } else {
    ""
  }

  str_wrap(str_c(scales::comma(counted$covered), " of ",
                 scales::comma(counted$total), " ", counted$unit,
                 " carry an amino acid at this position (",
                 round(100 * counted$covered / counted$total), "%).",
                 reconciliation), width = 95)
}

#' Host categories outside the top 90%, which the plot pools into "Other"
#'
#' Host order has a long tail of categories holding one or two sequences, each of which
#' still draws a full height bar. Only host order is pooled; host group has five categories
#' at most.
#'
#' @param frequencies as adaptation_column_frequencies() returns them
#' @param grouping "host_group" or "host_order"
#' @return the pooled category names, empty when nothing is pooled
#' @export
adaptation_pooled_hosts <- function(frequencies, grouping) {
  if(!identical(grouping, "host_order") || nrow(frequencies) == 0) {
    return(character())
  }

  ranked <-
    frequencies %>%
    dplyr::count(host_category, wt = n, name = "total") %>%
    arrange(desc(total)) %>%
    mutate(cumulative_share = cumsum(total) / sum(total))

  # keep categories up to and including the one that crosses 90%
  kept <- ranked %>% filter(lag(cumulative_share, default = 0) < 0.9) %>% pull(host_category)

  setdiff(ranked$host_category, kept)
}

#' The bars, as one row per host category and amino acid
#'
#' @param frequencies as adaptation_column_frequencies() returns them
#' @param pooled the categories to collapse into "Other", from adaptation_pooled_hosts()
#' @return tibble(host_category, amino_acid, n, sum_n, percentage, host_label)
#' @export
adaptation_plot_data <- function(frequencies, pooled = character()) {
  if(length(pooled) > 0) {
    frequencies <- frequencies %>%
      mutate(host_category = if_else(host_category %in% pooled, "Other", host_category))
  }

  plot_data <-
    frequencies %>%
    group_by(host_category, amino_acid) %>%
    summarise(n = sum(n)) %>% # total per amino acid in host category
    mutate(sum_n = sum(n)) %>% # total per host category
    mutate(percentage = n / sum_n) %>%
    ungroup %>%
    # the bars are proportions, so name the denominator: host order in particular produces
    # categories holding only a handful of representatives
    mutate(host_label = str_c(host_category, "\n(n = ", scales::comma(sum_n), ")")) %>%
    mutate(host_label = fct_reorder(host_label, sum_n, .desc = TRUE))

  if(length(pooled) > 0) { # the pooled bar reads better last than ranked by size
    plot_data <- plot_data %>%
      mutate(host_label = fct_relevel(host_label,
                                      str_subset(levels(host_label), "^Other\\n"),
                                      after = Inf))
  }

  plot_data
}

#' What was pooled, so a collapsed bar is not an unexplained category
#'
#' @param pooled the category names, from adaptation_pooled_hosts()
#' @return the caption line, or NULL when nothing was pooled
#' @export
adaptation_pooled_caption <- function(pooled) {
  if(length(pooled) == 0) {
    return(NULL)
  }

  str_wrap(str_c("Other pools ", length(pooled), " host order",
                 if(length(pooled) == 1) "" else "s", " outside the top 90% of sequences: ",
                 str_flatten_comma(sort(pooled)), "."), width = 95)
}

#' The amino acid frequency chart
#'
#' @param plot_data as adaptation_plot_data() returns it
#' @param subtitle the product, position and numbering reference on screen
#' @param caption the coverage, pooling and taxonomy lines, or NULL
#' @param x_label the host axis label
#' @param unit what a bar counts, for the tooltip - "sequences" or "clusters"
#' @return a ggplot, for girafe() to render
#' @export
adaptation_plot <- function(plot_data, subtitle, caption, x_label, unit) {
  # host order can produce a dozen or more categories, which will not fit horizontally at
  # the axis text size that suits the three host groups
  crowded_axis <- n_distinct(plot_data$host_label) > 4
  label_angle  <- if(crowded_axis) 45 else 0
  label_size   <- if(crowded_axis) rel(0.7) else rel(1) # relative to axis.text below

  plot_data %>%
    ggplot(aes(fill = amino_acid,
               y = n,
               x = host_label,
               tooltip = str_c(amino_acid, "<br>",
                               scales::comma(n), "/", scales::comma(sum_n), " ",
                               unit, "<br>",
                               scales::percent(percentage, accuracy = 0.1)),
               data_id = interaction(host_label, amino_acid)
               )) +
    geom_bar_interactive(position="fill", stat="identity") +
    labs(title = "Frequency of amino acids at consensus position",
         subtitle = subtitle,
         caption = caption,
         x = x_label,
         y = "Frequency",
         fill = "Amino Acid") +
    scale_fill_manual(values = AA_PALETTE, na.value = "grey50") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal(base_family='Open Sans') +
    theme(axis.text=element_text(size=rel(1.25),
                                 colour = "black"),
          axis.text.x = element_text(angle = label_angle,
                                     hjust = if(label_angle == 0) 0.5 else 1,
                                     size = label_size),
          axis.title=element_text(size=rel(1.5)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title =element_text(size=rel(1.75))
    )
}

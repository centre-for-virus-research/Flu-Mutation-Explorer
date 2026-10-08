#' Adaptation Mutations tool - catalogue table and amino acid frequency plot
#'
#' Self-contained: everything it needs is in `adaptationValues` and its own reactives. It
#' reads no tree state; requests from the other tools arrive through the `show` reactive
#' the parent passes in.
#'
#' The UI comes in two pieces because the layout interleaves them: a sidebar panel in the
#' left column and a tab body in the right. The parent decides where each goes and owns
#' the `input.tabselected` conditions that show or hide them, since the tab id lives in
#' the parent's namespace.

#' Rows the catalogue table shows at once
#'
#' Named rather than left as a literal in the table options because the hand-over from
#' Batch Screening has to work out which page a row is on, and a page length the two
#' disagreed about would land the user next to the site rather than on it.
ADAPTATION_PAGE_LENGTH <- 10L

#' Sidebar controls for the Adaptation Mutations tab
adaptation_sidebar_ui <- function(id) {
  ns <- NS(id)

                               wellPanel(
                               radioButtons(ns("view_query"), 
                                             label = tooltip(
                                               trigger = "Adaptations",
                                               "View all adaptation mutations for a given segment or query a protein or nucleotide sequence and auto-detect segment"
                                             ),
                                            choices = list("View" = "view", "Query" = "query"), 
                                            selected = "view",
                                            inline = FALSE),
                               
                               # conditionally display segment selection input
                               conditionalPanel(
                                 condition = paste0('input[\'', ns('view_query'), "\'] == \'view\'"),
                                 selectInput(ns("protein"),
                                             label = tooltip(
                                               trigger = "Segment",
                                               "Influenza genome segment and protein product. Segments 7 and 8 each encode a second product - M2 and NEP - whose catalogued positions are numbered against that protein, not against M1 or NS1."
                                             ),
                                             choices = product_choices,
                                             selected = "seg1",
                                             selectize = FALSE
                                 )
                               ), 
                               
                               conditionalPanel(
                                 condition = paste0('input[\'', ns('view_query'), "\'] == \'query\'"),
                                 textAreaInput(ns("query_adaptation"), 
                                               tooltip(
                                                 trigger = "Query sequence",
                                                 "Protein or nucleotide sequence to match against the reference. Nucleotide input is translated by BLAST."
                                               ),
                                               rows = 10, 
                                               width = NULL, 
                                               placeholder = NULL),
                                 actionButton(ns("submit_adaptation"), "Submit"),
                                 actionButton(ns("clear_adaptation"), "Clear")
                               ), # end conditionalPanel for query sequence

                               radioButtons(ns("host_grouping"),
                                            label = tooltip(
                                              trigger = "Group hosts by",
                                              "Taxonomic level for the host axis of the amino acid frequency chart"
                                            ),
                                            choices = list("Host group" = "host_group",
                                                           "Host order" = "host_order"),
                                            selected = "host_group",
                                            inline = FALSE),

                               # Offered only where the counts were built, so a deployment
                               # without that pipeline stage shows one plot rather than a
                               # control that cannot do anything
                               if(!is.null(alignment_counts)) {
                                 radioButtons(ns("sequence_scope"),
                                              label = tooltip(
                                                trigger = "Count",
                                                "Cluster representatives are the sequences on the tree, one per cluster. All sequences counts every sequence in the database at this position."
                                              ),
                                              choices = list("Cluster representatives" = "clusters",
                                                             "All sequences" = "sequences"),
                                              selected = "clusters",
                                              inline = FALSE)
                               }
                             ) # end wellPanel for adaptation mutations
}

#' Body of the Adaptation Mutations tab
adaptation_tab_ui <- function(id) {
  ns <- NS(id)

  tagList(
                            fluidRow(
                              column(width = 6, girafeOutput(ns("adaptation_plot")),
                                     markdown(
                                       "
                              Amino acid positions for HA mutations are based on mature peptide numbering.
                              Amino acid positions for mutations in all other segments are based on full-length protein numbering.
                              **These mammalian adaptation mutations are largely derived from _in vitro_ studies performed in specific viral strain backgrounds. Consequently, their adaptive significance may be context-dependent and may not apply to the submitted query protein sequence.**
                                       "
                                     )
                              ),
                              column(width = 6, dataTableOutput(ns("adaptationTable")))
                            )
  )
}

#' Server for the Adaptation Mutations tool
#'
#' @return a reactive carrying a confirmed request to show a position on the tree -
#'   list(product, reference, position) - or NULL. The parent hands it to tree_server();
#'   this tool cannot reach into the other, it can only ask.
adaptation_server <- function(id, show = reactive(NULL), reset = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
        adaptationValues <- reactiveValues(numbering = "H3") # adaptation mutations

    # Bumped by this tab's own Clear button and read by nobody here - it is returned to
    # the parent, merged with the other two tools' own copies, and handed back to every
    # tool as `reset`, so clicking Clear on any one tab resets all three.
    clearedVal <- reactiveVal(0L)

    # Sliced to zero rows rather than set to NULL outright, so the table stays whatever
    # shape displayed_mutations()/first_mappable_row()/etc. already expect - the same
    # empty-table path a query with no catalogued hits already exercises. Left alone if
    # it is already NULL, which is itself already "empty" everywhere that reads it.
    clear_mutations_table <- function() {
      if(!is.null(adaptationValues$mutations)) {
        adaptationValues$mutations <- adaptationValues$mutations[0, ]
      }
    }

    ### Incoming requests to show a position's plot ###

    # Batch Screening can ask for a position's bar plot the way it can ask the Tree tab
    # for a tip. Same shape of hand-over: one reactive, written by the tool that asks and
    # read by the tool that acts, so neither reaches into the other's state.
    #
    # Two steps, because the plot is driven by a selected row of the table rather than by
    # a position directly. Switching the product rebuilds that table, so the row cannot be
    # selected until it exists - the request is parked here and acted on below.
    pendingShow <- reactiveVal(NULL)

    # The row the table is to open on, once a request has been resolved to one.
    #
    # Read by the render rather than pushed through a dataTableProxy, which cannot win the
    # race it is in: a request that changes the product invalidates the table in the same
    # flush that resolves the request, and Shiny runs observers ahead of outputs, so the
    # proxy's selection lands on the table that is about to be replaced and the redraw wipes
    # the highlight. Handing the row to the render means the table is built with it
    # selected, rebuilt or not.
    #
    # The serial rides along so that clicking the same site twice is still a new value, and
    # still re-opens the table there after the user has paged away.
    openAt <- reactiveVal(NULL)

    observeEvent(show(), {
      request <- show()
      req(request$product, request$position)

      updateRadioButtons(session = session, inputId = "view_query", selected = "view")
      updateSelectInput(session = session, inputId = "protein", selected = request$product)

      pendingShow(request)
    })

    # observe() rather than observeEvent() on the table, because it has to depend on both
    # the parked request and the table: when the request names the product already on
    # screen, nothing rebuilds, and an observer waiting on the table alone would never
    # fire - the click would do nothing at all.
    observe({
      request <- pendingShow()
      req(request)

      mutations <- adaptationValues$mutations
      req(!is.null(mutations), nrow(mutations) > 0)

      # only act once the table is the requested product's, not the one it replaced
      req(identical(adaptationValues$product, request$product))

      row <- which(mutations$Position == request$position)[1]

      if(is.na(row)) {
        log_warn("Position {request$position} is not in the {request$product} table - not selected")
      } else {
        openAt(list(product = request$product, row = row, serial = request$serial))
        log_info("Showing {request$product} position {request$position} - row {row}")
      }

      pendingShow(NULL)
    })

    # Display all adaptation mutations for selected segment
    observe({
      req(input$view_query == "view") # require view_query input to be view
      
      adaptationValues$numbering <- "H3" # reset numbering
      adaptation_mutations <- adaptation_mutations_all[[input$protein]]

      adaptation_mutations %<>%
      select(-c(new_amino_acid,
                wt, how_was_the_mutation_discovered, notes, paper_s, doi)) %>% # columns not displayed
        dplyr::rename(Position = position,
                      Mutation = mutation,
                      Subtype = which_virus,
                      Paper = paper_html) # the links are built once in global.R
      
      if(input$protein == "seg4") {
        adaptation_mutations %<>%
          dplyr::rename(`Position(H5)` = position_h5,
                        `Mutation(H5)` = mutation_h5_numbering) # rename columns for display
      }

      # View mode is one product by construction - the dropdown picks it - but every row
      # says which one anyway, so that the table has the same shape as a query's and the
      # reactives below need not ask which mode produced it.
      adaptation_mutations %<>% mutate(.product = input$protein, .scheme = NA_character_)

      adaptationValues$mutations <- adaptation_mutations # store adaptation mutations for selected segment
      adaptationValues$products <- input$protein # every product the table covers
      adaptationValues$product <- input$protein # the key into transposed, discarded and numbering_refs
      adaptationValues$segment <- product_segment(input$protein) # the segment it belongs to, for display
    })

    # Selecting Query leaves whatever was last on screen - the View table, or an earlier
    # query's results - showing until a new query resolves, which reads as results for a
    # sequence that was never submitted. Emptied as soon as Query mode is chosen instead,
    # not only once Submit is pressed.
    observeEvent(input$view_query, {
      req(identical(input$view_query, "query"))
      clear_mutations_table()
    })

  ### Adaptation mutations ###

  # A proxy so a row selection left over from before the table is about to be redrawn
  # empty is not still showing for the one flush in between.
  adaptationTableProxy <- dataTableProxy("adaptationTable")

  # Puts this tab back to what it shows before anything has been submitted or selected.
  # Run only from `reset` below, so every tool resets itself the same way regardless of
  # which tab's Clear button caused it.
  #
  # clear_mutations_table() only where Query mode is what is showing: View mode's table
  # is never empty by default (it is the segment's own catalogue), so calling it while
  # View is on screen - reachable now that another tab's Clear can trigger this too -
  # would empty a table that was never a "result" to begin with.
  reset_tab <- function() {
    updateTextAreaInput(session, "query_adaptation", value = "")

    updateRadioButtons(session = session, inputId = "host_grouping", selected = "host_group")
    if(!is.null(alignment_counts)) {
      updateRadioButtons(session = session, inputId = "sequence_scope", selected = "clusters")
    }

    DT::selectRows(adaptationTableProxy, NULL)

    if(identical(input$view_query, "query")) {
      clear_mutations_table()
    }
  }

  # Only tells the parent - this tab resets when `reset` comes back, like the other two.
  observeEvent(input$clear_adaptation, clearedVal(clearedVal() + 1L))

  # Any tab's Clear button, this one's included, arrived through the parent - see
  # clearedVal above.
  observeEvent(reset(), {
    reset_tab()
  }, ignoreInit = TRUE)

  # User submits AA sequence for adaptation mutations
  # Detect segment
  # Detect H3 or H5 for HA
  # Perform pairwise alignment with reference sequence
  observeEvent(input$submit_adaptation, {
    req(input$view_query == "query") # require view_query option to be query
    req(input$query_adaptation) # require adaptation query sequence to be available

    # A resubmission otherwise leaves the previous query's results on screen for the
    # whole BLAST/alignment run - selecting Query mode already empties the table, but
    # this covers running a second query without leaving that mode in between.
    clear_mutations_table()

    # Reading the query and detecting its segment is shared with the Tree tool - see
    # R/query_segment.R. NULL means the query was unusable and the user has already been
    # told why.
    detected <- resolve_query_segment(input$query_adaptation)
    req(detected)

    segment <- detected$segment
    protein <- detected$product

    # The recognizer names the product, not only the segment, because its database is the
    # reference proteins themselves - one record per product. So no second scoring pass is
    # needed here to tell M2 from M1, whose length M2's positions all fall inside.
    if(!protein %in% available_products) {
      showNotification(str_c("No position data for ", product_label(protein),
                             " - it cannot be numbered"), type = "error")
      req(FALSE) # exit
    }

    log_info("Query resolved to {protein} ({product_label(protein)}) on segment {segment}")

    # Named to the user only where the segment carries more than one product; for the rest
    # "Detected segment 5" has already said it.
    if(sum(map_lgl(available_products, ~ product_segment(.x) == segment)) > 1) {
      showNotification(str_c("Detected ", product_label(protein)), type = "message")
    }

    updateSelectInput(session = session,
                      inputId = "protein",
                      selected = protein) # update product selection

    # Which of the two HA numbering schemes applies. seg4.fasta carries the mature H3
    # reference first and the mature H5 second, and the recognizer's database holds both as
    # seg4 and seg4_h5, so the hit that named the product named the scheme with it - no
    # second BLAST against Ref_H3_H5/ is needed to ask again.
    if(segment == 4) {
      showNotification(str_c("Using ", detected$scheme, " numbering"), type = "message")
    }

    adaptationValues$numbering <- if(identical(detected$scheme, "H5")) "H5" else "H3"

    in_play <- adaptation_in_play(detected$identified, segment)

    if(nrow(in_play) == 0) {
      showNotification(str_c("No product of segment ", segment, " could be numbered"),
                       type = "error")
      req(FALSE) # exit
    }

    # The residues the query carries against each product's reference, and the catalogue
    # rows they match - see adaptation_query_residues() for which of the two alignments
    # answers for which product.
    per_product <- pmap(in_play, function(sseqid, product, scheme) {
      read <- adaptation_query_residues(detected, product, scheme, sseqid, primary = protein)

      if(!is.null(read$error)) {
        log_warn("Pairwise alignment failed: {read$error}")
        showNotification(str_c("Pairwise alignment error: ", read$error,
                               ". Ensure the sequence contains only valid amino acid characters."),
                         type = "error")
      }

      residues <- read$residues
      if(is.null(residues)) return(NULL)

      log_info("Query covers {sum(!is.na(residues$query_amino_acid))} of ",
               "{nrow(residues)} positions of the {product} reference")

      adaptation_product_hits(product, scheme, residues)
    })

    names(per_product) <- in_play$product
    per_product <- compact(per_product)
    req(length(per_product) > 0)

    # A product with nothing catalogued can never contribute a row, so listing it would put
    # a Product column on the table to tell apart something that is not there. PB1-F2 and
    # PA-X are that case - their positions are searchable on the Tree tab, but neither has a
    # catalogued mutation yet - so a submitted segment 2 or 3 reads exactly as it did.
    #
    # Unless that is all there is: a query that is itself a PB1-F2 keeps its empty table,
    # which the plot explains, rather than being left with no product at all.
    catalogued <- keep(per_product, ~ .x$catalogued > 0)
    if(length(catalogued) > 0) per_product <- catalogued

    for(product in names(per_product)) {
      result <- per_product[[product]]
      if(result$assessable < result$catalogued) {
        log_info("{result$catalogued - result$assessable} catalogue entr(ies) for {product} ",
                 "have no position in the reference and were not checked")
      }
    }

    # A position the query does not reach is not an absent mutation - it is one
    # that could not be assessed, and saying so avoids a false negative
    unassessed <- sum(map_int(per_product, "unassessed"))

    if(unassessed > 0) {
      showNotification(str_c(unassessed, " catalogued position", if(unassessed == 1) "" else "s",
                             if(length(per_product) > 1) {
                               str_c(" across ", str_flatten_comma(map_chr(names(per_product), product_label)))
                             } else "",
                             " lie outside the region the query covers and could not be assessed"),
                       type = "message")
    }

    found <- map(per_product, "found") %>% purrr::list_rbind()

    if(nrow(found) == 0) {
      showNotification("No adaptation mutations found", type = "warning")
      found %<>%
        dplyr::rename(Paper = paper_html) %>% # rename column for display
        select(Mutation, Position, Subtype, Paper, .product, .scheme) # reorder/remove columns
    } else {
      found %<>%
        select(-c(paper_s, doi)) %>%
        dplyr::rename(Paper = paper_html) %>% # the links are built once in global.R
        relocate(Paper, .after = dplyr::last_col()) # where .keep = "unused" left it

      showNotification(str_c("Found ", nrow(found), " adaptation mutation",
                             if(nrow(found) == 1) "" else "s", " in the query",
                             if(length(per_product) > 1) {
                               str_c(" (", str_flatten_comma(
                                 found %>% dplyr::count(.product) %>%
                                   dplyr::mutate(text = str_c(map_chr(.product, product_label), " ", n)) %>%
                                   pull(text)), ")")
                             } else ""),
                       type = "warning")
    }

    # Product is named in the table only when there is more than one to tell apart. The
    # column is placed first: it is what a row is of, and every other column is read in
    # its light - position 31 means one site in M1 and another in M2.
    if(length(per_product) > 1) {
      found %<>%
        mutate(Product = map_chr(.product, product_label)) %>%
        relocate(Product)
    }

    adaptationValues$mutations <- found # the adaptation mutations found in the query
    adaptationValues$products <- names(per_product) # every product the table covers
    # The product the table is built from, and what selected_product() falls back to before
    # a row is chosen. The one the query is mainly of, unless that is the one with nothing
    # catalogued and another survived.
    adaptationValues$product <-
      if(protein %in% names(per_product)) protein else names(per_product)[1]
    adaptationValues$segment <- segment # the segment it belongs to, for display

  }) # end observeEvent for adaptation mutations query

  # Catalogue position of the adaptation mutation selected in the table
  selected_position <- reactive({
    req(adaptationValues$mutations %>% nrow > 0,
        input$adaptationTable_rows_selected) # require adaptation hits and user to select a row

    adaptationValues$mutations %>%
      dplyr::slice(input$adaptationTable_rows_selected) %>%
      pull(Position)
  })

  # Product the selected row belongs to, which is not the table's product any more: a
  # query over a whole genomic segment lists both of that segment's products, and position
  # 31 is one site in M1 and a different one in M2. Everything below that answers "about
  # what" - the numbering reference, the alignment column, the bar plot, the jump to the
  # tree - keys off this rather than off adaptationValues$product.
  #
  # adaptationValues$product stays the product the query is mainly of, and the table is
  # built from it. Were the table to depend on the selection instead, selecting a row would
  # rebuild the table and drop the selection that caused it.
  selected_product <- reactive({
    mutations <- adaptationValues$mutations
    row       <- input$adaptationTable_rows_selected

    if(is.null(mutations) || nrow(mutations) == 0 || !".product" %in% names(mutations) ||
       is.null(row) || length(row) != 1 || row > nrow(mutations)) {
      return(adaptationValues$product)
    }

    mutations$.product[[row]]
  })

  # Reference sequence whose numbering the selected position is expressed in. Resolved and
  # verified against the canonical protein at startup by product_reference()'s neighbours
  #
  # HA positions are reported in mature H3 or mature H5 numbering; View always uses H3.
  # in R/numbering_reference.R, rather than read straight from config, so a reference
  # dropping out of the clustered set cannot silently empty the module.
  #
  # HA positions are reported in mature H3 or mature H5 numbering; View always uses H3.
  numbering_reference <- reactive({
    req(selected_product())

    product_reference(selected_product(),
                      if(input$view_query == "view") "H3" else adaptationValues$numbering)
  })

  # Taxonomic level for the host axis of the frequency chart. Host order separates, for
  # example, Anseriformes from Galliformes and Primates from Artiodactyla, which host
  # group collapses into Birds / Human / Other Mammals.
  host_grouping <- reactive({
    if(is.null(input$host_grouping)) "host_group" else input$host_grouping
  })

  host_grouping_label <- reactive({
    if(host_grouping() == "host_order") "Host order" else "Host group"
  })

  # Segments 7 and 8 each offer two products, so a caption naming the segment alone
  # would not say which protein the positions are numbered against
  product_title <- reactive({
    req(selected_product())
    str_c("Segment ", adaptationValues$segment,
          " (", product_label(selected_product()), ")")
  })

  product_description <- reactive(str_to_lower(str_sub(product_title(), 1, 7)) %>%
                                    str_c(str_sub(product_title(), 8)))

  # The table's own caption names every product in it, where the plot below names the one
  # the selected row belongs to. A query over a whole genomic segment lists both.
  table_description <- reactive({
    products <- adaptationValues$products %||% adaptationValues$product
    req(length(products) > 0)

    str_c("segment ", adaptationValues$segment,
          " (", str_flatten_comma(map_chr(products, product_label)), ")")
  })

  # What each bar counts. "clusters" is one vote per cluster representative, the sequences
  # the tree is built from; "sequences" counts every sequence in the database at the
  # position. They differ sharply - 25 human HA representatives stand for two thirds of the
  # HA sequences - so the plot names the one on screen rather than leaving it implicit.
  sequence_scope <- reactive({
    if(is.null(alignment_counts) || is.null(input$sequence_scope)) "clusters" else input$sequence_scope
  })

  counting_sequences <- reactive(sequence_scope() == "sequences")

  scope_label <- reactive({
    if(counting_sequences()) "all sequences" else "cluster representatives"
  })

  scope_unit <- reactive({
    if(counting_sequences()) "sequences" else "clusters"
  })

  # First row whose position can be located in the alignment, so each segment opens on a
  # bar plot rather than on an explanation. HA row 1 is the signal peptide position -1,
  # which has no column in the mature HA alignment.
  first_mappable_row <- reactive({
    mutations <- adaptationValues$mutations

    if(is.null(adaptationValues$product) || is.null(mutations) || nrow(mutations) == 0) {
      return(1L)
    }

    # Each row is resolved against its own product's reference. The table can hold two
    # products at once, and it cannot go through numbering_reference() in any case: that
    # follows the selection, and this is what the selection starts at.
    mappable <- detect_index(seq_len(nrow(mutations)), function(row) {
      position <- mutations$Position[[row]]
      product  <- mutations$.product[[row]]
      reference <- product_reference(product, adaptationValues$numbering)

      !is.na(position) &&
        position > 0 &&
        !is.na(reference) &&
        !is.na(gapless_to_consensus(product, reference, position))
    })

    if(mappable == 0L) 1L else mappable # fall back to the first row if none can be mapped
  })

  # The row a hand-over from Batch Screening asked for, when there is one and it belongs
  # to the product being shown. NULL otherwise: the request outlives the render that
  # served it, so leaving one product and coming back must not select another product's
  # site, and a catalogue that has since shrunk must not select a row past its end.
  requested_row <- reactive({
    opening   <- openAt()
    mutations <- adaptationValues$mutations

    if(is.null(opening) || is.null(mutations) || nrow(mutations) == 0) return(NULL)
    if(!identical(opening$product, adaptationValues$product)) return(NULL)
    if(opening$row > nrow(mutations)) return(NULL)

    opening$row
  })

  # The row the table opens selected on, and the page that row is on. The catalogues run
  # to 268 rows for PB2 and 546 for HA, so a handed-over site is usually tens of pages in
  # - PB2 627 is row 212, page 22 - and selecting it without paging to it left the user on
  # page 1 with nothing on screen saying which site the plot was of.
  opening_row <- reactive(requested_row() %||% first_mappable_row())

  # 0-based index of the first row of that page, which is what DataTables starts on. Only
  # a hand-over moves it: opening on the first mappable row has always shown page 1.
  opening_display_start <- reactive({
    row <- requested_row()
    if(is.null(row)) 0L else ((row - 1L) %/% ADAPTATION_PAGE_LENGTH) * ADAPTATION_PAGE_LENGTH
  })

  # Explain an empty frequency result rather than rendering an empty plot
  frequency_message <- reactive({
    position <- selected_position()
    reference <- numbering_reference()

    if(is.na(position)) { # catalogue entries without a position
      "This mutation has no position recorded in the catalogue, so its amino acid frequencies cannot be shown."
    } else if(selected_product() == "seg4" && position < 0) {
      str_c("Position ", position, " lies in the HA signal peptide, which is not represented ",
            "in the mature HA alignment, so amino acid frequencies cannot be shown for it.")
    } else if(is.na(reference)) {
      str_c("No verified numbering reference is available for ", product_description(),
            " in this database version, so positions in it cannot be located in the alignment.")
    } else if(counting_sequences() &&
              !is.na(gapless_to_consensus(selected_product(), reference, position))) {
      # the column resolves, so this is an empty column rather than a numbering failure
      str_c("Position ", position, " maps to the alignment, but no sequence with a placed ",
            "host carries a residue there, so there is nothing to count.")
    } else {
      str_c("Position ", position, " could not be located in the alignment using reference ",
            reference, ", so no amino acid frequencies are available for it.")
    }
  })

  # Can this position be shown on the tree?
  #
  # Both tabs resolve a position through gapless_to_consensus(), but they choose the
  # reference differently: here it is the product's numbering reference, and on the Tree
  # tab it is whatever the Subtype dropdown holds. The jump is only correct when that
  # dropdown can be set to this reference.
  #
  # It always can today, both coming from subtype_references.rds, but the check stays
  # because the failure would not show: updateSelectInput() given a selection outside its
  # choices falls back silently to the first one, so the tree would be coloured by this
  # position in a different sequence, with nothing on screen to say so.
  jumpable <- function(position, reference, product) {
    !is.null(position) && !is.na(position) && position > 0 &&
      !is.null(reference) && !is.na(reference) &&
      reference %in% subtype_references(product) &&
      !is.na(gapless_to_consensus(product, reference, position))
  }

  # Which numbering scheme each position column of the current table is expressed in.
  #
  # Everything but HA has one, keyed by the product itself. HA in View mode has two side
  # by side - the catalogue records both - and they are numbered against different
  # references, so each column links against its own. In Query mode the BLAST hit has
  # already fixed the scheme and only one column survives.
  #
  # The names are the data frame's own column names; the values are keys into
  # numbering_refs, which holds product keys and the two HA schemes alike. Keyed on the
  # table's product, so the value is a default: where a query has listed two products,
  # displayed_mutations() takes each row's scheme from its own product instead.
  position_schemes <- reactive({
    product   <- adaptationValues$product
    mutations <- adaptationValues$mutations

    if(is.null(product) || is.null(mutations)) return(list())

    if(product != "seg4") {
      return(if("Position" %in% names(mutations)) list(Position = product) else list())
    }

    schemes <- list()

    if("Position" %in% names(mutations)) {
      schemes$Position <-
        if(input$view_query == "view" || adaptationValues$numbering == "H3") "ha_h3" else "ha_h5"
    }

    if("Position(H5)" %in% names(mutations)) {
      schemes[["Position(H5)"]] <- "ha_h5"
    }

    schemes
  })

  # What hovering a position says, and whether clicking it does anything.
  #
  # A position that cannot be shown says why rather than being silently inert: with the
  # numbers themselves as the links, an unlinked column would otherwise look like a bug.
  jump_tooltip <- function(position, scheme, product) {
    reference <- numbering_refs[[scheme]]
    label     <- product_label(product)

    if(is.null(position) || is.na(position)) {
      return(list(ok = FALSE, text = "No position recorded for this entry."))
    }

    if(position < 0) {
      return(list(ok = FALSE,
                  text = str_c("Position ", position, " is in the HA signal peptide, which has ",
                               "no column in the mature HA alignment.")))
    }

    if(is.null(reference) || is.na(reference)) {
      return(list(ok = FALSE,
                  text = str_c("No verified numbering reference for ", label,
                               " in this database version.")))
    }

    if(!reference %in% subtype_references(product)) {
      return(list(ok = FALSE,
                  text = str_c(label, " positions are numbered against ", reference,
                               ", which is not one of the Tree tab's subtype references, ",
                               "so the tree cannot be pointed at it.")))
    }

    if(is.na(gapless_to_consensus(product, reference, position))) {
      return(list(ok = FALSE,
                  text = str_c("Position ", position, " does not map to a column of the ",
                               label, " alignment using reference ", reference, ".")))
    }

    list(ok = TRUE,
         text = str_c("Show on the tree: segment ", product_segment(product), " coloured by the ",
                      "amino acid at ", label, " position ", position,
                      ", numbered against ", reference, ".",
                      if(product %in% secondary_products) {
                        str_c(" ", label, " shares segment ", product_segment(product),
                              "'s phylogeny, so the tree is that segment's.")
                      } else ""))
  }

  # The alignment column the selected catalogue position occupies, or NA where it has
  # none. Separate from adaptation_frequencies() because the coverage caption needs the
  # column itself rather than the counts drawn from it, and gapless_to_consensus() is a
  # lookup in an in-memory list - cheaper than threading it out through the frequencies.
  selected_column <- reactive({
    position  <- selected_position()
    reference <- numbering_reference()

    if(is.null(selected_product()) || is.na(position) ||
       is.na(reference) || position <= 0) {
      return(NA_integer_)
    }

    gapless_to_consensus(selected_product(), reference, position)
  })

  # How much of the data the bars rest on - see adaptation_coverage()
  coverage <- reactive({
    adaptation_coverage(selected_product(), selected_column(), numbering_reference(),
                        host_grouping(), counting_sequences())
  })

  # The coverage line of the plot caption, given the total the bars actually draw.
  #
  # A function of that total rather than a reactive, because it is only known once the
  # bars are built.
  coverage_caption <- function(plotted) adaptation_coverage_caption(coverage(), plotted)

  # The bars for the selected row - see adaptation_column_frequencies(). The guards are
  # here rather than there because each is a reason the position has no column at all,
  # which frequency_message() explains on screen.
  adaptation_frequencies <- reactive({
    position <- selected_position()
    reference <- numbering_reference()

    if(is.na(position)) { # catalogue entries without a position
      return(ADAPTATION_NO_FREQUENCIES)
    }

    # The HA alignment holds mature peptides only, so signal peptide positions
    # (negative numbering) have no alignment column to map to
    if(selected_product() == "seg4" && position < 0) {
      return(ADAPTATION_NO_FREQUENCIES)
    }

    if(is.na(reference)) {
      return(ADAPTATION_NO_FREQUENCIES)
    }

    req(position > 0) # require a position within the numbered protein

    log_info("Numbering position {position} in {selected_product()} against reference {reference}")

    consensus_position <- gapless_to_consensus(selected_product(), reference, position)

    if(is.na(consensus_position)) {
      log_warn("Position {position} did not map to an alignment column for reference {reference}")
      return(ADAPTATION_NO_FREQUENCIES)
    }

    adaptation_column_frequencies(selected_product(), consensus_position, reference,
                                  host_grouping(), counting_sequences())
  }) %>%
    bindCache(selected_product(), selected_position(), numbering_reference(), # the reference determines the answer, so it belongs in the key
              host_grouping(), sequence_scope()) # without the scope the other mode's plot comes back

  # The catalogue rows as the table shows them: a Tree column carrying a link on each row
  # whose position can actually be shown on the tree, and nothing on the rest.
  #
  # A reactive rather than a step inside renderDT, so what the table will contain can be
  # checked without rendering it - DT sends its rows over ajax, so they are not in the
  # rendered output at all.
  #
  # Row selection already drives the bar plot, so the link is a separate column rather
  # than a modal on every click of a row. It sends the position back as an input instead
  # of navigating. Safe to interpolate: the position is an integer from the catalogue,
  # and rows that fail the check get an empty cell rather than a link built from
  # anything a user supplied.
  displayed_mutations <- reactive({
    mutations <- adaptationValues$mutations
    schemes   <- position_schemes()

    if(is.null(mutations) || nrow(mutations) == 0 || length(schemes) == 0) {
      return(mutations)
    }

    input_id <- session$ns("tree_jump")

    for(column in names(schemes)) {
      # The table's scheme for this column. It settles HA, whose two numbering schemes sit
      # in two columns; for everything else a row's scheme is its own product, and a table
      # holding two products has two of those.
      column_scheme <- schemes[[column]]

      mutations[[column]] <- map2_chr(mutations[[column]], mutations$.product,
        function(position, product) {
          scheme <- if(product == "seg4") column_scheme else product
          hint   <- jump_tooltip(position, scheme, product)
          shown  <- if(is.na(position)) "" else as.character(position)

          # A tooltip either way. title= is what Bootstrap reads, and what the browser
          # falls back to showing if its JS is not there, so the text stays plain.
          if(!hint$ok) {
            return(str_c("<span data-bs-toggle='tooltip' title='", html_attribute(hint$text),
                         "' class='text-muted'>", shown, "</span>"))
          }

          # The scheme travels with the click: HA lists both its numbering schemes side by
          # side, and the position alone would not say which of them was clicked. So does
          # the product, because a click can land on a row other than the selected one and
          # the table may hold two products - position 31 is a site in M1 and another in M2.
          str_c("<a href='#' data-bs-toggle='tooltip' title='", html_attribute(hint$text), "' ",
                "onclick=\"Shiny.setInputValue('", input_id, "', ",
                "{position: ", as.integer(position), ", scheme: '", scheme, "', ",
                "product: '", product, "'}, ",
                "{priority: 'event'}); return false;\">", shown, "</a>")
        })
    }

    # .product and .scheme are the table's own bookkeeping, not columns to show. Dropped
    # here rather than hidden through DataTables, which would still count them when the
    # render below addresses the position columns by index.
    mutations %>% select(-dplyr::any_of(c(".product", ".scheme")))
  })

  # 0-based indices of the linked columns, for the DataTables render below. rownames are
  # off, so the first data column is 0. Taken from the displayed frame, which is what
  # DataTables is handed: it drops two columns and a multi-product table gains a Product
  # column ahead of the rest, so the indices are not those of adaptationValues$mutations.
  position_column_targets <- reactive({
    displayed <- displayed_mutations()
    if(is.null(displayed)) return(integer())
    which(names(displayed) %in% names(position_schemes())) - 1L
  })

  ### The jump to the tree ###

  # What the parent hands to the Tree tab. NULL until a jump is confirmed; every
  # confirmation writes a fresh list, so a second jump to the same position still fires.
  jumpVal <- reactiveVal(NULL)

  # A position clicked in the table. What the jump will do was on the tooltip the user
  # hovered to find the link, so this acts rather than asking again.
  observeEvent(input$tree_jump, {
    request   <- input$tree_jump
    position  <- request$position
    reference <- numbering_refs[[request$scheme]]

    # The product comes with the click rather than from the selection: the table can hold
    # two products, and a link is clicked where it sits, not where the selection is. A
    # payload naming something the app does not know is dropped.
    product <- request$product
    req(!is.null(product), product %in% names(products))

    # The link is rendered only on positions that pass, but the table is HTML and a
    # click can arrive from a page rendered before the product changed, so the check is
    # repeated here rather than trusted. A stale click is dropped, not followed.
    req(jumpable(position, reference, product))

    # HA is the reason the scheme is sent rather than inferred: its two position columns
    # sit side by side, so which one was clicked decides the reference, and with it every
    # residue the tree is about to show.
    req(request$scheme %in% names(numbering_refs))

    jumpVal(list(product = product, reference = reference, position = position))
  })

    ### Rendering ###

    # Adaptation mutations table
    output$adaptationTable <- renderDT({
        caption_text <- str_c("Adaptation mutations for ", table_description())

      mutations <- displayed_mutations()

      table_options <- list(
        pageLength = ADAPTATION_PAGE_LENGTH,
        displayStart = opening_display_start(), # the page the opening row is on
        dom = "frtip", # f adds the search box, so a site can be found without paging
        searchHighlight = TRUE,
        language = list(zeroRecords = "No adaptation mutations detected",
                        search = "Find:",
                        searchPlaceholder = "position or mutation"),

        # The positions are anchors, and a column of HTML would sort as text - "10" before
        # "9" - and would match the markup in the search box. So sorting and filtering see
        # the tag-stripped value and only the display sees the link. One render, applied to
        # whichever position columns this table has.
        columnDefs = if(length(position_column_targets()) > 0) {
          list(list(
            targets = position_column_targets(),
            render = JS(
              "function(data, type, row) {",
              "  if (type === 'display') return data;",
              "  var text = String(data === null ? '' : data).replace(/<[^>]*>/g, '');",
              "  if (type === 'sort' || type === 'type') {",
              "    var value = parseFloat(text);",
              "    return isNaN(value) ? null : value;",
              "  }",
              "  return text;",
              "}")
          ))
        },

        # Bootstrap tooltips have to be attached after every draw, because paging and
        # sorting replace the cells. Skipped rather than erroring if the Bootstrap
        # bundle is not on the page: the same text is in title=, so the browser shows
        # its own tooltip instead.
        drawCallback = JS(
          "function(settings) {",
          "  if (typeof bootstrap === 'undefined' || !bootstrap.Tooltip) return;",
          "  settings.nTable.querySelectorAll('[data-bs-toggle=tooltip]')",
          "    .forEach(function(node) {",
          "      if (!bootstrap.Tooltip.getInstance(node)) { new bootstrap.Tooltip(node); }",
          "    });",
          "}")
      )

      arguments <- list(mutations,
                        caption = caption_text,
                        options = table_options,
                        selection = list(mode = 'single', selected = opening_row()), # the handed-over site, else the first row that yields a plot
                        rownames = FALSE,
                        plugins = "searchHighlight", # mark the matched text
                        escape = FALSE)

      # HA's bare Position and Mutation columns are labelled with the scheme they are
      # expressed in, since the table can carry the other scheme's pair beside them. The
      # argument is added rather than passed as NULL, which datatable() would take as a
      # table with no column names at all.
      if(adaptationValues$product == "seg4") {
        scheme <- adaptationValues$numbering
        arguments$colnames <- set_names(c("Mutation", "Position"),
                                        str_c(c("Mutation(", "Position("), scheme, ")"))
      }

      do.call(datatable, arguments)
    })

    # Adaptation mutations plot
    output$adaptation_plot <- renderGirafe({
      # A product the catalogue has no entry for has no row to select, so every reactive
      # below would req() out and leave the panel blank with nothing said. PB1-F2 and
      # PA-X are that case: their position data is loaded and searchable on the Tree tab,
      # but neither has a catalogued mutation yet.
      # qualified because jsonlite::validate masks shiny::validate here
      shiny::validate(shiny::need(
        nrow(adaptationValues$mutations %||% data.frame()) > 0,
        str_c("No adaptation mutations are catalogued for ", product_description(),
              " yet, so there is no position to plot. Its positions can still be ",
              "searched on the Tree tab.")))

      frequencies <- adaptation_frequencies()

      # say why there is no plot rather than leaving the panel blank
      shiny::validate(shiny::need(nrow(frequencies) > 0, frequency_message()))

      subtitle <- str_c(product_title(),
                        " Position ", selected_position(),
                        " (numbering reference ", numbering_reference(), ")",
                        "\nCounting ", scope_label(),
                        sep = "")

      pooled    <- adaptation_pooled_hosts(frequencies, host_grouping())
      plot_data <- adaptation_plot_data(frequencies, pooled)

      # The counts reach only as far as the host taxonomy they were built from. Say so
      # rather than letting a partial tally read as the whole database.
      taxonomy_caption <- if(counting_sequences() && alignment_counts_partial) {
        str_wrap(str_c("Sequences whose host could not be placed taxonomically are not ",
                       "counted; rerun the taxonomy stage for complete coverage."), width = 95)
      } else {
        NULL
      }

      # Coverage first: it qualifies everything else in the plot, because the bars are
      # proportions and one carried by a fifth of the alignment draws at exactly the
      # same height as one carried by all of it.
      plot_caption <- str_flatten(c(coverage_caption(sum(plot_data$n)),
                                    adaptation_pooled_caption(pooled), taxonomy_caption),
                                  collapse = "\n")
      if(!nzchar(plot_caption)) plot_caption <- NULL

      gg <- adaptation_plot(plot_data,
                            subtitle = subtitle,
                            caption  = plot_caption,
                            x_label  = host_grouping_label(),
                            unit     = scope_unit())

        girafe(ggobj = gg, width_svg = 8, height_svg = 6,
               options = list(
                 opts_hover(css = "fill:cyan;"),
                 opts_tooltip(css = "background-color:white;color:black;padding:5px;border-radius:3px;")
               ))
    })

    # What this tool exports: a confirmed request to show a position on the tree, which
    # the parent passes to tree_server(), and a counter that bumps on this tab's own
    # Clear button, which the parent merges with the other two tools' copies into
    # `reset`. Nothing else crosses between them.
    list(jump = reactive(jumpVal()),
         cleared = reactive(clearedVal()))
  })
}

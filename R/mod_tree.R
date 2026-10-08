#' Tree tool - phylogeny, BLAST search and position search
#'
#' Self-contained: it holds the tree, metadata, BLAST hits and position state and shares
#' no reactive value with the other tools. What reaches it from them arrives through the
#' `jump` reactive the parent passes in.
#'
#' The UI comes in two pieces because the layout interleaves them: a sidebar panel in the
#' left column and a tab body in the right. The parent decides where each goes and owns
#' the `input.tabselected` conditions that show or hide them, since the tab id lives in
#' the parent's namespace.

#' Sidebar controls for the Tree tab
tree_sidebar_ui <- function(id) {
  ns <- NS(id)

                               wellPanel(
                               radioButtons(ns("tree_radio"), 
                                            label = tooltip(
                                              trigger = "Tree",
                                              "View all metadata for given segment or query a protein or nucleotide sequence and auto-detect segment or search amino acids at consensus position"
                                            ),
                                            choices = list("View" = "view", 
                                                           "Query" = "query",
                                                           "Position" = "position"), 
                                            selected = "view",
                                            inline = FALSE),
                               
                               # Segment select, for the modes where the user picks the segment
                               conditionalPanel( 
                                 condition = str_c("['view', 'position'].includes(input['", ns("tree_radio"), "'])"),
                                 
                                 selectInput(ns("tree"),
                                             label = tooltip(
                                               trigger = "Segment",
                                               "Influenza genome segment and protein product. Segments 7 and 8 each encode a second product - M2 and NEP - which share their segment's phylogeny but are numbered against their own protein, so a Position search on M2 reports M2 residues."
                                             ),
                                             choices = product_choices,
                                             selected = "seg1",
                                             selectize = FALSE),
                               ),
                               
                               conditionalPanel(
                                 condition = str_c("input['", ns("tree_radio"), "'] == 'query'"),
                                 textAreaInput(ns("query"), 
                                               tooltip(
                                                 trigger = "Query sequence",
                                                 "Protein or nucleotide sequence for the BLAST search. Nucleotide input is translated by BLAST."
                                               ),
                                               rows = 10
                                 ),
                                 actionButton(ns("submit"), "Submit"),
                                 actionButton(ns("clear"), "Clear")
                               ),
                               
                               # Query mode reports the segment BLAST found but does not offer to
                               # change it. The sequence stays in the box above, so a segment picked
                               # by hand would contradict it. It follows the sequence: submit another
                               # one, or use View or Position to look at a segment directly.
                               conditionalPanel(
                                 condition = str_c("input['", ns("tree_radio"), "'] == 'query' && ",
                                                   "output['", ns("query_segment_ready"), "']"),
                                 div(class = "form-group",
                                     tags$label(class = "control-label",
                                                tooltip(
                                                  trigger = "Segment",
                                                  "Segment BLAST detected for the query sequence. Submit another sequence to change it, or use View or Position to choose one."
                                                )),
                                     div(class = "form-control",
                                         style = "background-color: #e9ecef; cursor: not-allowed;",
                                         textOutput(ns("detected_segment"), inline = TRUE)))
                               ),
                               
                               # Position controls for position mode, and for query once a segment is
                               # known, so a hit can be followed straight to a site without losing the
                               # query highlight by switching mode
                               conditionalPanel(
                                 condition = str_c("input['", ns("tree_radio"), "'] == 'position' || ",
                                                   "(input['", ns("tree_radio"), "'] == 'query' && ",
                                                   "output['", ns("query_segment_ready"), "'])"),
                                 selectInput(ns("sequence_id"), 
                                             label = tooltip(
                                               trigger = "Subtype",
                                               "Reference subtype for amino acid positions"
                                             ), 
                                             choices = NULL,
                                             selectize = TRUE),
                                 
                                 numericInput(ns("position"), 
                                              label = tooltip("Position", 
                                                              "Position in numbering of the reference amino acid sequence. Changing the reference after a search shows the equivalent position in the new reference, and the tree colouring stays the same."),
                                              value = 1,
                                              min = 1),
                                 actionButton(ns("search"), "Search")
                               ),
                             )
}

#' Body of the Tree tab
tree_tab_ui <- function(id) {
  ns <- NS(id)

  tagList(
                            div(reactOutput(ns("taxoniumComponent")), id = "taxonium-container"),
                            
                            markdown(
                              "Maximum likelihood phylogeny of the clustered IAV dataset.
                            The curated database for this segment was clustered using *MMSeq2* on a 0.95 identity threshold.
                            Representative sequences were aligned using *MAFFT*.
                            Maximum likelihood trees were obtained with *IQ-TREE* and midpoint rooted for visualisation.

                              For *Position* search, amino acids are numbered based on the selected reference sequence.
                              Amino acid positions for HA mutations are based on mature peptide numbering. Amino acid positions for mutations in all other segments are based on full-length protein numbering.

                              Sequences with a missing amino acid at the queried position are shown as white nodes in the tree."
                            ),

                            # Clustering and the phylogenies are per segment, so selecting M2 shows
                            # segment 7's tree coloured by M2 residues rather than a phylogeny of M2.
                            # Said here because the paragraph above would otherwise read as the latter.
                            conditionalPanel(
                              condition = str_c("['", str_c(secondary_products, collapse = "','"),
                                                "'].includes(input['", ns("tree"), "'])"),
                              div(class = "alert alert-info",
                                  markdown(
                                    "**Some proteins (e.g., M2 and NEP) share their segment's phylogeny with other proteins.** Clustering and the
                                    trees are per segment, so this is the tree of the parent segment -
                                    the same representatives, in the same clusters - with tips coloured
                                    by the amino acid at the position you searched **in this reading
                                    frame**. It is not a phylogeny of the protein."
                                  ))
                            ),
                            hr(),
                            
                            dataTableOutput(ns("hitsTable"))
  )
}

#' Server for the Tree tool
#'
#' @param jump a reactive carrying a request from another tab to show a position on the
#'   tree - list(product, reference, position), optionally with query and nucleotide -
#'   or NULL. The tools still share no state: this is one reactive, passed in by the
#'   parent, read only here and written only by the requesting tab.
tree_server <- function(id, jump = reactive(NULL), reset = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    treeVal <- reactiveVal() # Newick string to pass to Taxonium component
    metaVal <- reactiveVal() # metadata to pass to Taxonium component
    hitsVal <- reactiveVal() # BLAST hits
    hitsTree <- reactiveVal(NULL) # tree those hits were searched against

    # The dropdown holds a product - "seg7", "seg7_M2" - because the position search is
    # per product: position 11 of M2 is a different residue from position 11 of M1. The
    # tree and its metadata are per segment, since clustering and the phylogenies are,
    # and M2 carries exactly segment 7's representatives in the same order. This is the
    # one place that difference is resolved.
    parent_tree <- reactive({
      req(input$tree)
      str_c("seg", product_segment(input$tree))
    })

    # Hits belong to the segment they were searched against, and are tagged with it
    # rather than cleared when the segment changes: updateSelectInput() in the submit
    # observer changes input$tree itself, so clearing on that would wipe the results
    # the same submission had just produced. Everything downstream reads this, so a
    # segment the hits were not searched against sees none of them.
    #
    # Compared by segment, not by product: BLAST searches a segment's database, and the
    # hits are the same sequences on the same tree whichever of that segment's products
    # is selected. Comparing the product would blank the table on switching M1 to M2.
    current_hits <- reactive({
      if(is.null(hitsVal()) || !identical(hitsTree(), parent_tree())) NULL else hitsVal()
    })
    
    # Strains of the selected BLAST hits, fed to Taxonium's own name search so the matching
    # tips are found and zoomed to rather than left for the user to hunt down. Tagged with
    # the tree they were found in and applied only to that tree: the hits table keeps its
    # contents when the segment changes, and hundreds of strains represent a cluster in
    # more than one segment, so an untagged list could highlight tips in a tree the search
    # was never run against. Tagging rather than clearing also keeps the search alive
    # through the metadata reset that follows an autodetected segment.
    searchVal <- reactiveVal(list(tree = NULL, names = character()))
    
    # Segment BLAST resolved for the submitted query, or NULL before one has been
    # submitted. Drives the segment and position boxes in Query mode: they are only
    # worth showing once there is a real segment behind them.
    querySegment <- reactiveVal(NULL)
    
    # The site last looked up, and the tree it was looked up in. Query mode restores
    # the segment it searched, which reloads that tree and would otherwise put the
    # plain segment metadata back - dropping the amino acid column, and with it the
    # colouring by residue and the entry for it in the Colour by menu.
    positionVal <- reactiveVal(list(tree = NULL, metadata = NULL, search = NULL))

    # The Position search the tree is coloured by right now - product, reference, the
    # alignment column and the position in that reference's numbering - or NULL when it
    # is coloured by nothing in particular. The column is what the colouring is made
    # from; the reference and position are how the user reads it. Changing the Subtype
    # dropdown re-expresses the same column in the new reference's numbering, so the tree
    # stays as it is. See the sequence_id observer below.
    activeSearch <- reactiveVal(NULL)

    # A jump waiting for the product change it asked for to come back from the browser.
    # Held rather than acted on immediately because changing the product re-runs the
    # observer that resets the metadata and the reference, which would undo the jump.
    pendingJump <- reactiveVal(NULL)

    # Bumped by this tab's own Clear button and read by nobody here - it is returned to
    # the parent, merged with the other two tools' own copies, and handed back to every
    # tool as `reset`, so clicking Clear on any one tab resets all three. See
    # reset_tab() below for what "reset" means on this tab.
    clearedVal <- reactiveVal(0L)

    # Label for the read-only segment shown in Query mode, in the dropdown's own
    # segment(product) form. Read out of product_choices rather than a parallel table:
    # querySegment() is the key the submit observer selects in the dropdown, so this is
    # by construction the label sitting in it - "7(M1)" for a query BLAST put on segment
    # 7, which is what the tab really did select.
    output$detected_segment <- renderText({
      req(querySegment())
      names(product_choices)[match(querySegment(), product_choices)]
    })
    
    output$query_segment_ready <- reactive(!is.null(querySegment()))
    outputOptions(output, "query_segment_ready", suspendWhenHidden = FALSE) # read while hidden

    # Metadata for one tree. Precomputed per segment in global.R, because cluster
    # composition is segment specific and a combined table keyed on strain could match the
    # wrong segment for the strains that represent a cluster in more than one of them.
    segment_metadata <- function(tree_name) {
      data.frame(status = "loaded",
                 filename = "test.csv",
                 data = initial_metadata_csv[[tree_name]], # metadata as CSV string
                 filetype = "meta_csv",
                 rows = "") # unique key to rerender component when metadata changes - empty for initial render
    }

  # Returning to Query restores the segment its sequence was searched against. The
    # segment select in View and Position moves input$tree freely, but the query box
    # still holds the sequence that was submitted, so Query goes back to what it found
    # rather than reporting whichever segment was last looked at. Only a new submission
    # changes it.
    observeEvent(input$tree_radio, {
      req(input$tree_radio == "query")
      req(querySegment())
      
      if(!identical(input$tree, querySegment())) {
        updateSelectInput(session = session, inputId = "tree", selected = querySegment())
      }
    })

  # Update tree data and reference select when user selects a tree
    observeEvent(input$tree, {
      # Keyed on the parent segment, not the product, so that switching between the two
      # products of one segment leaves the Newick untouched. Taxonium re-renders on a
      # change of filename, and reloading an identical tree would throw away the user's
      # pan and zoom to show them the same thing.
      selected_data <- data.frame(status = "loaded",
                                filename = parent_tree(),
                                data = tree_data[[parent_tree()]],
                                filetype = "nwk")

      treeVal(selected_data) # Newick string to pass to Taxonium component

      # Query mode returns to the site it was left on, so coming back from View or
      # Position finds the tree coloured by residue as before. Anywhere else the tree
      # loads with its own metadata, clearing search results.
      #
      # Compared by product: a site is a position in one reading frame, so the residues
      # remembered for M2 are not the ones to restore when M1 is selected.
      last_position <- positionVal()

      # A jump from the Adaptation Mutations table arrives here when it had to change the
      # product to get to its position. It carries its own reference and position, and
      # is applied below instead of the metadata reset - which would otherwise wipe the
      # very site the jump exists to show.
      jump <- pendingJump()
      jumping <- !is.null(jump) && identical(jump$product, input$tree)

      if(jumping) {
        activeSearch(NULL) # apply_position() sets it again at the end of this observer
      } else if(input$tree_radio == "query" &&
                identical(last_position$tree, input$tree) &&
                !is.null(last_position$metadata)) {
        metaVal(last_position$metadata)
        activeSearch(last_position$search)
      } else {
        metaVal(segment_metadata(parent_tree()))
        activeSearch(NULL) # the tree is plain again, so there is no position to follow
      }

      cluster_glue_segment <- subtype_references(input$tree) # subtype references for this segment

      # Default reference accessions for A/New_York/392/2004
      default_refs <- c(
        "1" = "NC_007373", # PB2
        "2" = "NC_007372", # PB1
        "3" = "NC_007371", # PA
        "4" = "NC_007366", # HA
        "5" = "NC_007369", # NP
        "6" = "NC_007368", # NA
        "7" = "NC_007367", # M
        "8" = "NC_007370"  # NS
      )
      
      target_ref <- default_refs[as.character(product_segment(input$tree))]
      
      # A jump names the reference its position is numbered against, and that is the only
      # one that gives the right residues. Otherwise the segment's usual default.
      if (jumping) {
        selected_ref <- jump$reference
      } else if (!is.na(target_ref) && target_ref %in% cluster_glue_segment) {
        selected_ref <- target_ref
      } else {
        selected_ref <- cluster_glue_segment[1]
      }

      updateSelectInput(session = session, inputId = "sequence_id",
                        choices = cluster_glue_segment,
                        selected = selected_ref)

      if(jumping) {
        apply_position(jump$product, jump$reference, jump$position)
        pendingJump(NULL)
      }
    })

  #3. BLAST
  #
  # notify = FALSE for a search that sets the tree up behind a click going somewhere
  # else: the Plot link in Batch Screening primes the tree without leaving the Adaptation
  # Mutations tab, and a "Found 500 hits" toast over a bar plot says nothing the user can
  # act on there. The log records it either way, because a background search that found
  # nothing is otherwise invisible.
  run_cluster_blast <- function(query_set, segment, notify = TRUE) {
    blast_db <- blast(db = file.path(conf$paths$data$blast_db, str_c("sgt_", segment)),
                      type = if(inherits(query_set, "DNAStringSet")) "blastx" else "blastp")

    hits_max <- predict(blast_db, query_set, verbose = TRUE,
                        custom_format = str_c("qseqid sseqid pident length mismatch gapopen ",
                                              "gaps qstart qend sstart send evalue bitscore"))

    if(nrow(hits_max) == 0) {
      hitsVal(NULL)
      hitsTree(NULL)
      searchVal(list(tree = NULL, names = character()))
      log_warn("Cluster BLAST against segment {segment} found nothing")
      if(notify) showNotification("No hits found in database", type = "warning")
      return(FALSE)
    }

    if(notify) {
      showNotification(if(nrow(hits_max) > max_blast_hits) {
                         str_c("Found ", nrow(hits_max), " hits, showing the top ", max_blast_hits)
                       } else {
                         str_c("Found ", nrow(hits_max), " hits")
                       }, type = "message")
    }

    displayed <- hits_max %>%
      inner_join(full_metadata, by = join_by("sseqid" == "primary_accession")) %>%
      arrange(evalue) %>%
      dplyr::slice_head(n = max_blast_hits) %>%
      select(-c(qseqid, gapopen, bitscore, qstart, qend, sstart, send,
                length, evalue, segment)) %>%
      relocate(strain) %>%
      dplyr::rename(Strain = strain,
                    `% Identity` = pident,
                    Mismatches = mismatch,
                    Gaps = gaps,
                    H = h_subtype,
                    N = n_subtype,
                    Host = host_group,
                    `Host Order` = host_order,
                    Cluster = cluster_members,
                    `Cluster Hosts` = cluster_hosts,
                    Accession = sseqid)

    # The tree zooms to the name it is handed, so the closest match is what the user
    # arrives on - displayed is sorted by E-value, so that is its first row.
    hitsVal(displayed)
    hitsTree(str_c("seg", segment))
    searchVal(list(tree = str_c("seg", segment), names = displayed$Strain[1]))

    log_info("Cluster BLAST segment {segment}: {nrow(hits_max)} hit(s), ",
             "closest {displayed$Strain[1]}")
    TRUE
  }

  # Puts this tab back to what it shows before anything has been searched. Run only from
  # `reset` below, so every tool resets itself the same way regardless of which tab's
  # Clear button caused it.
  reset_tab <- function() {
    # Back to View - the tab's own default mode - rather than left on Query or Position
    # with nothing behind it. Reachable now from outside this tab too: a reset arriving
    # while Query or Position was selected would otherwise leave those boxes on screen,
    # emptied but still there, which reads as broken rather than as reset.
    updateRadioButtons(session = session, inputId = "tree_radio", selected = "view")

    updateTextAreaInput(session, "query", value = "")

    # Segment box and the Subtype/Position/Search box are only shown once a query has
    # resolved to a segment - see query_segment_ready above. Clearing that hides both,
    # the same as if no query had ever been submitted.
    querySegment(NULL)

    # The rest puts the tree itself back to its pre-search state, not only the sidebar:
    # no BLAST hits, no highlighted tips, no site remembered to restore. Same reset
    # run_cluster_blast() uses for "no hits found" - see above - plus metaVal, which
    # nothing else recomputes without a change of input$tree.
    hitsVal(NULL)
    hitsTree(NULL)
    searchVal(list(tree = NULL, names = character()))
    positionVal(list(tree = NULL, metadata = NULL, search = NULL))
    activeSearch(NULL)
    metaVal(segment_metadata(parent_tree()))
  }

  # Only tells the parent - this tab resets when `reset` comes back, like the other two.
  observeEvent(input$clear, clearedVal(clearedVal() + 1L))

  # Any tab's Clear button, this one's included, arrived through the parent - see
  # clearedVal above.
  observeEvent(reset(), {
    reset_tab()
  }, ignoreInit = TRUE)

  # Run BLAST when user submits a query sequence.
  observeEvent(input$submit,{
    req(input$tree_radio == "query") # require view_query input to be view
    req(input$query) # require query sequence to be available

    # Reading the query and detecting its segment is shared with the Adaptation
    # Mutations tool - see R/query_segment.R. NULL means the query was unusable and the
    # user has already been told why.
    detected <- resolve_query_segment(input$query)
    req(detected)

    # A DNAStringSet for a nucleotide query, an AAStringSet for a protein one.
    # run_cluster_blast() picks blastx or blastp off the class, so a nucleotide query needs
    # nothing further here.
    seqs    <- detected$seqs
    segment <- detected$segment

    updateSelectInput(session = session,
                      inputId = "tree",
                      selected = str_c("seg", segment)) # update tree selection

    # Reveal the segment and position boxes now that the segment is known, and
    # clear any position left over from an earlier search so the box appears
    # empty for the user to fill in against their own hit
    querySegment(str_c("seg", segment))
    updateNumericInput(session = session, inputId = "position", value = NA)
    positionVal(list(tree = NULL, metadata = NULL, search = NULL)) # the site belonged to the old sequence
    activeSearch(NULL)

    run_cluster_blast(seqs, segment)
  }) # end observeEvent

  # Colour the tree by the residue at one position of one product, numbered against one
  # reference.
  #
  # A function rather than the body of the Search observer, because a jump from the
  # Adaptation Mutations table needs the same work done with values it supplies rather
  # than with whatever is in the inputs. Setting the inputs and firing the button would
  # not do: updateSelectInput() reaches the browser and comes back asynchronously, so
  # the search would run against the previous reference.
  #
  # Returns TRUE if the tree was recoloured.
  apply_position <- function(product, sequence_id, position) {
    consensus_aas <- consensus(product, sequence_id, position)

    if(is.null(consensus_aas) || length(consensus_aas) == 0) {
      showNotification(str_c("Position ", position, " could not be located in ",
                             product_label(product), " using reference ", sequence_id),
                       type = "error")
      return(FALSE)
    }

    # What the cluster behind each tip carries at this site, not only its representative -
    # a tip stands for up to tens of thousands of members, which is why cluster_hosts sits
    # beside it in the pop-up.
    #
    # Joined on the accession before it is dropped. Left join and a stated fallback: a
    # representative that was never clustered has no members to break down, and a blank
    # field in the pop-up would read as a missing value rather than an absent cluster.
    column    <- gapless_to_consensus(product, sequence_id, position)
    breakdown <- cluster_residue_summary(product, column)

    if(!is.null(breakdown)) {
      consensus_aas <-
        consensus_aas %>%
        left_join(breakdown, by = join_by(Node == representative)) %>%
        mutate(cluster_residues = replace_na(cluster_residues, "no cluster members"))
    }

    consensus_aas <-
      consensus_aas %>%
      inner_join(full_metadata, by = join_by(Node == primary_accession)) %>%
      relocate(strain) %>%  # first column, to match the tree tip labels
      select(-Node, -segment)

    metadata_text <- consensus_aas %>% format_csv

    selected_metadata <- data.frame(status = "loaded",
                                    filename = "test.csv",
                                    data = metadata_text,
                                    filetype = "meta_csv",
                                    rows = paste(c(product, sequence_id, position), collapse="-")) # the key that remounts the component
    search <- list(product = product, reference = sequence_id, column = column, position = position)

    metaVal(selected_metadata)
    positionVal(list(tree = product, metadata = selected_metadata, search = search)) # the site to return to
    activeSearch(search)
    TRUE
  }

  # Switching reference while a Position search is on the tree: the same alignment column
  # in the new reference's numbering. The tree is not touched - it is coloured by column,
  # and the column has not changed - so only the Position box moves.
  #
  # A reference with a gap at that column has no equivalent position, and a stale number
  # in the box would read as the answer, so the box is emptied and the user told. The
  # column is kept, so the next reference that does have a residue there converts from it.
  observeEvent(input$sequence_id, {
    search <- activeSearch()

    # Not following a search: nothing on the tree, a different product, or the dropdown
    # being refilled for one (its change arrives after the tree observer cleared the search)
    req(search, identical(search$product, input$tree), nzchar(input$sequence_id %||% ""))

    if(identical(search$reference, input$sequence_id)) return() # an echo of the search itself

    if(is.null(reference_numbering)) return() # no mapping built: the number stays as typed

    converted <- convert_reference_position(search$product, search$reference, input$sequence_id,
                                            position = search$position, column = search$column)

    new_label <- names(subtype_references(search$product))[match(input$sequence_id, subtype_references(search$product))]

    if(is.na(converted$column)) { # the reference the search came from is not in the table
      return()
    }

    activeSearch(list(product = search$product, reference = input$sequence_id,
                      column = converted$column, position = converted$position))

    if(is.na(converted$position)) {
      updateNumericInput(session = session, inputId = "position", value = NA)
      showNotification(str_c("No equivalent position: ", new_label, " has no ",
                             product_label(search$product), " residue at the alignment column of ",
                             "position ", search$position, " (a gap in that reference). The tree ",
                             "colouring is unchanged."),
                       type = "warning")
    } else {
      updateNumericInput(session = session, inputId = "position", value = converted$position)
    }
  })

  # Search for amino acids at the consensus position in the sequence data
  observeEvent(input$search, {
    # Query mode reaches here too, once BLAST has resolved a segment: the position
    # box is offered there so a hit can be followed to a site without switching mode
    req(input$tree_radio %in% c("position", "query"))

    if(is.null(input$position) || is.na(input$position)) { # the box starts empty in Query mode
      showNotification("Enter a position to search", type = "error")
      req(FALSE) # exit
    }

    req(input$position > 0)
    message("Consensus search at position ", input$position, " for sequence ", input$sequence_id, " in tree ", input$tree)

    # input$tree is already the product key consensus() wants ("seg1"), so it is passed
    # through rather than reduced to a segment number and looked up by position
    apply_position(input$tree, input$sequence_id, input$position)
  })

  # A position sent over from the Adaptation Mutations table, or from Batch Screening.
  #
  # It arrives already checked: the sending tab only offers the link where the position
  # maps to a column and its numbering reference is one this tab can be pointed at. What
  # is left is to show it - set the controls so the tab reads as if the user had driven
  # it, then colour the tree.
  #
  # A request carrying a query sequence also gets the cluster search, so the tab is
  # waiting on the closest matching tip rather than on the whole tree. focus = FALSE says
  # the user is not being brought here - the Plot link sets this tab up on its way to the
  # Adaptation Mutations tab - so that search runs without announcing itself.
  observeEvent(jump(), {
    request <- jump()
    req(request)

    focused   <- !identical(request$focus, FALSE)
    has_query <- !is.null(request$query) && nzchar(request$query)
    updateRadioButtons(session = session, inputId = "tree_radio",
                       selected = if(has_query) "query" else "position")
    updateNumericInput(session = session, inputId = "position", value = request$position)

    if(has_query) {
      updateTextAreaInput(session = session, inputId = "query", value = request$query)
      querySegment(str_c("seg", product_segment(request$product)))
      positionVal(list(tree = NULL, metadata = NULL, search = NULL))
      activeSearch(NULL)

      query_set <- if(isTRUE(request$nucleotide)) {
        DNAStringSet(request$query)
      } else {
        AAStringSet(request$query)
      }
      run_cluster_blast(query_set, product_segment(request$product), notify = focused)
    }

    if(identical(input$tree, request$product)) {
      # Already on that product, so the observer above will not fire and nothing else
      # will set the reference or the metadata
      updateSelectInput(session = session, inputId = "sequence_id",
                        choices = subtype_references(request$product),
                        selected = request$reference)
      apply_position(request$product, request$reference, request$position)
    } else {
      # Changing the product fires the observer above, which resets the metadata and the
      # reference. Hand it the request so it applies this instead of undoing it.
      pendingJump(request)
      updateSelectInput(session = session, inputId = "tree", selected = request$product)
    }
  })

  # Circle the selected hits on the tree
  #
  # Through Taxonium's name search rather than by replacing the metadata, so the tree keeps
  # the segment's own: every tip stays annotated and the colour scheme stays whatever the
  # user chose.
  observeEvent(input$hitsTable_rows_selected, {
    req(current_hits() %>% nrow > 0) # require hits for the displayed tree
    
    strains <-
      current_hits() %>% 
      dplyr::slice(input$hitsTable_rows_selected) %>% 
      pull(Strain) # tree tip labels are the strain names
    
    # Tagged with the segment, matching current_hits(): the highlighted tips are the same
    # tips of the same tree whichever product of that segment is selected
    searchVal(list(tree = parent_tree(), names = strains)) # feed the hit names to the tree search
  })

    ### Rendering ###

    output$taxoniumComponent <- renderReact({
      hit_search <- searchVal()
      
      # only search the tree the hits were actually found in
      hit_names <- if(identical(hit_search$tree, parent_tree())) hit_search$names else character()
      
      TaxoniumComponent(treeData = treeVal(), meta = metaVal(), search = hit_names)
    })

    # BLAST hits table
    output$hitsTable <- renderDT({
                                 req(current_hits()) # no table on a segment with no hits of its own
                                 current_hits()
                                 },
                                 caption = "BLAST hits",
                                 options = list(pageLength = 100),
                                 selection = list(mode = 'multiple', selected = c(1)), # select first row
                                 rownames = FALSE)

    list(cleared = reactive(clearedVal()))
  })
}

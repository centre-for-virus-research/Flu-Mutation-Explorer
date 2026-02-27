source("global.R")

ui<- function(id) {
  
  
  ns <- NS(id) # create namespace function
  
  segment_choices <- c("1(PB2)" = "seg1",
                       "2(PB1)" = "seg2",
                       "3(PA)" = "seg3",
                       "4(HA)" = "seg4",
                       "5(NP)" = "seg5",
                       "6(NA)" = "seg6",
                       "7(M1)" = "seg7",
                       "8(NS1)" = "seg8"
  )
  
  page_fluid(
    theme = bs_theme(bootswatch = "lumen"),
    
    tagList(
      tags$head(
        # local npm installation of Taxonium React component version 2.0.227
        tags$script(src="node_modules/taxonium-component/dist/taxonium-component.umd.js"), 
        
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
        tags$link(rel = "shortcut icon", href = "favicon-32x32.png"),
        includeHTML("google-analytics.html") # Google Analytics HTML
      ),
      CustomComponents # React and Taxonium
    ),
    
    page_navbar(
      title = "Flu-Mutation Explorer",
      addGFontHtmlDependency(family = c("Open Sans")),
      footer = tags$footer(
        fluidRow(
          column(6, tags$a(img(src = "mrc_uog_cvr_logo.png", title = "MRC-University of Glasgow Centre for Virus Research", alt = "MRC-University of Glasgow Centre for Virus Research logo", class = "img-fluid mx-auto footer-logo"), href="https://www.gla.ac.uk/researchinstitutes/iii/cvr/", target = "_blank")),
          column(3, tags$a(img(src = "rvc_logo.png", title = "Royal Veterinary College", alt = "Royal Veterinary College logo", class = "img-fluid mx-auto footer-logo rvc-logo"), href="https://www.rvc.ac.uk", target = "_blank")),
          column(3, tags$a(img(src = "cropped_logo.png", title = "FLU:TrailMap-One Health", alt = "FLU:TrailMap-One Health logo", class = "img-fluid mx-auto footer-logo"))),
        )
      ),      
      
      div(
        class = "container-fluid",
        fluidRow(
          column(
            width = 2,
            # insert conditionalPanel  image
            conditionalPanel(condition = "input.tabselected==0 || input.tabselected==3",
                             img(src = "Warhol Influenza.png", 
                                 title = "Flu-Mutation Explorer",
                                 alt = "Flu-Mutation Explorer pop art viruses logo",
                                 class = "img-fluid mx-auto d-block logo-image"
                             )),
            
            conditionalPanel(condition = "input.tabselected==1",
                             wellPanel(
                               radioButtons(ns("tree_radio"), 
                                            label = tooltip(
                                              trigger = "Tree",
                                              "View all metadata for given segment or query a sequence and auto-detect segment or search amino acids at consensus position"
                                            ),
                                            choices = list("View" = "view", 
                                                           "Query" = "query",
                                                           "Position" = "position"), 
                                            selected = "view",
                                            inline = FALSE),
                               
                               # conditionally display segment select for view or position options
                               conditionalPanel( 
                                 condition = str_c("['view', 'position'].includes(input['", ns("tree_radio"), "'])"),
                                 
                                 selectInput(ns("tree"), 
                                             label = tooltip(
                                               trigger = "Segment",
                                               "Influenza genome segment"
                                             ), 
                                             choices = segment_choices,
                                             selected = "seg1",
                                             selectize = FALSE),
                               ),
                               
                               conditionalPanel(
                                 condition = str_c("input['", ns("tree_radio"), "'] == 'query'"),
                                 textAreaInput(ns("query"), 
                                               tooltip(
                                                 trigger = "Query sequence",
                                                 "Amino acid sequence for BLAST search"
                                               ),
                                               rows = 10
                                 ),
                                 actionButton(ns("submit"), "Submit")
                               ),
                               
                               conditionalPanel(
                                 condition = paste0('input[\'', ns('tree_radio'), "\'] == \'position\'"),
                                 selectInput(ns("sequence_id"), 
                                             label = tooltip(
                                               trigger = "Subtype",
                                               "Reference subtype for amino acid positions"
                                             ), 
                                             choices = NULL,
                                             selectize = TRUE),
                                 
                                 numericInput(ns("position"), 
                                              label = tooltip("Position", 
                                                              "Position in numbering of reference amino acid sequence"), 
                                              value = 1,
                                              min = 1),
                                 actionButton(ns("search"), "Search")
                               ),
                             )
            ),
            conditionalPanel(condition = "input.tabselected==2",
                             wellPanel(
                               radioButtons(ns("view_query"), 
                                             label = tooltip(
                                               trigger = "Adaptations",
                                               "View all adaptation mutations for a given segment or query a sequence and auto-detect segment"
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
                                               "Influenza genome segment"
                                             ),
                                             choices = segment_choices,
                                             selected = "seg1",
                                             selectize = FALSE
                                 )
                               ), 
                               
                               conditionalPanel(
                                 condition = paste0('input[\'', ns('view_query'), "\'] == \'query\'"),
                                 textAreaInput(ns("query_adaptation"), 
                                               tooltip(
                                                 trigger = "Query sequence",
                                                 "Amino acid sequence to match with reference"
                                               ), 
                                               rows = 10, 
                                               width = NULL, 
                                               placeholder = NULL),
                                 actionButton(ns("submit_adaptation"), "Submit")
                               ) # end conditionalPanel for query sequence
                             ) # end wellPanel for adaptation mutations
            ), # end conditionalPanel adaptation mutations
          ),
          column(width = 10,
                 tabsetPanel(
                   id = "tabselected",
                   type = "tabs",
                   
                   tabPanel("Home", 
                            icon = icon("home"), 
                            value = 0,
                           fluidRow(
                              column(12, 
                                     card(
                                       card_header(
                                         class = "bg-dark",
                                         "Refreshing the Application"
                                       ),
                                       tags$p("If the application becomes unresponsive please refresh the page to continue.")
                                     )),
                              
                            ),  # end of fluidRow

                            fluidRow(
                              column(12, 
                                     card(
                                       card_header(
                                         class = "bg-dark",
                                         "Quick Start Guide"
                                       ),
                                       includeMarkdown("directions.md")
                                     )),
                              
                            ),  # end of fluidRow
                            fluidRow(
                              column(8,
                                     card(
                                       card_header(
                                         class = "bg-dark",
                                         "Status"
                                       ),
                                       tableOutput(ns("statusTable")),
                                       tableOutput(ns("statusSegmentsTable"))
                                     )),
                              
                              column(4, 
                                     card(
                                       card_header(
                                         class = "bg-dark",
                                         "Support"
                                       ),
                                       HTML("<a href='mailto:cvr-webresource-support@lists.cent.gla.ac.uk?Subject=Flu-Mutation-Explorer%20Support' target='_top'>Contact us</a>")
                                     ),
                                     
                                     card(
                                       card_header(
                                         class = "bg-dark",
                                         "Licence"
                                       ),
                                       markdown("
This work is licensed under the [GNU General Public Licence v3.0][gplv3].  
See the [GNU GPL v3][gplv3] for details.

[gplv3]: https://www.gnu.org/licenses/gpl-3.0.html
                                                ")
                                     )
                              ), # end column
                              
                              # column(6, 
                              #        
                              #        ),
                              
                              # #### FAVICON TAGS SECTION ####
                              # tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
                              
                              # bsModal("modalExample", "Instructional Video", "tabBut", size = "large" ,
                              #         p("Additional text and widgets can be added in these modal boxes. Video plays in chrome browser"),
                              #         iframe(width = "560", height = "315", url_link = "https://www.youtube.com/embed/0fKg7e37bQE")
                              # )
                              
                            ) # end fluidRow
                          ), # end tabPanel Home
                   
                   tabPanel("Tree", value=1,         
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
                            hr(),
                            
                            dataTableOutput(ns("hitsTable"))
                   ),
                   tabPanel("Adaptation Mutations", value=2,
                            fluidRow(
                              column(width = 6, girafeOutput(ns("adaptation_plot")),
                                     markdown(
                                       "
                              Amino acid positions for HA mutations are based on mature peptide numbering.
                              Amino acid positions for mutations in all other segments are based on full-length protein numbering.
                              **These mammalian adaptation mutations are largely derived from _in vitro_ studies performed in specific viral strain backgrounds. Consequently, their adaptive significance may be context-dependent and may not apply to the submitted query sequence.**
                                       "
                                     )
                              ),
                              column(width = 6, dataTableOutput(ns("adaptationTable")))
                            )
                   ),
                   tabPanel("About", value=3,
                            includeMarkdown("about.md")
                   )
                 )
          )
        )
      )      
    )
  )
}

server <- function(id) {
  moduleServer(id, function(input, output, session) {
    treeVal <- reactiveVal() # Newick string to pass to Taxonium component
    metaVal <- reactiveVal() # metadata to pass to Taxonium component
    hitsVal <- reactiveVal() # BLAST hits
    adaptationValues <- reactiveValues(numbering = "H3") # adaptation mutations
    
    initial_metadata <- data.frame(status = "loaded",
                                    filename = "test.csv",
                                    data = initial_metadata_csv, # metadata as CSV string
                                    filetype = "meta_csv",
                                    rows = "") # unique key to rerender component when metadata changes - empty for initial render
    
  # Update tree data and reference select when user selects a tree
    observeEvent(input$tree, {
      selected_data <- data.frame(status = "loaded",
                                filename = input$tree,
                                data = tree_data[[input$tree]],
                                filetype = "nwk")
      
      treeVal(selected_data) # Newick string to pass to Taxonium component
      metaVal(initial_metadata) # reset metadata to clear search results
      
      cluster_glue_segment <- # GLUE reference sequences for selected segment
        cluster_glue %>% 
        filter(segment == input$tree %>% str_extract("\\d+") %>% as.integer) %>% # extract segment number from tree name
        select(strain_display, representative) %>% 
        deframe # convert to named vector for updateSelectInput
      
      cluster_glue_segment

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
      
      segment_num <- input$tree %>% str_extract("\\d+")
      target_ref <- default_refs[segment_num]
      
      # Only use target reference if it already exists in cluster_glue_segment
      # Otherwise use the first available option
      if (!is.na(target_ref) && target_ref %in% cluster_glue_segment) {
        selected_ref <- target_ref
      } else {
        selected_ref <- cluster_glue_segment[1]
      }

      updateSelectInput(session = session, inputId = "sequence_id",
                        choices = cluster_glue_segment,
                        selected = selected_ref)
    })

    
    # Display all adaptation mutations for selected segment
    observe({
      req(input$view_query == "view") # require view_query input to be view
      
      adaptationValues$numbering <- "H3" # reset numbering
      adaptation_mutations <- adaptation_mutations_all[[input$protein]]

      adaptation_mutations %<>%
      select(-c(new_amino_acid,
                wt, how_was_the_mutation_discovered, notes)) %>%  # remove columns for display
        dplyr::rename(Position = position,
                      Mutation = mutation,
                      Subtype = which_virus) %>%   # rename columns for display
        rowwise %>%
        mutate(Paper = papers_dois(paper_s, doi), .keep = "unused") %>%  # TODO precompute this column
        ungroup # remove rowwise grouping
      
      if(input$protein == "seg4") {
        adaptation_mutations %<>%
          dplyr::rename(`Position(H5)` = position_h5,
                        `Mutation(H5)` = mutation_h5_numbering) # rename columns for display
      } 
      
      segment_number <- input$protein %>% str_extract("\\d+") %>% as.integer # extract segment number
      
      adaptationValues$mutations <- adaptation_mutations # store adaptation mutations for selected segment
      adaptationValues$segment <- segment_number # store segment name for selected segment
    })

  #3. BLAST
  # TODO: validate input
  # TODO: handle errors
  # TODO: handle long running tasks
  # Run BLAST when user submits a query sequence
  observeEvent(input$submit,{
    req(input$tree_radio == "query") # require view_query input to be view
    req(input$query) # require query sequence to be available
    
    trimmed <- 
      input$query %>% 
      str_trim %>% # trim whitespace
      str_replace_all("[\r\n-]" , "") # remove line breaks and gaps (-)
    
    req(trimmed) # require trimmed string to not be empty
    
    aas <- NULL # amino acids
    # try and catch error if query sequence is invalid
    tryCatch({
      aas <- AAStringSet(trimmed)
    }, error = function(e) {
      message(e$message)
      showNotification("Invalid query sequence", type = "error")
    })
    req(aas) # exit if AAStringSet not created
    
    # Autodetect BLAST search against reference database then get top hit
    blast_db_ref <- blast(db = "BLAST_segment_recognizer/ref_set_AA", type = "blastp")
    hits_max_ref <- 
      predict(blast_db_ref, aas, verbose = TRUE) %>%  # default args returns 500 hits
      arrange(evalue) # sort by E-value
    
    # check if hits_max_ref is empty
    # if no hits found, show error message
    if(nrow(hits_max_ref) == 0) {
      showNotification("Unable to identify segment", type = "error")
      req(FALSE) # exit
    } else {
      # get sseqid from top hit
      sseqid <- 
        hits_max_ref %>%
        dplyr::slice(1) %>% # get top hit
        pull(sseqid)  # get sseqid
      print(sseqid) # DEBUG
        
      # look up segment number in ref_set
      segment <- 
        ref_set %>% 
        filter(accession_version == sseqid) %>% # get segment number
        pull(segment)
      print(segment) # DEBUG
      
      showNotification(str_c("Detected segment ", segment), type = "message")
      
      updateSelectInput(session = session,
                        inputId = "tree",
                        selected = str_c("seg", segment)) # update tree selection
    
      # use autodetected segment
      blast_db <- blast(db = file.path(conf$paths$data$blast_db, str_c("sgt_", segment)), type = "blastp")
    } # end else

    hits_max <<- predict(blast_db, aas, verbose = TRUE # default args returns 500 hits
                       # , BLAST_args = "-max_target_seqs 200000 -evalue 10e-38"
                       # , BLAST_args = "-max_target_seqs 200000"
                       )
    
    if(nrow(hits_max) == 0) {
      showNotification("No hits found in database", type = "warning")
    } else {
      showNotification(str_c("Found ", nrow(hits_max), " hits"), type = "message")
      
      # PB1 TFEFTSFFYRYGFVANFSMELPSFGVSGI
      # PB2 ATKRLTVLGKDAGALTEDPDEGTAGVESAVLRGFLILGKEDKRYGPALSINELSNLAKGE
      # store BLAST result for table display
      hitsVal(hits_max %>%
                inner_join(full_metadata, by = join_by("sseqid" == "primary_accession")) %>% # join strain to BLAST hits
                arrange(evalue) %>% # sort by E-value
                select(-c(
                  qseqid, gapopen, bitscore, 
                  qstart, qend, sstart, send,
                  length, evalue, segment)) %>% # remove columns for display
                relocate(strain) %>% # relocate strain to first column to match tree node IDs
                dplyr::rename(Strain = strain, # rename columns for display
                              `% Identity` = pident, 
                              Mismatches = mismatch,
                              H = h_subtype, 
                              N = n_subtype,
                              Host = host_group,
                              `Host Order` = host_order,
                              Accession = sseqid
                )) # display names for BLAST hits table and Taxonium component
    } # end else
  }) # end observeEvent
  
  # Search for amino acids at the consensus position in the sequence data
  observeEvent(input$search, {
    req(input$tree_radio == "position" && input$position > 0) # require view_query input to be view and positive position
    message("Consensus search at position ", input$position, " for sequence ", input$sequence_id, " in tree ", input$tree) # DEBUG"
    segment_num <- input$tree %>% str_extract("\\d+") %>% as.integer # extract segment number from tree name
    
    # try and catch error if query sequence is invalid
      consensus_aas <- 
        consensus(segment_num, # extract segment number from tree name
                  input$sequence_id, 
                  input$position) 
      
      if(is.null(consensus_aas) || length(consensus_aas) == 0) {
        showNotification("No amino acids found at consensus position", type = "error")
        req(FALSE) # exit
      }
    
    consensus_aas <- 
      consensus_aas %>% 
      inner_join(full_metadata, by = join_by(Node == primary_accession)) %>%  # join other metadata to amino acids
      relocate(strain) %>%  # relocate strain to first column to match tree node IDs
      select(-Node, -segment) # remove primary accession and segment number from tree metadata display
    
    metadata_text <- consensus_aas %>% format_csv  # metadata as CSV string to pass to React component
    
    selected_metadata <- data.frame(status = "loaded",
                                    filename = "test.csv",
                                    data = metadata_text, # metadata as CSV string
                                    filetype = "meta_csv",
                                    rows = paste(c(input$tree, input$sequence_id, input$position), collapse="-") ) # unique key to rerender component
    metaVal(selected_metadata)
  })
  
  # Filter hits_max according to user selected hit rows
  observeEvent(input$hitsTable_rows_selected, {
    req(hitsVal() %>% nrow > 0) # require hits data to be available
    
    # metadata for Taxonium
    hits_meta <-
      hitsVal() %>% 
      dplyr::slice(input$hitsTable_rows_selected) %>% 
      # mutate(`Query Sequence` = input$query, .after = 1) %>% # add query sequence to metadata in first column
      mutate(`BLAST Hit` = "TRUE", .after = 1) %>% # add query hit to metadata in first column
      select(-Mismatches, -`% Identity`, -Accession) # remove columns for Taxonium metadata
    
    glimpse(hits_meta) # DEBUG
    
    metadata_text <- hits_meta %>% format_csv # metadata as CSV string to pass to React component

    selected_metadata <- data.frame(status = "loaded",
                                    filename = "test.csv",
                                    data = metadata_text, # metadata as CSV string
                                    filetype = "meta_csv",
                                    rows = paste(input$hitsTable_rows_selected, collapse="-") ) # selected rows to use as key to rerender component
    metaVal(selected_metadata)
  })
  
  ### Adaptation mutations ###
  
  # User submits AA sequence for adaptation mutations
  # Detect segment
  # Detect H3 or H5 for HA
  # Perform pairwise alignment with reference sequence
  observeEvent(input$submit_adaptation, {
    req(input$view_query == "query") # require view_query option to be query
    req(input$query_adaptation) # require adaptation query sequence to be available
      
    trimmed <-
      input$query_adaptation %>%
      str_trim %>% 
      str_replace_all("[\r\n-]" , "") # remove line breaks and gaps (-)
    
    req(trimmed) # require trimmed string to not be empty
      
    # TODO refactor common code into function
    aas <- NULL # amino acids
    # try and catch error if query sequence is invalid
    tryCatch({
      aas <- AAStringSet(trimmed)
    }, error = function(e) {
      message(e$message) # print error message to console
      showNotification("Invalid query sequence", type = "error")
    })
    req(aas) # exit if AAStringSet not created
      
    # Autodetect BLAST search against reference database then get top hit
    blast_db_ref <- blast(db = conf$paths$data$ref_set_aa, type = "blastp")
    hits_max_ref <-
      predict(blast_db_ref, aas, verbose = TRUE) %>%  # default args returns 500 hits
      arrange(evalue) # sort by E-value
    
    # check if hits_max_ref is empty
    # if no hits found, show error message
    if(nrow(hits_max_ref) == 0) {
      showNotification("Unable to identify segment", type = "error")
      req(FALSE) # exit
    } else {
      # get sseqid from top hit
      sseqid <-
        hits_max_ref %>%
        dplyr::slice(1) %>% # get top hit
        pull(sseqid)  # get sseqid
      print(sseqid) # DEBUG
        
      # look up segment number in ref_set
      segment <-
        ref_set %>%
        filter(accession_version == sseqid) %>% # get segment number
        pull(segment)
      print(segment) # DEBUG
      
      showNotification(str_c("Detected segment ", segment), type = "message")
      protein <- str_c("seg", segment) # protein name
        
      updateSelectInput(session = session,
                        inputId = "protein",
                        selected = protein) # update segment selection
      
      # If segment 4 (HA) then BLAST to determine H3 or H5
      if(segment == 4) {
        blast_db_h3_h5 <- blast(db = conf$paths$data$ref_h3_h5, type = "blastp")
        hits_max_h3_h5 <-
          predict(blast_db_h3_h5, aas, verbose = TRUE) %>%
          arrange(evalue) %T>% print # sort by E-value
        
        # get sseqid from top hit
        sseqid <-
          hits_max_h3_h5 %>%
          dplyr::slice(1) %>% # get top hit
          pull(sseqid)  # get sseqid
        print(sseqid) # DEBUG
        
        # if YP_308839.1 then H3, else H5
        if(sseqid == "YP_308839.1") {
          index <- 1 # AAStringSet index for H3
          showNotification("Using H3 numbering", type = "message")
        } else {
          index <- 2 # AAStringSet index for H5
          showNotification("Using H5 numbering", type = "message")
        }
      } else {
        index <- 1 # AAStringSet index for non-HA references - only one sequence for each
      }
    }
    
  # MKTIIALSYILCLAFGQNLPGNDSSTATLCLGHHAVLNGTIVKT              
    # H3 test       QKLPGNDNSTATLCLGWCNVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQ
    # H5 test DQICIGWQNNNSTEQVDTIMEKNVTVTHAQDILEKTHNGKLCDLNGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIV

    message("Protein: ", protein) # DEBUG    
    # tryCatch block to handle errors for pairwiseAlignment in case of invalid amino acid sequence
    pairwise <- NULL # pairwise alignment
    tryCatch({
      # Perform pairwise alignment with reference sequence TODO display alignment TODO refine alignment parameters gap penalty?
      pairwise <- pairwiseAlignment(aas, ref_seqs[[protein]][[index]], substitutionMatrix = "BLOSUM80")
      print(pairwise)
    }, error = function(e) {
      message(e$message) # print error message to console
      showNotification(str_c("Pairwise alignment error: ", e$message, ". Ensure the sequence contains only valid amino acid characters."), type = "error")
    })
    req(pairwise) # exit observer if alignment not created
    
    # Extract mismatches from alignment
    mismatches <- 
      pairwise %>% 
      mismatchTable %>%
      clean_names 
      
    # if mismatches has no rows
    if(nrow(mismatches) == 0) {
      print(mismatches)
      showNotification("No pairwise mismatches found", type = "default")
      # create empty table to display no results
      mismatches <- tibble(Mutation = character(), Position = integer(), 
                               Subtype = character(), Paper = character())
    } else {
      showNotification(str_c("Found", nrow(mismatches), "pairwise mismatches", sep = " "), 
                       type = "warning")
      
      # adaptationValues$numbering <- "H3"
      adaptation_mutations <- adaptation_mutations_all[[protein]]
      if(protein == "seg4") {
        if(sseqid == "YP_308669.1"){ # if H5
          adaptation_mutations %<>% 
            filter(position_h5 > 0) %>%  # filter out negative positions
            mutate(position = position_h5, 
                   mutation = mutation_h5_numbering, 
                   .keep = "unused") # use H5 position and mutation name
          adaptationValues$numbering <- "H5" # store numbering type
        } else { # H3
          adaptation_mutations %<>% 
            select(-position_h5, 
                   -mutation_h5_numbering) # remove H5 columns
          adaptationValues$numbering <- "H3" # store numbering type
        }
      }
      
      ## Check if mismatches are known adaptation mutations
      # get pattern amino acid and subject positions from mismatches
      # join with known adaptation mutations on position and amino acid mutation
      mismatches %<>% 
        inner_join(adaptation_mutations, 
                   by = c("subject_start" = "position", "pattern_substring" = "new_amino_acid")) %>% # join mismatches with known adaptation mutations
        select(-c(starts_with("pattern"), subject_end, subject_substring, wt, how_was_the_mutation_discovered, notes)) %>%  # remove columns for display
        dplyr::rename(Position = subject_start, 
                      Mutation = mutation,
                      Subtype = which_virus) # rename columns for display
      
      # if mismatches has no rows
      if(nrow(mismatches) == 0) {
        showNotification("No adaptation mutations found", type = "warning")
        mismatches %<>% 
          dplyr::rename(Paper = paper_s) %>% # rename column for display
          select(Mutation, Position, Subtype, Paper, -doi) # reorder/remove columns for display
      } else {
        mismatches %<>% 
          rowwise %>%
          mutate(Paper = papers_dois(paper_s, doi), .keep = "unused") %>%  # TODO precompute this column
          ungroup # remove rowwise grouping
        
        showNotification(str_c("Found", nrow(mismatches), "adaptation mutations", sep = " "), 
                         type = "warning")
      }
    } # end if mismatches
    
    # display results in table
    adaptationValues$mutations <- mismatches # store adaptation mutations for query sequence
    adaptationValues$segment <- segment # store segment name for query sequence
    
  }) # end observeEvent for adaptation mutations query
  
  # Calculate mutation frequencies according to user selected row
  adaptation_frequencies <- reactive({
    req(adaptationValues$mutations %>% nrow > 0, 
        input$adaptationTable_rows_selected) # require adaptation hits and user to select a row
    # get reference position of selected adaptation mutation
    position <- 
      adaptationValues$mutations %>% 
      dplyr::slice(input$adaptationTable_rows_selected) %>% 
      pull(Position)
    
    # for segment 4 (HA) get the consensus position from the mature peptide position
    # if view selected use H3 numbering 
    # H3 mature peptide accession NC_007366
    print(adaptationValues$segment) # DEBUG
    
    # segment <- input$protein %>% str_extract("\\d+")
    
    message("Position: ", position) # DEBUG
    message("Numbering: ", adaptationValues$numbering) # DEBUG
    
    if(adaptationValues$segment == 4){ # HA protein
      message("Segment 4 (HA)") # DEBUG  
      if(input$view_query == "view" | adaptationValues$numbering == "H3"){ # use H3 numbering for View
        if(position < 0) { # negative position relative to mature peptide numbering
          # Use KX121222 until NC_007366 full sequence added. Mature peptide starts at position 17 in gapless sequence.
          consensus_position <- gapless_to_consensus(adaptationValues$segment, "KX121222", 17 + position) # add negative position to first amino acid position
        } else { # positive numbering 
          consensus_position <- gapless_to_consensus(adaptationValues$segment, "NC_007366", position) # get consensus position for H3 mature HA peptide
        }
      } else { # H5 numbering
        consensus_position <- gapless_to_consensus(adaptationValues$segment, "NC_007362", position) # get consensus position for H3 mature HA peptide
      }
    } else { # for all other segments use the position as is
      # TODO consensus position for turkey England references
      req(position > 0) # require consensus position to be a number greater than 0
      
      # A/turkey/England/50-92/1991 accessions for segments 1-8
      # Loaded from config
      turkey_ref_id <- conf$references$turkey[[paste0("seg", adaptationValues$segment)]]
      log_info("Turkey Reference selected for segment {adaptationValues$segment}: {turkey_ref_id}")
      
      consensus_position <- gapless_to_consensus(adaptationValues$segment, 
                                                 turkey_ref_id, # accession for segment reference
                                                 position) 
    }
    
    adaptationValues$position <- position # store position for selected mutation
    message("Consensus position: ", consensus_position) # DEBUG
    
    if(is.na(consensus_position)) {
      return(data.frame(amino_acid = character(), host_group = character(), n = integer()))
    }

    # get amino acid at consensus position for all sequences in protein
    amino_acids <- 
      transposed %>% 
      pluck(adaptationValues$segment, consensus_position) %>% 
      unlist # convert to vector
    
    if(is.null(amino_acids) || length(amino_acids) <= 1) {
       return(data.frame(amino_acid = character(), host_group = character(), n = integer()))
    }
    
    aa_count <- 
      amino_acids[-1] %>% # remove reference amino acid
      enframe(name = "primary_accession", value = "amino_acid") %>% 
      inner_join(full_metadata) %>% 
      filter(!is.na(host_group)) %>% # remove sequences with missing host group
      filter(!host_group %in% c("Unknown", "Environment")) %>% # remove Unknown or Environment
      select(primary_accession, amino_acid, host_group) %>%
      mutate(amino_acid = factor(amino_acid, levels = get_aa_levels())) %>%
      dplyr::count(amino_acid, host_group)
    
    aa_count
  }) %>%
    bindCache(adaptationValues$segment, adaptationValues$mutations$Position[input$adaptationTable_rows_selected], adaptationValues$numbering)
  
  
  ### Rendering ### 
    output$taxoniumComponent <- renderReact({
      TaxoniumComponent(treeData = treeVal(), meta = metaVal())
    })
    
    # display system path DEBUG
    output$system_path <- renderText({
      Sys.getenv(c("PATH"))
    })
    
    # BLAST hits table
    output$hitsTable <- renderDT(hitsVal(),
                                 caption = "BLAST hits", 
                                 options = list(pageLength = 100),
                                 selection = list(mode = 'multiple', selected = c(1)), # select first row
                                 rownames = FALSE) 
    
    # Adaptation mutations table
    output$adaptationTable <- renderDT({
        caption_text <- str_c("Adaptation mutations for segment ", 
                              adaptationValues$segment,
                              sep = "")
      
      if(adaptationValues$segment == 4 && adaptationValues$numbering == "H5") {
        datatable(adaptationValues$mutations,
                  caption = caption_text,
                  options = list(pageLength = 10, 
                                 dom = "rtip",
                                 language = list(zeroRecords = "No adaptation mutations detected")),
                  colnames = c("Mutation(H5)" = "Mutation", "Position(H5)" = "Position"),
                  selection = list(mode = 'single', selected = c(1)), # select first row
                  rownames = FALSE,
                  escape = FALSE) # don't escape HTML
      } else if(adaptationValues$segment == 4 && adaptationValues$numbering == "H3") {
        datatable(adaptationValues$mutations,
                  caption = caption_text,
                  colnames = c("Mutation(H3)" = "Mutation", "Position(H3)" = "Position"),
                  options = list(pageLength = 10, 
                                 dom = "rtip",
                                 language = list(zeroRecords = "No adaptation mutations detected")),
                  selection = list(mode = 'single', selected = c(1)), # select first row
                  rownames = FALSE,
                  escape = FALSE) # don't escape HTML
      } else { # for all other segments
        datatable(adaptationValues$mutations,
                  caption = caption_text,
                  options = list(pageLength = 10, 
                                 dom = "rtip",
                                 language = list(zeroRecords = "No adaptation mutations detected")),
                  selection = list(mode = 'single', selected = c(1)), # select first row
                  rownames = FALSE,
                  escape = FALSE) # don't escape HTML
      }
    }) 
    
    output$statusTable <- renderTable(status, colnames = FALSE)
    output$statusSegmentsTable <- renderTable(status_segments, na = "") # override transparency for "NA" display
    
    # Adaptation mutations plot
    output$adaptation_plot <- renderGirafe({
      subtitle <- str_c("Segment ", adaptationValues$segment, 
                        " Position ", adaptationValues$position, 
                        sep = "")
      
      gg <-
        adaptation_frequencies() %>%
        group_by(host_group, amino_acid) %>%
        summarise(n = sum(n)) %>% # total per amino acid in host group
        mutate(sum_n = sum(n)) %>% # total per host group
        mutate(percentage = n / sum_n) %>% 
        ungroup %>% 
        ggplot(aes(fill = amino_acid,
                   y = n,
                   x = host_group,
                   tooltip = str_c(amino_acid, "<br>",
                                   n, "/", sum_n, " clusters<br>", 
                                   scales::percent(percentage, accuracy = 0.1)),
                   data_id = interaction(host_group, amino_acid)
                   )) +
        geom_bar_interactive(position="fill", stat="identity") +
        labs(title = "Frequency of amino acids at consensus position",
             subtitle = subtitle,
             x = "Amino Acid",
             y = "Frequency",
             fill = "Amino Acid") +
        scale_fill_manual(values = AA_PALETTE, na.value = "grey50") +
        scale_y_continuous(labels = scales::percent) +
        theme_minimal(base_family='Open Sans') +
        theme(axis.text=element_text(size=rel(1.25), 
                                     colour = "black"),
              axis.title=element_text(size=rel(1.5)),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title =element_text(size=rel(1.75))
        )
      
        girafe(ggobj = gg, width_svg = 8, height_svg = 6,
               options = list(
                 opts_hover(css = "fill:cyan;"),
                 opts_tooltip(css = "background-color:white;color:black;padding:5px;border-radius:3px;")
               ))
    })
    
  })
}

shinyApp(ui("app"), function(input, output) server("app"))

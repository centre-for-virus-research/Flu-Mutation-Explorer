#' Batch Screening - screen many sequences against the adaptation mutation catalogue
#'
#' The n-many counterpart of the Adaptation Mutations query box. Takes up to
#' MAX_SCREEN_SEQUENCES protein or nucleotide sequences, resolves each to a product by
#' BLAST against the reference proteins, and reports which catalogued positions each one
#' carries.
#'
#' Two tables rather than one. The long table is the record, but it is the wrong thing to
#' land on: measured on 100 whole segments it runs to nearly two thousand rows, and 57% of
#' them sit at the residue most of the database already carries. So the summary is the
#' entry point - one row per query, the outcome named - and selecting a row filters the
#' detail below it. That is what carries the limit at 2,000: the same set produces some
#' 37,000 detail rows, which is a table to filter rather than a table to read.

#' @rdname screening_ui
#' @export
screening_sidebar_ui <- function(id) {
  ns <- NS(id)

  wellPanel(
    radioButtons(ns("source"),
                 label = tooltip(
                   trigger = "Sequences",
                   "Paste FASTA or a bare sequence, or upload a FASTA file. Amino acid or nucleotide, but not both in one submission."
                 ),
                 choices = list("Paste" = "paste", "Upload FASTA" = "upload"),
                 selected = "paste",
                 inline = TRUE),

    conditionalPanel(
      condition = paste0('input[\'', ns('source'), "\'] == \'paste\'"),
      textAreaInput(ns("pasted"),
                    label = NULL,
                    rows = 10,
                    width = "100%",
                    placeholder = str_c(">sequence_1\nMERIKELRDLMSQSRTREILTKTTVDHMAIIKKYTSGRQEKNPSLRMKWMMAMKYPITADKRITEM...\n",
                                        ">sequence_2\n..."))
    ),

    conditionalPanel(
      condition = paste0('input[\'', ns('source'), "\'] == \'upload\'"),
      fileInput(ns("file"), label = NULL, accept = c(".fasta", ".fa", ".fas", ".txt"),
                buttonLabel = "Choose", placeholder = "No file")
    ),

    actionButton(ns("submit"), "Screen"),
    actionButton(ns("clear"), "Clear"),

    tags$hr(),

    downloadButton(ns("download"), "Download results", class = "btn-sm"),

    tags$p(class = "text-muted small mt-3",
           str_c("Up to ", format(MAX_SCREEN_SEQUENCES, big.mark = ","), " sequences, all ",
                 "nucleotide or all amino acid. Every sequence submitted gets a row, ",
                 "including those carrying no catalogued position. A large set takes a ",
                 "while: ", format(MAX_SCREEN_SEQUENCES, big.mark = ","),
                 " whole nucleotide segments is about a minute."))
  )
}

#' @rdname screening_ui
#' @export
screening_tab_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(width = 12,
             uiOutput(ns("status")),
             div(class = "mb-2",
                 uiOutput(ns("tableHeader"), inline = TRUE)),
               DT::dataTableOutput(ns("resultsTable")),
               markdown(
                 "
HA mutation positions use mature peptide numbering for H3 (A/New_York/392/2004) or H5 (A/goose/Guangdong/1/1996), in whichever of the H3
and H5 schemes the sequence matched, while mutations in all other segments use full-length protein numbering based on A/New_York/392/2004.
**Assessed** indicates the number of catalogued mutation positions covered by the submitted sequence. Partial sequences are screened only across the region they span; therefore, none found in a fragment does not indicate that no mutations are present in the complete protein.

**These mammalian adaptation mutations are largely derived from _in vitro_ studies performed in specific viral strain backgrounds. Consequently, their adaptive significance may be context-dependent and may not apply to the submitted sequences.**
                 "
               )
             )
    )
  )
}

#' Server for the Batch Screening tool
#'
#' @return a reactive carrying a confirmed request to show a position on the tree -
#'   list(product, reference, position, query, nucleotide) - or NULL, exactly as
#'   adaptation_server() does.
#' @export
screening_server <- function(id, reset = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {

    results <- reactiveVal(NULL)
    submitted_records <- reactiveVal(NULL)

    # Bumped by this tab's own Clear button and read by nobody here - it is returned to
    # the parent, merged with the other two tools' own copies, and handed back to every
    # tool as `reset`, so clicking Clear on any one tab resets all three.
    clearedVal <- reactiveVal(0L)

    ### Intake ###

    submitted_text <- reactive({
      if(identical(input$source, "upload")) {
        req(input$file)
        str_flatten(readLines(input$file$datapath, warn = FALSE), "\n")
      } else {
        req(input$pasted)
        input$pasted
      }
    })

    # Puts this tab back to what it shows before anything has been screened. Run only
    # from `reset` below, so every tool resets itself the same way regardless of which
    # tab's Clear button caused it.
    reset_tab <- function() {
      updateTextAreaInput(session, "pasted", value = "")
      shinyjs::reset("file")

      # Same reset the submit observer starts with, so a stale screen is not left showing
      # once the input it was run on is gone.
      results(NULL)
      submitted_records(NULL)
    }

    # Only tells the parent - this tab resets when `reset` comes back, like the other two.
    observeEvent(input$clear, clearedVal(clearedVal() + 1L))

    # Any tab's Clear button, this one's included, arrived through the parent - see
    # clearedVal above.
    observeEvent(reset(), {
      reset_tab()
    }, ignoreInit = TRUE)

    observeEvent(input$submit, {
      results(NULL)
      submitted_records(NULL)

      records <- parse_sequence_records(submitted_text())

      if(nrow(records) == 0) {
        showNotification("No sequences found in that input", type = "error")
        return()
      }

      if(nrow(records) > MAX_SCREEN_SEQUENCES) {
        showNotification(str_c(format(nrow(records), big.mark = ","),
                               " sequences submitted - the limit is ",
                               format(MAX_SCREEN_SEQUENCES, big.mark = ","),
                               ". Split the set and screen it in parts."),
                         type = "error", duration = 10)
        return()
      }

      # All one type or none. A mixed submission is refused rather than guessed at,
      # because the two are searched with different programs.
      nucleotide <- is_nucleotide(records$sequence)

      if(any(nucleotide) && any(!nucleotide)) {
        mixed <- records$name[!nucleotide]
        showNotification(
          str_c("That set mixes nucleotide and amino acid sequences. ",
                sum(nucleotide), " read as nucleotide, ", sum(!nucleotide), " as protein",
                if(length(mixed) <= 3) str_c(" (", str_flatten_comma(mixed), ")") else "",
                ". Submit one type at a time."),
          type = "error", duration = 12)
        return()
      }

      all_nucleotide <- all(nucleotide)

      if(is.null(conf$paths$data$reference_protein_db) ||
         length(Sys.glob(str_c(conf$paths$data$reference_protein_db, ".p*"))) == 0) {
        showNotification("The reference protein database is missing - expected at the reference_protein_db path in config.yml",
                         type = "error", duration = 12)
        return()
      }

      # start/pid logged here, done or failed/pid logged below once withProgress()
      # returns - since this process is shared by every concurrent session (no scheduler
      # configured in shiny-server.conf), a second request's start line falling inside
      # this window says it waited for this one's blastx call to return before its own
      # could run.
      pid <- Sys.getpid()
      started_at <- Sys.time()
      log_info("Screening {nrow(records)} {if(all_nucleotide) 'nucleotide' else 'protein'} sequence(s) [pid {pid}]")

      screened <- withProgress(
        message = str_c("Screening ", format(nrow(records), big.mark = ","), " sequences"),
        value = 0.3,
        {
          out <- tryCatch(screen_sequences(records, nucleotide = all_nucleotide),
                          error = function(e) {
                            log_error("Screening failed: {e$message}")
                            NULL
                          })
          incProgress(0.7)
          out
        })

      elapsed <- round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 2)

      if(is.null(screened)) {
        log_info("Screening failed after {elapsed}s [pid {pid}]")
        showNotification("The screen could not be completed - see the server log", type = "error")
        return()
      }

      log_info("Screening done for {nrow(records)} sequence(s) in {elapsed}s [pid {pid}]")

      submitted_records(records)
      results(screened)

      unresolved <- sum(is.na(screened$summary$product))
      showNotification(
        str_c(format(nrow(records), big.mark = ","), " screened - ",
              format(sum(screened$summary$hits), big.mark = ","),
              " catalogued position(s) found",
              if(unresolved > 0) {
                str_c("; ", format(unresolved, big.mark = ","),
                      " sequence(s) could not be identified")
              } else ""),
        type = "message", duration = 8)
    })

    ### Results ###

    output$status <- renderUI({
      res <- results()
      if(is.null(res)) {
        return(div(class = "alert alert-light border",
                   "Paste or upload sequences, then press Screen. Nucleotide sequences are translated by alignment, so whole genomic segments can be submitted; secondary reading frames (PB1-F2, PA-X, M2 and NEP) are also screened."))
      }
      NULL
    })

    shown_rows <- reactive({
      res <- results()
      req(res)
      screen_table(res)
    })

    output$tableHeader <- renderUI({
      res <- results()
      if(is.null(res)) return("Results")

      queries <- nrow(res$summary)
      found   <- nrow(res$hits)

      tagList(
        "Results",
        tags$span(class = "ms-2 fw-normal small",
                  str_c("(", format(queries, big.mark = ","), " sequence",
                        if(queries != 1) "s" else "",
                        ", ", format(found, big.mark = ","), " catalogued position",
                        if(found != 1) "s" else "", " found)")))
    })

    # The View column carries the two places a position can be looked at: the phylogeny,
    # coloured by the residue at that site, and the bar plot of what each host carries
    # there. The product and reference travel with the click, because an HA hit numbered in
    # H5 has to arrive at the tree with the H5 reference or every residue shown is the
    # wrong one. The position itself is plain text, saying which scheme it is expressed in.
    linked_positions <- reactive({
      rows <- shown_rows()
      if(nrow(rows) == 0) return(rows)
      records <- submitted_records()
      req(records)

      tree_id <- session$ns("tree_jump")
      plot_id <- session$ns("plot_jump")

      # Whether a query has a sequence behind it, as a named lookup rather than a filter()
      # over the whole submission per row: 2,000 whole segments produce some 37,000 hit
      # rows, so scanning every record for each of them costs about 30 s of the wait.
      usable_query <- records %>%
        dplyr::slice_head(n = 1, by = name) %>%
        dplyr::transmute(name, usable = !is.na(sequence) & nzchar(sequence))
      usable_query <- set_names(usable_query$usable, usable_query$name)

      # What a row's links need from its product, looked up once per product rather than
      # once per row. Per row, these lookups and the escaping below cost about 27 s of a
      # 2,000-segment screen's 37,000 rows; a screen only ever has a dozen products.
      by_product <- tibble(product = unique(stats::na.omit(rows$product))) %>%
        mutate(tree_product   = purrr::map_chr(product, screen_catalogue_key),
               tree_reference = purrr::map_chr(product, ~ screen_reference(.x) %||% NA_character_),
               tree_label     = purrr::map_chr(tree_product, products_label_safe))

      rows %>%
        left_join(by_product, by = "product") %>%
        mutate(
          # HA is catalogued in two numbering schemes and both read as "HA", so the
          # position has to say which one it is expressed in or it means nothing
          position_shown = case_when(is.na(position) ~ "—",
                                     is.na(scheme)   ~ as.character(position),
                                     .default        = str_c(position, " on ", scheme)),

          # The sequence's name, not the sequence: carrying the query in both onclicks
          # put 6.3 MB of sequence text into the DOM for a screen of 100. The observers
          # look it up in submitted_records() when the click arrives instead.
          query_js = js_string(qseqid),

          tree = str_c(
            "<a href='#' data-bs-toggle='tooltip' title='",
            html_attribute(str_c("Colour the tree by the residue at ",
                                 tree_label, " position ", position,
                                 ", numbered against ", tree_reference)),
            "' onclick=\"Shiny.setInputValue('", tree_id, "', ",
            "{position: ", position, ", product: '", tree_product,
            "', reference: '", tree_reference, "', query_id: ", query_js,
            "}, {priority: 'event'}); return false;\">Tree</a>"),

          # The Plot link sets the tree up on the way to the Adaptation Mutations tab,
          # so a user who looks at the bar plot and then opens the Tree tab finds it
          # already coloured at this site rather than blank.
          #
          # The tree's position travels separately, because the two tabs number this
          # row differently: the bar plot is keyed on plot_position, which for HA is
          # always H3, while a row the screen resolved to H5 belongs on the tree at its
          # H5 position against the H5 reference. Same residue and same alignment
          # column either way, but showing an H5 sequence at an H3 position says the
          # screen found something it did not.
          plot = if_else(is.na(plot_position), "", str_c(
            "<a href='#' data-bs-toggle='tooltip' title='",
            html_attribute(str_c("Amino acids by host at ",
                                 tree_label, " position ", plot_position)),
            "' onclick=\"Shiny.setInputValue('", plot_id, "', ",
            "{position: ", plot_position, ", product: '", tree_product,
            "', tree_position: ", position,
            ", tree_reference: '", tree_reference, "'",
            ", query_id: ", query_js,
            "}, {priority: 'event'}); return false;\">Plot</a>")),

          view = case_when(
            is.na(position) | is.na(product) ~ "",
            # unname() and coalesce(): a name that is not in the lookup comes back NA
            !coalesce(unname(usable_query[qseqid]), FALSE) ~ "",
            is.na(tree_reference) | position < 1 ~ "<span class='text-muted small'>not on tree</span>",
            .default = str_c(tree, if_else(nzchar(plot), " <span class='text-muted'>|</span> ", ""), plot))) %>%
        select(-c(tree_product, tree_reference, tree_label, query_js, tree, plot))
    })

    output$resultsTable <- DT::renderDT({
      rows <- linked_positions()

      if(nrow(rows) == 0) {
        return(DT::datatable(
          tibble(` ` = "No sequences screened yet."),
          rownames = FALSE, options = list(dom = "t")))
      }

      # Once per distinct paper and DOI rather than once per row: 37,000 rows carry some
      # 150 of them, and building the links per row cost about 9 s
      papers <- rows %>%
        distinct(reference, doi) %>%
        mutate(Paper = purrr::map2_chr(reference, doi, papers_dois))

      rows %>%
        left_join(papers, by = c("reference", "doi")) %>%
        transmute(
          Query = qseqid,
          Product = coalesce(screen_product_label(product), "—"),
          Mutation = coalesce(mutation, "—"),
          Position = position_shown,
          View = view,
          Outcome = outcome,
          Assessed = if_else(catalogued > 0, str_c(assessed, " of ", catalogued), "—"),
          Paper)
    },
    escape = FALSE,
    rownames = FALSE,
    selection = "none",
    filter = "top",
    options = list(pageLength = 25, scrollX = TRUE,
                   columnDefs = list(list(className = "dt-nowrap", targets = 0:6),
                                     list(orderable = FALSE, targets = 4))),
    callback = DT::JS(
      "table.on('draw.dt', function() {",
      "  var t = $('[data-bs-toggle=\"tooltip\"]');",
      "  if (window.bootstrap && t.length) { t.each(function(){ new bootstrap.Tooltip(this); }); }",
      "});"))

    ### Download ###

    # screen_download() shapes the rows - see R/screening.R. Before any screen it returns
    # the same columns with no rows, so the file is a header rather than an empty download
    # the browser cannot make sense of.
    output$download <- downloadHandler(
      filename = function() str_c("batch_screening_", format(Sys.Date(), "%Y%m%d"), ".csv"),
      contentType = "text/csv",
      content = function(file) {
        readr::write_csv(screen_download(results()), file, na = "")
      })

    ### Handing a position to the other two tools ###

    # Two requests, two reactives, both written here and read by the parent. Neither tool
    # is reached into; each is asked.
    jumpVal <- reactiveVal(NULL)
    showVal <- reactiveVal(NULL)

    # Every click is a fresh request, even when it names the site already asked for:
    # reactiveVal does not notify when set to a value identical to the one it holds, so a
    # second click on the same link would reach nobody. The serial makes each request its
    # own value.
    #
    # Read inside an observeEvent handler, which runs isolated, so this creates no
    # dependency and no loop.
    requests <- reactiveVal(0L)
    next_request <- function() {
      serial <- requests() + 1L
      requests(serial)
      serial
    }

    # What the Tree tab is asked for. The query sequence travels with it: without one the
    # tab colours the whole tree at that site, with one it also searches the segment's
    # clustered sequences and zooms to the closest match, which is the point of arriving
    # there from a screened sequence rather than from the catalogue.
    #
    # focus says whether the user is being taken to the tree or whether it is only being
    # set up behind a click going to the Adaptation Mutations tab.
    tree_request <- function(request, focus) {
      query <- clicked_query(request$query_id)

      list(product    = request$product,
           reference  = request$reference,
           position   = as.integer(request$position),
           query      = query,
           nucleotide = !is.null(query) && is_nucleotide(query),
           focus      = focus,
           serial     = next_request())
    }

    # The sequence a click names. Looked up rather than carried, so the table does not
    # hold a copy of every screened sequence per row - see the View column above.
    #
    # NULL when the id is not among the sequences on screen, which is a click from a
    # table rendered before the results changed. The tree then colours the site without a
    # cluster search rather than searching for a sequence that is no longer displayed.
    clicked_query <- function(id) {
      records <- submitted_records()
      if(is.null(id) || is.null(records) || nrow(records) == 0) return(NULL)

      sequence <- records$sequence[records$name == id]

      if(length(sequence) == 0 || is.na(sequence[1]) || !nzchar(sequence[1])) return(NULL)
      sequence[1]
    }

    # Both links can set the tree up, so both are held to what the Tree tab needs: a
    # product it holds position data for, and a reference to number it against.
    primes_tree <- function(request) {
      !is.null(request$reference) && nzchar(request$reference) &&
        request$product %in% names(discarded)
    }

    observeEvent(input$tree_jump, {
      request <- input$tree_jump
      req(request$product, request$reference, request$position)

      # The link is only drawn where the reference resolved, but a click can arrive from a
      # table rendered before the results changed, so this is checked again rather than
      # trusted - a stale click is dropped, not followed.
      req(primes_tree(request))

      jumpVal(tree_request(request, focus = TRUE))
    })

    observeEvent(input$plot_jump, {
      request <- input$plot_jump
      req(request$product, request$position)

      # The Adaptation Mutations tab offers only the products it has a catalogue for, so
      # a request for one it does not list would select nothing and leave the tab sitting
      # on whatever it showed before.
      req(request$product %in% names(adaptation_mutations_all))

      showVal(list(product = request$product,
                   position = as.integer(request$position),
                   serial = next_request()))

      # The tree is prepared at the same site on the way past. Silently, and without
      # changing tab - the click asked for the bar plot. A product the tree cannot show
      # is simply not prepared; the bar plot above does not depend on it.
      #
      # In the tree's own numbering, not the bar plot's: for an HA row the two differ,
      # and the tree is the one that has to agree with the product the screen reported.
      prime <- list(product   = request$product,
                    reference = request$tree_reference,
                    position  = request$tree_position,
                    query_id  = request$query_id)

      if(!is.null(prime$position) && primes_tree(prime)) {
        jumpVal(tree_request(prime, focus = FALSE))
      }
    })

    list(jump = reactive(jumpVal()),
         show = reactive(showVal()),
         cleared = reactive(clearedVal()))
  })
}

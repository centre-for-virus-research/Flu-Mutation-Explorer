source("global.R")

ui<- function(id) {
  
  
  ns <- NS(id) # create namespace function
  
  page_fluid(
    theme = bs_theme(bootswatch = "lumen"),
    
    tagList(
      useShinyjs(),
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
      title = "Flu Mutation Explorer",
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
                                 title = "Flu Mutation Explorer",
                                 alt = "Flu Mutation Explorer pop art viruses logo",
                                 class = "img-fluid mx-auto d-block logo-image"
                             )),
            
            conditionalPanel(condition = "input.tabselected==1",
                             tree_sidebar_ui(ns("tree"))
            ),
            conditionalPanel(condition = "input.tabselected==2",
                             adaptation_sidebar_ui(ns("adaptation"))
            ), # end conditionalPanel adaptation mutations
            conditionalPanel(condition = "input.tabselected==4",
                             screening_sidebar_ui(ns("screening"))
            ), # end conditionalPanel batch screening
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
                                       tags$p("If the application becomes unresponsive please refresh the page to continue. In R Shiny, timeouts are standard practice designed to save money/energy, reduce our carbon footprint, free up resources and keep the server running smoothly.")
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
                            tree_tab_ui(ns("tree"))
                   ),
                   tabPanel("Adaptation Mutations", value=2,
                            adaptation_tab_ui(ns("adaptation"))
                   ),
                   tabPanel("Batch Screening", value=4,
                            screening_tab_ui(ns("screening"))
                   ),
                   tabPanel("About", value=3,
                            # Wrapped so custom.css can style about.md's markdown tables
                            # without touching the Shiny tables on the other tabs.
                            div(class = "about-md", includeMarkdown("about.md"))
                   )
                 )
          )
        )
      )      
    )
  )
}

#' Parent server
#'
#' Owns only what belongs to the page itself - the Home tab's status tables and the tab
#' switching - and starts the three tools. They have separate namespaces and share no
#' reactive state; what passes between them passes through here.
server <- function(id) {
  moduleServer(id, function(input, output, session) {

    output$statusTable <- renderTable(status, colnames = FALSE)
    output$statusSegmentsTable <- renderTable(status_segments, na = "") # override transparency for "NA" display

    # Clicking Clear on any one tab returns every tab to its defaults, not only the one
    # it was clicked on. Each tool bumps its own "cleared" counter on its own Clear
    # button - see clearedVal in each R/mod_*.R - merged here into one and handed back
    # to all three as `reset`, so every tool resets itself the same way whether its own
    # button was clicked or another tab's was.
    resetAll <- reactiveVal(0L)

    # Batch Screening asks Adaptation Mutations for a bar plot through one reactive,
    # written by the tool that asks and read by the tool that acts, rather than shared.
    screening  <- screening_server("screening", reset = resetAll)
    adaptation <- adaptation_server("adaptation", show = screening$show, reset = resetAll)

    # Both tools can ask for a position on the tree, which reads one request, so they are
    # merged here. Whichever fired last wins.
    jump <- reactiveVal(NULL)
    observeEvent(adaptation$jump(), jump(adaptation$jump()), ignoreInit = TRUE)
    observeEvent(screening$jump(),  jump(screening$jump()),  ignoreInit = TRUE)

    tree <- tree_server("tree", jump = jump, reset = resetAll)

    observeEvent(screening$cleared(),  resetAll(resetAll() + 1L), ignoreInit = TRUE)
    observeEvent(adaptation$cleared(), resetAll(resetAll() + 1L), ignoreInit = TRUE)
    observeEvent(tree$cleared(),       resetAll(resetAll() + 1L), ignoreInit = TRUE)

    # The tabsetPanel id is unnamespaced, because the sidebar's conditionalPanel
    # conditions refer to it, so it has to be addressed through the root scope: on this
    # module's session updateTabsetPanel() would target "app-tabselected", which is
    # nothing.
    observeEvent(jump(), {
      request <- jump()
      req(request)

      # Batch Screening's Plot link also asks the tree for its site, so the tree is ready
      # if the user goes on to it, but marks the request focus = FALSE because the
      # observer below is taking them to the bar plot instead. Requests with no focus at
      # all come from Adaptation Mutations, whose only link is the tree one.
      req(!identical(request$focus, FALSE))

      updateTabsetPanel(session$rootScope(), "tabselected", selected = "1")
    })

    # The View column's other link asks for the bar plot rather than the tree, so it
    # lands on the Adaptation Mutations tab instead.
    observeEvent(screening$show(), {
      req(screening$show())
      updateTabsetPanel(session$rootScope(), "tabselected", selected = "2")
    })
  })
}

shinyApp(ui("app"), function(input, output) server("app"))

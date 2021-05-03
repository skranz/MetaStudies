#' Create and run the MetaStudies shiny app
#'
#' @param ... additional parameters passed to \code{shinyApp}
MetaStudiesApp = function(csv.file=NULL,...) {
  ui <- fluidPage(
    # titlePanel(h1("Optimal treatment assignment given covariates", align = "center")),
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose CSV File of estimates and standard errors",
          multiple = FALSE,
          accept = ".csv"
        ),
#          accept = c("text/csv",
#                     "text/comma-separated-values,text/plain",
#                     ".csv")),
        hr(),
        h3("Model parameters"),
        fluidRow(
          column(6, checkboxInput("symmetric", "Symmetric p(.)", value = FALSE),
                    checkboxGroupInput("cutoffs", "Cutoffs for p(.)",
                                          choiceNames = c("1.65", "1.96","2.33"),
                                          choiceValues = c(1.645, 1.960, 2.326),
                                    selected = 1.960)),
          column(6, radioButtons("modelmu", "Model for the distribution of effects",
                                 choices = c("Normal" = "normal",
                                             "Student-t" = "t" #,
                                                #"nonparametric" = 3
                                             ),
                                             selected = "normal"))
        ),
        #hr(),
        actionButton(inputId = "estimatebutton", label = "Estimate model")
      ),

      mainPanel(
        h2("Funnel plot, histogram of z-stats", align = "center"),
        fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                      plotOutput("funnel"),
                      plotOutput("hist")))
      )
    ),
    hr(),
    h2("Model estimates", align = "center"),
    h4("Distribution of true effects, conditional publication proabilities", align = "center"),
    column(12, align="center",
           tableOutput("estimatestable"),
           plotOutput("estplot", width = "70%")
    )
  )

  server <- function(input, output, session) {
    #object to store data and estimation results
    v = reactiveValues()

    load.data = function(file) {
      metadata=read.table(file,sep=",")
      v$X=metadata[,1]
      v$sigma=metadata[,2]
    }

    observe({
      cat("\nObserve input$file1 change")
      req(input$file1)
      file = input$file1$datapath
      restore.point("observe input$file1")
      load.data(file)
    })


    if (!is.null(csv.file)) {
      #v$file = csv.file
      load.data(csv.file)
    }

    # read data, generate funnel plot
    output$funnel=renderPlot({
      req(v$X)
      metastudies_funnel_plot(v$X,v$sigma)
    })

    # generate histogram
    output$hist=renderPlot({
      req(v$X)
      z_histogram(v$X, v$sigma)
    })

    #estimation
    observeEvent(input$estimatebutton,{
      v$cutoffs=as.numeric(unlist(input$cutoffs))
      v$symmetric=input$symmetric
      v$modelmu=input$modelmu
      if (!v$symmetric) v$cutoffs= c(-rev(v$cutoffs), 0 ,v$cutoffs)

      v$ms=metastudies_estimation(v$X,v$sigma,v$cutoffs, v$symmetric, model= v$modelmu)

      output$estplot = renderPlot({
          estimates_plot(v$ms)
      })
    })

    #render estimates to table
    output$estimatestable =  renderTable(rownames=TRUE,hover=TRUE,digits=3,
     {req(v$ms)
      v$ms$est_tab
     })
  }
  shinyApp(ui = ui, server = server,...)
}

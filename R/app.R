#' Create and run the MetaStudies shiny app
#'
#' @param csv.file if not NULL a file that will be initially used by the app
#' @param show.cor if TRUE add a tabPanel with correlations between estimates and standard
#'                 errors to get better insights whether the independence assumption
#'                 of Andrews & Kasy (2019) may be violated in the given data set.
#' @param ... additional parameters passed to \code{shinyApp}
MetaStudiesApp = function(csv.file=NULL,show.cor = TRUE,...) {
  res.ui = tagList(
    h2("Model estimates", align = "center"),
    h4("Distribution of true effects, conditional publication proabilities", align = "center"),
    column(12, align="center",
           tableOutput("estimatestable"),
           plotOutput("estplot", width = "70%")
    )
  )

  if (show.cor) {
    lower.ui = tabsetPanel(
      tabPanel("Results",res.ui),
      tabPanel("Specification tests",
        h3("Correlations between logs of estimate and standard deviation"),
        tableOutput("cortable"),
        helpText("
Explanation: One crucial assumption of the Andrews and Kasy (2019) approach is that in the latent distribution without publication bias the estimate and its standard error are statistically independent across tests.

While the latent distribution cannot be observed, the table above shows
some correlations for which absolute values not close to zero may indicate problems with respect to this assumption.

The first correlation uses an inverse probability weighting approach.
More precisely, we compute a weighted correlation, weighting each observation inversely with its estimated publication probability. Assuming all
assumptions are satisfied, this should yield a consistent estimator of the correlation in the unobserved latent distribution.

Later rows show the correlations between the estimate and standard
errors separately for each interval inside which a constant publication
probability is assumed using directly the observed data without inverse probability weighting. While those are not formal tests, it would be reassuring if these correlations are close to zero.")
      )
    )
  } else {
    lower.ui = tagList(hr(), res.ui)
  }


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
    lower.ui
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
      req(input$file1)
      file = input$file1$datapath
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

      if (show.cor) {
        v$cors = metastudy_X_sigma_cors(v$ms)
      }

      output$estplot = renderPlot({
          estimates_plot(v$ms)
      })
    })

    #render estimates to table
    output$estimatestable =  renderTable(rownames=TRUE,hover=TRUE,digits=3,
     {req(v$ms)
      v$ms$est_tab
     })

    if (show.cor) {
      # show correlations for specification tests
      output$cortable =  renderTable(rownames=TRUE,hover=TRUE,digits=3, striped=TRUE,{
        req(v$cors)
        tab = v$cors %>%
          filter(trans=="log") %>%
          mutate(
            cor = round(cor,3),
            label = case_when(
              mode=="ipv" ~ "Estimated latent distribution",
              TRUE ~ paste0("Observed distribution in range ",mode)
            ),
            ci = case_when(
              mode == "ipv" ~ "Not computed",
              TRUE ~ paste0("[",round(conf.cor.low,3),", ",round(conf.cor.up,3),"]")
            )
          ) %>%
          select(`Type`=label, `Correlation`=cor,`95% CI`=ci)
        tab
      })
    }

  }
  shinyApp(ui = ui, server = server,...)
}

#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(R6)
library(cubelyr)
library(stringr)

tidyMod <- read.csv("ModuleEigenproteins_tidy.csv")

StatsTable <- R6Class("StatsTable", list(

    ExperimentalGroups = NULL,
    RevGroups = NULL,
    EigenModules = NULL,
    StatsResults = NULL,
    StatsResultsNames = NULL,
    
    # Formats group names for CheckBox, removing underscores and Converting to Title Case
    SetExperimentalGroups = function(groups) {
        stopifnot(is.character(groups))
        groups <- sapply(groups, FUN = function(x) tools::toTitleCase( gsub('_', ' ', x) ) )
        self$ExperimentalGroups <- groups
        invisible(self)
    },
    
    #TODO: Exception catching
    GetGroupKey = function(name) {
        names(self$ExperimentalGroups)[which(self$ExperimentalGroups == name)]
    },
    
    SetEigenModules = function(modules) {
        stopifnot(is.character(modules))
        self$EigenModules <- modules
        invisible(self)
    },
    
    ConstructResultsTable = function() {
        stopifnot(length(self$ExperimentalGroups) > 0)
        groups <- names(self$ExperimentalGroups)
        modules <- self$EigenModules
        
        self$StatsResults <- array(
            data = NA,
            dim = c( length(groups), length(groups), 5,  length(modules)),
            dimnames = list(
                Compare1 = groups,
                Compare2 = groups,
                Test = c("ANOVA", "P1TT", "P2TT", "U1TT", "U2TT"),
                EigenModule = modules
            )
        )
        invisible(self)
    },
    
    ConstructResultsNames = function(){
        self$StatsResultsNames <- as.data.frame(as.tbl_cube(self$StatsResults))
        invisible(self)
    },
    
    GetResultsNames = function(idx){
        resultsFrame = self$StatsResultsNames[idx, ]
        names(resultsFrame)[5] <- "Results"
        resultsFrame
    },
    
    DisplayGroupNames = function(){
        unname(self$ExperimentalGroups)
    },
    
    
    PerformTTest = function(test = "P2TT", compare1, compare2, module) {
        stopifnot(is.character(compare1), is.character(compare2), is.character(module), is.character(test))
        mep1 = tidyMod[tidyMod$ExperimentalGroup == self$GetGroupKey(compare1) & tidyMod$Module == module, "Eigenprotein"]
        mep2 = tidyMod[tidyMod$ExperimentalGroup == self$GetGroupKey(compare2) & tidyMod$Module == module, "Eigenprotein"]
        #TODO: Add error checking to prevent paired testing for unpaired groups
        #      Honestly probably need to do something else to ensure that paired testing isn't ran erroneously
        #      Maybe Only allow one selection for unpaired T Tests
        if(mean(mep2) >= mean(mep1)) { one.side.hypothesis = "greater" } else { one.side.hypothesis = "less" }
        TTestArgumentKey = list(
            P1TT = list(x = mep2, y = mep1, alternative = one.side.hypothesis, paired = TRUE),
            P2TT = list(x = mep2, y = mep1, alternative = "t", paired = TRUE),
            U1TT = list(x = mep2, y = mep1, alternative = one.side.hypothesis, paired = FALSE),
            U2TT = list(x = mep2, y = mep1, alternative = "t", paired = FALSE)
        )
        self$StatsResults[ self$GetGroupKey(compare1), self$GetGroupKey(compare1), test, module] <- do.call(t.test, TTestArgumentKey[[test]])$p.value
        base::print(paste(compare1, compare2, test, module))
        self$ConstructResultsNames() #This is an inefficient way to update the names table
        invisible(self)
    },
    
    initialize = function(groups, modules) {
        self$SetExperimentalGroups(groups)
        self$SetEigenModules(modules)
        self$ConstructResultsTable()
        self$ConstructResultsNames()
    },

    print = function() {
        idx <- which(!is.na(self$StatsResults))
        if(length(idx) < 1) {
            cat("No statistical tests have been run", "\n")
        } else {
            results <- self$GetResultsNames(idx)
            results$Results <- round(results$Results, 4)
            results$Compare1 <- self$ExperimentalGroups[results$Compare1]
            results$Compare2 <- self$ExperimentalGroups[results$Compare2]
            res.names <- as.character(names(results))
            cat(sapply(names(results), FUN = function(x) str_pad(x, width = 20)), "\n")
            # TODO: Convert this nested loop to use sapply
            for ( i in seq_len( nrow(results) ) ) {
                for ( y in seq_len( ncol(results) ) ) {
                    cat(
                        str_pad(
                            results[i, y],
                            width = 20,
                            side = "left" 
                            )
                        )
                }
                cat('\n')
            }
        }
        invisible(self)
    },
    
    printPoorly = function() {
        idx <- which(!is.na(self$StatsResults))
        if(length(idx) < 1) {
            output <- paste("No statistical tests have been run", "\n")
        } else {
            results <- self$GetResultsNames(idx)
            results$Results <- round(results$Results, 4)
            results$Compare1 <- self$ExperimentalGroups[results$Compare1]
            results$Compare2 <- self$ExperimentalGroups[results$Compare2]
            res.names <- as.character(names(results))
            output <- paste(sapply(names(results), FUN = function(x) str_pad(x, width = 20)), "\n")
            # TODO: Convert this nested loop to use sapply
            for ( i in seq_len( nrow(results) ) ) {
                for ( y in seq_len( ncol(results) ) ) {
                    output <- paste(output, 
                        str_pad(
                            results[i, y],
                            width = 20,
                            side = "left" 
                        )
                    )
                }
                output <- paste(output, '\n')
            }
        }
        output
    }
    
    )
)

groups <- unique(tidyMod$ExperimentalGroup)
module <- unique(tidyMod$Module)

STT <- StatsTable$new(groups = unique(tidyMod$ExperimentalGroup),
                             modules = unique(tidyMod$Module))

# STT$PerformTTest("U2TT", compare1 = "Dia Low 5min", compare2 = "Nondia Low 5min", module = "MEblack")
# STT$GetGroupKey("Dia Low 5min")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30),
            
            # The list names are displayed in the UI, the list entries are internal
            #TODO: Need content-aware column sizing
            fluidRow(
                column(5,
                       checkboxGroupInput("compare1",
                                          "Comparison Conditions",
                       )
                        
                ), 
                column(5, offset = 1,
                       checkboxGroupInput("compare2",
                                          "Comparator Conditions"
                       )
                )
            ),
            
            selectInput("statsTest",
                        "Statistical Test to Perform",
                        choices = list(
                            ANOVA = "anova",
                            `T-Test: Paired, One-Tailed` = "P1TT",
                            `T-Test: Paired, Two-Tailed` = "P2TT",
                            `T-Test: Unpaired, One-Tailed` = "U1TT",
                            `T-Test: Unpaired, Two-Tailed` = "U2TT"
                        )
            ),

            textInput("tinput",
                      "Change checkbox options",
                      value = NULL),
            
            actionButton("runTest", "Run Test")
        ),
        
        

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           textOutput("testText")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
    
    observe({
        
        STT$SetExperimentalGroups(groups = unique(tidyMod$ExperimentalGroup))
        
        updateCheckboxGroupInput(session = getDefaultReactiveDomain(),
                                 inputId = "compare1",
                                 label = "Comparison Conditions",
                                 choices = unname(STT$ExperimentalGroups))
        
        updateCheckboxGroupInput(session = getDefaultReactiveDomain(),
                                 inputId = "compare2",
                                 label = "Comparator Conditions",
                                 choices = unname(STT$ExperimentalGroups))
        
    })
    
    observe({
        input$runTest
        req(input$runTest > 0)
        req(!is.null(input$statsTest))
        req(!is.null(input$compare1))
        req(!is.null(input$compare2))
        
        STT$PerformTTest(
            test = input$statsTest,
            compare1 = input$compare1,
            compare2 = input$compare2,
            module = "MEblack"
            )
    })
    
    StatsTable.reactivePrint <- eventReactive(input$runTest, {
        string_out <- STT$printPoorly() 
    })
                 
    
    output$testText <- renderText({
        StatsTable.reactivePrint()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

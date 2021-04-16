#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)

library(BiocManager)
options(repos = BiocManager::repositories())

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    # Application title
    titlePanel("GPICS120 prediction by Ju-Seog's Lab"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            fileInput(
                'geneExprfile',
                h4('Choose Gene expression file(txt/csv)'),
                accept = c('text/csv',
                           'text/comma-separated-values,text/plain',
                           '.csv')
            ),
            radioButtons(
                'standardizationType',
                'Standardization',
                c(
                    'Mendian-centering only' = 'medianCenter',
                    'Median-centering and dividing by SD' = 'devidedBySD'
                ),
                'medianCenter'
            ),
            actionButton("doPrediction", "Prediction", class = "btn-primary"),
            br(),
            br(),
            downloadButton('downloadResults', 'Download result table')
        ),
        
        # Show a plot of the generated distribution
        mainPanel(tabsetPanel(
            type = "tabs",
            tabPanel("How to use",
                     HTML("&nbsp; <p>Hi. This is the prediction tool for the analysis of Gastric cancer subtype using mRNA expression data.</p><p>Just upload your dataset. And press the prediction button. That's all."),
                     HTML("</p><p>&nbsp;</p><p><b>Step 1. Prepare of dataset.</b></p>"
                     ),
                     img(
                         src = "testDatasetExample.png",
                         width = 500,
                         height = 150
                     ),
                     HTML(
                         "<p> &nbsp;</p><p> The first line contains the labels Name(<em>HUGO Gene Nomenclature</em>) followed by the identifiers for each sample in the dataset.The dataset is the gene-level transcription estimates, as in log2(x+1) transformed normalized count. &nbsp;</p><p>&nbsp;</p><p><b>Step 2. Standardization. </b> &nbsp;</p><p> Select the data standardization method. &nbsp;</p><p><b>Step 3. Prediction.</b> &nbsp;</p><p> Press the predictiop. &nbsp;</p><p><b>Step 4. Check out the results.</b> &nbsp;</p><p>After analysis, You can find the results at the result tab. The results of dataset could be downloaded using the download button.</p>"
                     )
                     
            ),
            tabPanel("Analysis Graph",
                     plotlyOutput("resultPiePlot"),
                     plotlyOutput("resultSummaryPlot")
            ),
            tabPanel("Analysis table",
                     tableOutput("tablesTemp"),
            )
            
        ))
    )
))

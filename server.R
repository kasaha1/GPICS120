

library(tidyverse)
library(data.table)
library(shiny)
library(plotly)
library(heatmaply)
library(shinycssloaders)
library(HGNChelper)

library(BiocManager)
options(repos = BiocManager::repositories())
library(pamr)
library(impute)

source("functions.R")

trainingGeneMatrix <-
  fread("./www/trainingMatrix.txt") %>% as.data.frame()
trainingClass <-
  fread("./www/trainingClass.txt") %>% as.data.frame()
trainingGenenames <-
  fread("./www/geneNames.txt") %>% as.data.frame()

trainingX <- trainingGeneMatrix %>% data.matrix() %>% impute.knn()
trainingX <- trainingX$data

Genenames_precheck <-
  trainingGenenames %>% t %>% as.vector() %>% checkGeneSymbols
Genenames <- Genenames_precheck[3] %>% t %>% as.vector()

## transform 'class name' to numeric vector
levOfClass <- trainingClass %>% t %>% as.factor() %>% levels()
levOfClass_numeric <- levOfClass %>% as.factor() %>% as.numeric()
table.levels <-
  cbind(levOfClass, levOfClass_numeric) %>% as.data.frame()
trainingClass.m <-
  left_join(trainingClass, table.levels, by = c("class" = "levOfClass"))
trainingY <-
  trainingClass.m$levOfClass_numeric %>% t %>% as.character() %>% as.numeric()

## merging dataset for training
mydata <-
  list(
    x = trainingX,
    y = trainingY,
    genenames = Genenames,
    geneid = c(1:length(Genenames))
  )

## analysis start : Training
model <- pamr.train(mydata)

## analysis start : CrossValidation

model.cv <- pamr.cv(fit = model, data = mydata)

## Threshold 1.628 by analysis
Delta <- 0

## analysis start : centroid
centroid_gene <-
  pamr.listgenes(model, mydata, Delta, genenames = TRUE) %>% as.data.frame
centroid_gene$id <- centroid_gene$id %>% as.numeric()
centroid_gene <- centroid_gene %>% arrange(id)
colnames(centroid_gene)[c(3:(2 + length(levOfClass)))] <- levOfClass

# Modeling finish

# server functions----
shinyServer(function(input, output) {
  preparationText <- "Upload your dataset into the box on the left."
  
  # reactive function ----
  
  
  
  
  geneExprDataIn <- reactive({
    inFile <- input$geneExprfile
    req(inFile)
    f <-
      fread(inFile$datapath) %>% as.data.frame() %>% kasa.duplicationRemovalBySD()
    
    data.raw <- f
    
    HGNCcheck <- data.raw[1] %>% t() %>% as.vector()
    checkedSymbols <- checkGeneSymbols(HGNCcheck)
    data.raw[1] <- checkedSymbols[3]
    HGNCcheckdata_FALSE <-
      checkedSymbols %>% filter(Approved == "FALSE")
    HGNCcheckdata_FALSE <- HGNCcheckdata_FALSE[c(1, 3)]
    colnames(HGNCcheckdata_FALSE) <-
      c("YourGene.Symbols", "HGNCsymbols.Converted")
    # print(HGNCcheckdata_FALSE)
    output$tablesConvertedGeneSymbols <-
      renderTable(HGNCcheckdata_FALSE)
    colnames(data.raw)[1] <- c("genenames")
    
    
    testDataset.modi <-
      left_join(x = trainingGenenames,
                y = data.raw,
                by = c("genenames"))
    MissingValuseGenes <-
      testDataset.modi[!complete.cases(testDataset.modi), ]
    colnames(MissingValuseGenes)[1] <-
      c("Your symbols including Missing values")
    output$MissingValuesSymbols <-
      renderTable(MissingValuseGenes[1])
    
    
    return(testDataset.modi)
  })
  
  reactiveDataStandardization <- reactive({
    
    # standardization
    if (input$standardizationType == "NoneS") {
      data.raw <- geneExprDataIn() 
      # } else if (input$standardizationType == "medianCenter") {
      #     data.raw <- geneExprDataIn() %>% kasa.geneMedianCentering()
    } else if (input$standardizationType == "devidedBySD")
    {
      data.raw <-
        geneExprDataIn() %>% kasa.geneMedianCentering() %>% kasa.geneStandardization()
      
    }
    
    # impute
    if (ncol(data.raw) < 3) {
      rowmax.value <- 1
    } else {
      rowmax.value <- 0.5
    }  # rowmax decision by the number of samples
    
    testX.p <- data.raw[-1] %>% data.matrix() %>% impute.knn(rowmax = rowmax.value)
    testX <- testX.p$data
    
    return(testX)
  })
  
  reavticResultSummaryPlot <- reactive({
    dataPlotly <- DoPIC100prediction()
    fig <-
      plot_ly(
        dataPlotly,
        y = ~ posterior,
        x =  ~ Class,
        color =  ~ Class,
        text =  ~ paste(Sample),
        type = 'box',
        jitter = 0.3,
        boxpoints = 'all'
      )
    return(fig)
    
  })
  reavticResultPiePlot <- reactive({
    dataPlotly <- DoPIC100prediction()
    fig <- plot_ly(dataPlotly, labels = ~ Class, type = 'pie')
    return(fig)
    
  })
  reavticResultHeatmapPlot <- reactive({
    dataPlotly <- DoPIC100prediction()
    dataRaw <- reactiveDataStandardization()
    
    data.sample.order <-
      dataPlotly %>% arrange(Class, posterior)
    data.sample.order.vector <-
      data.sample.order[, 1] %>% as.character()
    data.sample.order.labels <-
      data.sample.order %>% mutate(labels = paste(Class, ":", Sample))
    data.sample.order.labels <-
      data.sample.order.labels$labels %>% as.vector()
    
    data.heatmap <- dataRaw
    
    data.heatmap <-
      data.heatmap[rev(1:nrow(data.heatmap)), data.sample.order.vector]
    
    
    fig <-
      plot_ly(
        z = data.heatmap,
        type = "heatmap",
        colors = colorRamp(c(
          "#009900", "#00FF00", "black", "#FF0000", "#990000"
        )),
        y = rev(Genenames),
        x = data.sample.order.labels,
        zauto = FALSE,
        zmin = -4,
        zmax = 4
      ) %>% layout(
        title = "The Heatmap of your dataset",
        font = list(family = "Arial"),
        xaxis = list(tickangle = 45)
      )
    
    return(fig)
  })
  # event reactive ----
  DoPIC100prediction <- eventReactive(input$doPrediction, {
    testX <- reactiveDataStandardization()
    
    preparationText <<- ""
    output$preparation <- renderText(preparationText)
    output$preparation2 <- renderText(preparationText)
    
    print("prediction start")
    
    res.class <-
      pamr.predict(
        fit = model,
        newx = testX,
        threshold = Delta,
        type = "class"
      ) %>% as.character() %>% as.numeric()
    res.class.t <- levOfClass[res.class]
    
    res.probability <-
      pamr.predict(
        fit = model,
        newx = testX,
        threshold = Delta,
        type = "posterior"
      )
    res.probability.m <- res.probability %>% round(digits = 3)
    res.probability.t <-
      apply(res.probability, 1, max) %>% round(digits = 3)
    
    res.ID <- colnames(testX)
    
    res.table <-
      cbind(res.ID,
            res.class.t,
            res.probability.t,
            res.probability.m) %>% as.data.frame()
    colnames(res.table) <-
      c(
        "Sample",
        "Class",
        "posterior",
        "post.CGS1",
        "post.CGS2",
        "post.CGS3",
        "post.CGS4",
        "post.CGS5",
        "post.CGS6"
      )
    return(res.table)
  })
  
  # output
  output$downloadResults <- downloadHandler(
    filename = function() {
      "Result.txt"
    },
    content = function(file) {
      contents.table <- DoPIC100prediction()
      write_delim(contents.table, file, delim = "\t", na = "")
    },
    contentType = "text/plain"
  )
  output$downloadExample <- downloadHandler(
    filename = function() {
      "testMatrix.txt"
    },
    content = function(file) {
      contents.table <- fread("./www/testDataset/testMatrix.txt")
      write_delim(contents.table, file, delim = "\t", na = "")
    },
    contentType = "text/plain"
  )
  # ouput rendering
  output$tablesTemp <- renderTable(DoPIC100prediction())
  output$resultSummaryPlot <-
    renderPlotly(reavticResultSummaryPlot())
  output$resultPiePlot <- renderPlotly(reavticResultPiePlot())
  output$resultHeatmapPlot <-
    renderPlotly(reavticResultHeatmapPlot())
  output$preparation <- renderText(preparationText)
  output$preparation2 <- renderText(preparationText)
  
})

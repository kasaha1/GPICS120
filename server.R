
library(tidyverse)
library(data.table)
library(shiny)
library(plotly)
library(heatmaply)

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

Genenames <- trainingGenenames %>% t %>% as.vector()

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

## Threshold 0 by analysis
Delta <- 0

## analysis start : centroid
centroid_gene <-
    pamr.listgenes(model, mydata, Delta, genenames = TRUE) %>% as.data.frame
centroid_gene$id <- centroid_gene$id %>% as.numeric()
centroid_gene <- centroid_gene %>% arrange(id)
colnames(centroid_gene)[c(3:(2 + length(levOfClass)))] <- levOfClass

# Modeling finish

shinyServer(function(input, output) {
    # reactive function
    geneExprDataIn <- reactive({
        inFile <- input$geneExprfile
        req(inFile)
        f <-
            fread(inFile$datapath) %>% as.data.frame() %>% kasa.duplicationRemovalBySD()
        
        data.raw <- f
        colnames(data.raw)[1] <- c("genenames")
        testDataset.modi <-
            left_join(x = trainingGenenames,
                      y = data.raw,
                      by = c("genenames"))
       
        # status <<- "File upload complete. The next step for analysis is ready."
        # output$status <- renderText(status)
        return(testDataset.modi)
    })
    
    reactiveDataStandardization <- reactive({
        if (input$standardizationType == "medianCenter") {
            data.raw <- geneExprDataIn() %>% kasa.geneMedianCentering()
        } else {
            
            data.raw <-
                geneExprDataIn() %>% kasa.geneMedianCentering() %>% kasa.geneStandardization()
           
        }
        testX.p <- data.raw[-1] %>% data.matrix() %>% impute.knn()
        testX <- testX.p$data
        
        return(testX)
    })
    
    reavticResultSummaryPlot <- reactive({
        dataPlotly <- DoPIC100prediction()
        fig <- plot_ly(dataPlotly,y = ~posterior, x=~Class, color=~Class, text=~paste(Sample), type = 'box', jitter=0.3,boxpoints = 'all')
        return(fig)
        
    })
    reavticResultPiePlot <- reactive({
        dataPlotly <- DoPIC100prediction()
        fig <- plot_ly(dataPlotly, labels = ~Class, type = 'pie')
        return(fig)
        
    })
    # event reactive
    DoPIC100prediction <- eventReactive(input$doPrediction, {
        output$preparation <- renderText("")
        output$preparation2 <- renderText("")
        
        testX <- reactiveDataStandardization()
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
                "post.A",
                "post.B",
                "post.C",
                "post.D",
                "post.E",
                "post.F"
            )
        return(res.table)
    })
    reavticResultHeatmapPlot <- reactive({
        dataPlotly <- DoPIC100prediction()
        dataRaw <- reactiveDataStandardization()
        
        data.sample.order <- dataPlotly %>% arrange(Class,posterior) 
        data.sample.order.vector <- data.sample.order[,1] %>% as.character()
        data.sample.order.labels <- data.sample.order %>% mutate(labels=paste(Class,":",Sample))
        data.sample.order.labels <- data.sample.order.labels$labels %>% as.vector()
        
        data.heatmap <- dataRaw 
        
        data.heatmap <- data.heatmap[rev(1:nrow(data.heatmap)),data.sample.order.vector]
        
        
        fig <- plot_ly(z = data.heatmap, type = "heatmap",colors = colorRamp(c("#009900","#00FF00","black","#FF0000","#990000")),y=rev(Genenames),x=data.sample.order.labels, zauto = FALSE, zmin = -4, zmax = 4) %>% layout(title="The Heatmap of your dataset",font=list(family="Arial"),xaxis=list(tickangle=45))
        
        return(fig)
    })
    # output
    output$downloadResults <- downloadHandler(
        filename = function() {
            "Result.txt"
        },
        content = function(file) {
            contents.table <- DoPIC100prediction()
            write_delim(contents.table, file, delim = "\t",na = "")
        },
        contentType = "text/plain"
    )
    # ouput rendering
    output$tablesTemp <- renderTable(DoPIC100prediction())
    output$resultSummaryPlot <- renderPlotly(reavticResultSummaryPlot())
    output$resultPiePlot <- renderPlotly(reavticResultPiePlot())
    output$resultHeatmapPlot <- renderPlotly(reavticResultHeatmapPlot())
    output$preparation <- renderText("Insert Your Dataset")
    output$preparation2 <- renderText("Insert Your Dataset")
})

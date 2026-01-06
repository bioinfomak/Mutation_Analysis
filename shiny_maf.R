# ALS TGS – Shiny App for MAF Exploration
# ----------------------------------
# Save this as app.R and run shiny::runApp()

suppressPackageStartupMessages({
  library(shiny)
  library(maftools)
  library(data.table)
  library(ggplot2)
})

# ---- Load merged MAF ----
maf <- read.maf("FILE_NAME")

als_genes <- c(
  "SOD1","C9orf72","TARDBP","FUS","OPTN","TBK1",
  "NEK1","SETX","VCP","UBQLN2","ATXN2","ALS2"
)

# ---- UI ----
ui <- fluidPage(
  titlePanel("ALS Targeted Sequencing – MAF Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "view",
        "Select Analysis",
        choices = c(
          "Overview",
          "Oncoplot",
          "ALS Gene Panel",
          "Single Gene",
          "TMB",
          "Ti/Tv",
          "Data Table"
        )
      ),
      
      conditionalPanel(
        condition = "input.view == 'Single Gene'",
        selectInput("gene", "Select Gene", choices = sort(unique(maf@gene.summary$Hugo_Symbol)))
      )
    ),
    
    mainPanel(
      plotOutput("mainPlot", height = "600px"),
      tableOutput("summaryTable")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  output$mainPlot <- renderPlot({
    
    if (input$view == "Overview") {
      plotmafSummary(
        maf = maf,
        rmOutlier = TRUE,
        addStat = "median",
        dashboard = TRUE
      )
      
    } else if (input$view == "Oncoplot") {
      oncoplot(
        maf = maf,
        top = 25,
        showTumorSampleBarcodes = TRUE,
        removeNonMutated = TRUE
      )
      
    } else if (input$view == "ALS Gene Panel") {
      maf_als <- subsetMaf(maf = maf, genes = als_genes, mafObj = TRUE)
      oncoplot(maf_als, top = length(als_genes))
      
    } else if (input$view == "Single Gene") {
      lollipopPlot(
        maf = maf,
        gene = input$gene,
        AACol = "HGVSp_Short"
      )
      
    } else if (input$view == "TMB") {
      tmb_res <- tmb(maf)
      boxplot(
        tmb_res$total_perMB,
        main = "Tumor Mutation Burden",
        ylab = "Mutations / Mb"
      )
      
    } else if (input$view == "Ti/Tv") {
      titv_res <- titv(maf)
      plotTiTv(titv_res)
    }
  })
  
  output$summaryTable <- renderTable({
    if (input$view == "Data Table") {
      head(maf@data, 50)
    } else if (input$view == "Overview") {
      getSampleSummary(maf)
    } else {
      NULL
    }
  })
}

# ---- Run app ----
shinyApp(ui = ui, server = server)

############################################################
# ALS Shiny MAF Explorer – FINAL CLEAN VERSION
# Targeted + WES | SALS vs FALS | With p-values
############################################################

suppressPackageStartupMessages({
  library(shiny)
  library(maftools)
  library(data.table)
  library(ggplot2)
})

# ===============================
# LOAD MAF FILES
# ===============================
maf_T_M_SALS <- read.maf("PDXRL_M_SALS_merged.maf", isTCGA = FALSE)
maf_T_F_SALS <- read.maf("PDXRL_F_SALS_merged.maf", isTCGA = FALSE)
maf_T_FALS   <- read.maf("PDXRL_FALS_merged.maf",   isTCGA = FALSE)
maf_NORMAL   <- read.maf("PDXRL_NORMAL_merged.maf", isTCGA = FALSE)

maf_W_M_SALS <- read.maf("WES_M_SALS_merged.maf", isTCGA = FALSE)
maf_W_F_SALS <- read.maf("WES_F_SALS_merged.maf", isTCGA = FALSE)
maf_W_FALS   <- read.maf("WES_FALS_merged.maf",   isTCGA = FALSE)

# ===============================
# ALS GENE PANEL
# ===============================
als_genes <- c(
  "SOD1","C9orf72","TARDBP","FUS","OPTN","TBK1",
  "NEK1","SETX","VCP","UBQLN2","ATXN2","ALS2"
)

# ===============================
# UI
# ===============================
ui <- fluidPage(
  
  titlePanel("ALS MAF Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput(
        "cohort", "Select Cohort",
        c(
          "Targeted M-SALS" = "T_M_SALS",
          "Targeted F-SALS" = "T_F_SALS",
          "Targeted FALS"   = "T_FALS",
          "Normal"          = "NORMAL",
          "WES M-SALS"      = "W_M_SALS",
          "WES F-SALS"      = "W_F_SALS",
          "WES FALS"        = "W_FALS"
        )
      ),
      
      selectInput(
        "view", "Select Analysis",
        c(
          "Overview",
          "Oncoplot",
          "ALS Gene Panel",
          "Single Gene",
          "Somatic Interactions",
          
          "Targeted M-SALS vs FALS",
          "Targeted F-SALS vs FALS",
          
          "Targeted M-SALS vs Normal",
          "Targeted F-SALS vs Normal",
          "Targeted FALS vs Normal",
          
          "WES M-SALS vs Normal",
          "WES F-SALS vs Normal",
          "WES FALS vs Normal",
          
          "Targeted vs WES M-SALS",
          "Targeted vs WES F-SALS",
          "Targeted vs WES FALS",
          
          "TMB",
          "Ti/Tv",
          
          "Gene Mutation Frequency",
          "Data Table"
        )
      ),
      
      conditionalPanel(
        condition = "input.view == 'Single Gene'",
        uiOutput("gene_ui")
      ),
      
      hr(),
      downloadButton("download_plot", "Download Plot"),
      downloadButton("download_table", "Download Table")
    ),
    
    mainPanel(
      plotOutput("mainPlot", height = "650px"),
      tableOutput("summaryTable")
    )
  )
)

# ===============================
# SERVER
# ===============================
server <- function(input, output, session) {
  
  maf_active <- reactive({
    switch(input$cohort,
           "T_M_SALS" = maf_T_M_SALS,
           "T_F_SALS" = maf_T_F_SALS,
           "T_FALS"   = maf_T_FALS,
           "NORMAL"   = maf_NORMAL,
           "W_M_SALS" = maf_W_M_SALS,
           "W_F_SALS" = maf_W_F_SALS,
           "W_FALS"   = maf_W_FALS
    )
  })
  
  output$gene_ui <- renderUI({
    selectInput(
      "gene", "Select Gene",
      sort(unique(maf_active()@gene.summary$Hugo_Symbol))
    )
  })
  
  # ===============================
  # SAFE PLOT FUNCTION
  # ===============================
  make_plot <- function(view, maf, gene = NULL) {
    
    if (view == "Overview") {
      plotmafSummary(maf, dashboard = TRUE)
      
    } else if (view == "Oncoplot") {
      oncoplot(maf, top = 25)
      
    } else if (view == "ALS Gene Panel") {
      g <- intersect(als_genes, maf@gene.summary$Hugo_Symbol)
      if (length(g) == 0) {
        plot.new()
        text(0.5, 0.5, "No ALS genes detected")
        return()
      }
      oncoplot(subsetMaf(maf, genes = g), top = length(g))
      
    } else if (view == "Single Gene") {
      lollipopPlot(maf, gene = gene)
      
    } else if (view == "Somatic Interactions") {
      somaticInteractions(maf, top = 50)
      
      # ===============================
      # mafCompare with p-values
      # ===============================
    } else if (grepl("vs", view)) {
      
      cmp <- switch(view,
                    "Targeted M-SALS vs FALS" =
                      mafCompare(maf_T_M_SALS, maf_T_FALS, "T_M_SALS", "T_FALS", minMut = 5),
                    "Targeted F-SALS vs FALS" =
                      mafCompare(maf_T_F_SALS, maf_T_FALS, "T_F_SALS", "T_FALS", minMut = 5),
                    "Targeted M-SALS vs Normal" =
                      mafCompare(maf_T_M_SALS, maf_NORMAL, "T_M_SALS", "Normal", minMut = 5),
                    "Targeted F-SALS vs Normal" =
                      mafCompare(maf_T_F_SALS, maf_NORMAL, "T_F_SALS", "Normal", minMut = 5),
                    "Targeted FALS vs Normal" =
                      mafCompare(maf_T_FALS, maf_NORMAL, "T_FALS", "Normal", minMut = 5),
                    "WES M-SALS vs Normal" =
                      mafCompare(maf_W_M_SALS, maf_NORMAL, "W_M_SALS", "Normal", minMut = 5),
                    "WES F-SALS vs Normal" =
                      mafCompare(maf_W_F_SALS, maf_NORMAL, "W_F_SALS", "Normal", minMut = 5),
                    "WES FALS vs Normal" =
                      mafCompare(maf_W_FALS, maf_NORMAL, "W_FALS", "Normal", minMut = 5),
                    "Targeted vs WES M-SALS" =
                      mafCompare(maf_T_M_SALS, maf_W_M_SALS, "Targeted", "WES", minMut = 5),
                    "Targeted vs WES F-SALS" =
                      mafCompare(maf_T_F_SALS, maf_W_F_SALS, "Targeted", "WES", minMut = 5),
                    "Targeted vs WES FALS" =
                      mafCompare(maf_T_FALS, maf_W_FALS, "Targeted", "WES", minMut = 5)
      )
      
      forestPlot(cmp, pVal = 0.05)
      
    } else if (view == "TMB") {
      ggplot(tmb(maf), aes(y = total_perMB)) +
        geom_boxplot() +
        ylab("Mutations / Mb") +
        theme_minimal()
      
    } else if (view == "Ti/Tv") {
      plotTiTv(titv(maf))
      
    } else {
      plot.new()
    }
  }
  
  output$mainPlot <- renderPlot({
    make_plot(input$view, maf_active(), input$gene)
  })
  
  output$summaryTable <- renderTable({
    if (input$view == "Gene Mutation Frequency")
      maf_active()@gene.summary[order(-MutatedSamples)]
    else if (input$view == "Data Table")
      head(maf_active()@data, 50)
  })
  
  output$download_plot <- downloadHandler(
    filename = function() paste0(gsub(" ", "_", input$view), ".png"),
    content = function(file) {
      png(file, 2400, 1800, res = 300)
      make_plot(input$view, maf_active(), input$gene)
      dev.off()
    }
  )
  
  output$download_table <- downloadHandler(
    filename = function() paste0(gsub(" ", "_", input$view), ".csv"),
    content = function(file) {
      fwrite(
        if (input$view == "Gene Mutation Frequency")
          maf_active()@gene.summary
        else maf_active()@data,
        file
      )
    }
  )
}

shinyApp(ui, server)

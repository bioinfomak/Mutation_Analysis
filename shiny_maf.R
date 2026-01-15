############################################################
# ALS Shiny MAF Explorer
# Targeted + WES | SALS (M/F), FALS, Normal
# NO COHORT MERGING
############################################################

suppressPackageStartupMessages({
  library(shiny)
  library(maftools)
  library(data.table)
  library(ggplot2)
})

# ===============================
# LOAD TARGETED MAF FILES
# ===============================

maf_T_M_SALS <- read.maf("PDXRL_M_SALS_merged.maf", isTCGA = FALSE)
maf_T_F_SALS <- read.maf("PDXRL_F_SALS_merged.maf", isTCGA = FALSE)
maf_T_FALS   <- read.maf("PDXRL_FALS_merged.maf",   isTCGA = FALSE)
maf_NORMAL   <- read.maf("PDXRL_NORMAL_merged.maf", isTCGA = FALSE)

# ===============================
# LOAD WES MAF FILES
# ===============================

maf_W_M_SALS <- read.maf("WES_M_SALS_merged.maf", isTCGA = FALSE)
maf_W_F_SALS <- read.maf("WES_F_SALS_merged.maf", isTCGA = FALSE)
maf_W_FALS   <- read.maf("WES_FALS_merged.maf",   isTCGA = FALSE)

# ===============================
# ADD SEX ANNOTATION
# ===============================

maf_T_M_SALS@clinical.data$Sex <- "Male"
maf_T_F_SALS@clinical.data$Sex <- "Female"
maf_W_M_SALS@clinical.data$Sex <- "Male"
maf_W_F_SALS@clinical.data$Sex <- "Female"

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
  
  titlePanel("ALS MAF Explorer (Targeted + WES | No SALS Merging)"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput(
        "cohort",
        "Select Cohort",
        choices = c(
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
        "view",
        "Select Analysis",
        choices = c(
          "Overview",
          "Oncoplot",
          "ALS Gene Panel",
          "Single Gene",
          "Somatic Interactions",
          
          "Targeted M-SALS vs Normal",
          "Targeted F-SALS vs Normal",
          "Targeted FALS vs Normal",
          
          "WES M-SALS vs Normal",
          "WES F-SALS vs Normal",
          "WES FALS vs Normal",
          
          "Targeted M-SALS vs WES M-SALS",
          "Targeted F-SALS vs WES F-SALS",
          "Targeted FALS vs WES FALS",
          
          "TMB",
          "Ti/Tv",
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
    switch(
      input$cohort,
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
      "gene",
      "Select Gene",
      choices = sort(unique(maf_active()@gene.summary$Hugo_Symbol))
    )
  })
  
  make_plot <- function(view, maf, gene = NULL) {
    
    if (view == "Overview") {
      plotmafSummary(maf, dashboard = TRUE)
      
    } else if (view == "Oncoplot") {
      oncoplot(maf, top = 25)
      
    } else if (view == "ALS Gene Panel") {
      oncoplot(subsetMaf(maf, als_genes, mafObj = TRUE),
               top = length(als_genes))
      
    } else if (view == "Single Gene") {
      lollipopPlot(maf, gene = gene)
      
    } else if (view == "Somatic Interactions") {
      somaticInteractions(maf, top = 50)
      
    } else if (view == "Targeted M-SALS vs Normal") {
      forestPlot(mafCompare(maf_T_M_SALS, maf_NORMAL,
                            "Targeted M-SALS", "Normal"))
      
    } else if (view == "Targeted F-SALS vs Normal") {
      forestPlot(mafCompare(maf_T_F_SALS, maf_NORMAL,
                            "Targeted F-SALS", "Normal"))
      
    } else if (view == "Targeted FALS vs Normal") {
      forestPlot(mafCompare(maf_T_FALS, maf_NORMAL,
                            "Targeted FALS", "Normal"))
      
    } else if (view == "WES M-SALS vs Normal") {
      forestPlot(mafCompare(maf_W_M_SALS, maf_NORMAL,
                            "WES M-SALS", "Normal"))
      
    } else if (view == "WES F-SALS vs Normal") {
      forestPlot(mafCompare(maf_W_F_SALS, maf_NORMAL,
                            "WES F-SALS", "Normal"))
      
    } else if (view == "WES FALS vs Normal") {
      forestPlot(mafCompare(maf_W_FALS, maf_NORMAL,
                            "WES FALS", "Normal"))
      
    } else if (view == "Targeted M-SALS vs WES M-SALS") {
      forestPlot(mafCompare(maf_T_M_SALS, maf_W_M_SALS,
                            "Targeted M-SALS", "WES M-SALS"))
      
    } else if (view == "Targeted F-SALS vs WES F-SALS") {
      forestPlot(mafCompare(maf_T_F_SALS, maf_W_F_SALS,
                            "Targeted F-SALS", "WES F-SALS"))
      
    } else if (view == "Targeted FALS vs WES FALS") {
      forestPlot(mafCompare(maf_T_FALS, maf_W_FALS,
                            "Targeted FALS", "WES FALS"))
      
    } else if (view == "TMB") {
      ggplot(tmb(maf), aes(x = "", y = total_perMB)) +
        geom_boxplot() +
        ylab("Mutations / Mb") +
        theme_minimal()
      
    } else if (view == "Ti/Tv") {
      plotTiTv(titv(maf))
    }
  }
  
  output$mainPlot <- renderPlot({
    make_plot(input$view, maf_active(), input$gene)
  })
  
  output$summaryTable <- renderTable({
    if (input$view == "Data Table") {
      head(maf_active()@data, 50)
    }
  })
  
  # ===============================
  # DOWNLOAD PLOT
  # ===============================
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(gsub(" ", "_", input$view), "_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 2400, height = 1800, res = 300)
      make_plot(input$view, maf_active(), input$gene)
      dev.off()
    }
  )
  
  # ===============================
  # DOWNLOAD TABLE
  # ===============================
  
  output$download_table <- downloadHandler(
    filename = function() {
      paste0(gsub(" ", "_", input$view), "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      
      df <- switch(
        input$view,
        "Data Table" = maf_active()@data,
        "Overview"   = getSampleSummary(maf_active()),
        NULL
      )
      
      if (!is.null(df)) {
        fwrite(df, file)
      }
    }
  )
  
}

# ===============================
# RUN APP
# ===============================

shinyApp(ui, server)


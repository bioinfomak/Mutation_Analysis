############################################################
# ALS Targeted Sequencing – Shiny MAF Explorer
# Includes Sporadic ALS vs Normal & Familial ALS comparisons
############################################################

suppressPackageStartupMessages({
  library(shiny)
  library(maftools)
  library(data.table)
  library(ggplot2)
})

# ===============================
# LOAD MAF FILES (in app dir)
# ===============================

maf_male <- read.maf(
  "PDXRL_M_SALS_merged.maf",
  removeDuplicatedVariants = TRUE,
  isTCGA = FALSE
)

maf_female <- read.maf(
  "PDXRL_F_SALS_merged.maf",
  removeDuplicatedVariants = TRUE,
  isTCGA = FALSE
)

maf_normal <- read.maf(
  "PDXRL_NORMAL_merged.maf",
  removeDuplicatedVariants = TRUE,
  isTCGA = FALSE
)

maf_FALS <- read.maf(
  "PDXRL_FALS_merged.maf",
  removeDuplicatedVariants = TRUE,
  isTCGA = FALSE
)

# ===============================
# Add Sex annotation to Male/Female MAFs
# ===============================

maf_male@clinical.data$Sex   <- "Male"
maf_female@clinical.data$Sex <- "Female"

# ===============================
# Create merged Sporadic ALS MAF (Male + Female)
# ===============================

maf_SALS <- merge_mafs(list(maf_male, maf_female))

# ===============================
# ALS gene panel
# ===============================

als_genes <- c(
  "SOD1","C9orf72","TARDBP","FUS","OPTN","TBK1",
  "NEK1","SETX","VCP","UBQLN2","ATXN2","ALS2"
)

# ===============================
# UI
# ===============================

ui <- fluidPage(
  
  titlePanel("ALS Targeted Sequencing – MAF Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "cohort",
        "Select Cohort",
        choices = c(
          "Male ALS"   = "M",
          "Female ALS" = "F",
          "Normal"     = "N",
          "Sporadic ALS (M+F)" = "SALS",
          "Familial ALS" = "FALS"
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
          "Male vs Female Oncoplot",
          "Male vs Female Gene Comparison",
          "Sporadic ALS vs Normal Comparison",
          "Familial ALS vs Normal Comparison",
          "Sporadic ALS vs Familial ALS Comparison",
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
      
      downloadButton("download_plot", "Download Figure"),
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
  
  # ---- Select active MAF by cohort input ----
  maf_active <- reactive({
    switch(
      input$cohort,
      "M" = maf_male,
      "F" = maf_female,
      "N" = maf_normal,
      "SALS" = maf_SALS,
      "FALS" = maf_FALS
    )
  })
  
  # ---- Dynamic gene selector for Single Gene view ----
  output$gene_ui <- renderUI({
    selectInput(
      "gene",
      "Select Gene",
      choices = sort(unique(maf_active()@gene.summary$Hugo_Symbol))
    )
  })
  
  # ---- Plotting function for all views ----
  make_plot <- function(view, maf, gene = NULL) {
    
    if (view == "Overview") {
      
      plotmafSummary(
        maf,
        rmOutlier = TRUE,
        addStat = "median",
        dashboard = TRUE
      )
      
    } else if (view == "Oncoplot") {
      
      oncoplot(
        maf,
        top = 25,
        removeNonMutated = TRUE
      )
      
    } else if (view == "ALS Gene Panel") {
      
      maf_als <- subsetMaf(maf, genes = als_genes, mafObj = TRUE)
      
      if (nrow(maf_als@data) == 0) {
        plot.new()
        text(0.5, 0.5, "No ALS genes found", cex = 1.3)
      } else {
        oncoplot(
          maf_als,
          top = length(als_genes),
          removeNonMutated = TRUE
        )
      }
      
    } else if (view == "Single Gene") {
      
      lollipopPlot(
        maf,
        gene = gene,
        AACol = "HGVSp_Short"
      )
      
    } else if (view == "Somatic Interactions") {
      
      somaticInteractions(
        maf,
        top = 50,
        pvalue = c(0.05, 0.01)
      )
      
    } else if (view == "Male vs Female Oncoplot") {
      
      maf_MF <- merge_mafs(list(maf_male, maf_female))
      maf_MF@clinical.data$Sex <- factor(maf_MF@clinical.data$Sex)
      
      oncoplot(
        maf = maf_MF,
        top = 30,
        clinicalFeatures = "Sex",
        sortByAnnotation = TRUE,
        annotationColor = list(Sex = c(Male = "#377EB8", Female = "#E41A1C")),
        removeNonMutated = TRUE
      )
      
    } else if (view == "Male vs Female Gene Comparison") {
      
      comp <- mafCompare(
        m1 = maf_male,
        m2 = maf_female,
        m1Name = "Male",
        m2Name = "Female"
      )
      
      forestPlot(comp)
      
    } else if (view == "Sporadic ALS vs Normal Comparison") {
      
      comp <- mafCompare(
        m1 = maf_SALS,
        m2 = maf_normal,
        m1Name = "Sporadic ALS",
        m2Name = "Normal"
      )
      
      forestPlot(comp)
      
    } else if (view == "Familial ALS vs Normal Comparison") {
      
      comp <- mafCompare(
        m1 = maf_FALS,
        m2 = maf_normal,
        m1Name = "Familial ALS",
        m2Name = "Normal"
      )
      
      forestPlot(comp)
      
    } else if (view == "Sporadic ALS vs Familial ALS Comparison") {
      
      comp <- mafCompare(
        m1 = maf_SALS,
        m2 = maf_FALS,
        m1Name = "Sporadic ALS",
        m2Name = "Familial ALS"
      )
      
      forestPlot(comp)
      
    } else if (view == "TMB") {
      
      tmb_res <- tmb(maf)
      
      ggplot(tmb_res, aes(x = "", y = total_perMB)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.6) +
        ylab("Mutations per Mb") +
        theme_minimal()
      
    } else if (view == "Ti/Tv") {
      
      plotTiTv(titv(maf))
      
    }
  }
  
  # ---- Render main plot ----
  output$mainPlot <- renderPlot({
    make_plot(
      view = input$view,
      maf  = maf_active(),
      gene = input$gene
    )
  })
  
  # ---- Render data table ----
  output$summaryTable <- renderTable({
    
    maf <- maf_active()
    
    if (input$view == "Data Table") {
      
      head(maf@data, 50)
      
    } else if (input$view == "Overview") {
      
      getSampleSummary(maf)
      
    } else if (input$view == "Male vs Female Gene Comparison") {
      
      comp <- mafCompare(
        m1 = maf_male,
        m2 = maf_female,
        m1Name = "Male",
        m2Name = "Female"
      )
      
      comp$results[order(comp$results$pval), ][1:25, ]
      
    } else if (input$view == "Sporadic ALS vs Normal Comparison") {
      
      comp <- mafCompare(
        m1 = maf_SALS,
        m2 = maf_normal,
        m1Name = "Sporadic ALS",
        m2Name = "Normal"
      )
      
      comp$results[order(comp$results$pval), ][1:25, ]
      
    } else if (input$view == "Familial ALS vs Normal Comparison") {
      
      comp <- mafCompare(
        m1 = maf_FALS,
        m2 = maf_normal,
        m1Name = "Familial ALS",
        m2Name = "Normal"
      )
      
      comp$results[order(comp$results$pval), ][1:25, ]
      
    } else if (input$view == "Sporadic ALS vs Familial ALS Comparison") {
      
      comp <- mafCompare(
        m1 = maf_SALS,
        m2 = maf_FALS,
        m1Name = "Sporadic ALS",
        m2Name = "Familial ALS"
      )
      
      comp$results[order(comp$results$pval), ][1:25, ]
      
    } else {
      NULL
    }
  })
  
  # ---- Download plot handler ----
  output$download_plot <- downloadHandler(
    
    filename = function() {
      paste0(
        gsub(" ", "_", input$view),
        "_",
        Sys.Date(),
        ".png"
      )
    },
    
    content = function(file) {
      png(file, width = 2400, height = 1800, res = 300)
      
      make_plot(
        view = input$view,
        maf  = maf_active(),
        gene = input$gene
      )
      
      dev.off()
    }
  )
  
  # ---- Download table handler ----
  output$download_table <- downloadHandler(
    
    filename = function() {
      paste0(
        gsub(" ", "_", input$view),
        "_",
        Sys.Date(),
        ".csv"
      )
    },
    
    content = function(file) {
      
      maf <- maf_active()
      
      df <- switch(
        input$view,
        
        "Data Table" = maf@data,
        
        "Overview" = getSampleSummary(maf),
        
        "Male vs Female Gene Comparison" = {
          comp <- mafCompare(
            m1 = maf_male,
            m2 = maf_female,
            m1Name = "Male",
            m2Name = "Female"
          )
          comp$results
        },
        
        "Sporadic ALS vs Normal Comparison" = {
          comp <- mafCompare(
            m1 = maf_SALS,
            m2 = maf_normal,
            m1Name = "Sporadic ALS",
            m2Name = "Normal"
          )
          comp$results
        },
        
        "Familial ALS vs Normal Comparison" = {
          comp <- mafCompare(
            m1 = maf_FALS,
            m2 = maf_normal,
            m1Name = "Familial ALS",
            m2Name = "Normal"
          )
          comp$results
        },
        
        "Sporadic ALS vs Familial ALS Comparison" = {
          comp <- mafCompare(
            m1 = maf_SALS,
            m2 = maf_FALS,
            m1Name = "Sporadic ALS",
            m2Name = "Familial ALS"
          )
          comp$results
        },
        
        NULL
      )
      
      if (!is.null(df)) {
        fwrite(df, file)
      }
    }
  )
  
}

# ===============================
# Run the Shiny app
# ===============================

shinyApp(ui = ui, server = server)


# ---- Run app ----
shinyApp(ui = ui, server = server)

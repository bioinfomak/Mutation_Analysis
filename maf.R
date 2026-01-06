############################################################
## ALS TGS – Merge & Basic Analysis of MAF files
############################################################

## 1. Load libraries
suppressPackageStartupMessages({
  library(maftools)
  library(data.table)
})

## 2. Set working directory
setwd("PATH")

## 3. List all MAF files
maf_files <- list.files(
  pattern = "\\.maf$",
  full.names = TRUE
)

stopifnot(length(maf_files) > 0)

message("Found ", length(maf_files), " MAF files")

## 4. Read each MAF
maf_list <- lapply(maf_files, function(f) {
  message("Reading: ", basename(f))
  read.maf(f)
})

## Merge all MAFs
merged_maf <- merge_mafs(
  mafs = maf_list
)

## Save merged MAF
write.table(
  merged_maf@data,
  file = "FILE_NAME",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

maf <- read.maf("FILE_NAME")

############################################################
## 7. BASIC QC & SUMMARY
############################################################

## Number of samples
cat("Total samples:", length(getSampleSummary(maf)$Tumor_Sample_Barcode), "\n")

## Top mutated genes
print(getGeneSummary(maf)[1:20, ])

## Mutation classification summary
getSampleSummary(maf)
############################################################
## 8. PLOTS (optional but recommended)
############################################################

## Mutation load per sample
# Overall summary plot
plotmafSummary(
  maf = maf,
  rmOutlier = TRUE,
  addStat = "median",
  dashboard = TRUE
)

oncoplot(
  maf = maf,
  top = 25,
  showTumorSampleBarcodes = TRUE,
  removeNonMutated = TRUE
)

plotmafSummary(
  maf = maf,
  rmOutlier = TRUE,
  addStat = "median"
)

# Variant classification
maf@variant.classification.summary

# Variant types
maf@variant.type.summary


als_genes <- c(
  "SOD1","C9orf72","TARDBP","FUS","OPTN","TBK1",
  "NEK1","SETX","VCP","UBQLN2","ATXN2", "ALS2"
)

subsetMaf(
  maf = maf,
  genes = als_genes,
  mafObj = TRUE
) -> maf_als

oncoplot(maf_als, top = length(als_genes))


lollipopPlot(
  maf = maf,
  gene = "ALS2",
  AACol = "HGVSp_Short"
)

tmb <- tmb(maf)

boxplot(
  tmb$total_perMB,
  main = "Tumor Mutation Burden (per Mb)",
  ylab = "Mutations / Mb"
)

titv <- titv(maf)

plotTiTv(titv)


somaticInteractions(
  maf = maf,
  top = 25,
  pvalue = c(0.05)
)
############################################################
## 9. Save R object for later use
############################################################
saveRDS(merged_maf, file = "FILE_NAME.rds")


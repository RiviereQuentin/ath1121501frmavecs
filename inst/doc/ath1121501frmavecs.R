## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  if (!require("remotes", quietly = TRUE)) {
#    install.packages("remotes")
#  }
#  if (!require("BiocManager", quietly = TRUE)) {
#    install.packages("BiocManager")
#  }
#  options(repos = BiocManager::repositories())
#  getOption("repos")
#  BiocManager::install(
#    "RiviereQuentin/ath1121501frmavecs",
#    dependencies = TRUE,
#    build_vignettes = TRUE,
#    force = TRUE
#  )

## -----------------------------------------------------------------------------
library(ath1121501frmavecs)
library(GEOquery)
library(affy)
library(frma)

# 1. Download the CEL files
GSM_considered <- c("GSM433634", "GSM1179807", "GSM433644")
for (sample in GSM_considered) {
  getGEOSuppFiles(
    sample,
    makeDirectory = FALSE,
    baseDir = getwd(),
    ## Or change to another directory
    filter_regex = "cel.gz || CEL.gz",
    fetch_files = TRUE
  )
}

# 2. Obtain the CEL file paths
celpaths <-
  grepl(
    pattern = paste0(c(
      ".*(", paste0(GSM_considered, collapse = "|"), ").*"
    ), collapse = ""),
    x = list.files(getwd()),
    ## Replace getwd() by baseDir
    ignore.case = TRUE
  )
celpaths <- list.files(getwd())[celpaths]
celpaths <-
  celpaths[grepl(pattern = "cel.gz", celpaths, ignore.case = TRUE)]

# 3. Load and preprocess the data
cel_batch <- read.affybatch(celpaths)
data(ath1121501frmavecs)
cel_frma <-
  frma(
    object = cel_batch,
    target = "full",
    input.vec = ath1121501frmavecs,
    verbose = TRUE
  )

# 4. Extract the fRMA values
cel_e <- exprs(cel_frma)
cel_e <- as.data.frame(cel_e)
head(cel_e)

## -----------------------------------------------------------------------------
library(ath1121501.db)
probeset2gene <- mapIds(
  ath1121501.db,
  keys = rownames(cel_e),
  column = "TAIR",
  keytype = "PROBEID"
)
#cel_e <- cel_e[rownames(cel_e) %in% probeset2gene$array_element_name,]
cel_e$AGI <-
  probeset2gene[match(x = rownames(cel_e), table = names(probeset2gene))]
multiple_ps <- table(cel_e$AGI)
multiple_ps <- names(multiple_ps[multiple_ps > 1])
todel_ps <- list()
for (gene in multiple_ps) {
  considered <- cel_e[cel_e$AGI == gene, ]
  mu_considered <-
    apply(as.matrix(considered[, seq(1, ncol(considered) - 1)]), 1, mean)
  mu_considered <-
    mu_considered[order(mu_considered, decreasing = TRUE)]
  todel_ps[[gene]] <- names(mu_considered)[2:length(mu_considered)]
}
todel_ps <- do.call(c, todel_ps)
cel_e <- cel_e[!(rownames(cel_e) %in% todel_ps), ]


# ath1121501frmavecs

## Introduction

The *ath1121501frmavecs* data package allows the implementation of the frozen Robust Multiarray Analysis (fRMA) procedure for *Arabidopsis thaliana* (for microarray data compatible with the GPL198 Gene Expression Omnibus platform). Such a pre-processing method of gene expression arrays was developed by McCall et al. (2010) to ensure the comparability of datasets preprocessed separately. This not only allows to obtain reliable results from small batches of samples but also to reduce the computational burden to perform meta-analysis and save time when including new samples to a database.

The fRMA procedure is a platform-specific method, based on the pre-computation of parameter vectors holding the reference normalization distribution, the probe effect estimates, the within batch residual variance, the between batch residual variance, and the within probeset average standard deviation. GPL198 corresponds to the most popular platform for *Arabidopsis thaliana*, hosting data obtained with the Affymetrix Arabidopsis ATH1 Genome Array and annotated with the ATH1-121501 CDF file. Upon recommendation of the authors, *ath1121501frmavecs* was built based on 100 triplicates originating from 100 different series, using the *frmaTools* package (McCall and Irizarry, 2011). Those 100 triplicates were related to a large variety of study subjects, encompassing biotic treatments, abiotic stresses, control conditions (in seedlings), development (different tissues, organs, in control condition) and chemical treatments (hormones, growth regulators, ...). Details about these training samples are provided in the manual of the package (make a link).

## Installation

*ath1121501frmavecs*  is an annotation package for the *fRMA* Bioconductor package, allowing to preprocess micro-array data for *Arabidopsis thaliana* given the **GPL198** platform (Array: Affymetrix Arabidopsis ATH1 Genome Array; CDF: ATH1-121501). It require a version of R >= 3.5.0. It can be installed from the github repository as follows:

```
    if(!require("remotes", quietly = TRUE)){  
        install.packages("remotes")
        }
    if(!require("BiocManager", quietly = TRUE)){  
        install.packages("BiocManager")
    }
   options(repos = BiocManager::repositories())
   getOption("repos")
   BiocManager::install("RiviereQuentin/ath1121501frmavecs",                     
     dependencies = TRUE,                     
     build_vignettes = TRUE,
     force = TRUE)
```

## Usage

In the following example, we show how to compute fRMA values for three samples selected on the **GPL198** platform, namely: GSM433634, GSM1179807, and GSM433644 (GEO sample IDs). First, we download the corresponding CEL files with the *GEOquery* package and retrieve the paths of those files on our local computer; then we load the data with *affy*, and preprocess the micro-array data and extract the expression values with *frma*.

```
	library(ath1121501frmavecs)
	library(GEOquery)
	library(affy)
	library(frma)

	# 1. Download the CEL files
	GSM_considered <- c("GSM433634", "GSM1179807", "GSM433644")
	for (sample in GSM_considered){
	  getGEOSuppFiles(sample, 
		            makeDirectory = FALSE,
		            baseDir = getwd(), ## Or change to another directory
		            filter_regex = "cel.gz || CEL.gz",
		            fetch_files = TRUE)
	  }

	# 2. Obtain the CEL file paths
	celpaths <- grepl(pattern = paste0(c(".*(", paste0(GSM_considered, collapse = "|"), ").*"), collapse = ""), 
		          x = list.files(getwd()), ## Replace getwd() by baseDir
		          ignore.case = TRUE)
	celpaths <- list.files(getwd())[celpaths]
	celpaths <- celpaths[grepl(pattern = "cel.gz", celpaths, ignore.case = TRUE)]

	# 3. Load and preprocess the data
	cel_batch <- read.affybatch(celpaths)
	data(ath1121501frmavecs)
	cel_frma <- frma(object = cel_batch, target = "full", input.vec = ath1121501frmavecs, verbose = TRUE)

	# 4. Extract the fRMA values
	cel_e <- exprs(cel_frma)
	cel_e <- as.data.frame(cel_e)
	head(cel_e)
```

Finally, we can annotate the probesets with the gene names and select, for each gene model, only one probeset, with the highest levels of signal on average. The annotations of the probesets are available within the *ath1121501.db* package. 

```
	library(ath1121501.db)
	probeset2gene <- mapIds(ath1121501.db,
		                keys = rownames(cel_e),
		                column = "TAIR", keytype = "PROBEID")
	#cel_e <- cel_e[rownames(cel_e) %in% probeset2gene$array_element_name,]
	cel_e$AGI <- probeset2gene[match(x = rownames(cel_e), table = names(probeset2gene))]
	multiple_ps <- table(cel_e$AGI)
	multiple_ps <- names(multiple_ps[multiple_ps > 1])
	todel_ps <- list()
	for (gene in multiple_ps){
	  considered <- cel_e[cel_e$AGI == gene,]
	  mu_considered <- apply(as.matrix(considered[,seq(1, ncol(considered)-1)]), 1, mean)
	  mu_considered <- mu_considered[order(mu_considered, decreasing = TRUE)]
	  todel_ps[[gene]] <- names(mu_considered)[2:length(mu_considered)]
	}
	todel_ps <- do.call(c, todel_ps)
	cel_e <- cel_e[!(rownames(cel_e) %in% todel_ps),]
```


## References

McCall MN, Bolstad BM, Irizarry RA. Frozen robust multiarray analysis (fRMA). Biostatistics. 2010 Apr;11(2):242-53. doi: 10.1093/biostatistics/kxp059. Epub 2010 Jan 22. PMID: 20097884; PMCID: PMC2830579.

McCall MN, Irizarry RA. Thawing Frozen Robust Multi-array Analysis (fRMA). BMC Bioinformatics. 2011 Sep 16;12:369. doi: 10.1186/1471-2105-12-369. PMID: 21923903; PMCID: PMC3180392.

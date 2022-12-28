# Benchmarking HLA genotype prediction from scRNA data

This repo contains code corresponding to the manuscript ["Prediction of HLA genotypes from single-cell transcription data"](https://www.biorxiv.org/content/10.1101/2022.06.09.495569v1)

### Abstract 

- The human leukocyte antigen (HLA) locus plays a central role in adaptive immune function and has significant clinical implications for tissue transplant compatibility and allelic disease associations. Studies using bulk-cell RNA sequencing have demonstrated that HLA transcription may be regulated in an allele-specific manner and single-cell RNA sequencing (scRNA-seq) has the potential to better characterize these expression patterns. However, quantification of allele-specific expression (ASE) for HLA loci requires sample-specific reference genotyping due to extensive polymorphism. While genotype prediction from bulk RNA sequencing is well described, the feasibility of predicting HLA genotypes directly from single-cell data is unknown. Here we evaluate and expand upon several computational HLA genotyping tools by comparing predictions from human single-cell data to gold-standard, molecular genotyping. The highest 2-field accuracy averaged across all loci was 76% by arcasHLA and increased to 86% using a composite model of multiple genotyping tools. We also developed a highly accurate model (AUC 0.93) for predicting HLA-DRB345 copy number in order to improve genotyping accuracy of the HLA-DRB locus. Genotyping accuracy improved with read depth and was reproducible at repeat sampling. Using a metanalytic approach, we also show that HLA genotypes from PHLAT and OptiType can generate ASE ratios that are highly correlated (R2 = 0.8 and 0.94, respectively) with those derived from gold-standard genotyping.

### Notes on organization

- Each figure in the manuscript corresponds to a numbered directory. 
  - Each figure directory contains a similarly numbered R markdown notebook (`.Rmd`) with the code relevant to reproducing the analysis
  - Each notebook has a corresponding `.nb.html` file, which can be downloaded and opened in a local browser to view code an inline figures, without having to run the notebook in R
- Data relevant to reproducing the analysis can be found in the `data/` directory 
- A summarized guide to implement the models developed in this manuscript can be found in the `model_.../` directories

### Notes on models

- We include two models in this project:
  - `model_HLAD_DRB_kNN` - A kNN model that predicts copy numbers of the HLA-DRB345 locus
  - `model_genotype_composite` - A decision tree model that determines the most accurate HLA genotype when presented outputs from arcasHLA, OptiType, PHLAT, and HLAminer.
- Both directories contain example data and an R notebook that shows how each model can be applied to this data. 

### Notes on reproducibility

- The code in this repo is capable of reproducing the majority of the manuscript's analyses as is.
- Using the `here` package in R, nearly all file paths in the associated code will recognize this repo as the base path, regardless of what local directory the repo is forked to.
- However, due to the size of certain raw files (`fastq`, `bam`, etc.) that are too large for standard github repositories, some code that reflects sequencing, rather than analytic workflows, has not been refactored to reflect relative file paths. These will need to be changed manually.
  - However, such pipelines represent standard implementations of published tools (e.g. `HISAT`, `arcasHLA`), so our specific code need not be strictly followed.
- For analysis downstream of these raw sequencing pipelines, the relevant intermediate output files have been included in the `data\` folder for use with analysis pipelines. 
- It is best to start a new R session when replicating different figure's analysis notebooks. Some libraries across different notebooks have conflicting name spaces. 



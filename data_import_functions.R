library(tidyverse)

# # For tests
# isb_path <- "/labs/khatrilab/solomonb/covid/isb"
# isb_samples <- read_tsv("/labs/khatrilab/solomonb/covid/isb/logs/210217_232725/parallel.log") %>% 
#   separate(Command, into = c(NA, "sample"), sep = " ") %>% 
#   pull(sample) %>% unique()

### Import ISB molecular genotyping
invitro_import <- function(path = "/local-scratch/datasets/ISBOtherCovid/hla-haplotyping/scisco-genetics/data/hla_2020-10-23_1457.tsv"){
  read_tsv(path) %>% 
    select(`Sample ID`, `Allele 1`, `Allele 2`) %>%
    pivot_longer(cols = !`Sample ID`, names_to = "allele_id", values_to = "allele") %>% 
    filter(!grepl("Not_Present", allele)) %>% 
    mutate(allele_id = str_extract(allele_id, "[1-2]"),
           genotyper = "invitro") %>% 
    rename("sample" = "Sample ID")
}
invitro_import()


### Run arcasHLA merge
# Refreshes the genotypes.tsv file created by arcasHLA merge
arcasHLA_merge <- function(path="/labs/khatrilab/solomonb/covid/isb/arcasHLA"){
  command <- sprintf(
    "source /labs/khatrilab/solomonb/miniconda3/etc/profile.d/conda.sh;
    conda activate samtools;
    ARCAS_DIR=%s;
    arcasHLA merge --i $ARCAS_DIR --o $ARCAS_DIR;
    conda deactivate",
  path)
  system(command)
}
# arcasHLA_merge()


### arcasHLA import
arcas_import <- function(path, call_merge=F){
  if(call_merge==T){
    arcasHLA_merge()
  }
  read_tsv(sprintf("%s/genotypes.tsv", path)) %>% 
    pivot_longer(-subject, names_to = "allele_id", values_to = "allele") %>% 
    mutate(allele_id = str_sub(allele_id, -1),
           genotyper = "arcasHLA") %>% 
    dplyr::rename(sample = subject)
}
# arcas_import(path = sprintf("%s/arcasHLA", isb_path))


### PHLAT import
phlat_import <- function(path, sample){
  sample_path <- sprintf("%s/%s_HLA.sum", path, sample)
  if (file.exists(sample_path)){
    read_tsv(sample_path) %>% 
      select(contains("Allele")) %>% 
      pivot_longer(cols = everything(), names_to = "allele_id", values_to = "allele")%>% 
      mutate(allele_id = str_sub(allele_id, -1),
             genotyper = "phlat")
  } else {
    NULL
  }
}
# phlat_import(path = sprintf("%s/phlat", isb_path), sample = "INCOV001-CV")

### OptiType import
optitype_import <- function(path, sample){
  sample_path <- sprintf("%s/%s_result.tsv", path, sample)
  if (file.exists(sample_path)){
    suppressWarnings(read_tsv(sample_path)) %>% 
      select(A1:C2) %>% 
      pivot_longer(cols = everything(), names_to = "allele_id", values_to = "allele")%>% 
      mutate(allele_id = str_sub(allele_id, -1),
             genotyper = "optitype")
  } else {
    NULL
  }
}
# optitype_import(path = sprintf("%s/optitype", isb_path), sample = "INCOV001-CV")


### HLAMiner import
hlaminer_import <- function(path, sample){
  sample_path <- sprintf("%s/HLAminer_HPRA_%s.csv", path, sample)
  if (file.exists(sample_path)){
    tibble(lines = read_lines(sample_path)) %>% 
      filter(grepl("^HLA|\t", lines)) %>% 
      separate(lines, into = c("hla", "prediction", "allele"), sep = "\t", 
               fill = "right", extra = "drop") %>% 
      mutate_all(function(x) ifelse(x=="",NA,x)) %>% 
      fill(hla, prediction, .direction = "down") %>% 
      drop_na() %>% 
      separate(allele, into = c("allele", "score", "expected", "confidence"), 
               sep = ",", fill = "right", extra = ) %>% 
      group_by(hla, prediction) %>% 
      mutate_at(vars(c(score, expected, confidence)), as.numeric) %>% 
      # Clean up names
      mutate(allele_id = map_chr(prediction, str_extract, "[1-9]")) %>% 
      mutate(allele = gsub("[A-Z]$", "", allele)) %>% 
      # Select highest score
      group_by(hla, allele_id) %>% 
      top_n(1, score) %>% 
      ungroup() %>% 
      select(allele_id, allele) %>% 
      mutate(genotyper = "hlaminer")
  } else {
    NULL
  }
}
# hlaminer_import(path = sprintf("%s/hla_miner", isb_path), sample = "INCOV003-BL")


### Combined import function
combine_HLA_import <- function(path, samples){
  suppressMessages({
    invitro <- invitro_import()
    
    arcas <- arcas_import(path = sprintf("%s/arcasHLA", path))
    
    phlat <- tibble(sample = samples) %>% 
      mutate(data = map(sample, function(x) {
        phlat_import(path = sprintf("%s/phlat", path), sample = x)
      })) %>% 
      unnest(data)
    
    opti <- tibble(sample = samples) %>% 
      mutate(data = map(sample, function(x) {
        optitype_import(path = sprintf("%s/optitype", path), sample = x)
      })) %>% 
      unnest(data)
    
    miner <- tibble(sample = samples) %>% 
      mutate(data = map(sample, function(x) {
        hlaminer_import(path = sprintf("%s/hla_miner", path), sample = x)
      })) %>% 
      unnest(data)
    
    bind_rows(arcas, phlat, opti, miner, invitro) %>% 
      drop_na()
  })
}

x <- combine_HLA_import(path = isb_path, samples = isb_samples)

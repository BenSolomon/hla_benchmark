require(tidyverse)

# # For tests
# isb_path <- "/labs/khatrilab/solomonb/covid/isb"
# isb_samples <- read_tsv("/labs/khatrilab/solomonb/covid/isb/logs/210217_232725/parallel.log") %>% 
#   separate(Command, into = c(NA, "sample"), sep = " ") %>% 
#   pull(sample) %>% unique()

### Import ISB molecular genotyping
invitro_import <- function(path = "hla_2020-10-23_1457.tsv", exclude_comment_samples=F){
  df <- read_tsv(path)
  # Remove samples with comments (e.g. Sample contamination)
  if (exclude_comment_samples == T){
    exclude_samples <- df %>% filter(!is.na(Comments)) %>% pull(`Sample ID`) %>% unique()
    print(
      sprintf("The following samples were removed from invitro typing due to comments: %s", 
              paste(exclude_samples, collapse = " "))
      )
    df <- df %>% filter(is.na(Comments))
  }
  # Format data
  df %>% 
    select(`Sample ID`, `Allele 1`, `Allele 2`) %>%
    pivot_longer(cols = !`Sample ID`, names_to = "allele_id", values_to = "allele") %>% 
    filter(!grepl("Not_Present", allele)) %>% 
    mutate(allele_id = str_extract(allele_id, "[1-2]"),
           genotyper = "invitro") %>% 
    dplyr::rename("sample" = "Sample ID")
}
# invitro_import()


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


### scHLAcount genotype import
scHLA_genotype_import <- function(path, sample){
  sample_path <- sprintf("%s/%s_results/labels.tsv", path, sample, sample)
  if (file.exists(sample_path)){
    suppressMessages(read_tsv(sample_path, col_names = "allele")) %>% 
      filter(grepl("\\*", allele)) %>% 
      separate(allele, into = "locus", sep = "\\*", extra = "drop", remove = F) %>% 
      group_by(locus) %>% 
      mutate(allele_id = 1:n()) %>% 
      ungroup() %>% 
      mutate(genotyper = "scHLAcount") %>% 
      select(-locus)
  } else {
    NULL
  }
}
# schla_path <- "/labs/khatrilab/solomonb/covid/isb/scHLAcount/scHLAcount_genotyping/consolidated_output"
# scHLA_genotype_import(schla_path, "INCOV003-BL_S5")



### Combined import function
# path - base path
# *_subdir - subdirectory of each genotyper within path
combine_HLA_import <- function(
    path, 
    samples, 
    invitro_path="hla_2020-10-23_1457.tsv", 
    filter_invitro = F, 
    expand_invitro = T,
    arcas_subdir = "arcasHLA",
    phlat_subdir = "phlat",
    opti_subdir = "optitype",
    miner_subdir = "hla_miner",
    schla_subdir = "scHLAcount/results/211024_180304"
  ){
  suppressMessages({
    if (is.null(invitro_path)){
      invitro <- NULL
      expand_invitro <- F
    } else {
      invitro <- invitro_import(path = invitro_path, exclude_comment_samples = filter_invitro)
    }
    
    ### Every subject has molecular HLA typing, but may occur in the BL, AC, or CV timepoint
    ### This expands the molecular HLA typing from the timepoint it occurs to all other time points 
    ### This allows for scRNA genotyping from a given subject and timepoint to be assessed even if that
    ### subjects molecular typing happened at a different time point. Since HLA is genetic, the molecular
    ### typing should not change
    if (expand_invitro == T){
      subjects <- invitro %>% pull(sample) %>% str_replace("-.*$","") %>% unique()
      timepoints <- invitro %>% pull(sample) %>% str_replace("^.*-","") %>% unique()
      invitro <- expand_grid(sample = subjects, time = timepoints) %>% 
        left_join(invitro %>% 
                    separate(sample, into = c("sample", NA), sep = "-"),
                  by = "sample") %>% 
        unite("sample", sample, time, sep = "-")
    }
    
    arcas <- arcas_import(path = sprintf("%s/%s", path, arcas_subdir))

    phlat <- tibble(sample = samples) %>%
      mutate(data = map(sample, function(x) {
        phlat_import(path = sprintf("%s/%s", path, phlat_subdir), sample = x)
      })) %>%
      unnest(data)

    opti <- tibble(sample = samples) %>%
      mutate(data = map(sample, function(x) {
        optitype_import(path = sprintf("%s/%s", path, opti_subdir), sample = x)
      })) %>%
      unnest(data)

    miner <- tibble(sample = samples) %>%
      mutate(data = map(sample, function(x) {
        hlaminer_import(path = sprintf("%s/%s", path, miner_subdir), sample = x)
      })) %>%
      unnest(data)

    schla <- tibble(sample = samples) %>%
      mutate(data = map(sample, function(x) {
        scHLA_genotype_import(path = sprintf("%s/%s", path, schla_subdir), sample = x)
      })) %>%
      unnest(data) %>% 
      mutate_all(as.character)
    
    bind_rows(arcas, phlat, opti, miner, schla, invitro) %>%
      drop_na() %>%
      filter(sample %in% samples)
  })
}

# x <- combine_HLA_import(path = isb_path, samples = isb_samples[1:10], filter_invitro = T)
# x <- combine_HLA_import(path = isb_path, samples = isb_samples)


# Format HLA table
format_hla_table <- function(hla_table){
  hla_table %>% 
    separate(allele, into = c("locus", "fields"), sep = "\\*", 
             fill = "right", extra = "drop") %>% 
    separate(fields, into = c("field_1", "field_2", "field_3"), sep = ":", 
             fill = "right", extra = "drop") %>% 
    mutate(field_2 = paste(field_1, field_2, sep = "_"), 
           field_3 = paste(field_2, field_3, sep = "_")) %>% 
    mutate_at(vars(field_1:field_3), function(x) ifelse(grepl("NA",x),NA,x))
}
# format_hla_table(combine_HLA_import(path = isb_path, samples = isb_samples))


# Import HLA loci alignment stats from arcas logs
hla_mapping_stats_import <- function(samples, log_dir){
  tibble(sample = samples) %>% 
    mutate(data = map(sample, function(x){
      log_path <- sprintf("%s/%s.genotype.log",log_dir,x)
      df <- tibble(lines = read_lines(log_path))
      if (any(grepl("error", df$lines, ignore.case = T))){
        NA
      } else {
        df %>% 
          mutate(lines = gsub("\t", "", lines)) %>% 
          filter(grepl("^HLA", lines)) %>% 
          separate(lines, into = c("locus", "abundance", "reads", "classes"), sep = " +")
      }
    })) %>%
    unnest(data) %>% 
    select(!(contains("data"))) %>% 
    separate(locus, into = c(NA, "locus"), sep = "-") %>% 
    mutate_at(c("reads", "classes"), as.numeric)
}
# arcas_log_dir <- sprintf("%s/arcasHLA", isb_path)
# alignment_stats_df <- hla_mapping_stats_import(isb_samples, arcas_log_dir)

# Given path to log file, extracts times and duration of each step 
# in sequencing pipeline as a data frame
parse_log_time <- function(path){
  df <- tryCatch({suppressMessages(
    read_tsv(path, col_names = "lines") %>% 
      filter(grepl("^### ", lines)) %>% 
      mutate(lines = gsub("### |\\[|\\]", "", lines)) %>% 
      separate(lines, into = c("status", "component", "time"), sep = "___") %>% 
      mutate(time = parse_date_time(time, "mdy HMS p", tz = "America/Los_Angeles")) %>% 
      pivot_wider(names_from = status, values_from = time) %>% 
      mutate(process_time = difftime(COMPLETE, START, unit = "hours"))
  )}, error = function(c) NA)
  return(df)
}
# sample_path <- "/labs/khatrilab/solomonb/covid/isb/logs/210217_232725/INCOV019-AC/INCOV019-AC_pipeline.log"
# parse_log_time(sample_path)

# Gets bp length of set of alleles via IMGTHLA reference sequences
get_allele_length <- function(alleles, IMGTHLA_path = "/labs/khatrilab/solomonb/references/IMGTHLA"){
  path <- sprintf("%s/hla_gen.fasta", IMGTHLA_path)
  command_string <- sprintf("cat %s | awk '$1 ~ /^>/ {print}'", path)
  tibble(header = system(command_string, intern = T)) %>% 
    separate(header, into = c("id", "allele", "bp"), sep = " ", extra = "drop") %>% 
    filter(allele %in% alleles) %>% 
    separate(allele, into = "locus", sep = "\\*", extra = "drop") %>% 
    select(-id) %>% 
    mutate_at(c("bp"), as.numeric)
}
# get_allele_length(c("B*07:02:01:01"))


# Basic scHLAcount count matrix import
scHLA_import <- function(sample, result_path, label_path, barcode_path){
  pool <- str_split(sample, "_",simplify = T)[2]
  mm <- Matrix::readMM(result_path)
  dimnames(mm) <- list(
    allele = read_tsv(label_path, col_names = "allele") %>% pull(allele), 
    cell = read_tsv(barcode_path, col_names = "barcode") %>% pull(barcode))
  data.frame(allele=rownames(mm)[mm@i + 1], cell=colnames(mm)[mm@j + 1], count=mm@x) %>% 
    mutate(cell = sprintf("%s.%s", pool, cell))
}

# Applies a fixed order to paired alleles within a given sample
# Necessary for ratios to ensure which allele is numerator/denominator
apply_allele_order <- function(label_path, data){
  # Key to fix order of alleles
  allele_order_key <- read_tsv(label_path, col_names = "allele") %>% 
    filter(grepl("\\*", allele)) %>% 
    separate(allele, into = c("gene", "allele_id"), sep = "\\*", remove = F) %>% 
    group_by(gene) %>% 
    mutate(allele_order = 1:n(),
           n_alleles_expected = n()) %>% 
    ungroup() %>% 
    select(allele, allele_order, n_alleles_expected)
  # Merge key to data
  data %>% 
    left_join(allele_order_key, by = "allele")
}

# Records number of expected alleles based on genotype provided to scHLAcount,
# the number of observed alleles, the counts for each allele as well as counts 
# mapped to an HLA-gene but not a specific allele. Calculates sums and allele ratio
gene_sums_and_ratio <- function(data){
  data %>% 
    separate(allele, into = c("gene", "fields"), sep = "\\*", remove = F, fill = "right") %>%
    group_by(cell, gene) %>% 
    # Count total reads per gene, including those not assigned at the allele level
    mutate(gene_sum_all = sum(count)) %>% 
    # Remove rows for counts not assigned at allele level
    filter(!is.na(fields)) %>% 
    # Count total reads per gene, only those that were typed at the allele level
    mutate(gene_sum_typed = sum(count)) %>% 
    # Calculate allele ratio
    mutate(allele_ratio = count / gene_sum_typed) %>% 
    # Count observed number of alleles
    mutate(n_alleles_observed = n()) %>% 
    drop_na() %>% 
    select(cell, allele, allele_order, 
           gene, fields, 
           n_alleles_expected, n_alleles_observed, 
           everything())
}

# Controller function for:
### scHLA_import()
### apply_allele_order() 
### gene_sums_and_ratio()
scHLA_data_processing <- function(sample, result_dir, barcode_dir){
  # Check for all file paths
  result_path <- sprintf("%s/%s_results/count_matrix.mtx", result_dir, sample)
  label_path <- sprintf("%s/%s_results/labels.tsv", result_dir, sample)
  barcode_path <- sprintf("%s/%s_barcode.tsv", barcode_dir, sample)
  check_files <- file.exists(c(result_path, label_path, barcode_path))
  if (any(!check_files)){return("Invalid file")}
  
  suppressMessages(tryCatch({
    df <- scHLA_import(sample = sample, 
                       result_path = result_path,
                       label_path = label_path,
                       barcode_path = barcode_path)
    
    df <- apply_allele_order(label_path = label_path,
                             data = df)
    
    gene_sums_and_ratio(df)
  }, error = function(c) NA))
}



# scHLA_data_processing(
#   sample="INCOV005-BL_S7",
#   result_dir=sprintf("%s/scHLAcount/output/invitro", isb_path),
#   barcode_dir=sprintf("%s/scHLAcount/barcodes", isb_path)
# )

library(tidyverse)
library(checkmate)


### Allele matching function
allele_match <- function(ref,comp,ground_truth=T, verbose=F, count_NA_match = F){
  display_message <-  ifelse(verbose == F, "message", "none")
  
  suppressMessages({
    # If allele ends with an NA string, convert entire object to NA
    ref[grepl("NA$",ref)] <- NA
    comp[grepl("NA$",comp)] <- NA
    
    # Count non-NA alleles
    n_ref <- sum(!is.na(ref))
    n_comp <- sum(!is.na(comp))
    
    # How to handle NAs in both reference and comparison
    if(n_ref==0 & n_comp==0){
      # Default: Return NA
      if (count_NA_match == F){message("No alleles for comparison");return(NA)}
      # Alternative: Consider perfect match. Relevant for DRB345 where lacking a
      # prediction when there is no reference is actually informative
      else {return(1)}
    }
    
    # If ref is ground_truth, return NA if all ref alleles are NA
    if(n_ref==0 & ground_truth==T){message("No ground truth alleles");return(NA)}
    # Cannot have more than two alleles
    if(n_ref > 2 | n_comp > 2){message("Invalid number of alleles");return(NA)}
    
    
    # If comparing to a ground_truth, only penalize NAs in comp, not those in ref
    if (n_ref == 1 & ground_truth == T){return(sum(any(ref[!is.na(ref)] %in% comp, na.rm = T)))}
    
    # If unequal number of NAs between ref and comp, fill smaller with random placeholder
    ## If both genotypes only have one allele, don't want to fill with random placeholder
    ## Because this would automatically create a non-match
    if (n_ref != n_comp){
      comp <-  c(comp[!is.na(comp)], sample(c("fill_1","fill_2"), 2-n_comp, replace = F))
      ref <-  c(ref[!is.na(ref)], sample(c("fill_3","fill_4"), 2-n_ref, replace = F))
    }
    
    # Calculate bray-curtis similarity between alleles
    accuracy <- bind_rows(
      tibble(label = "ref", allele = ref),
      tibble(label = "comp", allele = comp)
    ) %>%
      table() %>% 
      vegan::vegdist(method = "bray") %>% 
      as.numeric() %>% 
      (function(x) 1-x)
    return(accuracy)
  },
  classes = display_message) # Displays messages based on verbose
}

# ### Allele matching function ===========================================OLD
# allele_match <- function(ref,comp){
#   n_ref <- sum(!is.na(ref))
#   n_comp <- sum(!is.na(comp))
#   if (n_comp == 1){
#     comp <- c(comp, NA)
#   }
#   # No meaningful comparison if reference is NA
#   if (all(is.na(ref))){
#     NA 
#     # Cannot have more than two alleles
#   } else if (n_ref >2 | n_comp > 2){
#     return(NA)
#     print("Invalid number of alleles")
#     # If only one allele typed in reference, only need one match in comparison
#   } else if (n_ref == 1){
#     sum(any(ref == comp, na.rm = T))
#     # Matches when all alleles present
#     # Forward and reverse match search since order does not matter
#   } else if (n_ref == 2){
#     max(sum(ref == comp, na.rm = T)/2, sum(ref == rev(comp), na.rm = T)/2)
#   }
# }

# unit_test <- tribble(
#   ~reference, ~comparison, ~expected_standard, ~expected_agreement,
#   c("03", "03"), c("03", "03"), 1, 1,
#   c("03", "01"), c("03", "01"), 1, 1,
#   "03", c("03", "03"), 1, 0.5,
#   c("03", "01"), c("03", "03"), 0.5, 0.5,
#   c("03", "03"), c("03", "01"), 0.5, 0.5,
#   c("03", "01"), c(NA, "01"), 0.5, 0.5,
#   c("03", "03"), "03", 0.5, 0.5,
#   c("01", NA), c(NA, NA), 0, 0,
#   "01", c(NA, NA), 0, 0,
#   c("01", NA), c("01", NA), 1, 1,
#   "01", c("01", NA), 1, 1,
#   c("03", "03"), c("01", "01"), 0, 0,
#   c("03", "01"), NA, 0, 0,
#   c(NA, NA), c(NA, NA), NA, NA,
#   c(NA, NA), NA, NA, NA,
#   NA, NA, NA, NA,
#   c("03","03","03"), c("01", "01"), NA, NA
# )
# unit_test
# 
# unit_test_output <- unit_test %>%
#   mutate(
#     calculated_standard = map2_dbl(reference, comparison, allele_match, ground_truth = T),
#     calculated_agreement = map2_dbl(reference, comparison, allele_match, ground_truth = F), 
#     match_standard = expected_standard == calculated_standard,
#     match_agreement = expected_agreement == calculated_agreement)
# unit_test_output


compare_hla<- function(hla_df, 
                       reference = "invitro", 
                       exclude_missing=T, 
                       penalize_extra_drb345 = T, 
                       penalize_extra_classic = F,
                       match_drb345_na = T){
  # Input checks
  arg_col <- makeAssertCollection()
  hla_df_columns <- c("sample", "allele_id", "locus", "field_1", "field_2", "field_3", "genotyper")
  assertNames(names(hla_df), permutation.of = hla_df_columns, .var.name = "hla_df column names", add = arg_col)
  assertChoice(reference, unique(hla_df$genotyper), add = arg_col)
  assertLogical(exclude_missing, add = arg_col)
  assertLogical(penalize_extra_drb345, add = arg_col)
  assertLogical(penalize_extra_classic, add = arg_col)
  if (arg_col$isEmpty()==F) {map(arg_col$getMessages(),print);reportAssertions(arg_col)}
  
  # Function
  df <- hla_df %>%
    select(-allele_id) %>%
    group_by(sample, locus, genotyper) %>%
    summarise_all(list) %>%
    pivot_longer(cols = !c(sample, locus, genotyper), names_to = "field", values_to = "allele") %>%
    pivot_wider(names_from = "genotyper", names_prefix = "genotyper_", values_from = "allele", values_fill = list(NA)) %>%
    dplyr::rename("reference" = contains(reference)) %>%
    pivot_longer(contains("genotyper"),  names_to = "genotyper", values_to = "allele",names_prefix = "genotyper_")
  if (exclude_missing == T){
    df <- df %>%
      filter_at(vars(c(reference, allele)), function(x) !is.na(x))
      # mutate(n_reference = map_dbl(reference, function(x) sum(!is.na(x) & !grepl("NA", x)))) %>%
      # mutate(n_allele = map_dbl(allele, function(x) sum(!is.na(x) & !grepl("NA", x))))
      # filter_at(vars(contains("n_")), function(x) x ==2) # Removed b/c excluding samples where invitro n = 1 and sample n = 2, which is unwanted
      # filter(n_allele == 2)
  }
  df_drb345 <- df %>% 
    filter(grepl("^DRB[345]",locus)) %>% 
    mutate(accuracy = map2_dbl(reference, allele, allele_match, 
                               ground_truth = !penalize_extra_drb345, 
                               count_NA_match = match_drb345_na))
  df_other <- df %>% 
    filter(grepl("^[ABC]",locus) | grepl("^D[PQR][AB]1", locus)) %>% 
    mutate(accuracy = map2_dbl(reference, allele, allele_match, 
                               ground_truth = !penalize_extra_classic, 
                               count_NA_match = F))
  df <- bind_rows(
    df_other,
    df_drb345
  )
  # if (na_drop == T){
  #   df <- df %>% drop_na()
  # }
  df
}
# compare_hla(hla_df = all_hla, reference = "invitro", exclude_missing = T, compare_function = "accuracy")



### Expects format from accuracy_functions::compare_hla
### Summarizes accuracy across locus, field, and genotyper
calculate_summary_df <- function(df){
  suppressMessages({df %>% 
      group_by(locus, field, genotyper) %>% 
      summarise(mean_accuracy = mean(accuracy, na.rm = T),
                sd = sd(accuracy, na.rm=T),
                se = sd(accuracy, na.rm=T)/sqrt(n()),
                n =) %>% 
      ungroup() %>% 
      mutate(mean_accuracy = ifelse(mean_accuracy == 0, NA, mean_accuracy), 
             sd = ifelse(mean_accuracy == 0, NA, sd),
             se = ifelse(mean_accuracy == 0, NA, se),
             locus = factor(locus, levels = sort(unique(locus), decreasing = T)))
  })
}



### Expects format from data_import_functions::format_hla_table
### Tallys number of alleles predicted for each position {0,1,2}
allele_tally <- function(df){
  df %>% 
    select(-allele_id) %>% 
    pivot_longer(cols = field_1:field_3, names_to = "field", values_to = "allele") %>% 
    pivot_wider(names_from = "genotyper", values_from = "allele", values_fn = function(x) sum(!is.na(x)), values_fill = 0) %>% 
    filter(hlaminer <= 2) %>%     ### NEED TO DEAL WITH HLAMINER AMBIGOUS ALLELES ###
    ungroup() %>% 
    pivot_longer(cols = !sample:field, names_to = "genotyper", values_to = "count")
}




# Replaces SRA sample names with experiment sample names for PMID samples
# Only requires input df has a single column names "sample" with SRA IDs
pmid_sample_rename <- function(df){
  # Import naming key
  pmid_key <- suppressMessages(read_csv("/labs/khatrilab/solomonb/covid/pmid30518681/pmid30518681_samples.csv") %>% 
                                 select(sample_name = sample, run_name = Run))
  
  # Check if names already properly formatted
  name_check <- any(df$sample %in% pmid_key$sample_name)
  if (name_check){
    message("Names already formatted")
    return(df)
  }
  
  # Format names if needed
  df %>% 
    left_join(pmid_key, c("sample" = "run_name")) %>% 
    select(-sample) %>% 
    rename("sample" = sample_name)
}

# Calculate ratio of DRB345/DRB1 reads for each sample
# Expects output of hla_mapping_stats_import
get_drb345_ratios <- function(df){
  df %>% 
    filter(grepl("^DRB[1345]", locus)) %>% 
    select(sample, locus, reads) %>% 
    mutate(reads = as.numeric(reads)) %>% 
    pivot_wider(names_from = "locus", values_from = 'reads', values_fill = 0) %>% 
    mutate_at(vars(DRB3:DRB5), funs(./DRB1)) %>% 
    select(-DRB1) %>% 
    pivot_longer(!sample, names_to = 'locus', values_to ='frequency')
}

# Converts invitro HLA typing into counts of predicted alleles for each DRB345 gene
get_drb345_copy_number <- function(invitro_path){
  invitro <- invitro_import(path = invitro_path, exclude_comment_samples = T)
  invitro %>% 
    separate(allele, into=c('locus','allele'), sep = '\\*') %>% 
    filter(grepl('^DRB[345]',locus)) %>%
    select(-allele, -genotyper) %>% 
    pivot_wider(names_from ='locus', values_from = 'allele_id', values_fn = length, values_fill = 0) %>% 
    pivot_longer(!sample, names_to = 'locus', values_to = 'copy_number') 
}

# Combined outputs of get_drb345_ratios and get_drb345_copy_number into single DF
create_drb345_df <- function(ratio_df, copy_number_df){
  ratio_df %>% 
    pivot_wider(names_from = 'locus', values_from ='frequency', values_fill=0) %>% 
    left_join(copy_number_df, by = c('sample')) %>% 
    drop_na()
}

# Creates an expanded df of all valid loci-field-genotyper combinations for drb345
# Used to join data to include entries for absent predictions
# Takes a list of samples create all combinations along
# Valid combinations determined from isb data since many samples
make_drb_scaffold <- function(samples){
  isb_hla <- readRDS("all_hla_expanded.RDS")
  valid_combinations <- isb_hla %>% 
    filter(grepl("^DRB[345]", locus)) %>% 
    select(locus, contains("field"), genotyper) %>% 
    pivot_longer(contains("field"), names_to = "field", values_to = "allele") %>% 
    drop_na() %>% 
    select(-allele) %>% 
    distinct() %>% 
    unite(valid_combination, sep = "__") %>% 
    pull(valid_combination)
  df <- expand.grid(
    sample = unique(samples),
    allele_id = c("1","2"),
    valid_combination = valid_combinations
  ) %>% 
    separate(valid_combination, into = c("locus", "genotyper", "field"), sep = "__")
  return(df)
}

# Calculates drb345 accuracy based at all valid locus-field-genotyper combinations 
# (e.g. no field_3 for HLAminer), even if not seen in data.
# score_absent_predictions determines the absence of a prediction is itself considered
# informative and given a pseudovalue "X" which is then scored by allele_match.
# Expects df in the same format as from format_hla_table  
calculate_drb345_accuracy <- function(df, score_absent_predictions=F){
  # Input checks
  arg_col <- makeAssertCollection()
  df_columns <- c("sample", "allele_id", "locus", "field_1", "field_2", "field_3", "genotyper")
  assertNames(names(df), permutation.of = df_columns, .var.name = "df column names", add = arg_col)
  assertLogical(score_absent_predictions, add = arg_col)
  if (arg_col$isEmpty()==F) {map(arg_col$getMessages(),print);reportAssertions(arg_col)}
  
  #Function
  drb_scaffold <- make_drb_scaffold(unique(df$sample))
  
  df <- df %>% 
    pivot_longer(contains("field"), names_to = "field", values_to = "allele") %>% 
    right_join(drb_scaffold, by = names(drb_scaffold)) 
  
  if (score_absent_predictions==T){
    df <- df %>% mutate(allele = ifelse(is.na(allele), "X", allele))
  }
  
  df %>% 
    pivot_wider(names_from = "field", values_from = "allele", values_fn = function(x) ifelse(length(x) > 1, NA, x)) %>% 
    compare_hla(
      hla_df = .,
      reference = "invitro",
      exclude_missing = F,
      match_drb345_na = F)
}

# Wrapper to combine accuracy calculation for drb345 and non-drb345 loci
# Expects a df in the form of format_hla_table
calculate_all_hla_accuracy <- function(df){
  non_drb345_df <- df %>% 
    filter(!(grepl("^DRB[345]", locus))) %>% 
    compare_hla(
      hla_df = ., 
      reference = "invitro", 
      exclude_missing = F)
  
  drb345_df <- all_hla_drb345_filtered %>% 
    calculate_drb345_accuracy()
  
  df <- bind_rows(non_drb345_df, drb345_df)
  
  return(df)
}



# Returns best genotyper as string based on a set heirachry 
# arcas, opti, phlat are all length-2 allele vectors or NA
# field and locus are strings
composite_picker <- function(locus, arcas, opti, phlat, field, include_phlat=T){
  # Convert input allele vector to a list, so case_when isn't vectorized
  for (g in c("arcas", "opti", "phlat")){
    assign(g, list(get(g)))
  }
  # Pick composite
  if (locus %in% c("A","B","C")){
    output <- case_when(
      !is.na(opti) & field != "field_3" ~ "optitype",
      !is.na(phlat) & field != "field_3" & include_phlat == T ~ "phlat",
      !is.na(arcas) ~ "arcasHLA",
      TRUE ~ NA_character_
    )
  } else if (locus %in% c("DQB1")){
    output <- case_when(
      !is.na(phlat) & include_phlat == T ~ "phlat",
      !is.na(arcas) ~ "arcasHLA",
      TRUE ~ NA_character_
    )
  } else {
    output <- case_when(
      !is.na(arcas) ~ "arcasHLA",
      !is.na(phlat) & include_phlat == T ~ "phlat",
      TRUE ~ NA_character_
    )
  }
  return(output)
}
# composite_picker(locus = "A", arcas = c("01","01"), opti = c("02","02"), phlat = c("03","03"), field = "field_1")


# Returns allele values  in arcas, phlat, or opti 
# based on genotyper selection specified by selection
# Used in combination with composite_picker
composite_selector <- function(arcas, opti, phlat, selection){
  output <- case_when(
    is.na(selection) ~ list(NA),
    selection == "arcasHLA" ~ list(arcas),
    selection == "optitype" ~ list(opti),
    selection == "phlat" ~ list(phlat),
    TRUE ~ list(NA)
  )
  return(unlist(output))
}
# composite_selector(arcas = 1, opti = 2, phlat = 3, selection = "optitype") 


# Takes an accuracy df from compare_hla or calculate_all_hla_accuracy
# and applies composite_picker and composite_selector to create a composite genotyper
create_composite_df <- function(df, include_phlat = T, composite_name = "composite"){
  df <- df %>% 
    pivot_wider(names_from = "genotyper", values_from = c("accuracy", "allele")) %>% 
    # Pick genotyper with composite_picker
    mutate(composite_source = try(pmap_chr(
      list(
        locus = locus, 
        arcas = allele_arcasHLA,
        opti = allele_optitype, 
        phlat = allele_phlat, 
        field = field,
        include_phlat = include_phlat), 
      composite_picker))) %>% 
    # Isolate composite accuracy values with composite_selector
    mutate(accuracy_comp = try(pmap_dbl(
      list(
        arcas = accuracy_arcasHLA, 
        opti = accuracy_optitype, 
        phlat = accuracy_phlat, 
        selection = composite_source), 
      composite_selector))) %>% 
    # Isolate composite alleles with composite_selector
    mutate(allele_comp = try(pmap(
      list(
        arcas = allele_arcasHLA, 
        opti = allele_optitype, 
        phlat = allele_phlat, 
        selection = composite_source), 
      composite_selector))) 
  # Pivot df back to original configuration
  df <- df %>% 
    ungroup() %>% 
    pivot_longer(
      contains(c("allele", "accuracy")),
      names_to = c(".value", "genotyper"),
      names_sep = "_"
    ) %>% 
    mutate(genotyper = ifelse(genotyper == "comp", composite_name, genotyper))
  return(df)
}


# Averages accuracy values across MHC1, MHC2, and all MHC loci 
# Expects output of compare_hla or calculate_all_hla_accuracy
locus_average_accuracy <- function(df){
  # First average DRB3, DRB4, DRB5 accuracy into single DRB345 term
  drb345_average <- df %>% 
    select(-composite_source) %>% 
    filter(grepl("^DRB[345]", locus)) %>% 
    group_by(sample, field, genotyper) %>% 
    summarise(accuracy = mean(accuracy, na.rm = T)) %>% 
    mutate(locus = "DRB345")
  
  # Calculate averages across MHC1 and MHC2 independently
  average_byClass <- df %>% 
    filter(!grepl("^DRB[345]", locus)) %>% 
    bind_rows(drb345_average) %>% 
    ungroup() %>% 
    mutate(class = ifelse(grepl("^D", locus), "MHC II", "MHC I")) 
  
  # Calculate averages all MHC
  average_all <- df %>% 
    filter(!grepl("^DRB[345]", locus)) %>% 
    bind_rows(drb345_average) %>% 
    ungroup() %>% 
    mutate(class = "MHC All")
  
  # Combine MHC-split and combined averages into single df
  average_combined <- bind_rows(average_byClass, average_all) %>% 
    group_by(sample, class, field, genotyper) %>% 
    summarise(accuracy = mean(accuracy, na.rm = T)) %>% 
    rename("locus" = "class") 
  
  return(average_combined)
}

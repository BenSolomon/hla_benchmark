require(tidyverse)
require(tidymodels)
require(checkmate)
require(lubridate)

# Bray-Curtis similarity is ideal metric for matching when ref and comp both have 2 alleles
# Will be used inside of allele_match, which also has specifications for when n allele != 2
bray_match <- function(ref, comp){
  bc_output <- bind_rows(
    tibble(label = "ref", allele = ref),
    tibble(label = "comp", allele = comp)
  ) %>%
    table() %>% 
    vegan::vegdist(method = "bray") %>% 
    as.numeric() %>% 
    (function(x) 1-x)
  return(bc_output)
}


### Allele matching function
allele_match <- function(ref,comp,ground_truth=T, verbose=F, count_NA_match = F, mode = "accuracy"){
  # Input checks
  arg_col <- makeAssertCollection()
  assertChoice(mode, c("accuracy", "success"), add = arg_col)
  assertLogical(ground_truth, add = arg_col)
  assertLogical(verbose, add = arg_col)
  assertLogical(count_NA_match, add = arg_col)
  if (arg_col$isEmpty()==F) {map(arg_col$getMessages(),print);reportAssertions(arg_col)}
  
  display_message <-  ifelse(verbose == F, "message", "none") # Effect of verbose on suppressMessages
  
  suppressMessages({
    # If allele ends with an NA string, convert entire object to NA 
    # e.g. 01:01:NA should be NA for field 3
    ref[grepl("NA$",ref)] <- NA
    comp[grepl("NA$",comp)] <- NA
    
    # Count non-NA alleles
    n_ref <- sum(!is.na(ref))
    n_comp <- sum(!is.na(comp))
    
    # How to handle complete NAs in both reference and comparison
    if(n_ref==0 & n_comp==0){
      # Default: Return NA
      if (count_NA_match == F){message("No alleles for comparison");return(NA)}
      # Alternative: Consider perfect match. Relevant for DRB345 where lacking a
      # prediction when there is no reference is actually informative
      else {return(1)}
    }
    
    # Cannot have more than two alleles
    if(n_ref > 2 | n_comp > 2){message("Invalid number of alleles");return(NA)}
    # Cannot have 0 ref alleles when ref is a ground truth
    if(n_ref==0 & ground_truth==T){message("No ground truth alleles");return(NA)}
    
    # Success = if comparison prediction made, how well matched (i.e. no NA penalty)
    # Conditions where bc_match cannot be used to calculate success (i.e. n_allele != 2)
    if (mode == "success"){
      if(n_comp == 0){return(NA)} # If no prediction, discard
      if (n_ref != n_comp){
        ref <- unique(ref[!is.na(ref)])  # Only 
        success <- sum(comp %in% ref)/n_comp
        return(success)}
    }
    
    # Accuracy = how well matched even if missing predications (i.e. NAs penalized)
    # Conditions where bc_match cannot be used to calculate accuracy (i.e. n_allele != 2)
    if(mode == "accuracy"){
      ## If ref is a ground truth, only penalize NAs in comp, not those in ref
      if (n_ref == 1 & ground_truth == T){
        ref <- ref[!is.na(ref)]
        accuracy <- sum(any(ref %in% comp, na.rm = T))
        return(accuracy)}
      # If unequal number of NAs between ref and comp, fill smaller with random placeholder
      ## If both genotypes only have one allele, don't want to fill with random placeholder
      ## Because this would automatically create a non-match
      if (n_ref != n_comp & mode == "accuracy"){
        comp <-  c(comp[!is.na(comp)], sample(c("fill_1","fill_2"), 2-n_comp, replace = F))
        ref <-  c(ref[!is.na(ref)], sample(c("fill_3","fill_4"), 2-n_ref, replace = F))}     
    }
    
    # Calculate bray-curtis similarity between alleles 
    # Used for success and accuracy when n_alleles == 2
    # or expanded to 2 by accuracy placeholders
    return(bray_match(ref = ref, comp = comp))
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
#   ~reference, ~comparison, ~expected_success, ~expected_accuracy, ~expected_agreement,
#   c("03", "03"), c("03", "03"), 1, 1, 1,
#   c("03", "01"), c("03", "01"), 1, 1, 1,
#   "03", c("03", "03"), 1, 1, 0.5,
#   c("03", "01"), c("03", "03"), 0.5, 0.5, 0.5,
#   c("03", "03"), c("03", "01"), 0.5, 0.5, 0.5,
#   c("03", "01"), c(NA, "01"), 1, 0.5, 0.5,
#   c("03", "03"), "03", 1, 0.5, 0.5,
#   c("01", NA), c(NA, NA), NA, 0, 0,
#   "01", c(NA, NA), NA, 0, 0,
#   c("01", NA), c("01", NA), 1, 1, 1,
#   "01", c("01", NA), 1, 1, 1,
#   c("03", "03"), c("01", "01"), 0, 0, 0,
#   c("03", "01"), NA, NA, 0, 0,
#   c(NA, NA), c(NA, NA), NA, NA, NA,
#   c(NA, NA), NA, NA, NA, NA,
#   NA, NA, NA, NA, NA,
#   c("03","03","03"), c("01", "01"), NA, NA, NA
# )
# # unit_test
# 
# unit_test_output <- unit_test %>%
#   mutate(
#     calculated_success = map2_dbl(reference, comparison, allele_match, ground_truth = T, mode = "success"),
#     match_success = expected_success == calculated_success,
#     calculated_accuracy = map2_dbl(reference, comparison, allele_match, ground_truth = T),
#     match_accuracy = expected_accuracy == calculated_accuracy,
#     calculated_agreement = map2_dbl(reference, comparison, allele_match, ground_truth = F),
#     match_agreement = expected_agreement == calculated_agreement) %>% 
#   select(reference, comparison, contains("success"), contains("accuracy"), contains("agreement"))
# unit_test_output


compare_hla<- function(hla_df, 
                       reference = "invitro", 
                       method = "accuracy", 
                       penalize_extra_drb345 = T, 
                       penalize_extra_classic = F,
                       match_drb345_na = T){
  # Input checks
  arg_col <- makeAssertCollection()
  hla_df_columns <- c("sample", "allele_id", "locus", "field_1", "field_2", "field_3", "genotyper")
  assertNames(names(hla_df), permutation.of = hla_df_columns, .var.name = "hla_df column names", add = arg_col)
  assertChoice(reference, unique(hla_df$genotyper), add = arg_col)
  assertChoice(method, c("accuracy", "success"), add = arg_col)
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
  df_drb345 <- df %>% 
    filter(grepl("^DRB[345]",locus)) %>% 
    mutate(accuracy = map2_dbl(reference, allele, allele_match, 
                               ground_truth = !penalize_extra_drb345, 
                               count_NA_match = match_drb345_na,
                               mode = method))
  df_other <- df %>% 
    filter(grepl("^[ABC]",locus) | grepl("^D[PQR][AB]1", locus)) %>% 
    mutate(accuracy = map2_dbl(reference, allele, allele_match, 
                               ground_truth = !penalize_extra_classic, 
                               count_NA_match = F,
                               mode = method))
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
calculate_summary_df <- function(df, var = "accuracy", remove_invalid=T){
  df <- suppressMessages({df %>% 
      group_by(locus, field, genotyper) %>% 
      summarise(mean_accuracy = mean(!!sym(var), na.rm = T),
                sd = sd(!!sym(var), na.rm=T),
                se = sd(!!sym(var), na.rm=T)/sqrt(n()),
                n =) %>% 
      ungroup() %>% 
      mutate(locus = factor(locus, levels = sort(unique(locus), decreasing = T)))
  })
  if (remove_invalid == T){
    df <- df %>% exclude_genotyper_fields()
  }
  return(df)
}

### Alternative to calculate_summary_df used to find mean of accuracy comparisons
### Groupss must be held in different columns specified by var_1 and var_2
### Summarizes accuracy across locus, field, and genotyper and calculates P-values across groups
### Expects format similar to output of accuracy_functions::compare_hla
accuracy_difference <- function(df, var_1, var_2, sig_function = "wilcox.test"){
  sig_function <- get(sig_function)
  
  df <- df %>% 
    ungroup() %>% 
    mutate(locus = factor(locus, levels = sort(unique(locus), decreasing = T)))
  
  # Calculate accuracy and agreement means
  # Then find difference between variables specified by var_1 and var_2
  mean_df <- df %>%
    group_by(locus, field, genotyper) %>%
    summarise_if(
      is.numeric, # This should capture "accuracy" and "agreement" columns
      mean, 
      na.rm = T) %>%
    mutate(mean_difference = !!sym(var_1) - !!sym(var_2))
  
  # Calculate p-values for paired test of difference between var_1 and var_2
  # P-values adjusted using FDR based on number of loci for given genotyper-field combo
  p_df <- df %>%
    group_by(locus, field, genotyper) %>%
    nest() %>%
    mutate(data = map(data, function(df) {
      tryCatch(
        suppressWarnings(broom::tidy(sig_function(df[[var_1]], df[[var_2]], paired=T)) %>% dplyr::select(p.value)),
        error = function(c) NA
      )
    })) %>%
    unnest(data) %>%
    group_by(field, genotyper) %>%
    mutate(p_val_adj = round(p.adjust(
      p.value, method = "fdr", n = n()
    ), digits = 3))
  
  # Join mean values with their respective p-values
  left_join(mean_df,
            p_df,
            by = c("locus", "field", "genotyper")) %>%
    exclude_genotyper_fields()
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

# Predict DRB345 loci based on KNN and DRB345/DRB1 read ratios
# Expects output of create_drb345_df
apply_drb345_knn <- function(df, model){
  df %>% 
    bind_cols(
      model %>% predict(new_data = df),
      model %>% predict(new_data = df, type = 'prob')
    ) %>% 
    mutate(copy_number = factor(copy_number))
}

# Filter DF of DRB345 alleles based on predicted number of DRB345 loci from KNN
# drb_predictions expects output of apply_drb345_knn
# all_hla_df expects output of combine_HLA_import
filter_drb_by_KNN <- function(drb_predictions, all_hla_df){
  drb_key <- drb_predictions %>% 
    select(sample, locus, copy_number, .pred_class)
  drb_filtered <- all_hla_df %>% 
    filter(genotyper != "invitro" & grepl("^DRB[345]", locus)) %>% 
    left_join(drb_key, by = c("sample", "locus")) %>% 
    mutate(.pred_class = as.character(.pred_class)) %>% 
    mutate_at(vars(c("copy_number", ".pred_class")), function(x) ifelse(is.na(x), "0", x)) %>% 
    filter(allele_id <= .pred_class) %>%  
    select(-copy_number, -.pred_class)
  return(drb_filtered)
}
# filter_drb_by_KNN(drb345_3p_predicted, pmid_3p_hla)

# Wrapper around filter_drb_by_KNN that takes filtered DRB345 alleles and applies them to DF of all HLA loci
# drb_predictions expects output of apply_drb345_knn
# all_hla_df expects output of combine_HLA_import
# samples is vector of unique samples to be included
filter_drb_in_all_hla <- function(drb_predictions, all_hla_df, samples){
  drb_filtered <- filter_drb_by_KNN(drb_predictions = drb_predictions, all_hla_df = all_hla_df)
  all_hla_df %>% 
    filter(!(grepl("^DRB[345]", locus) & genotyper != "invitro")) %>% 
    bind_rows(drb_filtered) %>% 
    filter(sample %in% samples)
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
calculate_drb345_accuracy <- function(df, method = "accuracy", match_drb345_na=F){
  # Input checks
  arg_col <- makeAssertCollection()
  df_columns <- c("sample", "allele_id", "locus", "field_1", "field_2", "field_3", "genotyper")
  assertNames(names(df), permutation.of = df_columns, .var.name = "df column names", add = arg_col)
  assertLogical(match_drb345_na, add = arg_col)
  if (arg_col$isEmpty()==F) {map(arg_col$getMessages(),print);reportAssertions(arg_col)}
  
  # Build a scaffold of valid sample, field, genotyper combinations
  drb_scaffold <- make_drb_scaffold(unique(df$sample))
  
  df <- df %>% 
    pivot_longer(contains("field"), names_to = "field", values_to = "allele") %>%
    # 1st scaffold bind: prevents "X" assignment if not a valid combination
    right_join(drb_scaffold, by = names(drb_scaffold)) 
  
  # Modify scaffold for second application
  drb_scaffold <- drb_scaffold %>% 
    select(-allele_id) %>% 
    filter(genotyper != "invitro") %>% 
    distinct()
  
  df %>% 
    pivot_wider(names_from = "field", values_from = "allele", values_fn = function(x) ifelse(length(x) > 1, NA, x)) %>% 
    compare_hla(
      hla_df = .,
      reference = "invitro",
      method = method,
      match_drb345_na = match_drb345_na) %>% 
    # 2nd scaffold bind: excludes non-valid combinations created by pivot 
    right_join(drb_scaffold, by = names(drb_scaffold)) 
}

# Wrapper to combine accuracy calculation for drb345 and non-drb345 loci
# Expects a df in the form of format_hla_table
calculate_all_hla_accuracy <- function(df, method = "accuracy", match_drb345_na=F){
  non_drb345_df <- df %>% 
    filter(!(grepl("^DRB[345]", locus))) %>% 
    compare_hla(
      hla_df = ., 
      reference = "invitro", 
      method = method)
  
  drb345_df <- df %>% 
    filter(grepl("^DRB[345]", locus)) %>% 
    calculate_drb345_accuracy(method = method, match_drb345_na = match_drb345_na)
  
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
    # select(-composite_source) %>% 
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

# Converts accuracy df to format to be used in composite decision tree modeling
# Accepts input from accuracy functions like calculate_all_hla_accuracy
# for_fit = T when input for building a model, for_fit=F when input for applying a model
make_tree_input <- function(df, genotypers = c("arcasHLA", "optitype", "phlat"), for_fit = F){
  output <- df %>% 
    ungroup() %>% 
    filter(genotyper %in% genotypers) %>% # Limit to only genotypers specified
    mutate(present = map_lgl(allele, function(x) !all(is.na(x)))) %>% # ID predictions that only generated NAs
    mutate(accuracy = ifelse(present == T, accuracy, NA)) %>% # Convert 0 accuracy to NA based on previous line
    select(sample, locus, field, genotyper, accuracy, present) %>%  
    # Pivot wider then longer to insert NAs at missing combinations
    pivot_wider(names_from = "genotyper", values_from = c("accuracy", "present")) %>% 
    pivot_longer(contains("accuracy"), names_to = "genotyper", values_to = "accuracy", names_prefix = "accuracy_") %>% 
    mutate_if(is.character, factor) %>% 
    mutate_if(is.logical, factor)
  
  # If using to fit the tree model, keep only the best predictions to train the model on
  # Find the maximum accuracy value for a sample-locus-field combo and keep only those genotypers acheiving that max
  if (for_fit == T){
    output <- output %>%      
      group_by(sample, locus, field) %>% 
      mutate(max_accuracy = ifelse(all(is.na(accuracy)), NA, max(accuracy, na.rm = T)))%>%
      ungroup() %>%
      filter(accuracy == max_accuracy & !is.na(max_accuracy)) %>% 
      select(-accuracy, -max_accuracy) 
  }
  
  # Convert sample to rownames, add index to make unique
  output <- output %>% 
    mutate(sample = paste(sample, 1:n(), sep = "_")) %>% 
    column_to_rownames("sample") %>% 
    
    return(output)
}
# accuracy_df %>% make_tree_input()

# Fits a classification decision tree to output of make_tree_input
# Outputs a list with several components
#     "top_models_stats" = best results of grid search
#     "best_model_id" = the actual best result of the grid search
#     "model" = this best model fit to the training set
#     "testing_output" = predictions of the model on the test set
#     "full_output" = predictions of the model on the training + test set
make_tree <- function(df,
                      train_split = 0.7, 
                      outcome = "genotyper",
                      genotypers = c("arcasHLA", "optitype", "phlat"),
                      grid_search_levels = 4,
                      n_cores = 1,
                      seed1=123,seed2=234,seed3=345){
  # browser()
  df <- make_tree_input(df, for_fit = T, genotypers = genotypers)
  
  formula_string <- sprintf("%s ~ .", outcome)
  
  # Create train/test
  set.seed(seed1)
  df_split <- initial_split(df, prop =train_split, strata = "genotyper")
  df_train <- training(df_split)
  df_test <- testing(df_split)
  
  # Create training cross-folds
  set.seed(seed2)
  df_folds <- vfold_cv(df_train)
  
  # Set tree parameters
  tree_spec <- decision_tree(
    cost_complexity = tune(),
    tree_depth = tune(),
    min_n = tune()
  ) %>%
    set_engine("rpart") %>%
    set_mode("classification")
  
  # Set tuning parameter search
  tree_grid <- grid_regular(cost_complexity(), tree_depth(), min_n(), levels = grid_search_levels)
  
  # Grid search
  doParallel::registerDoParallel(cores = n_cores)
  formula_string <- sprintf("%s ~ .", "genotyper")
  set.seed(seed3)
  tree_rs <- tune_grid(
    tree_spec,
    as.formula(formula_string),
    resamples = df_folds,
    grid = tree_grid,
    metrics = metric_set(accuracy, sens, spec, roc_auc)
  )
  
  # Select best model parameters
  final_tree <- finalize_model(tree_spec, select_best(tree_rs, "roc_auc"))
  final_fit <- fit(final_tree, as.formula(formula_string), df_train)
  
  # Predict outcomes on test set
  testing_output <- df_test %>% 
    bind_cols(
      final_fit %>% predict(new_data = df_test),
      final_fit %>% predict(new_data = df_test, type = 'prob') 
    )
  
  # Calculate AUC
  prediction_var <- setdiff(names(testing_output)[grepl(".pred", names(testing_output))], ".pred_class")
  if (length(prediction_var)==2){prediction_var <- prediction_var[1]}
  auc <- testing_output %>% roc_auc(truth = genotyper, contains(prediction_var))
  
  # Compile results in list
  final_list <- list(
    "top_models_stats" = show_best(tree_rs, "roc_auc"),
    "best_model_stats" = select_best(tree_rs, "roc_auc"),
    "model" = final_fit,
    "auc" = auc
  )
  return(final_list)
}


# # Filter accuracy DF using predictions from classification tree made by make_tree
# # Accepts input from accuracy functions like calculate_all_hla_accuracy
# filter_by_tree <- function(df, model_list, name = "composite"){
#   key <- model_list$full_predictions %>% 
#     rownames_to_column("sample") %>% 
#     separate(sample, into = c("sample", NA), sep = "_") %>% 
#     select(sample, locus, field, .pred_class) %>% 
#     distinct()
#   output <- df %>% 
#     left_join(key, by = c("sample", "locus", "field")) %>% 
#     filter(genotyper == .pred_class) %>% 
#     mutate(genotyper = name) %>% 
#     select(-.pred_class)
#   return(output)
# }

# Takes and accuracy DF an applies the tree model
# Filters best genotypers and renames as composite
apply_tree <- function(df, model_list, name = "Composite"){
  # Convert to model format
  df_model <- df %>% 
    make_tree_input(for_fit = F) 
  # return(df_model)
  # Generate predictions from model
  df_model <- df_model %>% 
    bind_cols(
      model_list$model %>% predict(new_data = df_model),
      model_list$model %>% predict(new_data = df_model, type = 'prob') 
    )
  # Pull best prediction based on model
  key <- df_model %>% 
    rownames_to_column("sample") %>% 
    separate(sample, into = c("sample", NA), sep = "_") %>% 
    select(sample, locus, field, .pred_class) %>% 
    distinct()
  # Filter original data based on best predictions and rename as composite
  df <- df %>%     
    left_join(key, by = c("sample", "locus", "field")) %>% 
    filter(genotyper == .pred_class) %>% 
    mutate(genotyper = name) %>%
    select(-.pred_class)
  return(df)
}

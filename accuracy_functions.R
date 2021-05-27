library(tidyverse)
library(checkmate)


### Allele matching function
allele_match <- function(ref,comp,ground_truth=T, verbose=F){
  display_message <-  ifelse(verbose == F, "message", "none")
  
  suppressMessages({
    # If allele ends with an NA string, convert entire object to NA
    ref[grepl("NA$",ref)] <- NA
    comp[grepl("NA$",comp)] <- NA
    
    # Count non-NA alleles
    n_ref <- sum(!is.na(ref))
    n_comp <- sum(!is.na(comp))
    
    # Return NA if all alleles are NA
    if(n_ref==0 & n_comp==0){message("No alleles for comparison");return(NA)}
    # If ref is ground_truth, return NA if all ref alleles are NA
    if(n_ref==0 & ground_truth==T){message("No ground truth alleles");return(NA)}
    # Cannot have more than two alleles
    if(n_ref > 2 | n_comp > 2){message("Invalid number of alleles");return(NA)}
    
    
    # If comparing to a ground_truth, only penalize NAs in comp, not those in ref
    if (n_ref == 1 & ground_truth == T){return(sum(any(ref == comp, na.rm = T)))}
    
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


compare_hla<- function(hla_df, reference = "invitro", exclude_missing=T, compare_function = "accuracy", na_drop = T){
  # Input checks
  arg_col <- makeAssertCollection()
  hla_df_columns <- c("sample", "allele_id", "locus", "field_1", "field_2", "field_3", "genotyper")
  assertNames(names(hla_df), permutation.of = hla_df_columns, .var.name = "hla_df column names", add = arg_col)
  assertChoice(reference, unique(hla_df$genotyper), add = arg_col)
  assertLogical(exclude_missing, add = arg_col)
  assertChoice(compare_function, c("accuracy", "agreement"), add = arg_col)
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
      filter_at(vars(c(reference, allele)), function(x) !is.na(x)) %>% 
      mutate(n_reference = map_dbl(reference, function(x) sum(!is.na(x) & !grepl("NA", x)))) %>% 
      mutate(n_allele = map_dbl(allele, function(x) sum(!is.na(x) & !grepl("NA", x)))) %>% 
      filter_at(vars(contains("n_")), function(x) x ==2)
  }
  if (compare_function == "accuracy"){
    df <- df %>% mutate(accuracy = map2_dbl(reference, allele, allele_match))
  } else if (compare_function == "agreement"){
    df %>% mutate(accuracy = map2_dbl(reference, allele, allele_agreement))
  } 
  if (na_drop == T){
    df <- df %>% drop_na()
  }
  df
}
# compare_hla(hla_df = all_hla, reference = "invitro", exclude_missing = T, compare_function = "accuracy")
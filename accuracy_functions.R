library(tidyverse)
library(checkmate)

### Allele matching function
allele_match <- function(ref,comp){
  n_ref <- sum(!is.na(ref))
  n_comp <- sum(!is.na(comp))
  if (n_comp == 1){
    comp <- c(comp, NA)
  }
  # No meaningful comparison if reference is NA
  if (all(is.na(ref))){
    NA 
    # Cannot have more than two alleles
  } else if (n_ref >2 | n_comp > 2){
    return(NA)
    print("Invalid number of alleles")
    # If only one allele typed in reference, only need one match in comparison
  } else if (n_ref == 1){
    sum(any(ref == comp, na.rm = T))
    # Matches when all alleles present
    # Forward and reverse match search since order does not matter
  } else if (n_ref == 2){
    max(sum(ref == comp, na.rm = T)/2, sum(ref == rev(comp), na.rm = T)/2)
  }
}

# unit_test <- bind_rows(
#   as_tibble(list(reference = list(c("03", "03")), comparison = list(c("03", "03")), expected_accuracy = 1.0)),
#   as_tibble(list(reference = list(c("03", "01")), comparison = list(c("03", "01")), expected_accuracy = 1)),
#   as_tibble(list(reference = list("03"), comparison = list(c("03", "03")), expected_accuracy = 1.0)),
#   as_tibble(list(reference = list(c("03", "01")), comparison = list(c("03", "03")), expected_accuracy = 0.5)),
#   as_tibble(list(reference = list(c("03", "03")), comparison = list(c("03", "01")), expected_accuracy = 0.5)),
#   as_tibble(list(reference = list(c("03", "01")), comparison = list(c(NA, "01")), expected_accuracy = 0.5)),
#   as_tibble(list(reference = list(c("03", "03")), comparison = list("03"), expected_accuracy = 0.5)),
#   as_tibble(list(reference = list(c("01", NA)), comparison = list(c(NA, NA)), expected_accuracy = 0.0)),
#   as_tibble(list(reference = list(c("01")), comparison = list(c(NA, NA)), expected_accuracy = 0.0)),
#   as_tibble(list(reference = list(c("01", NA)), comparison = list(c("01", NA)), expected_accuracy = 1.0)),
#   as_tibble(list(reference = list(c("01")), comparison = list(c("01", NA)), expected_accuracy = 1.0)),
#   as_tibble(list(reference = list(c("03", "03")), comparison = list(c("01", "01")), expected_accuracy = 0.0)),
#   as_tibble(list(reference = list(c("03", "01")), comparison = list(NA), expected_accuracy = 0.0)),
#   as_tibble(list(reference = list(c(NA, NA)), comparison = list(c(NA, NA)), expected_accuracy = NA)),
#   as_tibble(list(reference = list(c(NA, NA)), comparison = list(NA), expected_accuracy = NA)),
#   as_tibble(list(reference = list(NA), comparison = list(NA), expected_accuracy = NA)),
#   as_tibble(list(reference = list(c("03","03","03")), comparison = list(c("01", "01")), expected_accuracy = NA))
# )
# 
# unit_test_output <- unit_test %>% 
#   mutate(calculated_accuracy = map2_dbl(reference, comparison, allele_match))
# unit_test_output


compare_hla<- function(hla_df, reference = "invitro", exclude_missing=T, compare_function = "accuracy"){
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
    mutate(field_2 = paste(field_1, field_2, sep = "_"), field_3 = paste(field_2, field_3, sep = "_")) %>%
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
    df %>% 
      mutate(accuracy = map2_dbl(reference, allele, allele_match)) %>% 
      drop_na()
  } else if (compare_function == "agreement"){
    df %>%
      mutate(accuracy = map2_dbl(reference, allele, allele_agreement))%>%
      drop_na()
  } 
}
# compare_hla(hla_df = all_hla, reference = "invitro", exclude_missing = T, compare_function = "accuracy")
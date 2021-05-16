library(tidyverse)
library(flextable)


### Expects format from accuracy_functions::compare_hla
### Summarizes accuracy across locus, field, and genotyper
calculate_summary_df <- function(df){
  suppressMessages({df %>% 
      group_by(locus, field, genotyper) %>% 
      summarise(mean_accuracy = mean(accuracy),
                sd = sd(accuracy),
                se = sd(accuracy)/sqrt(n())) %>% 
      ungroup() %>% 
      mutate(mean_accuracy = ifelse(mean_accuracy == 0, NA, mean_accuracy), 
             sd = ifelse(mean_accuracy == 0, NA, sd),
             se = ifelse(mean_accuracy == 0, NA, se),
             locus = factor(locus, levels = sort(unique(locus), decreasing = T)))
  })
}

### Prettifies and orders HLA fields
reformat_hla_field <- function(field_vector){
  suppressWarnings(
    field_vector %>% 
      factor() %>% 
      fct_recode(
        "Field 1" = "field_1",
        "Field 2" = "field_2",
        "Field 3" = "field_3"
      ) %>% 
      fct_relevel(
        "Field 1",
        "Field 2",
        "Field 3"
      ))
}

### Prettifies and orders HLA genotypers
reformat_hla_genotyper <- function(genotyper_vector){
  suppressWarnings(
    genotyper_vector %>% 
      factor() %>% 
      fct_recode(
        "HLAminer" = "hlaminer",
        "OptiType" = "optitype",
        "PHLAT" = "phlat"
      ) %>% 
      fct_relevel(
        "arcasHLA",
        "OptiType",
        "PHLAT",
        "HLAminer"
      ))
}

gg_hla_prediction_frequency <- function(df, loci = NULL, color = "viridis"){
  if (is.null(loci)){
    loci <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1","DRB3","DRB4","DRB5")
  }
  df %>% 
    mutate(field = reformat_hla_field(field)) %>% 
    mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
    filter(frequency <= 1) %>% 
    filter(locus %in% loci) %>% 
    mutate(count = frequency*2) %>% 
    ggplot(aes(x = genotyper))+
    geom_bar(aes(fill = factor(count)), position = "fill", color = "black", size = 0.25) +
    facet_grid(locus ~ field)+
    theme_bw()+
    scale_fill_viridis_d(option = color)+
    scale_y_continuous(breaks = c(0,0.5,1.0))+
    labs(x="",y="Proportion of samples",fill="Number of \npredicted \nalleles") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background.x = element_blank())
}

### Expects format from figure_format_functions::calculate_summary_df
### Creates heatmap of mean accuracy
### Colors refer to those within the viridis package
gg_summary_hla_accuracy <- function(df, color_label = "Accuracy", color = "plasma"){
  df %>% 
    mutate(field = reformat_hla_field(field)) %>% 
    mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
    mutate(class = ifelse(grepl("^[ABC]",locus),"mhc1","mhc2")) %>% 
    ggplot(aes(x = genotyper, y = locus)) +
    geom_tile(aes(fill = mean_accuracy), color = "white")+
    # geom_text(aes(label = round(mean_accuracy, digits = 1)))+
    # ggrepel::geom_text_repel(aes(label = round(mean_accuracy, digits = 2)), 
    #                          size = 3.5, 
    #                          bg.color = "white",
    #                          force = 0,
    #                          force_pull = 0)+
    facet_grid(class~field, scales = "free_y", space = "free_y")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background.x = element_blank(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank()) +
    scale_fill_viridis_c(option = color, limits = c(0,1), na.value = "transparent") +
    labs(x = "", y = "", fill = color_label)
}

### Expects format from figure_format_functions::calculate_summary_df
### Creates a table of mean accuracy and standard error
flex_summary_hla_accuracy <- function(df){
  df <- df %>% 
    mutate(field = reformat_hla_field(field)) %>% 
    mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
    mutate(cell_value = sprintf("%s %s %s", round(mean_accuracy,2),"\u00B1",round(se,2) )) %>%
    mutate(cell_value = ifelse(grepl("NA", cell_value), NA, cell_value)) %>% 
    select(-sd,-se,-class, -mean_accuracy) %>% 
    pivot_wider(names_from = c("field","genotyper"), values_from = "cell_value", names_sep = "-") 
  df_key <- suppressWarnings({tibble(col_keys = names(df)) %>% 
    separate(col_keys, into = c("field", "genotyper"), sep = "-", remove = F) %>%
    drop_na()})
  df %>% 
    flextable() %>% 
    colformat_char(na_str = "---") %>% 
    set_header_df(mapping = df_key, key = "col_keys") %>% 
    merge_h(part = "header") %>% 
    theme_vanilla() %>% 
    vline(j=c(1,5,9), border = fp_border_default()) %>% 
    fix_border_issues()
}

### Expects format from accuracy_functions::compare_hla
### Plots accuracy by genotyper and allele, with information of total times allele was seen
### Only displays one HLA field, set by field_selection
gg_allele_hla_accuracy <- function(df, color_label = "Accuracy", field_selection = "field_2", color = "plasma"){
  suppressMessages({
    accuracy_df %>% 
      unnest(reference) %>% 
      filter(field == field_selection) %>%
      filter(allele != "NA") %>% 
      mutate(match = map2_dbl(reference, allele, function(x,y) sum(x %in% y, na.rm = T))) %>% 
      group_by(locus, reference, genotyper) %>% 
      summarise(accuracy = mean(match), n = n()) %>%
      ggplot(aes(x=reference, y = genotyper))+
      geom_point(aes(fill = accuracy, size = n), shape = 21, color = "black")+
      facet_wrap(~ locus, ncol = 1, scales = "free_x",  strip.position="right")+
      scale_fill_viridis_c(option = color) +
      scale_size_continuous(breaks = c(5,25,50,100), range = c(1,10)) +
      theme_bw()+
      labs(x = "", y = "", fill = color_label, size = "Count of allele \nin dataset")+
      theme(axis.text.x = element_text(angle= 45, hjust = 1))
  })
}

library(tidyverse)
library(flextable)


### Prettifies and orders HLA fields
reformat_hla_field <- function(field_vector, reverse = F){
  suppressWarnings(
    output <- field_vector %>% 
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
  if (reverse == T){output <- fct_rev(output)}
  return(output)
}

### Prettifies and orders HLA genotypers
reformat_hla_genotyper <- function(genotyper_vector, reverse = F){
  suppressWarnings(
    output <- genotyper_vector %>% 
      factor() %>% 
      fct_recode(
        "Ground truth" = "invitro",
        "HLAminer" = "hlaminer",
        "OptiType" = "optitype",
        "PHLAT" = "phlat",
        "Composite" = "composite"
      ) %>% 
      fct_relevel(
        "Composite",
        "Ground truth",
        "arcasHLA",
        "OptiType",
        "PHLAT",
        "HLAminer"
      ))
  if (reverse == T){output <- fct_rev(output)}
  return(output)
}

### Prettifies and orders HLA loci
reformat_hla_loci <- function(loci_vector, reverse = F){
  suppressWarnings(
    output <- loci_vector %>% 
      factor() %>% 
      fct_recode(
        "All MHC" = "MHC All"
      ) %>% 
      fct_relevel(
        "All MHC",
        "MHC II",
        "MHC I"
      ))
  if (reverse == T){output <- fct_rev(output)}
  return(output)
}


### Expects format from calculation_functions::allele_tally
### Plots relative proportion of 0,1,2 allele predictions across all samples as a bar chart
gg_hla_prediction_tally <- function(df, loci = NULL, color = "viridis", reverse = F){
  if (is.null(loci)){
    loci <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1","DRB3","DRB4","DRB5")
  }
  viridis_direction <- ifelse(reverse == T, -1, 1)
  df %>% 
    ungroup() %>% 
    mutate(field = reformat_hla_field(field)) %>% 
    mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
    # filter(frequency <= 1) %>%
    filter(locus %in% loci) %>% 
    mutate(count = count) %>% 
    ggplot(aes(x = genotyper))+
    geom_bar(aes(fill = factor(count)), position = "fill", color = "black", size = 0.25) +
    facet_grid(locus ~ field)+
    theme_bw()+
    scale_fill_viridis_d(option = color, direction = viridis_direction)+
    scale_y_continuous(breaks = c(0,0.5,1.0))+
    labs(x="",y="Proportion of samples",fill="Number \nof alleles") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background.x = element_blank())
}

### Expects format from calculation_functions::calculate_summary_df
### Creates heatmap of mean accuracy
### Colors refer to those within the viridis package
gg_summary_hla_accuracy <- function(df, color_label = "Accuracy", color = "plasma", reverse = F){
  viridis_direction <- ifelse(reverse == T, -1, 1)
  df %>% 
    mutate(field = reformat_hla_field(field)) %>% 
    mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
    mutate(locus = reformat_hla_loci(locus)) %>% 
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
    scale_fill_viridis_c(option = color, 
                         limits = c(0,1), 
                         na.value = "transparent", 
                         direction = viridis_direction) +
    labs(x = "", y = "", fill = color_label)
}

### Expects format from calculation_functions::calculate_summary_df
### Creates a table of mean accuracy and standard error
flex_summary_hla_accuracy <- function(df){
  df <- df %>% 
    mutate(field = reformat_hla_field(field)) %>% 
    mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
    mutate(cell_value = sprintf("%s %s %s", round(mean_accuracy,2),"\u00B1",round(se,2) )) %>%
    mutate(cell_value = ifelse(grepl("NA", cell_value), NA, cell_value)) %>% 
    select(-sd,-se, -mean_accuracy) %>% 
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

### Expects format from calculation::compare_hla
### Plots accuracy by genotyper and allele, with information of total times allele was seen
### Only displays one HLA field, set by field_selection
gg_allele_hla_accuracy <- function(df, color_label = "Accuracy", field_selection = "field_2", color = "plasma"){
  suppressMessages({
    df %>% 
      unnest(reference) %>% 
      filter(field == field_selection) %>%
      filter(allele != "NA") %>% 
      mutate(match = map2_dbl(reference, allele, function(x,y) sum(x %in% y, na.rm = T))) %>% 
      group_by(locus, reference, genotyper) %>% 
      summarise(accuracy = mean(match), n = n()) %>%
      mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
      mutate(locus = reformat_hla_loci(locus)) %>% 
      ggplot(aes(x=genotyper, y = reference))+
      geom_point(aes(fill = accuracy, size = n), shape = 21, color = "black")+
      facet_wrap(~ locus, nrow = 1, scales = "free_y",)+
      scale_fill_viridis_c(option = color) +
      scale_size_continuous(breaks = c(5,25,50,100), range = c(1,10)) +
      theme_bw()+
      labs(x = "", y = "", fill = color_label, size = "Count of allele \nin dataset")+
      theme(axis.text.x = element_text(angle= 45, hjust = 1))
  })
}

gg_individual_hla_accuracy <- function(df, color_label = "Accuracy", field_selection = "field_2", color = "plasma"){
  suppressMessages({
    df %>% 
      filter(field == field_selection) %>% 
      mutate(field = reformat_hla_field(field)) %>% 
      mutate(genotyper = reformat_hla_genotyper(genotyper, reverse = T)) %>% 
      mutate(locus = reformat_hla_loci(locus)) %>% 
      ggplot(aes(x = sample, y = genotyper, fill = accuracy)) + 
      geom_tile() +
      theme_classic() +
      facet_grid(locus ~.) +
      scale_fill_viridis_c(na.value = "transparent", option = color) +
      theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
      labs(x = "Sample", y = "", fill = color_label)
  })
}

gg_multilevel_roc <- function(df){
  roc_plot <- df %>% 
    roc_curve(truth = copy_number, .pred_0:.pred_2) %>% 
    autoplot()
  ggplot_build(roc_plot)$data[[1]] %>% 
    mutate(n_alleles = fct_recode(PANEL, "0"="1", "1"="2", "2"="3")) %>% 
    ggplot(aes(x=x,y=y,color=n_alleles))+
    geom_path()+
    geom_abline(slope=1,intercept=0, linetype="dotted")+
    theme_bw() +
    scale_color_brewer(palette = "Set1") +
    labs(x = "1 - Specificity", y = "Sensitivity", color = "Number\nof alleles")
}


# Expects output of create_drb345_df
# A baseline set can be provided that adds a grey reference layer of matching visualizations (used for PMID data)
gg_drb345_ratio <- function(sample_drb345_df, baseline_drb345_df = NULL){
  sample_drb345_df<- sample_drb345_df %>% 
    pivot_longer(DRB3:DRB5, names_to='locus_2', values_to='drb345_ratio') %>% 
    filter(locus == locus_2) %>% 
    mutate(copy_number = factor(copy_number))
  
  if (!is.null(baseline_drb345_df)){
    baseline_drb345_df <- baseline_drb345_df %>%
      pivot_longer(DRB3:DRB5, names_to='locus_2', values_to='drb345_ratio') %>%
      filter(locus == locus_2) %>%
      mutate(copy_number = factor(copy_number))
  }
  
  sample_drb345_layers <- list(
    geom_density(alpha = 0.1, aes(y = ..scaled.., color = copy_number, fill = copy_number)),
    geom_point(aes(y=0), color = "black", size = 1.5),
    geom_point(aes(y=0, color = copy_number), size = 1)
  )
  
  baseline_drb345_layers <- list(
    geom_density(data = baseline_drb345_df, aes(y = ..scaled..), fill = NA, color = "grey50"),
    geom_point(data = baseline_drb345_df, aes(y=1), color = "grey25", size = 0.75),
    geom_point(data = baseline_drb345_df, aes(y=1), color = "grey75", size = 0.25)
  )
  
  formatting_layers <- list(
    theme_bw(),
    scale_fill_brewer(palette = "Set1"),
    scale_color_brewer(palette = "Set1"),
    facet_grid(copy_number ~ locus, scales = "free"),
    theme(legend.position = "none"),
    labs(x = "DRB345 / DRB1 ratio", y = "Density")
  )
  
  plt <- ggplot(data = sample_drb345_df, aes(x=drb345_ratio))
  
  if (!is.null(baseline_drb345_df)){
    plt <- plt + baseline_drb345_layers
  }
  
  plt <- plt +
    sample_drb345_layers +
    formatting_layers
  
  
  print(plt)
}

# Plot violin-jitter plots of sequencing pipeline task runtimes
# Input is df with at least columns of "sample", "component", and "process_time"
# Can take output of parse_log_time
gg_runtime <- function(df, include_tasks = NULL){
  # Set default tasks
  if (is.null(include_tasks)){
    include_tasks <- c("HISAT-AND-SORT", "ARCAS-EXTRACT", "ARCAS-GENOTYPE", 
                       "PHLAT", "OPTITYPE", "HLAMINER")}
  
  # Input checks
  arg_col <- makeAssertCollection()
  possible_tasks <- c(
    "FASTQ-COMPILE", "PRE-FASTQC", "TRIM", "POST-FASTQC", "HISAT-AND-SORT",
    "INDEX", "ARCAS-EXTRACT", "ARCAS-GENOTYPE", "PHLAT", "OPTITYPE","HLAMINER"
  )
  df_columns <- c("sample", "component", "process_time")
  assertNames(include_tasks, subset.of = possible_tasks, .var.name = "possible tasks", add = arg_col)
  assertNames(colnames(df), must.include = df_columns, .var.name = "possible tasks", add = arg_col)
  if (arg_col$isEmpty()==F) {map(arg_col$getMessages(),print);reportAssertions(arg_col)}
  
  # Create plot
  plt <- df %>%
    select(sample, component, process_time) %>%
    mutate(process_time = as.numeric(process_time)) %>%
    filter(component %in% include_tasks)  %>%
    mutate(component = factor(component, levels = c(
      "FASTQ-COMPILE", "PRE-FASTQC", "TRIM", "POST-FASTQC", "HISAT-AND-SORT",
      "INDEX", "ARCAS-EXTRACT", "ARCAS-GENOTYPE", "PHLAT", "OPTITYPE","HLAMINER"
    ))) %>%
    drop_na() %>%
    ggplot(aes(x = component, y = process_time))+
    geom_violin(scale = "width")+
    geom_jitter(size = 0.1, alpha = 0.2, width = 0.25) +
    geom_boxplot(width = 0.1, outlier.size = 0)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Process time (hr)")
  
  return(plt)
}






### OLD WORK

### Replaced by gg_hla_predicition_tally
### Older version that plotted allele frequency rather than count
gg_hla_prediction_frequency <- function(df, loci = NULL, color = "viridis", reverse = F){
  if (is.null(loci)){
    loci <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1","DRB3","DRB4","DRB5")
  }
  viridis_direction <- ifelse(reverse == T, -1, 1)
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
    scale_fill_viridis_d(option = color, direction = viridis_direction)+
    scale_y_continuous(breaks = c(0,0.5,1.0))+
    labs(x="",y="Proportion of samples",fill="Number \nof alleles") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background.x = element_blank())
}


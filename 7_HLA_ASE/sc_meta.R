library(tidyverse)
library(meta)
library(here)

source(here("helper_functions/data_import_functions.R"))
scHLAcount_dir <- here("data/isb/scHLAcount/")

srt <- readRDS("/labs/khatrilab/hongzheng/webapps/isb/srt_isb260.meta.rds")
cells <- c("CD14 Monocyte", "CD4 T", "CD8 T", "NK", "B", "CD16 Monocyte", "cDC", "pDC")

# Function to perform meta-analysis on dataframe where
# each row is a cell and columns:
# `observed` (reads of dominant allele)
# `gene_sum_typed` (total cell reads)
# `expected` (reads of one allele expected if 50:50 allele_1 : allele_2)
sc_meta <- function(df){
  l <- nrow(df)
  m <- metabin(event.e = observed, n.e = gene_sum_typed, event.c = expected, n.c = gene_sum_typed,
               data = df,
               method = "Inverse",
               incr = 0.1,
               sm = "OR")
  data.frame(summary(m)$random) %>% 
    mutate(n_cells = l) %>% 
    select(TE, seTE, lower, upper, n_cells) 
}

# Import data based on sample and genotyper
cell_stats <- expand_grid(
  genotyper = c("invitro", "arcasHLA", "optitype", "phlat", "hlaminer"),
  sample = read_lines(sprintf("%s/BL_fastq_files.txt", scHLAcount_dir))
  ) %>%
  # head(5) %>% # Specify limit to number of meta-analyses
  mutate(data = map2(sample, genotyper, function(s,g){
    result_path = sprintf("%s/output/%s", scHLAcount_dir, g)
    barcode_path = sprintf("%s/barcodes", scHLAcount_dir)
    scHLA_data_processing(sample = s,
                          result_dir = result_path,
                          barcode_dir = barcode_path)
  })) %>% unnest(data) %>% 
  filter(!is.na(cell)) %>% 
  mutate(sample = gsub("_[A-Z][0-9]$","",sample)) %>% # Consolidate samples
  left_join(srt %>% select(celltype, cell), by = "cell") %>% # Add celltypes
  filter(celltype %in% cells) # Keep only standard cell types


meta_df <- cell_stats %>% 
  # Keep only most expressed allele (or random if 50:50)
  group_by(sample, genotyper, gene, cell) %>% 
  slice_max(order_by = allele_ratio, n = 1, with_ties = F) %>% 
  ungroup() %>% 
  # Fill out contingency table
  mutate(complement = gene_sum_typed - count, 
         expected = 0.5*gene_sum_typed) %>% 
  rename(observed = count) %>% 
  select(sample, genotyper, celltype, gene, observed, expected, gene_sum_typed) %>% 
  group_by(sample, genotyper, celltype, gene) %>%
  # Nest and run meta-analysis
  nest() %>%
  ungroup() %>% 
  mutate(data = map(data,function(x) {sc_meta(x)})) %>%
  unnest(data)

write_csv(meta_df, here("7_HLA_ASE/meta_analysis_results_2.csv"))




  
#   
# plan(multisession, workers = 3)
# tictoc::tic()
# nothing <- 1:48 %>% future_map(., function(x) Sys.sleep(1))
# tictoc::toc()
# 
# ggplot(z %>% unnest(n_cells), aes(x = celltype, y = TE)) + 
#   # geom_point() + 
#   geom_violin() + 
#   facet_grid(genotyper ~ gene) +
#   stat_summary(fun = mean, geom="point")+
#   stat_summary(fun.data = mean_se, geom="errorbar")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# plan(multicore, workers = 10);tictoc::tic();nothing <- 1:10 %>% future_map(., function(x) Sys.sleep(1));tictoc::toc()
# 
# plan(multicore, workers = 10);tictoc::tic();test_output <- test_meta %>% mutate(data = future_map(data,function(x) {sc_meta(x)})) %>% unnest(data) ;tictoc::toc()
# plan(sequential);tictoc::tic();test_output <- test_meta %>% mutate(data = map(data,function(x) {sc_meta(x)})) %>% unnest(data);tictoc::toc()
# plan(multicore, workers = 10);tictoc::tic();test_output <- test_meta %>% mutate(data = map(data,function(x) {nrow(x)})) %>% unnest(data) ;tictoc::toc()
# 
# plan(multicore, workers = 10);tictoc::tic();test_output <- test_meta %>% head(10) %>% mutate(data = future_map(data,function(x) {a <- nrow(x); b <- a*0; c <- b+1; Sys.sleep(c)})) %>% unnest(data) ;tictoc::toc()
# plan(sequential);tictoc::tic();test_output <- test_meta %>% head(10) %>% mutate(data = future_map(data,function(x) {a <- nrow(x); b <- a*0; c <- b+1; Sys.sleep(c)})) %>% unnest(data) ;tictoc::toc()
# 
# sc_meta <- function(df){
#   m <- metabin(event.e = observed, n.e = gene_sum_typed, event.c = expected, n.c = gene_sum_typed,
#                data = df,
#                method = "Inverse",
#                incr = 0.1,
#                sm = "OR")
#   data.frame(summary(m)$random) %>% 
#     select(TE, seTE, lower, upper) 
#   # %>% 
#   # mutate_all(exp) # Convert from log-odds to odds
# }
# 
# sc_meta_future <- function(df){
#   metabin(event.e = observed, n.e = gene_sum_typed, event.c = expected, n.c = gene_sum_typed,
#                data = df,
#                method = "Inverse",
#                incr = 0.1,
#                sm = "OR")
# }
# 
# plan(multicore, workers = 10);tictoc::tic();test_output <- future_map(test_meta$data,function(x) {sc_meta_future(x)});tictoc::toc()
# 
# plan(multicore, workers = 10);tictoc::tic();test_output <- future_map(test_meta$data,function(x) {metabin(event.e = x$observed, n.e = x$gene_sum_typed, event.c = x$expected, n.c = x$gene_sum_typed,
#                                                                                                           method = "Inverse",
#                                                                                                           incr = 0.1,
#                                                                                                           sm = "OR")});tictoc::toc()
# 
# 
# 
# a

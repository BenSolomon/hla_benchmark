library(tidyverse)
library(cowplot)
library(here)
source(here("helper_functions/data_import_functions.R"))
source(here("helper_functions/figure_format_functions.R"))
source(here("helper_functions/calculation_functions.R"))
isb_path <- here("data/isb")
isb_samples <- read_tsv(here("data/isb/logs/parallel.log")) %>% 
  separate(Command, into = c(NA, "sample"), sep = " ") %>% 
  pull(sample) %>% unique()

df <- combine_HLA_import(path = isb_path, samples = isb_samples, expand_invitro = T)

df <- df %>% 
  separate(allele, into = c("locus"), sep = "\\*", extra = "drop", remove = F) %>% 
  filter(
    grepl("(AC|BL)$", sample), # Only AC, BL samples
    grepl("^[ABCD]", locus), # Only HLA-A,B,C,D genes
    !grepl("^(DM|DO|DRA)|[67]$",locus) # Further limit HLA-D genes
  ) %>% 
  pivot_wider(names_from = c("locus", "allele_id"), values_from = "allele", values_fill = NA, values_fn = function(x) paste(x, collapse = "; ")) %>% 
  mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
  arrange(genotyper, sample) %>% 
  mutate_all(na_if, "NULL")

write_csv(df, here("data/isb/compiled_isb_genotypes.csv"))

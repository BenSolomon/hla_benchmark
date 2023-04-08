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

invitro_hla <- invitro_import(path = here("data/pmid30518681/pmid30518681_scisco_hla.tsv"), 
                              exclude_comment_samples = T) %>% 
  mutate(sample = gsub("STU-HD-","",sample))

invitro_hla_c <- invitro_hla %>% filter(sample == "C")

# Expand results to technical and time replicates
invitro_hla <- bind_rows(
  invitro_hla,
  invitro_hla_c %>% mutate(sample = "Ck"),
  invitro_hla_c %>% mutate(sample = "C1"),
  invitro_hla_c %>% mutate(sample = "C2")
)

sample_metadata <- read_csv(here("data/pmid30518681/pmid30518681_samples.csv"))
pmid_3p_srr <- sample_metadata %>% filter(platform == "scRNA") %>% pull(Run)
pmid_3p_samples <- sample_metadata %>% filter(platform == "scRNA") %>% pull(sample)
pmid_bulk_srr <- sample_metadata %>% filter(platform == "bulkRNA") %>% pull(Run)
pmid_bulk_samples <- sample_metadata %>% filter(platform == "bulkRNA") %>% pull(sample)

### 3'
pmid_3p_hla <- combine_HLA_import(
  path = here("data/pmid30518681/scRNA"),
  samples = pmid_3p_srr,
  invitro_path = NULL) %>%  # Will manually add invitro
  rename(Run = sample) %>% 
  left_join(sample_metadata, by = "Run") %>% 
  mutate(sample = ifelse(is.na(sample), original_name, sample)) %>% 
  select(sample, allele_id, allele, genotyper) %>% 
  bind_rows(invitro_hla) %>% 
  filter(sample %in% pmid_3p_samples)
pmid_3p_hla

pmid_3p_hla %>% 
  separate(allele, into = c("locus"), sep = "\\*", extra = "drop", remove = F) %>% 
  filter(
    grepl("^[ABCD]", locus), # Only HLA-A,B,C,D genes
    !grepl("^(DM|DO|DRA)|[67]$",locus) # Further limit HLA-D genes
  ) %>% 
  pivot_wider(names_from = c("locus", "allele_id"), values_from = "allele", values_fill = NA, values_fn = function(x) paste(x, collapse = "; ")) %>% 
  mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
  arrange(genotyper, sample) %>% 
  mutate_all(na_if, "NULL") %>% 
  write_csv(here("data/pmid30518681/compiled_nih_3p_genotypes.csv"))
### Bulk
pmid_bulk_hla <- combine_HLA_import(
  path = here("data/pmid30518681/bulkRNA"),
  samples = pmid_bulk_srr,
  invitro_path = NULL) %>%  # Will manually add invitro
  rename(Run = sample) %>% 
  left_join(sample_metadata, by = "Run") %>% 
  mutate(sample = ifelse(is.na(sample), original_name, sample)) %>% 
  select(sample, allele_id, allele, genotyper) %>% 
  bind_rows(invitro_hla) %>% 
  filter(sample %in% pmid_bulk_samples)
pmid_bulk_hla


pmid_bulk_hla %>% 
  separate(allele, into = c("locus"), sep = "\\*", extra = "drop", remove = F) %>% 
  filter(
    grepl("^[ABCD]", locus), # Only HLA-A,B,C,D genes
    !grepl("^(DM|DO|DRA)|[67]$",locus) # Further limit HLA-D genes
  ) %>% 
  pivot_wider(names_from = c("locus", "allele_id"), values_from = "allele", values_fill = NA, values_fn = function(x) paste(x, collapse = "; ")) %>% 
  mutate(genotyper = reformat_hla_genotyper(genotyper)) %>% 
  arrange(genotyper, sample) %>% 
  mutate_all(na_if, "NULL") %>% 
  write_csv(here("data/pmid30518681/compiled_nih_bulk_genotypes.csv"))

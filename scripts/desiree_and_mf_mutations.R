#!/share/comailab/kramundson/miniconda3/bin/Rscript
# mf_mutations.R
# Kirk Amundson
# 9 March 2021

# USAGE: Rscript mf_mutations.R <vcf>
# where <vcf> is a gzip raw or filtered VCF (4.2 format)

# Packages
library(tidyverse)
library(stringr)

# Functions

## Parse VCF, extracting column names by call to shell
read_vcf <- function(file, gzip = T, multi = F) {
    
  if (multi) {
    file_to_search = file[1]
  }
  else { 
    (file_to_search = file)
  }
  
  if (gzip) {
    vcf_header <- system(paste0("gunzip -c ", file_to_search, " | head -n 1000 | grep '^#C'"), intern = T) %>% 
      str_remove("^#") %>% 
      str_remove_all("[0-9]x_") %>% 
      str_split(pattern = "\t") %>% 
      pluck(1)
  } 
  else {
    vcf_header <- system(paste0("grep '^#C'", file), intern = T) %>% 
      str_remove("^#") %>% 
      str_remove_all("[0-9]x_") %>% 
      str_split(pattern = "\t") %>% 
      pluck(1)
  }
  
  vcf <- map_dfr(file, function(x) read_tsv(x, col_names = vcf_header, comment = "#", na = c(".", "./.", "././.", "./././.")))  
  
  return(vcf)
}

## Parse VCF-like flat TSV I typically make when filtering raw calls
read_flat_tsv <- function(file) {
  vcf <- file %>% 
    map_dfr(function(x) read_tsv(x, col_names = T, na = c(".", "./.", "././.", "./././.", "NA")))
}

## Separate sample-specific fields of a VCF into their own respective columns
sep <- function(...) {
  dots <- list(...)
  separate_(..., into = paste(dots[[2]], attributes, sep = "_"), convert = T, sep = ":")
}

# Parse command line arguments
args <- commandArgs(trailingOnly = T)
file <- args[1]

# Read VCF
vcf <- read_vcf(file)
attributes <- names(table(vcf$FORMAT)) %>% 
  str_split(":") %>% 
  pluck(1)

# Parse column names to work on
samples_to_open <- colnames(vcf)[-c(1:9)]

# Separate sample-specific attributes
opened <- vcf %>% 
  Reduce(f = sep, x = samples_to_open)

# Get MF mutations
mf_mutations <- opened %>% 
  filter(if_any(matches("MF.+GT"), ~ . != USDA_Desiree_GT))

# Write mf mutations to file
out <- paste0("mf-mutations-", file) %>% str_replace("vcf.gz", "tsv.gz") %>% gzfile()
write_tsv(mf_mutations, out, col_names = T)

# Get Desiree mutations, USDA Desiree included
des_mutations <- opened %>%
  select(-matches("KPT")) %>% 
  filter(if_any(matches("Cornell_Desiree_GT|CNRS_Desiree_GT|IPK_Desiree_GT|USDA_Desiree_GT"), ~ . != CIP_Desiree_GT))

# Write Desiree mutations to file
out2 <- paste0("des-mutations-", file) %>% str_replace("vcf.gz", "tsv.gz") %>% gzfile()
write_tsv(des_mutations, out2, col_names = T)

# Get Desiree mutations, USDA Desiree excluded
true_des_mutations <- opened %>% 
  select(-matches("USDA_Desiree")) %>% 
  filter(if_any(matches("Cornell_Desiree_GT|CNRS_Desiree_GT|IPK_Desiree_GT"), ~ . != CIP_Desiree_GT))

# Write true Desiree mutations to file
out3 <- paste0("des-mutations-usda-out-", file) %>% str_replace("vcf.gz", "tsv.gz") %>% gzfile()
write_tsv(true_des_mutations, out3, col_names = T)
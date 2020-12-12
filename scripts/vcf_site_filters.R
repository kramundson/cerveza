#!/share/comailab/kramundson/miniconda3/bin/Rscript
# 2020_1116
# Kirk Amundson
# vcf_site_filters.R

# USAGE: Rscript vcf_site_filters.R <input vcf>

# Applies site quality filters to per-interval VCF. For this potato run, I called variants by chromosome
# This script only takes in one chromosome at a time, to parallelize, run at the command line using GNU parallel as follows:
# NOTE: modidfy shell expansion as necessary
# parallel -j <number of jobs to run in parallel> Rscript 2019_0225_filter_cheat_MM.R ::: *.vcf.gz

library(tidyverse)
library(stringr)

# generic_qual_filter <- function(file, mom, dad) {
generic_qual_filter <- function(file) {
  # parse header, requires zgrep installed on system
  print("Parsing VCF header")
  vcf_header <- system(paste("zgrep '#C'", file), intern = T) %>%
    str_replace("^#C", "C") %>%
    str_replace_all("[0-9]x_", "") %>%
    str_split(pattern = "\t")
  
  # read in file
  print("Reading in VCF")
  vcf <- read_tsv(file, col_names = vcf_header[[1]], comment = "#", na = c("NA", ".", "./.", "././.", "./././."))
  
  # parse INFO column; this is the per-locus attributes
  # this assumes that the INFO column format of the first row is valid for all rows.
  # For the MM cheat dataset, I validated this at the command line:
  # zgrep -v "^#" all-calls.vcf.gz | cut -f 8 | sed -e 's/=[A-Za-z0-9.,-]\+//g' | uniq -c
  # returns 26135465 AB;ABP;AC;AF;AN;AO;CIGAR;DP;DPB;DPRA;EPP;EPPR;GTI;LEN;MEANALT;MQM;MQMR;NS;NUMALT;ODDS;PAIRED;PAIREDR;PAO;PQA;PQR;PRO;QA;QR;RO;RPL;RPP;RPPR;RPR;RUN;SAF;SAP;SAR;SRF;SRP;SRR;TYPE
  info <- str_split(vcf$INFO[1], ";")[[1]] %>% 
    str_replace("=.+", "")
  print(info)
  
  # parse FORMAT column, this is the per-sample attributes
  print("Parsing FORMAT")
  attributes <- str_split(names(table(vcf$FORMAT)), ":", simplify = F)
  print(attributes)[[1]]
  
  # apply filters, opening up per-sample columns as necessary
  print("Applying site quality filters")
  
  vcf_filt_1 <- vcf %>% 
    mutate(INFO = str_replace_all(INFO, "[A-Za-z]*=", "")) %>%
    separate(INFO, into = info, sep = ";", convert = T) %>% 
    filter(FORMAT != "GT:GQ:DP:AD:RO:QR:AO:QA:GL")
    
  vcf_filt_2 <- vcf_filt_1 %>% 
    filter(NUMALT == 1)
    
  vcf_filt_3 <- vcf_filt_2 %>% 
    filter(CIGAR == "1X")
    
  vcf_filt_4 <- vcf_filt_3 %>% 
    filter(QUAL >= 20)
  
  vcf_filt_5 <- vcf_filt_4 %>% 
    mutate(MQM = as.numeric(MQM)) %>% 
    filter(MQM >= 50)
  
  vcf_filt_6 <- vcf_filt_5 %>% 
    filter(MQMR >= 50)

  vcf_filt_7 <- vcf_filt_6 %>% 
    mutate(MQ.diff = abs(MQM - MQMR)) %>% 
    filter(MQ.diff < 10)
  
  vcf_filt_8 <- vcf_filt_7 %>% 
    mutate(RPPR = as.numeric(RPPR)) %>% 
    filter(RPPR <= 20)
  
  vcf_filt_9 <- vcf_filt_8 %>% 
    mutate(RPP = as.numeric(RPP)) %>% 
    filter(RPP <= 20)
  
  vcf_filt_10 <- vcf_filt_9 %>% 
    mutate(EPP = as.numeric(EPP)) %>% 
    filter(EPP <= 20)
  
  vcf_filt_11 <- vcf_filt_10 %>% 
    mutate(EPPR = as.numeric(EPPR)) %>% 
    filter(EPPR <= 20)

  vcf_filt_12 <- vcf_filt_11 %>% 
    mutate(SAP = as.numeric(SAP)) %>% 
    filter(SAP <= 20) 
  
  vcf_filt_13 <- vcf_filt_12 %>% 
    mutate(SRP = as.numeric(SRP)) %>% 
    filter(SRP <= 20)

  # write out filtered chromosome variant as tsv, then when all are done, read back in, open up sample-specific attributes, and filter on those.
  # want to generate a list of parental SNPs for each family comparison
  writeout <- gsub("calls", "filtered-calls", file)
  print(paste("Writing out filtered calls to:", writeout))
  write.table(arrange(vcf_filt_13, POS), gzfile(writeout),
              quote=F, sep = '\t', eol = '\n',
              na = "NA", row.names = F, col.names = T)
  
  chrom <- unique(vcf_filt_13$CHROM)
  
  summary_df <- data.frame("chrom" = chrom,
                           "step" = c("raw",
                                      "FORMAT != GT:GQ:DP:AD:RO:QR:AO:QA:GL",
                                      "NUMALT == 1",
                                      "CIGAR == 1X",
                                      "QUAL >= 20",
                                      "MQM >= 50",
                                      "MQMR >= 50",
                                      "MQ.diff < 10",
                                      "RPPR <= 20",
                                      "RPP <= 20",
                                      "EPP <= 20",
                                      "EPPR <= 20",
                                      "SAP <= 20",
                                      "SRP <= 20"),
                           "count" = c(nrow(vcf),
                                       nrow(vcf_filt_1),
                                       nrow(vcf_filt_2),
                                       nrow(vcf_filt_3),
                                       nrow(vcf_filt_4),
                                       nrow(vcf_filt_5),
                                       nrow(vcf_filt_6),
                                       nrow(vcf_filt_7),
                                       nrow(vcf_filt_8),
                                       nrow(vcf_filt_9),
                                       nrow(vcf_filt_10),
                                       nrow(vcf_filt_11),
                                       nrow(vcf_filt_12),
                                       nrow(vcf_filt_13)))
  
  writeout_2 <- gsub("calls", "filtered-calls-summary", file)
  write.table(summary_df, writeout_2, quote = T, sep = "\t",
              row.names = F, col.names = T, eol = "\n",
              na = "NA")
                           
}

args = commandArgs(trailingOnly = TRUE)
file <- args[1]
generic_qual_filter(file)

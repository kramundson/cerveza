#!/bin/bash

wget http://solanaceae.plantbiology.msu.edu/data/S_tuberosum_Group_Phureja_chloroplast_DM1-3-516-R44.fasta.zip
wget http://solanaceae.plantbiology.msu.edu/data/S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta.zip
wget http://solanaceae.plantbiology.msu.edu/data/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz

gunzip DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz
unzip S_tuberosum_Group_Phureja_chloroplast_DM1-3-516-R44.fasta.zip
unzip S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta.zip

sed -i.bak -e 's/ .\+//g' -e 's/-/_/g' -e '$a\' S_tuberosum_Group_Phureja_chloroplast_DM1-3-516-R44.fasta
sed -i.bak -e 's/ .\+//g' -e 's/-/_/g' -e 's/\?/qmark_/g' -e '$a\' S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta

cat DM_1-3_516_R44_potato_genome_assembly.v6.1.fa \
    S_tuberosum_Group_Phureja_chloroplast_DM1-3-516-R44.fasta \
    S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta \
    > DMv6_chloro_mito.fasta

samtools faidx DMv6_chloro_mito.fasta
picard CreateSequenceDictionary R=DMv6_chloro_mito.fasta
cut -f 1-2 DMv6_chloro_mito.fasta.fai > DMv6_chloro_mito.genome
bwa index DMv6_chloro_mito.fasta
bedtools makewindows -g DMv6_chloro_mito.genome -w 100000 > DMv6_chloro_mito_100k_windows.bed
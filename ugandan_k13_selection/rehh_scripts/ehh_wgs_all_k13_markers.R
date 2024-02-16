################################################################################
title: "EHH and IHH Analysis on Plasmodium falciparum Genomic Data"
################################################################################

knitr::opts_chunk$set(echo = TRUE)
library(vcfR)
library(rehh)
library(ape)
library(ggplot2)
library(dplyr)

## Introduction

# This document presents an analysis of Extended Haplotype Homozygosity (EHH) and Integrated Haplotype Homozygosity (IHH) in Plasmodium falciparum genomic data, focusing on the kelch-13 622I mutation and its potential co-selected alleles.

## Data Loading

# hap_file: This is the path to your VCF file. Make sure the path is correctly specified.
# 
# min_perc_geno.mrk: This parameter sets the minimum percentage of genotype information required for a marker (SNP) to be included in the analysis. Setting it to 70 means that any SNP with less than 70% genotype information across all 
samples will be excluded. This is a way to filter out SNPs with too much missing data.
# 
# vcf_reader: Specifies the method used to read the VCF file. Using "data.table" is a good choice for larger datasets as it is efficient and fast.
# 
# verbose: Setting this to TRUE will provide additional information during the function execution, which can be helpful for debugging or understanding the process.
# 
# polarize_vcf: If set to FALSE, the VCF file will not be polarized. Polarization is about aligning alleles with a reference. If your analysis doesn’t require this or if it’s already been done, you can set it to FALSE.
# 
# remove_multiple_markers: If this parameter is included and set to TRUE, it would remove positions with more than one marker (e.g., multiallelic sites). This can be important because rehh is designed to work with biallelic markers. If 
your VCF contains multiallelic sites, you might want to set this to TRUE.
# 
# Other parameters you might consider include:
# 
# chr.name: If your analysis is focused on a specific chromosome, you can specify it here. This is particularly useful if your VCF contains multiple chromosomes but you're only interested in one.
# 
# region: If you are only interested in a specific region of the chromosome, you can specify the start and end positions with this parameter.
# 
# max_perc_missing.ind: Similar to min_perc_geno.mrk, but for individuals. It specifies the maximum percentage of missing genotype data allowed for an individual to be included in the analysis.
# 
# 
# Extended params summary:
# 
# hap_file: The path to the VCF file.
# 
# map_file (optional): Path to a map file, if available, which contains information about the physical or genetic positions of the markers.
# 
# min_perc_geno.hap (optional): Minimum percentage of genotyping information required for a haplotype to be included in the analysis.
# 
# min_perc_geno.mrk: Minimum percentage of genotype information required for a marker (SNP) to be included in the analysis. Setting it to 100% means all markers must have complete genotyping across all samples.
# 
# min_maf (optional): Minimum minor allele frequency required for a marker to be included. This is useful for filtering out rare variants.
# 
# chr.name (optional): Specifies the chromosome to be analyzed if the VCF file contains multiple chromosomes.
# 
# popsel (optional): Used to select specific populations for the analysis if your dataset includes multiple populations.
# 
# recode.allele: Indicates whether alleles should be recoded. The default value of FALSE will retain the original allele coding.
# 
# allele_coding: Specifies the coding of the alleles. The default "12" coding means alleles will be coded as 1 and 2.
# 
# haplotype.in.columns (optional): Determines the format of the output data.
# 
# remove_multiple_markers: If TRUE, removes positions with more than one marker (e.g., multiallelic sites). Useful for ensuring that only biallelic markers are analyzed.
# 
# polarize_vcf: If TRUE, the alleles in the VCF file will be polarized with respect to a reference.
# 
# capitalize_AA: Capitalizes the alleles, which can be important for consistency in allele representation.
# 
# vcf_reader: The method used to read the VCF file. "data.table" is efficient for large datasets.
# 
# position_scaling_factor (optional): A factor to scale the positions of the markers, useful in certain genetic analyses.
# 
# verbose: If TRUE, provides additional information during the execution of the function, useful for monitoring progress and debugging.

knitr::opts_chunk$set(echo = TRUE)
library(vcfR)
library(rehh)
library(ape)
library(ggplot2)
library(dplyr)

# Load VCF data
vcf_path <- "/Users/george/Bailey_lab/ugandan_k13_selection_project/sWGS_files/vcfs/chr13_norm.COI2orbelow.update.norm.SNPs.merged.ploidyfix.missgeno5.maf1.vcf.gz" # Replace with your VCF file path

# question - are these phased - if not, do they need to be
# question - are these polaized to reference, if not do they need to be - try polarize_vcf=T
chr_h<-rehh::data2haplohh(hap_file=vcf_path,
                          min_perc_geno.mrk =90,
                          vcf_reader = "data.table", 
                          verbose = T , 
                          polarize_vcf = F,
                          remove_multiple_markers=T)

#pull position specific marker (arbitrary at this time)
chr.res <- scan_hh(chr_h,discard_integration_at_border = F)

#chr.res<-rownames_to_column(chr.res,"name")

chr.res <- chr.res %>%
  mutate(name=row.names(chr.res))

plotting_dir <- "/Users/george/Bailey_lab/ugandan_k13_selection_project/inital_exporatory_analysis/ehh_plots/"
site_list <- c("1724974", "1725133", "1725259", "1725277", "1725295", "1725316", "1725340", "1725370", "1725382", "1725385", "1725388", "1725389", "1725418", "1725454", "1725521", "1725556", "1725570", "1725592", "1725626", "1725652", 
"1725662", "1725676")
allele_list <- c("675V", "622I", "580Y", "574L", "568G", "561H", "553L", "543T", "539T", "538V", "537I", "537D", "527H", "515K", "493H", "481V", "476I", "469YorF", "458Y", "449A", "446I", "441L")

calc_plot_ehh_various_alleles <- function(site,allele){
  
  if (length(chr.res$name[chr.res$POSITION == site]) != 0){
    
    marker<-chr.res$name[chr.res$POSITION == site]
    
    #calculate ehh in order to pull positions ONLY 
    WGS_x<-calc_ehh(chr_h, 
                    mrk = marker, 
                    include_nhaplo = T,
                    include_zero_values =T,
                    polarized = F,
                    phased = T, 
                    discard_integration_at_border =T)
    
    png(paste0(plotting_dir,allele,"_ehh_WGS.png"))
    
    plot((WGS_x), xlim = c(1650000, 1800000), 
         main = paste0("Loci: ",allele))
    
    dev.off()
    
  }
}

for (each_site_number in 1:length(site_list)){
  
  site <- site_list[each_site_number]
  allele <- allele_list[each_site_number]
  print(paste("plotting", allele))
  
  calc_plot_ehh_various_alleles(site,allele)
}


# Microsatellite_Comparison_Scripts


This repository contains a set of tools for analyzing and comparing microsatellite regions from variant call format (VCF) files and PCR-derived allele data. These scripts assist in comparing haplotypes across tools, verifying genotype concordance, and standardizing PCR allele calls.


## Overview

The repository includes the following Python scripts:

1. automated_hapcaller_comparisons.py – Parses VCF files to analyze haplotypes in specified regions, counts indels, and normalizes allele lengths.

2. diff_refs_genotype_concordance_hap_downsampling.py – Compares genotypes between two VCFs (truth vs. called) across specified regions, accommodating reference mismatches.

3. pcr_data_excel_parsing.py – Parses PCR-based allele call data from a CSV file, normalizes allele lengths per region, and prints standardized genotypes.




## 1. automated_hapcaller_comparisons.py

### Description: 
Analyzes a VCF file to extract and compare indel-based haplotypes within user-defined genomic regions. It identifies complex regions with multiple indels and normalizes allele length differences per sample.

### Usage: 

python automated_hapcaller_comparisons.py -v <VCF_FILE> -r "\<REGIONS>\"

### Arguments:

-v, --vcf: Path to the input VCF file

-r, --regions: Semicolon-separated list of regions. Format: 'chr1:20000-20100; chr2:30000-31000'




## 2. diff_refs_genotype_concordance_hap_downsampling.py

### Description:
Performs a detailed comparison between a "truth" and a "called" VCF file, focusing on genotype concordance within specified regions. Handles differences in reference sequences by comparing allele length shifts or exact sequences.

### Usage:
python diff_refs_genotype_concordance_hap_downsampling.py --truth <TRUTH_VCF> --call <CALL_VCF> --regions "\<REGIONS>\"

### Arguments:

--truth: Path to the ground-truth/reference VCF file

--call: Path to the comparison/called VCF file

--regions: Semicolon-separated list of regions. Format: 'chr1:20000-20100; chr2:30000-31000'




## 3. pcr_data_excel_parsing.py

### Description:
Processes a CSV file containing PCR-derived allele lengths. For each microsatellite region, it normalizes allele values and converts them into standardized genotypes.

### Requirements:
Input CSV file where:
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;First column = sample name
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Every pair of subsequent columns = allele calls for one region

### Usage:
Update the following line in the script:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;filename = 'fill in'

with the actual CSV file path, then run:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;python pcr_data_excel_parsing.py




## Requirements:

Python 3.6+

vcfpy


## Install dependencies:

pip install vcfpy

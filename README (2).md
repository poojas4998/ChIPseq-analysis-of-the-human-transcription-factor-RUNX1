# Project Description

This project aims to investigate the role of the RUNX1 protein in breast cancer cells. The RUNX protein family plays a crucial role in regulating various biological processes to determine proper cell fate. The RUNX1 protein acts as an oncogene and tumor suppressor in breast cancer cells by participating in pathways that change chromatin structure. The authors of this study used ChIP-seq and genome-wide chromatin conformation capture techniques to investigate the effects of suppressing the expression of RUNX1 in MCF-7 breast cancer cells. They found that RUNX1 is involved in long-range chromatin interactions and regulates a set of genes that are involved in cellular processes such as cell cycle progression and DNA damage response.

The project aims to replicate the results of the paper using the replicate ChIP-seq data for analysis. The samples used in this study come from RUNX1 depleted and MCF-7 control breast cancer cells. The workflow for this project includes data acquisition, reproducible peak finding and filtering, peak annotation and motif finding, and generating bigWig files for correlation analysis and visualization.

# Contributors

Group 5

Manasa Rapuru - Data Curator (manasarapuru)

Pragya Rawat - Programmer (rpragya17)

Pooja Savla - Analyst (poojas4998)

Vrinda Jethalia - Biologist (vrindajethalia799)

# Repository Contents

data_curator_simple.snake - snakemake workflow to conduct data processing and check data quality 

Programmer_script.Rmd - qsub used to find reproducible peaks from provided data. Includes script used to plot pie chart of peaks annotated to nearest genomic feature. 

index_make.qsub - This code is used to generate bigWig files from BAM files using the bamCoverage tool. The --bam option specifies the input BAM file, and the --outFileName option specifies the name of the output bigWig file. Each command generates a separate bigWig file for each of the input BAM files. The generated bigWig files can be used for visualization and downstream analysis.

bigwigsummary.qsub - The first code calculates a matrix of read densities over specified genomic regions for two IP samples. It then plots the profiles for both samples.
The second code calculates a summary of read coverage for multiple bigWig files and generates a clustered heatmap of the Pearson correlation values between them.

compute_matrix.qsub - This is a bash script that loads required modules and uses deeptools to generate scaled matrices and plots for two different IP samples. The script computes scaled matrices using computeMatrix function with input bigWig files and gene bed file, and generates plots for each sample using plotProfile function with the corresponding matrix file as input. The script also prints some information about the job and its execution.

biologist_script.R - This is an R script to obtain differentially expressed genes, a bar chart that provides information on the number of DEGs (with +/- 5kb of TSS) interacting with RUNX1, and a heatmap of the HI-C data.

BF 528 Project 3 Report(1).pdf - Final written report containing all figures and results as well as discussion of these results 

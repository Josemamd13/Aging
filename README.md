# Metabolism and Ageing Study

This repository contains code and documentation related to a research project focusing on understanding the changes in metabolism at the tissue level with ageing and the differences between genders.

## Introduction

The ageing of populations due to increased life expectancy is a significant biomedical challenge in the 21st century. With longer life spans come challenges like increased multimorbidities, altered metabolism, mitochondrial dysfunction, and higher healthcare costs. Understanding the role of metabolism in ageing is crucial but remains an area with much to discover.

## Objectives

The main objectives of this study are:
- To analyze the changes in metabolism at the tissue level during ageing in both men and women.
- To build tissue-specific metabolic models using GTEx consortium gene expression data.
- To explore age-associated functional changes in metabolic tasks within these tissue-specific models.

## Methodology

### Data Source
We will utilize the GTEx consortium data, providing gene expression data from 46 human tissues.

### Approach
- Building tissue-specific metabolic models separately for men and women across different age ranges (young, adult, old).
- Integrating gene expression datasets into the Human Metabolic model to create context-specific metabolic models.
- Exploring age-associated functional changes using standardized metabolic tasks.


# R Script Overview
## 00_ExploratoryDataAnalysis.R

### The line if(args[1]=="PrepareTables"){ checks whether the first argument provided when running the script is equal to "PrepareTables". If this condition is met, the subsequent block of code enclosed within {} will be executed. If the condition isn't satisfied, the code inside the curly braces won't be executed, and the script will move on to the following instructions...

#### Description
The provided code appears to be an R script that performs several data preparation and processing operations for analyzing differential gene expression and exploring gene expression and phenotypic data.

#### Data Preparation:
- Reads sample annotation and phenotype files from the GTEx project.
- Filters and combines sample annotation and phenotype data.
- Creates and saves processed metadata tables in CSV format.
- Performs manipulations to extract and aggregate specific sample and phenotype information.

#### Visualizations:
- Generates bar charts for age distribution in different states (young, adult, old) based on gender and overall.
- Produces visualizations of observation counts per age state for males and females.
  
#### Processing Gene Expression Data:
- Reads gene read data.
- Filters and organizes sample metadata to match gene expression data.
- Divides and saves gene expression data by tissue type into separate CSV files.
- Creates summary tables of observation counts per tissue type, age state, and gender, saved as separate CSV files for further analysis.

### The line if(args[1]=="CreateRankedListOfGenes"){ checks whether the first argument provided when running the script is equal to "CreateRankedListOfGenes". If this condition is met, the subsequent block of code enclosed within {} will be executed. If the condition isn't satisfied, the code inside the curly braces won't be executed, and the script will move on to the following instructions...

#### Description
This script generates ranked lists of genes derived from differential expression results. It takes as input files with two columns: the first column contains ENSG IDs, and the second column contains Fold Change (FC) values. The script creates two types of ranked lists: one based on ENSG IDs and the other using Gene Symbols.

The script is divided into two main sections:

1. ENSG Rank

This section processes the differential expression results based on ENSG IDs:

- Reads the files, sorts data by FC values, and selects ENSG IDs and FC columns.
- Generates a new .rnk file containing sorted ENSG IDs and FC values in the specified ENSG directory.
  
2. GeneSymbol Rank

This section utilizes the biomaRt library to obtain Gene Symbols corresponding to ENSG IDs:

- Reads the files, sorts data by FC values, and extracts ENSG IDs.
- Retrieves corresponding Gene Symbols for ENSG IDs using the biomaRt package.
- Merges Gene Symbols with the original data based on ENSG IDs, selects Gene Symbols and FC columns, and generates a new .rnk file in the specified GeneSymbol directory.

### The line if(args[1]=="GSEA"){ checks whether the first argument provided when running the script is equal to "GSEA". If this condition is met, the subsequent block of code enclosed within {} will be executed. If the condition isn't satisfied, the code inside the curly braces won't be executed, and the script will move on to the following instructions...

#### Description
This script facilitates Gene Set Enrichment Analysis (GSEA) on pre-ranked gene lists using the GeneSymbol ranked files obtained from differential expression analysis. It executes GSEA with specified parameters for multiple tissues (The script is designed to be run with the command line and requires GSEA (version 4.3.2) and Java to be installed. Ensure that the GSEA command-line interface (gsea-cli.sh) is accessible.).


## 01_DifferentialExpresiion.R
### Description
The following R script (`01_....R`) contains a function named `diffexpress` aimed at performing differential expression analyses based on various arguments such as tissue type, sex, and comparison between groups. It involves various steps as outlined below:

### Libraries: Loads necessary R libraries
library(dplyr)
library(data.table)
library(edgeR)
library(limma)

### Directory Creation: Checks and creates directories for storing results if they don’t exist.

### Function Definition: Defines the diffexpress function to perform differential expression analyses based on tissue type, sex, and comparison parameters.

### Data Loading: Loads metadata and count data from respective files.

### Data Filtering: Filters samples based on specified conditions:

- Minimum sample counts
- Filtering by sex and tissue type
  
### Differential Expression Analysis

- Normalizes the data and creates design matrices for differential expression.
- Conducts differential expression analysis using linear models (lmFit) and statistical tests (eBayes).
- Writes the results of the analysis into output files in the "Results/DifferentialExpression/" directory.

In summary, the script seems to automate the process of conducting differential expression analysis for gene expression data based on various conditions, such as tissue type, sex, and age comparison.

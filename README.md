# IFN-MS interim project

*Guillermo Reales*
25/01/2023


This folder contains all the data and code for the IFN pQTL MS project.
The project aims to generate PGS on IFNAR2 pQTL GWAS files and use individual MS data to correlate PGS for given protein levels and disease status.

This README will be updated as I go.

## Data

The data comes from Lundtoft et al., 2022, who performed GWAS on a set of immune-related proteins in 6 PBMCs.

We started by IFNAR2, for which MS-protective tagged SNPs were associated with increased levels in B cells (CD20) and T cells (CD4 and CD8).

## Methods

### Data processing

We copied the files to `data/raw-summary`, updated their columns and computed SE from SE (see `code/02-process-summary-data-1.R`).
The datasets come in hg19, which we'll use for convenience in later steps.
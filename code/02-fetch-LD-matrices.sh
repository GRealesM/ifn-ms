#!/bin/bash
## Fetching UKBB LD matrices
## Guillermo Reales
## 2023-02-01

# We'll now fetch the UKBB LD matrices from Prive et al. 
# These are LD matrices filtered by the HapMap3 variants and we'll use them to generate our RápidoPGS models.

mkdir ../data/LD

# We used the code below to fetch the matrices, but since we already have them from a different project, we'll simply copy them
#mkdir ../data/LD
#cd    ../data/LD
#
#wget https://ndownloader.figshare.com/articles/13034123/versions/3 -O ukbbmat.zip
#unzip ukbbmat.zip
#rm ukbbmat.zip

# However, since we already have them from a different project, we'll simply copy them
cp /home/gr440/rds/rds-cew54-basis/05-PGS/v3/references/UKBB_mat/*.rds ../data/LD
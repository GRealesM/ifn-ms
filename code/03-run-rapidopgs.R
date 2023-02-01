## Running RápidoPGS
## Guillermo Reales
## 2023-02-01

# Now we have everything, we'll run RápidoPGS

setwd("/home/gr440/rds/rds-cew54-basis/Projects/ifn-ms/code")

# Load packages
library(RapidoPGS) 

# Create ds list, path to LD matrices and out dir
indir  <- "../data/proc-summary/"
dslist <- dir(indir, "-fp.gz")


outdir <- "../data/pgs-models/"
if(!dir.exists(outdir)) dir.create(outdir)


LDmatrices <- "../data/LD"
N <- "N" # Sample size will be in the N column
ncores <- 1 # In principle, no need for more than one

# Here we're using an array job to run the models, so we'll create an array in bash containing
# the indexes to select the traits.
# And we'll use those to generate the models
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args) # Array index will be a 1-6 integer.

args <- dslist[args] # This way we transform the index number into the real trait name
dsn  <- gsub("-fp.gz", "", args)

ds <- fread(paste0(indir, args))

# Remove MHC region (chr6:20M-40M) from ds before running RápidoPGS
ds <- ds[ !(CHR == 6 & BP > 20000000 & BP < 40000000)]

pgsmodel  <- rapidopgs_multi(ds, trait="quant", LDmatrices=LDmatrices,  N=N, ncores=ncores, alpha.block = 1e-3, alpha.snp = 0.1)
pgsmodel <- merge(pgsmodel[, .(CHR, BP, ppi_susie, weight)], ds, by=c("CHR", "BP"))

fwrite(pgsmodel, paste0(outdir, dsn, "-pgs.gz"), sep="\t")

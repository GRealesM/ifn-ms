# Processing Summary data - Step 1
 
# Here we'll 
# 1 - Import the raw GWAS summary stats 
# 2 - Get the other allele from .bim file
# 3 - Rename columns for pipeline
# 4 - Compute SE from BETA and STAT (aka Z-score)

setwd("/home/gr440/rds/rds-cew54-basis/Projects/ifn-ms/code")

# Load packages
library(data.table)
library(IMDtools)
setDTthreads(10)

# Get file names
fp  <- '../data/raw-summary/'
op  <- '../data/proc-summary/'
fn  <- dir(fp, pattern = ".gz")

# Import .bim file
bim  <- fread(paste0(fp, "BCU_filtered.bim"), header = FALSE)

# Import the frequency file to incorporate the allele frequencies, computed from 1000GP EUR individuals (503)
# We'll do a little processing using a sample file to align the ALT_FREQ in our files
frf  <- fread("../data/help-files/EUR_Freqs_hg19snps.txt")
frf[, c("BETA", "SE", "P"):=0] # Mock BETA/SE columns to fool g.align
setnames(frf, c("CHR19", "BP19"), c("CHR38", "BP38")) # to fool g.align, not really hg38
man  <- bim[, .(V1, V4, V6, V5)]
names(man) <- c("CHR38", "BP38", "REF", "ALT")
alf <- g.align(ds = frf, manifest = man) # Align frequencies from Europeans to the alleles in the files
alf  <- alf[, .(CHR38, BP38, REF,ALT, ALT_FREQ)]
names(alf)[1:2]  <- c("CHR", "BP")

# Process files
sapply(fn, function(i){
    message("Workin' on ", i,"...")
    f1 <- fread(paste0(fp, i))
    f1 <- merge(f1, bim[, .(V2, V5, V6)], by.x=c("SNP", "A1"), by.y = c("V2", "V5"))
    f1 <- merge(f1, alf, by.x=c("CHR", "BP", "V6", "A1"), by.y=c("CHR", "BP", "REF", "ALT"))
    f1[, SE:= BETA/STAT] # Assuming STAT is the Z-score...
    f1 <- f1[!is.na(BETA)] # Remove missing SNPs
    setnames(f1, c("CHR", "BP", "SNP", "A1", "V6"), c("CHR19", "BP19","SNPID", "ALT", "REF"))
    f1 <- f1[, .(CHR19, BP19, SNPID, REF, ALT, ALT_FREQ, BETA, SE, P)]
    setorder(f1, CHR19, BP19)
    nfn  <- gsub(".gz", "-plr.gz", i)
    fwrite(f1, paste0(op, nfn), sep="\t")
    message("File ", nfn, " ready!")
})


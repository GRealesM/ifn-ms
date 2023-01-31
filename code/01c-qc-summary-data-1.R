# QC Summary data - Version 1
 
# Here we'll 
# 1 - Import the processed datasets from the previous step
# 2 - Compute sdY and neff
# 3 - Filter variants outside established criteria. Generate QC figures from it.
# 4 - Filter variants by the HapMap3 variants to match the variants in the 
# UKBB LD matrices we'll use for RápidoPGS.

setwd("/home/gr440/rds/rds-cew54-basis/Projects/ifn-ms/code")

# Load packages
library(data.table)
library(IMDtools)
library(ggplot2)
setDTthreads(10)

# Define helper functions

sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}


# Get file names
fp  <- '../data/proc-summary/'
fn  <- dir(fp, pattern = ".gz")

# Load HapMap3 variants

# We obtained the original file from LDpred2 paper, and made a local copy
#hm3 <- readRDS(url("https://ndownloader.figshare.com/files/25503788"))
#saveRDS(hm3, "../data/help-files/hm3_variants.RDS")

hm3 <- readRDS("../data/help-files/hm3_variants.RDS")
hm3 <- as.data.table(hm3[,1:4])
names(hm3)  <- c("CHR38", "BP38", "REF", "ALT") # Not really hg38, but we'll need those names for g.align down the road



# Get dataset, process dataset

fn1  <- fn[2]
nfn  <- gsub("-plr.gz", "", fn1)

f1 <- fread(paste0(fp,fn1))



### Let's compute sdY and n_eff
f1[,sdy:=sdY.est(SE^2, ALT_FREQ, N)]
sd_ss <- with(f1, sdy / sqrt(N * SE^2))
sd_val <- with(f1, sqrt(2*ALT_FREQ*(1-ALT_FREQ)))

is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.2) | sd_ss < 0.1 | sd_val < 0.05 
table(is_bad)

dt <- data.table(ss = sd_ss, val = sd_val, is_bad = is_bad)

tp  <- ggplot(dt, aes(sd_val, sd_ss, color = is_bad, alpha = I(0.5))) +
  geom_point()+
  theme_minimal() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")+
  ggtitle(paste0(nfn, " QC"))
  
ggsave(tp, paste0("../figures/sd_", nfn, "_qc.png"), height = 7, width = 7, units = "in", bg = "white")
f1  <- f1[!is_bad]

# Align with HapMap3 variants
alr <- 




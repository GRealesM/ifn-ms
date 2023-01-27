## Here we'll get the summary data from the files we're interested in and process them for use with RápidoPGS

mkdir ../data/raw-summary
mkdir ../data/proc-summary
mkdir ../data/help-files

cp ~/rds/rds-cew54-wallace-share/Data/GWAS-summary/ifn-pqtls/*IFNAR2* ./data/raw-summary

cd ../data/raw-summary

# Simplify names
rename "P1F_transformed.I1_GM_" "" *.gz
rename "P1F_transformed.I2_GM_" "" *.gz
rename "of_" "" *.gz
rename "assoc.linear." "" *.gz

# Get .bim file to allocate SNPs, too
cp ~/rds/rds-cew54-wallace-share/Data/GWAS-summary/ifn-pqtls/BCU_filtered.bim ./

# Retrieve SNP allele frequencies from EUR individuals from 1000GP
grep "EUR" ~/rds/rds-cew54-basis/95-1000genomes/reference_hg19/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1, $1}' > ../data/help-files/EUR_1KG_samples.txt


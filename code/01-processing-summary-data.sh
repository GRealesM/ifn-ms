## Here we'll get the summary data from the files we're interested in and process them for use with RápidoPGS

mkdir ../data/raw-summary
mkdir ../data/proc-summary
mkdir ../data/help-files

hdir="../data/help-files/"

cp ~/rds/rds-cew54-wallace-share/Data/GWAS-summary/ifn-pqtls/*IFNAR2* ../data/raw-summary

cd ../data/raw-summary

# Simplify names
rename "P1F_transformed.I1_GM_" "" *.gz
rename "P1F_transformed.I2_GM_" "" *.gz
rename "of_" "" *.gz
rename "assoc.linear." "" *.gz

cd ../../code

# Get .bim file to allocate SNPs, too
cp ~/rds/rds-cew54-wallace-share/Data/GWAS-summary/ifn-pqtls/BCU_filtered.bim ../data/raw-summary

# Retrieve SNP allele frequencies from EUR individuals from 1000GP
grep "EUR" ~/rds/rds-cew54-basis/95-1000genomes/reference_hg19/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1, $1}' > "$hdir"EUR_1KG_samples.txt

# Prepare a SNP list from a sample file, since they all have the same SNPs
zcat ../data/raw-summary/IFNAR2_CD4p.gz | awk 'NR>1{print $1, $3, $3, $3}' > "$hdir"hg19snps.txt

echo -e "CHR19\tBP19\tSNPID\tREF\tALT\tALT_FREQ\tOBS_CT" > "$hdir"EUR_Freqs_hg19snps.txt

# We extract frequencies and append to our Freqs file
for chr in {1..22};
do
	echo "Extracting from chr$chr..."
	grep -P "^"$chr" " "$hdir"hg19snps.txt > "$hdir"tmpsnps.txt
	plink2 --bfile ~/rds/rds-cew54-basis/95-1000genomes/reference_hg19/chr$chr --keep "$hdir"EUR_1KG_samples.txt --extract range "$hdir"tmpsnps.txt --freq cols=+pos --out "$hdir"temp
	tail -n+2 "$hdir"temp.afreq >> "$hdir"EUR_Freqs_hg19snps.txt
done

# Cleanup
rm "$hdir"temp* "$hdir"tmpsnps.txt

# Then we can run the R script to process the files.
Rscript 01b-process-summary-data-1.R
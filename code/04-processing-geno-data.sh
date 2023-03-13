#!/bin/bash
## Processing genotype data
## Guillermo Reales
## 2023-03-11

## Now we'll process the genotype data from the IMSGC study, and apply QC to them.
#
# We have the following datasets to process:
# - CCC2_Cases_UKC
# - CCC2_Cases_UKN
# - CCC2_Cases_UKP
# - CCC2_Cases_UKW
# - CCC2_Controls_Illu58C
# - CCC2_Controls_IluNBS
#
# They come in Oxford genotype format (ie. .gen.gz) so we'll first convert them to plink (bed) format.

# The files have already been copied to the directory below. They contain the data organised by chromosomes
dirs=(CCC2_Cases_UKC CCC2_Cases_UKN CCC2_Cases_UKP CCC2_Cases_UKW CCC2_Controls_Illu58C CCC2_Controls_IlluNBS)
fstarts=(MS_UKC MS_UKN MS_UKP MS_UKW 58C NBS)
dpath="../data/raw-geno/"

# The first issue is that their sample file comes with the wrong header, so we'll need to replace the second line
samplefiles=($(ls "$dpath"*/*.sample))

for i in ${samplefiles[@]}; do  sed -i '2s/.*/0 0 0 D D B/g' $i; done

# Then we'll convert all files to plink format


for i in seq {0..5};
do
    for f in seq {01..22};
        do
        plink2 --gen $dpath${dirs[$i]}"/"${fstarts[$i]}_"$f"_illumina.gen.gz 'ref-unknown' --sample $dpath${dirs[$i]}"/"${fstarts[$i]}_illumina.sample --oxford-single-chr chr$f --missing-code '-1' --threads 25 --make-bed --out $dpath${dirs[$i]}"/"${fstarts[$i]}_$f
        done
done
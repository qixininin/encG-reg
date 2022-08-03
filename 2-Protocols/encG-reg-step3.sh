#!/bin/bash
# Input User ID
user=$1
coflag=$2

# Users are asked to generate an m-by-k matrix whose elements are sample from N(0,1/m)
# Obtain the number of markers (m)
awk 'END{print NR}' ${coflag}.snpA1 > ${coflag}.m

# This Rscript automatically reads parameters stored in Golden.m and Golden.k
# And generate Golden.key, an m-by-k matrix
Rscript GenerateRandMat.R ${coflag}

# Combine key with SNPID and A1 alleles by columns
paste -d "\t" ${coflag}.snpA1 ${coflag}.key > ${coflag}.snpA1key

# Extract SNPs in your bfiles by plink
awk '{print $1}' ${coflag}.snpA1 > ${coflag}.snp
plink --bfile ${user} --extract ${coflag}.snp --make-bed --out ${user}.${coflag}

# Then merge 1KG-CHN with your bfiles into new plink files named Golden.merged
plink --bfile 1KG-EUR.${coflag} --bmerge ${user}.${coflag} --make-bed --out ${user}.${coflag}.merged

# Use plink2 "--score" to generate encrypted genotype matrix.
# variance-standardize: genotypes are scaled by SNP
k=`cat ${coflag}.k`
plink2 --bfile ${user}.${coflag}.merged --score ${coflag}.snpA1key 1 2 variance-standardize --score-col-nums 3-$(($k+2)) --out ${user}.${coflag}.merged

# Return Golden.${user}.sscore to server agent

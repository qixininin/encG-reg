# Step 1: Within-cohort quality control (cohort: QC)
# Set environmental variables (!!!!NEED TO BE MODIFIED!!!!)
# Please input *YOUR_COHORT_ID* with the Cohort ID we provide
# Please input *YOUR_PREFIX* with your plink bfile prefix

user=$1
bfile=$2

# Inclusion criteria
# (1) autosome SNPs only
# (2) SNPs with minor allele frequency (MAF) > 0.05 only
# (3) SNPs with missing rate < 0.2 only
# Suggested plink command:

plink --bfile ${bfile} --autosome --snps-only --maf 0.05 --geno 0.2 --no-pheno --make-bed --out ${user}
plink --bfile ${user} --freq --out ${user}
awk 'END{print NR}' ${user}.fam > ${user}.n
# Return ${user}.frq to server agent.

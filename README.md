---
editor_options: 
  markdown: 
    wrap: sentence
---

# encG-reg

## 1-Simulations

This github contains simulation codes of GRM, encGRM and encG-reg.
- Fig 1 Resolution for varying relatedness using GRM, encGRM and encG-reg.
- Fig 3 Sampling variance of GRM, encGRM and encG-reg in simulations.
- Fig S1 Heatmap presenting the role random matrix played in matrix multiplication.
- Fig S2 Validation for the sampling variance for GRM (assumption: binomial distribution).

## 2-Protocols

We offered a user-friendly protocol for encG-reg.
This protocol can be automated, such as by a web server that coordinates the study.
There are four steps in total, where steps 1 and 3 are performed by each collaborator and steps 2 and 4 are performed by a central analyst.
We provide suggested commands for your possible reference, and your environment should have plink1.9, plink2.0 and R installed.

### Step 1 Within-cohort quality controls

``` shell
# Set environmental variables (!!!!NEED TO BE MODIFIED!!!!)
# Please replace *YOUR_COHORT_ID* with the Cohort ID we provide
# Please replace *YOUR_PREFIX* with your plink bfile prefix
user=*YOUR_COHORT_ID* bfile=*YOUR_PREFIX*

# Inclusion criteria
# (1) autosome SNPs only
# (2) SNPs with minor allele frequency (MAF) > 0.05 only
# (3) SNPs with missing rate < 0.2 only Suggested plink command: 
plink --bfile ${bfile} --autosome --snps-only --maf 0.05 --geno 0.2 --no-pheno --make-bed --out ${user} plink --bfile ${user} --freq --out ${user}
```

Return \${user}.frq to server agent.

To be a very good citizen in this collaboration, the suggested name of \*.frq file is "YOUR_COHORT_ID.frq"

-   Step 2 Determine m and k
-   Step 3 Encrypt genotype matrix
-   Step 4 Perform encG-reg across cohorts

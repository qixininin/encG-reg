---
editor_options: 
  markdown: 
    wrap: sentence
---

# encG-reg

## 1-Simulations

This github contains simulation codes of GRM, encGRM and encG-reg.

-   Fig 1 Resolution for varying relatedness using GRM, encGRM and encG-reg.

-   Fig 3 Sampling variance of GRM, encGRM and encG-reg in simulations.

-   Fig S1 Heatmap presenting the role random matrix played in matrix multiplication.

-   Fig S2 Validation for the sampling variance for GRM (assumption: binomial distribution).

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
plink --bfile ${bfile} --autosome --snps-only --maf 0.05 --geno 0.2 --no-pheno --make-bed --out ${user}
plink --bfile ${user} --freq --out ${user}
```

Return \${user}.frq to server agent.

To be a very good citizen in this collaboration, the suggested name of \*.frq file is "YOUR_COHORT_ID.frq"

[plink1.9-frq](https://www.cog-genomics.org/plink/1.9/formats#frq)

| Column  | Meaning                       |
|---------|-------------------------------|
| CHR     | Chromosome code               |
| SNP     | Variant identifier            |
| A1      | Allele 1 (usually minor)      |
| A2      | Allele 2 (usually major)      |
| MAF     | Allele 1 frequency            |
| NCHROBS | Number of allele observations |
|         |                               |

| CHR | SNP         | A1  | A2  | MAF     | NCHROBS |
|-----|-------------|-----|-----|---------|---------|
| 1   | 1:13273     | C   | G   | 0.0601  | 416     |
| 1   | 1:14599     | A   | T   | 0.0601  | 416     |
| 1   | 1:14604     | G   | A   | 0.0601  | 416     |
| 1   | rs75454623  | G   | A   | 0.3966  | 416     |
| 1   | rs78601809  | T   | G   | 0.5     | 416     |
| 1   | 1:15777     | G   | A   | 0.08413 | 416     |
| 1   | rs200482301 | G   | T   | 0.4736  | 416     |
| 1   | 1:54716     | T   | C   | 0.1298  | 416     |
| 1   | rs3107975   | C   | T   | 0.137   | 416     |
|     |             |     |     |         |         |

### Step 2 Determine m and k

Upon the ".frq" files received, server agent identifies the shared SNPs across cohorts and choose the optimal SNP set, which will be used for randomization algorithm.
As the genotypes are generated in their respective platforms, to make life easier server agent excludes: palindromic bi-allelic loci, say A-T, G-C; strand-flipped loci, say A-G in one cohort but T-C in another.

#### Step 2.1 QC examination

To examine across-cohort quality control, we used CONVERGE data set as the reference control to reveal any possible mistake made in Step 1.
This examination includes MAF density plot between CONVERGE and every data set from the collaborators, and plot special shift between major and minor alleles when MAF approaches 0.5.
Examination reports are given in Appendix.

#### Step 2.2 Shared SNPs

We took the intersection of all SNP lists among all cohorts based on their SNP ID in ".frq" files.
In total, 1462 SNPs were in common among 12 cohorts of 1KG-CHN, UKB-CHN, CONVERGE, MESA, ALS, SYSU, BIG1, BIG2, Fudan, Yikon1, Yikon2 and Westlake.
Intersection information and MAF density plots are also given in Appendix.

#### Step 2.3 Genetic background across-cohort

We conduct principal component analysis (PCA) based on reported allele frequencies and use the population genetics Fst statistic to verify the genetic origin of each cohort.
PCA plot and pseudo-structure plot are given in Appendix.

#### Step 2.4 Determine m and k

According to Eq 1 and Eq 2, server agent then determines m and k upon the survived SNPs.
The number of shared SNPs are enough for identifying 1-degree relatedness, we would offer a list of 500 shared SNPs, whose m_e is 477 and the corresponding minimal number of k is 757.

Server agent will send GenerateRandMat.R, random seed, k, an SNP list, and 1KG-CHN binary format files to each collaborator.

### Step 3 Encrypt genotype matrix

### Step 4 Perform encG-reg across cohorts

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

We provide suggested commands for your possible reference, and your environment should have **plink1.9**, **plink2.0** and **R** installed.

### Step 1 Within-cohort quality controls

Set environmental variables **(!!!!NEED TO BE MODIFIED!!!!)**

Please replace **YOUR_COHORT_ID** with the Cohort ID.

Please replace **YOUR_PREFIX** with your plink bfile prefix.

``` shell
user=*YOUR_COHORT_ID*
bfile=*YOUR_PREFIX*
```

Inclusion criteria:

(1) Autosome SNPs only;

(2) SNPs with minor allele frequency (MAF) \> 0.05 only;

(3) SNPs with missing rate \< 0.2 only.

Suggested plink command:

``` shell
plink --bfile ${bfile} --autosome --snps-only --maf 0.05 --geno 0.2 --no-pheno --make-bed --out ${user}
plink --bfile ${user} --freq --out ${user}
```

We adopt one of the plink1.9 formats ".frq" [plink1.9 formats:frq](https://www.cog-genomics.org/plink/1.9/formats#frq) as a standard format for sharing summary data. ".frq" files include the following contents.

| Header  | Contents                      |
|---------|-------------------------------|
| CHR     | Chromosome code               |
| SNP     | Variant identifier            |
| A1      | Allele 1 (usually minor)      |
| A2      | Allele 2 (usually major)      |
| MAF     | Allele 1 frequency            |
| NCHROBS | Number of allele observations |

An example of **1KG-CHN.frq** is:

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

Return **\*.frq** to central site. To be a good citizen in this collaboration, the suggested name of **\*.frq** file will be **"YOUR_COHORT_ID.frq"**.

### Step 2 Determine m and k

Upon the **\*.frq** files received, central site identifies the shared SNPs across cohorts and choose the optimal SNP set, which will be used for randomization algorithm. As the genotypes are generated in their respective platforms, to make life easier central site excludes: palindromic bi-allelic loci, say A-T, G-C; strand-flipped loci, say A-G in one cohort but T-C in another.

#### Step 2.1 QC examination

To examine across-cohort quality control, we used CONVERGE data set as the reference control to reveal any possible mistake made in [Step 1](###Step-1-Within-cohort-quality-controls). This examination includes MAF density plot between CONVERGE and every data set from the collaborators, and plot special shift between major and minor alleles when MAF approaches 0.5.

#### Step 2.2 Shared SNPs

We took the intersection of all SNP lists among all cohorts based on their SNP ID in **\*.frq** files. In total, 1462 SNPs were in common among 12 cohorts of 1KG-CHN, UKB-CHN, CONVERGE, MESA, SBWCH, CAS1, CAS2, Fudan, Yikon1, Yikon2 and WBBC.

#### Step 2.3 Genetic background across-cohort

We conduct principal component analysis based on reported allele frequencies (**fPCA**) and use the population genetics Fst statistic to verify the genetic origin of each cohort (**fStructure**).

#### Step 2.4 Determine m and k

Central site determines m and k upon the survived SNPs. The number of shared SNPs are enough for identifying 1st-degree relatedness, we would offer a list of shared SNPs.

Central site sends **GenerateRandMat.R**, **random seed**, **k**, **an SNP list**, and **1KG-CHN binary plink format files** to each collaborator.

### Step 3 Encrypt genotype matrix

This step is similar to generate risk profile score in genetic prediction. Although a routine profile scoring step is very unlikely misconducted alone, unfortunate systematic mistakes may creep in because of some discordant reference alleles across the cohorts. As foolproof verification, every collaborator will receive 1KG-CHN and merge into their own data [Step 3.2](####Step-3.2-Merge-with-1KG-CHN).

Set environmental variables **(!!!!NEED TO BE MODIFIED!!!!)**

Again, please replace **YOUR_COHORT_ID** with the Cohort ID.

``` shell
user=*YOUR_COHORT_ID*
```

#### Step 3.1 Implement randomization

Users are asked to generate an m-by-k matrix whose elements are sample from N(0,1/m)

Obtain the number of markers (m)

``` shell
awk 'END{print NR}' Golden.snpA1 > Golden.m
```

This Rscript automatically reads parameters stored in Golden.m and Golden.k and generate Golden.key, an m-by-k matrix.

``` shell
Rscript GenerateRandMat.R Golden
```

You are supposed to see a matrix like this

|         |         |         |      |
|---------|---------|---------|------|
| 0.0373  | -0.0250	| 0.0309  | \... | 
| -0.0123	| 0.0443  | -0.0060	| \... | 
| -0.0159	| -0.0019	| 0.0426  | \... | 


Combine key with SNPID and A1 alleles by columns.

``` shell
paste -d "\t" Golden.snpA1 Golden.key > Golden.snpA1key
```

You are supposed to see a data table like this

| SNP | A1  | K_1     | K_2     | K_3     | \... |
|-----|-----|---------|---------|---------|------|
| rs1 | G   | 0.0373  | -0.0250 | 0.0309  | \... |
| rs2 | A   | -0.0123 | 0.0443  | -0.0060 | \... |
| rs3 | G   | -0.0159 | -0.0019 | 0.0426  | \... |

#### Step 3.2 Merge with 1KG-CHN (foolproof verification)

First extract SNPs in your bfiles by plink.

``` shell
awk '{print $1}' Golden.snpA1 > Golden.snp
plink --bfile ${user} --extract Golden.snp --make-bed --out ${user}.extract
```

Then merge 1KG-CHN with your bfiles into new plink files named Golden.merged

These files will be used in [Step 3.3](####Step-3.3-Genotype-encryption) plink2 "--score".

``` shell
plink --bfile 1KG-CHN.extract --bmerge ${user}.extract --make-bed --out Golden.merged
```

#### Step 3.3 Genotype encryption

Users are asked to encrypt their genotype matrix with the random matrix by plink2.0 and return the encrypted genotype matrix to central site.

Use plink2 "--score" to generate encrypted genotype matrix. (variance-standardize: genotypes are scaled by SNP)

``` shell
k=`cat Golden.k`
plink2 --bfile Golden.merged --score Golden.snpA1key 1 2 variance-standardize --score-col-nums 3-$(($k+2)) --out Golden.${user}
```

Return **Golden.\${user}.sscore** to central site.

We adopt one of the plink2.0 formats ".sscore" [plink2.0 formats:sscore](https://www.cog-genomics.org/plink/2.0/formats#sscore) as a standard sharing format for sharing encrypted genotype data. ".sscore" files include the following contents. 

| FID   | IID   | ALLELE_CT | NAMED_ALLELE_DOSAGE_SUM   | SCORE1_AVG | SCORE2_AVG | SCORE3_AVG | \...  |
|-------|-------|-----------|---------------------------|------------|------------|------------|-------|
| FID1  | IID1  | 996       | 273                       | -4.151E-04 | 5.563E-04  | -3.861E-04 | \...  |
| FID2  | IID2  | 970       | 267                       | -8.676E-05 | 1.800E-04  | 1.612E-03  | \...  |
| FID3  | IID3  | 990       | 273                       | -5.005E-04 | 3.104E-05  | -6.440E-04 | \...  |

#### Step 3.2 Merge with 1KG-CHN (foolproof verification)

First extract SNPs in your bfiles by plink.

``` shell
awk '{print $1}' Golden.snpA1 > Golden.snp
plink --bfile ${user} --extract Golden.snp --make-bed --out ${user}.extract
```

Then merge 1KG-CHN with your bfiles into new plink files named Golden.merged.

These files will be used in [Step 3.3](####Step-3.3-Genotype-encryption) plink2 "--score".

``` shell
plink --bfile 1KG-CHN.extract --bmerge ${user}.extract --make-bed --out Golden.merged
```

#### Step 3.3 Genotype encryption

Users are asked to encrypt their genotype matrix with the same random matrix by plink2.0 and return the encrypted genotype matrix to central site.

Use plink2 "--score" to generate encrypted genotype matrix.

``` shell
k=`cat Golden.k`
plink2 --bfile Golden.merged --score Golden.snpA1key 1 2 variance-standardize --score-col-nums 3-$(($k+2)) --out Golden.${user}
```

(variance-standardize: genotypes are scaled by SNP)

Return **Golden.\${user}.sscore** to central site.

### Step 4 Perform encG-reg across cohorts

Cohort-wise comparison for overlapping relatives will be conducted by central site. A foolproof implementation in [Step 3.2](####Step-3.2-Merge-with-1KG-CHN) leads to at least 1KG-CHN samples consistently identified as "overlap" between every pair of cohorts in [Step 4](###-Step-4-Perform-encG-reg-across-cohorts). Looking forward other possible overlapping that may pop out as expected as unexpected.

Bingo!

## Citations
Zhang Q-X, Liu T, Guo X, Zhen J, Yang M-y, Khederzadeh S, et al. (2024) Searching across-cohort relatives in 54,092 GWAS samples via encrypted genotype regression. PLoS Genet 20(1):e1011037. https://doi.org/10.1371/journal.pgen.1011037    


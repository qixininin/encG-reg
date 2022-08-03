# encG-reg

## 1-Simulations
This github contains simulation codes of GRM, encGRM and encG-reg.
- Fig S1 Heatmap presenting the role random matrix played in matrix multiplication. 
- Fig S2 Sampling variance of GRM, encGRM and encG-reg in simulations. 
- Fig S3 Validation for the sampling variance for GRM (assumption: binomial distribution). 
- Fig S4 Resolution for varying relatedness using GRM, encGRM and encG-reg. 


## 2-Protocols
We offered a user-friendly protocol for encG-reg. This protocol can be automated, such as by a web server that coordinates the study. 
There are four steps in total, where steps 1 and 3 are performed by each collaborator and steps 2 and 4 are performed by a central analyst.
We provide suggested commands for your possible reference, and your environment should have plink1.9, plink2.0 and R installed.

- Step 1 Within-cohort quality controls
- Step 2 Determine m and k
- Step 3 Encrypt genotype matrix
- Step 4 Perform encGRM across cohorts

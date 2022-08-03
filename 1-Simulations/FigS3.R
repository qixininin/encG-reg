## Figure S3 Validation for the sampling variance for GRM (assumption: binomial distribution)
## Note: This may take you ~1.5h
## encGRMsource.R includes functions:
##                RandomMatrixEncryption()
##                GenerateGeno()
##                GenerateGeno_r()
##                DprimetoD()
##                GenerateLDGeno()

setwd("~/Desktop/Cryptography/github/1-Simulations/")
source("encGRMsource.R")
library(ggplot2)
library(ggpubr)
library(egg)

n1 = 1000                       # sample size for pop 1
n2 = 1000                       # sample size for pop 2
M = 2000                        # number of markers
count=10                        # replication
p = list()
for(r in 0:3)
{
  seta=(0.5)^r
  freqvec = seq(0.05, 0.5, 0.1)
  the = (1-2*seta+seta/(2*freqvec*(1-freqvec)))/M
  obs = matrix(NA, count, length(freqvec))
  
  ## a range of MAFs
  for(f in 1:length(freqvec))
  {
    freq = rep(freqvec[f], M)
    for(c in 1:count)
    {
      Gr = GenerateGeno_r(freq, n1, r, T)
      X1 = Gr[[1]]
      X2 = Gr[[2]]
      
      ## Scale by allele frequency
      mean = apply(X1, 2, mean)/2
      X1 = t(apply(X1, 1, function(x) {(x-2*mean)/sqrt(2*mean*(1-mean))}))
      mean = apply(X2, 2, mean)/2
      X2 = t(apply(X2, 1, function(x) {(x-2*mean)/sqrt(2*mean*(1-mean))}))
      
      ## GRM
      K12 = tcrossprod(X1, X2) / M
      obs[c,f] = var(diag(K12))
    }
  }
  
  rst = data.frame(grp=rep(c("observed","theoretical"), each=length(freqvec)), 
                   maf=as.factor(rep(freqvec, 2)), 
                   score=c(apply(obs, 2, mean), the),
                   sd=c(apply(obs, 2, sd), rep(0,length(freqvec))))
  p[[r+1]] = ggplot(rst, aes(x=maf, y=score, fill=grp, color=grp))+
    geom_bar(stat="identity", width = 0.6, position=position_dodge(width = 0.8))+
    geom_errorbar(aes(ymin=score-sd, ymax=score+sd), color="black", width=0.2, position=position_dodge(0.8))+
    geom_text(aes(label=round(score,4)),vjust=-1.2, color="black", size=2.5, position=position_dodge(0.8))+
    scale_color_manual(values=c("chartreuse3","grey80"), name="") +
    scale_fill_manual(values=c("chartreuse3","grey80"), name="") +
    labs(x="MAF", y="var(GRM)", title = c("Identical","1st-degree","2nd-degree","3rd-degree")[r+1])+
    theme_article()+
    theme(legend.position = "bottom", 
          plot.title = element_text(hjust = 0.5, size = 9), 
          axis.line = element_line(colour = "black"),
          panel.border = element_blank()) 
  print(r+1)
}

png("varGRM-binomial.png", width = 1200, height = 350)
ggpubr::ggarrange(p[[1]],p[[2]],p[[3]],p[[4]], ncol = 4, nrow = 1, common.legend = T, legend = "right")
dev.off()

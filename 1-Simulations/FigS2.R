## Figure S2 Validation for the sampling variance for GRM, encGRM, and encG-reg
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

title = c("Identical","1st-degree","2nd-degree","3rd-degree")
#### GRM ####
n1 = 1000                       # sample size for pop 1
n2 = 1000                       # sample size for pop 2
Mvec = seq(1000, 2000, 250)     # number of markers
Mnum = length(Mvec)
p = list()
dtp =list()
for(r in 0:3)
{
  rst = vector()
  for(i in 1:Mnum)
  {
    M = Mvec[i]
    freq = runif(M, 0.05, 0.5)
    ## r-degree relatives
    Gr = GenerateGeno_r(freq, n1, r, FALSE)
    X1 = Gr[[1]]
    X2 = Gr[[2]]
    X1 = scale(X1)
    X2 = scale(X2)
    ## GRM calculated from individual genotype
    K12 = tcrossprod(X1, X2) / M
    rst = c(rst, var(diag(K12)))
  }
  dtp[[r+1]] = data.frame(grp = rep(c("Theoretical", "Observed"), each = Mnum),
                          r = r,
                          M = as.factor(rep(Mvec, 2)),
                          var = c((1+(0.25)^r)/Mvec, rst))
  p[[r+1]] = ggplot(data = dtp[[r+1]], aes(x=M, y=var, color=grp, fill=grp))+
    geom_bar(stat="identity", width = 0.6, position=position_dodge(width = 0.8))+
    scale_color_manual(values=c("chartreuse3","grey80"), name="") +
    scale_fill_manual(values=c("chartreuse3","grey80"), name="") +
    scale_y_continuous(limits = c(0,0.0025))+
    labs(title = title[r+1], y="Var(GRM)")+
    theme_article()+
    theme(plot.title = element_text(hjust = 0.5, size=12),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.border = element_blank())
}

p[[1]] = p[[1]]
p[[2]] = p[[2]] + theme(axis.title.y = element_blank())
p[[3]] = p[[3]] + theme(axis.title.y = element_blank())
p[[4]] = p[[4]] + theme(axis.title.y = element_blank())

#### encGRM ####
n1 = 1000                       # sample size for pop 1
n2 = 1000                       # sample size for pop 2
Mvec = seq(1000, 2000, 250)     # number of markers
Kvec = ceiling(Mvec/(1/0.8-1))  # column numbers of random matrix
Mnum = length(Mvec)
q = list()
dtq = list()
for(r in 0:3)
{
  rst = vector()
  for(i in 1:Mnum)
  {
    M = Mvec[i]
    K = Kvec[i]
    freq = runif(M, 0.05, 0.5)
    ## r-degree relatives
    Gr = GenerateGeno_r(freq, n1, r, FALSE)
    X1 = Gr[[1]]
    X2 = Gr[[2]]
    X1 = scale(X1)
    X2 = scale(X2)
    ## encryption
    S = matrix(rnorm(M*K,sd=sqrt(1/M)), M, K)
    X1hat = tcrossprod(X1, t(S))
    X2hat = tcrossprod(X2, t(S))
    ## encGRM calculated from encrypted genotype
    K12hat = tcrossprod(X1hat, X2hat) / K
    rst = c(rst, var(diag(K12hat)))
  }
  dtq[[r+1]] = data.frame(grp = rep(c("Theoretical", "Observed"), each = Mnum),
                          r = r,
                          M = as.factor(rep(Mvec, 2)),
                          var = c((1+(0.25)^r)/Mvec+(1+(0.25)^r)/Kvec, rst))
  q[[r+1]] = ggplot(data = dtq[[r+1]], aes(x=M, y=var, color=grp, fill=grp))+
    geom_bar(stat="identity", width = 0.6, position=position_dodge(width = 0.8))+
    scale_color_manual(values=c("deepskyblue3","grey80"), name="") +
    scale_fill_manual(values=c("deepskyblue3","grey80"), name="") +
    scale_y_continuous(limits = c(0,0.003), breaks = c(0, 0.001, 0.002, 0.003), labels = c("0.0000","0.0010","0.0020","0.0030"))+
    labs(title = title[r+1], y="Var(encGRM)")+
    theme_article()+
    theme(plot.title = element_text(hjust = 0.5, size=12),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.border = element_blank())
}

q[[1]] = q[[1]]
q[[2]] = q[[2]] + theme(axis.title.y = element_blank())
q[[3]] = q[[3]] + theme(axis.title.y = element_blank())
q[[4]] = q[[4]] + theme(axis.title.y = element_blank())


#### encG-reg ####
n1 = 1000                       # sample size for pop 1
n2 = 1000                       # sample size for pop 2
Mvec = seq(1000, 2000, 250)     # number of markers
Kvec = ceiling(Mvec/(1/0.8-1))  # column numbers of random matrix
Mnum = length(Mvec)
s = list()
dts = list()
for(r in 0:3)
{
  rst = vector()
  for(i in 1:Mnum)
  {
    M = Mvec[i]
    K = Kvec[i]
    freq = runif(M, 0.05, 0.5)
    ## r-degree relatives
    Gr = GenerateGeno_r(freq, n1, r, FALSE)
    X1 = Gr[[1]]
    X2 = Gr[[2]]
    X1 = scale(X1)
    X2 = scale(X2)
    ## encryption
    S = matrix(rnorm(M*K,sd=sqrt(1/M)), M, K)
    X1hat = tcrossprod(X1, t(S))
    X2hat = tcrossprod(X2, t(S))
    X1hat = t(apply(X1hat, 1, scale))
    X2hat = t(apply(X2hat, 1, scale))
    ## regression with encrypted genotype
    B = c()
    for(j in 1:n1)
    {
      mod = lm(X1hat[j,]~X2hat[j,])
      B = c(B, summary(mod)$coefficients[2,1])
    }
    rst = c(rst, var(B))

  }
  dts[[r+1]] = data.frame(grp = rep(c("Theoretical", "Observed"), each = Mnum),
                          r = r,
                          M = as.factor(rep(Mvec, 2)),
                          var = c((1-(0.25)^r)/Mvec+(1-(0.25)^r)/Kvec, rst))
  s[[r+1]] = ggplot(data = dts[[r+1]], aes(x=M, y=var, color=grp, fill=grp))+
    geom_bar(stat="identity", width = 0.6, position=position_dodge(width = 0.8))+
    scale_color_manual(values=c("orangered","grey80"), name="") +
    scale_fill_manual(values=c("orangered","grey80"), name="") +
    scale_y_continuous(limits = c(0,0.0015))+
    labs(title = title[r+1], y="Var(encG-reg)")+
    theme_article()+
    theme(plot.title = element_text(hjust = 0.5, size=12),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.border = element_blank())
}

s[[1]] = s[[1]]
s[[2]] = s[[2]] + theme(axis.title.y = element_blank())
s[[3]] = s[[3]] + theme(axis.title.y = element_blank())
s[[4]] = s[[4]] + theme(axis.title.y = element_blank())

png("varGRM_encGRM_encGreg.png", width = 1200, height = 850)
egg::ggarrange(ggpubr::ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], nrow = 1, ncol = 4, labels = c("A1","A2","A3","A4"), common.legend = T, legend = "right"),
               ggpubr::ggarrange(q[[1]], q[[2]], q[[3]], q[[4]], nrow = 1, ncol = 4, labels = c("B1","B2","B3","B4"), common.legend = T, legend = "right"),
               ggpubr::ggarrange(s[[1]], s[[2]], s[[3]], s[[4]], nrow = 1, ncol = 4, labels = c("C1","C2","C3","C4"), common.legend = T, legend = "right"),nrow = 3)
dev.off()

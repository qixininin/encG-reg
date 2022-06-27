## Figure S2 Validation for the sampling variance for GRM and encGRM
## encGRMsource.R includes functions:
##                RandomMatrixEncryption()
##                GenerateGeno()
##                GenerateGeno_r()
##                DprimetoD()
##                GenerateLDGeno()

setwd("~/Desktop/Cryptography/github/")
source("encGRMsource.R")
library(ggplot2)
library(egg)

color = c("olivedrab3","lightyellow3")
#### FigS2 a-d ####
n1 = 1000                       # sample size for pop 1
n2 = 1000                       # sample size for pop 2
Mvec = seq(1000, 2000, 250)     # number of markers
Mnum = length(Mvec)
p = list()
for(r in 0:3)
{
  rst = vector()
  title = c("Identical","1st-degree","2nd-degree","3rd-degree")[r+1]
  for(i in 1:Mnum)
  {
    M = Mvec[i]
    freq = runif(M, 0.05, 0.5)
    ## r-degree relatives
    Gr = GenerateGeno_r(freq, n1, r, FALSE)
    X1 = Gr[[1]]
    X2 = Gr[[2]]
    X1 = t(apply(X1, 1, scale))
    X2 = t(apply(X2, 1, scale))  
    
    ## GRM calculated from individual genotype
    K12 = tcrossprod(X1, X2) / M
    rst = c(rst, var(diag(K12)))
  }
  dt = data.frame(grp = rep(c("Theoretical", "Observed"), each = Mnum),
                  r = r,
                  M = as.factor(rep(Mvec, 2)),
                  var = c((1-(0.25)^r)/Mvec, rst))
  p[[r+1]] = ggplot(data = dt, aes(x=M, y=var, color=grp, fill=grp))+
    geom_bar(stat="identity", width = 0.6, position=position_dodge(width = 0.8))+
    scale_color_manual(values=color, name="") +
    scale_fill_manual(values=color, name="") +
    scale_y_continuous(limits = c(0,0.0015))+
    labs(title = title,
         y="Var(GRM)")+
    theme_article()+
    theme(plot.title = element_text(hjust = 0.5, size=9),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank())
}


#### FigS2 e-h ####
n1 = 1000                       # sample size for pop 1
n2 = 1000                       # sample size for pop 2
Mvec = seq(1000, 2000, 250)     # number of markers
Kvec = ceiling(Mvec/(1/0.8-1))  # column numbers of random matrix
Mnum = length(Mvec)
q = list()
for(r in 0:3)
{
  rst = vector()
  title = c("Identical","1st-degree","2nd-degree","3rd-degree")[r+1]
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
    ## encGRM calculated from encrypted genotype
    K12hat = tcrossprod(X1hat, X2hat) / K
    rst = c(rst, var(diag(K12hat)))
  }
  dt = data.frame(grp = rep(c("Theoretical", "Observed"), each = Mnum),
                  r = r,
                  M = as.factor(rep(Mvec, 2)),
                  var = c((1-(0.25)^r)/Mvec+(1-(0.25)^r)/Kvec, rst))
  q[[r+1]] = ggplot(data = dt, aes(x=M, y=var, color=grp, fill=grp))+
    geom_bar(stat="identity", width = 0.6, position=position_dodge(width = 0.8))+
    scale_color_manual(values=color, name="") +
    scale_fill_manual(values=color, name="") +
    scale_y_continuous(limits = c(0,0.0015))+
    labs(title = title,
         y="Var(encGRM)")+
    theme_article()+
    theme(plot.title = element_text(hjust = 0.5, size=9),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank())
}

#### FigS2 i-l ####
n1 = 200                        # sample size for pop 1
n2 = 200                        # sample size for pop 2
M  = 500                        # number of markers
freq = runif(M, 0.05, 0.5)      # maf for markers
Kvec = seq(2*M, 10*M, 2*M)      # column numbers of random matrix
knum = length(Kvec)
count = 20                      # repeats

## genotype without LD
X1 = GenerateGeno(freq, n1)
X2 = GenerateGeno(freq, n2)
X1 = scale(X1)
X2 = scale(X2)
K21 = tcrossprod(X2, X1) / M
K11 = tcrossprod(X1, X1) / M
Me = 1/var(K11[lower.tri(K11)])
corrK21 = matrix(NA, knum, count)
time_LE = matrix(NA, knum, count)
for(c in 1:count)
{
  for(i in 1:knum)
  {
    t0 = proc.time()
    K = Kvec[i]
    A = RandomMatrixEncryption(X1, X2, M, K, FALSE)
    K21hat = tcrossprod(A[[2]], A[[1]])/ K
    t1 = proc.time()
    time_LE[i,c] = (t1-t0)[1]
    corrK21[i,c] = cor(as.vector(K21), as.vector(K21hat))^2
  }
}
corr2 = data.frame(grp=rep(c("observed","theoretical"), each=knum), K=as.factor(rep(Kvec, 2)), 
                   rsq=c(apply(corrK21, 1, mean), 1/(1+Me/Kvec)),
                   sd=c(apply(corrK21, 1, sd), rep(0,knum)))

## genotype with LD
Dprime = runif(M, 0.5, 0.9)
ld = DprimetoD(freq, Dprime)
X1 = GenerateLDGeno(freq, ld, n1)
X2 = GenerateLDGeno(freq, ld, n2)
X1 = scale(X1)
X2 = scale(X2)
K21 = tcrossprod(X2, X1) / M
K11 = tcrossprod(X1, X1) / M
Me = 1/var(K11[lower.tri(K11)])
corrK21 = matrix(NA, knum, count)
time_LD = matrix(NA, knum, count)
for(c in 1:count)
{
  for(i in 1:knum)
  {
    t0 = proc.time()
    K = Kvec[i]
    A = RandomMatrixEncryption(X1, X2, M, K, FALSE)
    K21hat = tcrossprod(A[[2]], A[[1]])/ K
    corrK21[i,c] = cor(as.vector(K21), as.vector(K21hat))^2
    t1 = proc.time()
    time_LD[i,c] = (t1-t0)[1]
  }
}
corr1 = data.frame(grp=rep(c("observed","theoretical"), each=knum), K=as.factor(rep(Kvec, 2)), 
                   rsq=c(apply(corrK21, 1, mean), 1/(1+Me/Kvec)),
                   sd=c(apply(corrK21, 1, sd), rep(0,knum)))

## plot rsq
corr1$rsq = round(corr1$rsq,3)
corr2$rsq = round(corr2$rsq,3)
p1 = ggplot(corr2, aes(x=K, y=rsq, fill=grp, color=grp))+
  geom_bar(stat="identity", width = 0.6, position=position_dodge(width = 0.8))+
  geom_errorbar(aes(ymin=rsq-sd, ymax=rsq+sd), color="black", width=0.2, position=position_dodge(0.8))+
  geom_text(aes(label=round(rsq,4)),vjust=-1.2, color="black", size=2.5, position=position_dodge(0.8))+
  scale_color_manual(values=color, name="") +
  scale_fill_manual(values=color, name="") +
  labs(title = "without LD",
       y=expression(R^2))+
  theme_article()+
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 9), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank()) 
p2 = ggplot(corr1, aes(x=K, y=rsq, fill=grp, color=grp))+
  geom_bar(stat="identity", width = 0.6, position=position_dodge(width = 0.8))+
  geom_errorbar(aes(ymin=rsq-sd, ymax=rsq+sd), color="black", width=0.2, position=position_dodge(0.8))+
  geom_text(aes(label=round(rsq,4)),vjust=-1.2, color="black", size=2.5, position=position_dodge(0.8))+
  scale_color_manual(values=color, name="") +
  scale_fill_manual(values=color, name="") +
  labs(title = "with LD",
       y=expression(R^2))+
  theme_article()+
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 9), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank()) 

## plot time
time1 = data.frame(K=as.factor(Kvec), t=apply(time_LE, 1, mean), sd=apply(time_LE, 1, sd))
time2 = data.frame(K=as.factor(Kvec), t=apply(time_LD, 1, mean), sd=apply(time_LD, 1, sd))
p3 = ggplot(time1, aes(x=K, y=t, group=1))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=t-sd, ymax=t+sd), width=.2,
                position=position_dodge(0.05))+
  labs(title = "without LD",
       y="CPU time(secs)")+
  theme_article()+
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size=9), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank()) 

p4 = ggplot(time2, aes(x=K, y=t, group=1))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=t-sd, ymax=t+sd), width=.2,
                position=position_dodge(0.05))+
  labs(title = "with LD",
       y="CPU time(secs)")+
  theme_article()+
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size = 9), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank()) 


png("FigS2.png", width = 1800, height = 1280)
ggarrange(p[[1]]+theme(legend.position = "none"),
          p[[2]]+theme(legend.position = "none"),
          p[[3]]+theme(legend.position = "none"),
          p[[4]]+theme(legend.position = "none"),
          q[[1]]+theme(legend.position = "none"),
          q[[2]]+theme(legend.position = "none"),
          q[[3]]+theme(legend.position = "none"),
          q[[4]]+theme(legend.position = "right"),
          p1, p2, p3, p4, ncol=4,
          labels = c("a","b","c","d","e","f","g","h","i","j","k","l"))
dev.off()


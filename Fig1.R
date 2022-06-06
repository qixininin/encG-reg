## Figure 1 Resolution for detecting relatives in simulated encrypted genotypes
## encGRMsource.R includes functions:
##                RandomMatrixEncryption()
##                GenerateGeno()
##                GenerateGeno_r()

setwd("~/Desktop/")
source("encGRMsource.R")
library(ggplot2)
library(egg)

color = c("indianred","sienna1","goldenrod2","olivedrab3","lightyellow3")
n1    = 200                           # sample size in cohort 1
n2    = 200                           # sample size in cohort 2
alpha = 0.05/(n1*n2)                  # type I error
beta  = 0.1                           # type II error two tail test
D     = 4                             # top degree of relatedness included in simulation
r     = c(0:(D-1))                    # a vector of various degree of relatedness included in simulation
n_cp  = rep(10, D)                    # numbers for each degree of relatedness
n1_0  = n1-sum(n_cp)                  # remaining sample size of unrelated individuals in cohort 1
n2_0  = n2-sum(n_cp)                  # remaining sample size of unrelated individuals in cohort 2
Mvec  = c(50,500,5000,10000)          # a vector of various number of markers
setavec=(0.5)^r                       # a vector of expected relatedness score for various degree of relatedness
setavec=setavec*0.9                   # a discount of expected relatedness score
p=1
dt = list()
for(m in 1:length(Mvec))
{
  for(k in 1:m)
  {
    M = Mvec[m]
    # Determine K
    seta = setavec[k]
    tmp = ( ( qnorm(1-beta)*sqrt(1-seta^2) + qnorm(1-alpha) ) / seta )^2 
    K = ceiling(1/(1/tmp-1/M))
    
    # Genotypes 
    freq = runif(M, 0.05, 0.5)
    X1 = matrix(NA, 0, M)
    X2 = matrix(NA, 0, M)
    for(i in 1:D)
    {
      Gr = GenerateGeno_r(freq, n_cp[i], r[i], FALSE)
      X1 = rbind(X1, Gr[[1]])
      X2 = rbind(X2, Gr[[2]])
    }
    X1 = rbind(X1, GenerateGeno(freq, n1_0))
    X2 = rbind(X2, GenerateGeno(freq, n2_0))
    X1 = scale(X1)
    X2 = scale(X2)

    ## Encryption
    A = RandomMatrixEncryption(X1, X2, M, K, TRUE)
    vG21hat = as.vector( tcrossprod(A[[2]], A[[1]])/ K )
    vG21    = as.vector( tcrossprod(X2, X1)        / M )
    
    ## Save data
    yin = qnorm(0.05/n1/n2)*sqrt(1/M)
    xin = qnorm(0.05/n1/n2)*sqrt(1/M+1/K)
    lab = lab_relations_true(n1, n2, n_cp)
    R = k-1
    dt[[p]] = data.frame(x=vG21hat, y=vG21, lab=as.factor(lab), R = R, M=M, K=K, yin=yin, xin=xin)
    
    p=p+1
  }
}

## plot
dt0 = do.call("rbind", dt)
dt0$K = as.factor(dt0$K)
dt0$M = as.factor(dt0$M)
dt0$R = as.factor(dt0$R)

png("Fig1.png", width = 780, height = 880)
ggplot(dt0, aes(x=x, y=y)) + 
  geom_point(aes(color=lab), size=0.3) +
  geom_hline(aes(yintercept = yin), color="blue", alpha=0.5, linetype = 5) +
  geom_hline(aes(yintercept = -yin), color="blue", alpha=0.5, linetype = 5) +
  geom_vline(aes(xintercept = xin), color="red", alpha=0.5, linetype = 5) +
  geom_vline(aes(xintercept = -xin), color="red", alpha=0.5, linetype = 5) +
  scale_x_continuous(limits=c(-1.15, 1.15)) +
  xlab("Estimated relatedness") +
  scale_y_continuous(limits=c(-1.15, 1.15)) +
  ylab("Relatedness") +
  scale_color_manual(values=color, name="Relations", 
                     labels=c("Identical","1st-degree","2nd-degree","3rd-degree","Unrelated")) +
  facet_grid(M~R, labeller = label_both) + 
  theme_article()+
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        legend.title = element_text(size = 9,face="bold"))
dev.off()

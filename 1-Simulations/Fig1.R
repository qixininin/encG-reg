## Figure 2 Resolution for varying relatedness using GRM, encGRM and encG-reg
## encGRMsource.R includes functions:
##                RandomMatrixEncryption()
##                GenerateGeno()
##                GenerateGeno_r()

setwd("~/Desktop/Cryptography/github/1-Simulations/")
source("encGRMsource.R")
library(ggplot2)
library(ggpubr)
library(egg)
library(svglite)
library(gridExtra)
library(grid)
library(cowplot)
library(jcolors)

color = c("grey80","#66CD00","#009ACD","#FF4500")
n1    = 200                           # sample size in cohort 1
n2    = 200                           # sample size in cohort 2
alpha = 0.05/(n1*n2)                  # type I error
beta  = 0.1                           # type II error two tail test
D     = 3                             # top degree of relatedness included in simulation
r     = c(0:(D-1))                    # a vector of various degree of relatedness included in simulation
n_cp  = c(10,10,10,0)                 # numbers for each degree of relatedness
n1_0  = n1-sum(n_cp)                  # remaining sample size of unrelated individuals in cohort 1
n2_0  = n2-sum(n_cp)                  # remaining sample size of unrelated individuals in cohort 2
setavec=(0.5)^r                       # a vector of expected relatedness score for various degree of relatedness
Mvec = ceiling(2 * ( ( qnorm(1-beta)*sqrt(1+setavec^2) + qnorm(1-alpha) ) / setavec )^2) 
# a vector of various number of markers based on Eq3 (2 times larger)

#### GRM vs encG-reg ####
count=1
dt = list()
for(m in 1:length(Mvec))
{
  M = Mvec[m]
  for(k in 1:m)
  {
    K = ceiling(1/(1/(( ( qnorm(1-beta)*sqrt(1-setavec[k]^2) + qnorm(1-alpha) ) / setavec[k] )^2 )-1/Mvec[m])) # based on Eq4
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
    S = matrix(rnorm(M*K,sd=sqrt(1/M)), M, K)
    X1hat <- tcrossprod(X1, t(S))
    X2hat <- tcrossprod(X2, t(S))
    X1hat <- t(apply(X1hat, 1, scale))
    X2hat <- t(apply(X2hat, 1, scale))
    
    
    ## GRM and encG-reg
    G21hat   = tcrossprod(X2hat, X1hat)/ K
    G21      = tcrossprod(X2, X1) / M
    
    ## Save data
    yin = qnorm(0.05/n1/n2)*sqrt(1/M)
    xin = qnorm(0.05/n1/n2)*sqrt(1/M+1/K)
    lab = lab_relations_true(n1, n2, n_cp)
    R = k-1
    dt[[count]] = data.frame(x=as.vector(G21hat), y=as.vector(G21), lab=as.factor(lab), R=R, M=M, K=K, yin=yin, xin=xin)
    dt[[count]] = dt[[count]] %>% dplyr::arrange(lab)
    count=count+1
  }
}

## main plot
p = list()
for (i in 1:(count-1))
{
  M = dt[[i]]$M[1]
  K = dt[[i]]$K[1]
  p[[i]] = ggplot(dt[[i]], aes(x=x,y=y,grp=lab))+
    geom_point(aes(color=lab), size=0.3) +
    geom_hline(aes(yintercept = yin), color="blue", alpha=0.5, linetype = 5) +
    geom_hline(aes(yintercept = -yin), color="blue", alpha=0.5, linetype = 5) +
    geom_vline(aes(xintercept = xin), color="red", alpha=0.5, linetype = 5) +
    geom_vline(aes(xintercept = -xin), color="red", alpha=0.5, linetype = 5) +
    scale_x_continuous(name = "", limits=c(-1.15, 1.15)) +
    scale_y_continuous(name = "",limits=c(-1.15, 1.15)) +
    scale_color_manual(values=color, name="Relatedness", 
                       labels=c("Unrelated","Identical","1st-degree","2nd-degree")) +  
    labs(title = bquote(italic(m[e])~"="~.(M)~","~italic(k)~"="~.(K)))+
    theme_article()+
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 9), 
          axis.line = element_line(colour = "black"),
          panel.border = element_blank())
}

## distribution plot
dpy = list()
for(m in 1:length(Mvec))
{
  M = Mvec[m]
  dpy[[m]] = ggplot()+
    stat_function(fun = dnorm, args = list(  0, sqrt(           1/M)), n = 1000, color = color[1])+
    stat_function(fun = dnorm, args = list(  1, sqrt((1+(0.25)^0)/M)), n = 1000, color = color[2])+
    stat_function(fun = dnorm, args = list(1/2, sqrt((1+(0.25)^1)/M)), n = 1000, color = color[3])+
    stat_function(fun = dnorm, args = list(1/4, sqrt((1+(0.25)^2)/M)), n = 1000, color = color[4])+
    scale_x_continuous(name = "", limits=c(-1.15, 1.15), position = "top") +
    scale_y_reverse(name = "") +
    coord_flip()+
    theme_article()+
    theme(legend.position = "none", 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
}
dpx = list()
for(k in 1:D)
{
  m = D
  M = Mvec[m]
  K = ceiling(1/(1/(( ( qnorm(1-beta)*sqrt(1-setavec[k]^2) + qnorm(1-alpha) ) / setavec[k] )^2 )-1/Mvec[m]))
  dpx[[k]] = ggplot()+
    stat_function(fun = dnorm, args = list(  0, sqrt(           1/M +            1/K)), n = 1000, color = color[1])+
    stat_function(fun = dnorm, args = list(  1, sqrt((1-(0.25)^0)/M + (1-(0.25)^0)/K)), n = 1000, color = color[2])+
    stat_function(fun = dnorm, args = list(1/2, sqrt((1-(0.25)^1)/M + (1-(0.25)^1)/K)), n = 1000, color = color[3])+
    stat_function(fun = dnorm, args = list(1/4, sqrt((1-(0.25)^2)/M + (1-(0.25)^2)/K)), n = 1000, color = color[4])+
    scale_x_continuous(name = "", limits=c(-1.15, 1.15)) +
    scale_y_continuous(name = "") +
    theme_article()+
    theme(legend.position = "none", 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
}

## blank plot
blankPlot = ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

## fig lengend
leg <- get_legend(p[[1]] + 
                    theme(legend.position = "right", legend.title = element_text(size = 12)) + 
                    guides(color = guide_legend(override.aes = list(size=2))))

## output
final = egg::ggarrange(dpy[[1]] ,   p[[1]]   , blankPlot, blankPlot, 
                       dpy[[2]] ,   p[[2]]   , p[[3]]   , blankPlot, 
                       dpy[[3]] ,   p[[4]]   , p[[5]]   , p[[6]]   , 
                       blankPlot,   dpx[[1]] , dpx[[2]] , dpx[[3]] , 
                       nrow = 4 , ncol = 4   , widths = c(1,2,2,2),heights = c(2,2,2,1))
final = annotate_figure(final, left = "Relatedness by GRM", bottom = "Relatedness by encG-reg")
png("simu-MvsK-encG-reg.png", width = 750, height = 700)
# pdf("simu-MvsK-encG-reg.pdf", width = 10, height = 9)
ggpubr::ggarrange(final, as_ggplot(leg), nrow = 1, widths = c(6,1))
dev.off()


#### GRM vs encGRM ####
count=1
dt = list()
for(m in 1:length(Mvec))
{
  M = Mvec[m]
  for(k in 1:m)
  {
    K = ceiling(1/(1/(( ( qnorm(1-beta)*sqrt(1+setavec[k]^2) + qnorm(1-alpha) ) / setavec[k] )^2 )-1/Mvec[m]))
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
    S = matrix(rnorm(M*K,sd=sqrt(1/M)), M, K)
    X1hat <- tcrossprod(X1, t(S))
    X2hat <- tcrossprod(X2, t(S))
    
    ## GRM and encG-reg
    G21hat   = tcrossprod(X2hat, X1hat)/ K
    G21      = tcrossprod(X2, X1) / M
    
    ## Save data
    yin = qnorm(0.05/n1/n2)*sqrt(1/M)
    xin = qnorm(0.05/n1/n2)*sqrt(1/M+1/K)
    lab = lab_relations_true(n1, n2, n_cp)
    R = k-1
    dt[[count]] = data.frame(x=as.vector(G21hat), y=as.vector(G21), lab=as.factor(lab), R=R, M=M, K=K, yin=yin, xin=xin)
    dt[[count]] = dt[[count]] %>% dplyr::arrange(lab)
    count=count+1
  }
}

## main plot
p = list()
for (i in 1:(count-1))
{
  M = dt[[i]]$M[1]
  K = dt[[i]]$K[1]
  p[[i]] = ggplot(dt[[i]], aes(x=x,y=y,grp=lab))+
    geom_point(aes(color=lab), size=0.3) +
    geom_hline(aes(yintercept =  yin), color="blue", alpha=0.5, linetype = 5) +
    geom_hline(aes(yintercept = -yin), color="blue", alpha=0.5, linetype = 5) +
    geom_vline(aes(xintercept =  xin), color="red" , alpha=0.5, linetype = 5) +
    geom_vline(aes(xintercept = -xin), color="red" , alpha=0.5, linetype = 5) +
    scale_x_continuous(name = "", limits=c(-1.15, 1.15)) +
    scale_y_continuous(name = "", limits=c(-1.15, 1.15)) +
    scale_color_manual(values=color, name="Relatedness", 
                       labels=c("Unrelated","Identical","1st-degree","2nd-degree")) +  
    labs(title = bquote(italic(m[e])~"="~.(M)~","~italic(k)~"="~.(K)))+
    theme_article()+
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 9), 
          axis.line = element_line(colour = "black"),
          panel.border = element_blank())
}

## distribution plot
dpy = list()
for(m in 1:length(Mvec))
{
  M = Mvec[m]
  dpy[[m]] = ggplot()+
    stat_function(fun = dnorm, args = list(  0, sqrt(           1/M)), n = 1000, color = color[1])+
    stat_function(fun = dnorm, args = list(  1, sqrt((1+(0.25)^0)/M)), n = 1000, color = color[2])+
    stat_function(fun = dnorm, args = list(1/2, sqrt((1+(0.25)^1)/M)), n = 1000, color = color[3])+
    stat_function(fun = dnorm, args = list(1/4, sqrt((1+(0.25)^2)/M)), n = 1000, color = color[4])+
    scale_x_continuous(name = "", limits=c(-1.15, 1.15), position = "top") +
    scale_y_reverse(name = "") +
    coord_flip()+
    theme_article()+
    theme(legend.position = "none", 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
}
dpx = list()
for(k in 1:D)
{
  m = D
  K = ceiling(1/(1/(( ( qnorm(1-beta)*sqrt(1+setavec[k]^2) + qnorm(1-alpha) ) / setavec[k] )^2 )-1/Mvec[m]))
  M = Mvec[m]
  dpx[[k]] = ggplot()+
    stat_function(fun = dnorm, args = list(  0, sqrt(           1/M +            1/K)), n = 1000, color = color[1])+
    stat_function(fun = dnorm, args = list(  1, sqrt((1+(0.25)^0)/M + (1+(0.25)^0)/K)), n = 1000, color = color[2])+
    stat_function(fun = dnorm, args = list(1/2, sqrt((1+(0.25)^1)/M + (1+(0.25)^1)/K)), n = 1000, color = color[3])+
    stat_function(fun = dnorm, args = list(1/4, sqrt((1+(0.25)^2)/M + (1+(0.25)^2)/K)), n = 1000, color = color[4])+
    scale_x_continuous(name = "", limits=c(-1.15, 1.15)) +
    scale_y_continuous(name = "") +
    theme_article()+
    theme(legend.position = "none", 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
}

## blank plot
blankPlot = ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

## fig lengend
leg <- get_legend(p[[1]] + 
                    theme(legend.position = "right", legend.title = element_text(size = 12)) + 
                    guides(color = guide_legend(override.aes = list(size=2))))

## output
final = egg::ggarrange(dpy[[1]] ,   p[[1]]   , blankPlot, blankPlot, 
                       dpy[[2]] ,   p[[2]]   , p[[3]]   , blankPlot, 
                       dpy[[3]] ,   p[[4]]   , p[[5]]   , p[[6]]   , 
                       blankPlot,   dpx[[1]] , dpx[[2]] , dpx[[3]] , 
                       nrow = 4 , ncol = 4   , widths = c(1,2,2,2),heights = c(2,2,2,1))
final = annotate_figure(final, left = "Relatedness by GRM", bottom = "Relatedness by encGRM")
png("simu-MvsK-encGRM.png", width = 750, height = 700)
# pdf("simu-MvsK-encGRM.pdf", width = 10, height = 9)
ggpubr::ggarrange(final, as_ggplot(leg), nrow = 1, widths = c(6,1))
dev.off()

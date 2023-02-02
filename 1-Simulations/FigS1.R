## Figure S1 Heatmap presenting the role random matrix played in matrix multiplication
## encGRMsource.R includes functions:
##                GenerateGenoMatrix()

library(ggplot2)
library(egg)
library(viridis)
source("encGRMsource.R")
M = 100
K = 500
S = matrix(rnorm(M*K,sd=sqrt(1/K)), M, K)
SS = tcrossprod(S, S)

data = expand.grid(Y=seq(1:M),X=rev(seq(1:M)))
data$Z = as.vector(diag(M))
A = ggplot(data, aes(X,Y,fill=Z))+
  geom_tile()+
  labs(title=expression(I))+
  scale_fill_viridis(discrete=FALSE,direction = -1)+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,face = "bold"))

data = expand.grid(Y=seq(1:M),X=rev(seq(1:M)))
data$Z = as.vector(SS)
B = ggplot(data, aes(X,Y,fill=Z))+
  geom_tile()+
  labs(title=expression(SS^T))+
  scale_fill_viridis(discrete=FALSE,direction = -1)+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,face = "bold"))

n1=20
n2=25
freq=runif(M, 0.05, 0.5)
X1 = GenerateGenoMatrix(n1, freq, F)
X2 = GenerateGenoMatrix(n2, freq, F)
X1SSX2 = X1 %*% SS %*% t(X2)
X1X2 = X1 %*% t(X2)


data = expand.grid(Y=seq(1:n1),X=rev(seq(1:n2)))
data$Z = as.vector(X1X2)
C = ggplot(data, aes(X,Y,fill=Z))+
  geom_tile()+
  labs(title=expression(X[1]*X[2]))+
  scale_fill_viridis(discrete=FALSE,direction = -1)+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,face = "bold")) 

data = expand.grid(Y=seq(1:n1),X=rev(seq(1:n2)))
data$Z = as.vector(X1SSX2)
D = ggplot(data, aes(X,Y,fill=Z))+
  geom_tile()+
  labs(title=expression(X[1]*S*S^T*X[2]))+
  scale_fill_viridis(discrete=FALSE,direction = -1)+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,face = "bold")) 

png("simu-RandomMatrix.png", width=700, height=750)
ggarrange(A,B,C,D,nrow = 2)
dev.off()

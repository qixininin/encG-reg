## Figure S1 Heatmap presenting the role random matrix played in matrix multiplication
## encGRMsource.R includes functions:
##                GenerateGenoMatrix()

setwd("~/Desktop/")
source("encGRMsource.R")
M = 100
K = 500
S = matrix(rnorm(M*K,sd=sqrt(1/K)), M, K)
SS = tcrossprod(S, S)
heatmap(diag(M), symm = T, revC = T, 
        Rowv = NA, Colv = NA, labRow = NA, labCol = NA, main = "Identity matrix",
        col = hcl.colors(250, "Blues", rev = T))
heatmap(SS, symm = T, revC = T, 
        Rowv = NA, Colv = NA, labRow = NA, labCol = NA, main = "S-by-t(S)",
        col = hcl.colors(250, "Blues", rev = T))

n1=20
n2=25
freq=runif(M, 0.05, 0.5)
X1 = GenerateGenoMatrix(n1, freq, F)
X2 = GenerateGenoMatrix(n2, freq, F)
X1SSX2 = X1 %*% SS %*% t(X2)
X1X2 = X1 %*% t(X2)
heatmap(X1X2, revC = T, 
        Rowv = NA, Colv = NA, labRow = NA, labCol = NA, main = "X1X2",
        col = hcl.colors(250, "Blues", rev = T))
heatmap(X1SSX2, revC = T, 
        Rowv = NA, Colv = NA, labRow = NA, labCol = NA, main = "X1SSX2",
        col = hcl.colors(250, "Blues", rev = T))

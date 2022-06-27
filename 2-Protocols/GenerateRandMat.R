# This R file is used for users to calculate the same random matrix.

arg = commandArgs(T)
coflag = arg[1]

set.seed(as.numeric(read.table(paste0(coflag, ".seed"))))
K = as.numeric(read.table(paste0(coflag, ".k")))
M = as.numeric(read.table(paste0(coflag, ".m")))
S = matrix(rnorm(M*K,sd=sqrt(1/M)), M, K)

write.table(format(S, digits = 3), file=paste0(coflag, ".key"), quote = F, col.names = F, row.names = F)
######### Genome function include LD information #########
GenerateLDAllele <- function(freq, ld) # return gMat( M*2 ) for one individual two chromosome
{
  gMat = matrix(0, nrow=length(freq), ncol=2)
  for(i in 1:length(freq)) 
  {
    for(j in 1:2)
    {
      idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1)
      if(i == 1) 
      {
        gMat[i,j] = idx 
      }
      else 
      {
        d = runif(1, 0, 1)
        a = gMat[i-1,j]
        f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1])
        f2 = ifelse(a == 0, freq[i], 1-freq[i]) 
        gMat[i,j] = ifelse(d < (f1 * f2 +ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j]) 
      }
    }
  }
  return(gMat)
}
GenerateLDAllele_r <- function(freq, ld, r, c) # return gMat( M*2 ) for one individual two chromosome
{
  ld = ld*(1-c)^r
  gMat = matrix(0, nrow=length(freq), ncol=2)
  for(i in 1:length(freq))
  {
    for(j in 1:2)
    {
      idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1) 
      if(i == 1) 
      {
        gMat[i,j] = idx 
      }
      else 
      {
        d = runif(1, 0, 1)
        a = gMat[i-1,j]
        f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1])
        f2 = ifelse(a == 0, freq[i], 1-freq[i])
        gMat[i,j] = ifelse(d < (f1 * f2 +ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j]) 
      }
    }
  }
  return(gMat)
}
GenerateLDGeno <- function(freq, ld, N) # return unrelated g( N*M )
{
  g = matrix(0, nrow=N, ncol=length(freq))
  if(N==0) return(g)
  for(h in 1:N)
  {
    gMat = matrix(0, nrow=length(freq), ncol=2) 
    for(i in 1:length(freq))
    {
      for(j in 1:2)
      {
        idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1) 
        if(i == 1) 
        {
          gMat[i,j] = idx 
        }
        else # for not first site
        {
          d = runif(1, 0, 1) 
          a = gMat[i-1,j]
          f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1]) 
          f2 = ifelse(a == 0, freq[i], 1-freq[i]) 
          gMat[i,j] = ifelse(d < (f1 * f2 + ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j]) 
        }
      }
    }
    g[h,] = gMat[,1] + gMat[,2]
  }
  return(g)
}
DprimetoD <- function(freq, Dprime)
{
  M=length(freq)
  f1=freq
  f2=1-freq
  f1r=f1[-1];f1r=c(f1r,0)
  f2r=f2[-1];f2r=c(f2r,0)
  D=apply(cbind(f1*f2r,f2*f1r), 1, min)*Dprime
  return(D)
}
######### Genome function exclude LD information #########
GenerateAllele <- function(freq) # return gMat( 2*M )
{
  M = length(freq)
  gMat = matrix(NA, 2, M)
  for(i in 1:M)
  {
    gMat[, i] = rbinom(2, 1, freq[i])
  }
  return(gMat)
}
GenerateGeno <- function(freq, N) # return unrelated g( N*M )
{
  M = length(freq)
  g = matrix(NA, N, M)
  if(N==0) return(g)
  for(i in 1:M)
  {
    g[,i] = rbinom(N, 2, freq[i])
  }
  return(g)
}
GenerateGeno_r <- function(freq, N, r, sibflag=FALSE) # return r-degree relation g = list( N*M , N*M ) 
{
  M = length(freq)
  g = list(matrix(NA, N, M), matrix(NA, N, M))
  if(N==0) return(g) 
  if(r==0)
  {
    g[[1]]=g[[2]]=GenerateGeno(freq, N)
    return(g)
  }
  ibdscore = 0.5^r
  for(h in 1:N)
  {
    gMat1 = GenerateAllele(freq) #(2*M)
    gMat2 = GenerateAllele(freq) #(2*M)
    # sib alike
    if(sibflag)
    {
      for(j in 1:2)
      {
        ibd = sample(1:M, ceiling(M*ibdscore))
        gMat2[j,ibd] =  gMat1[j,ibd]
      }
    }
    # father-son alike
    else
    {
      ibd = sample(1:M, ceiling(M*ibdscore*2))
      gMat2[1,ibd] =  gMat1[1,ibd]
    }
    g[[1]][h, ] = apply(gMat1, 2, sum)
    g[[2]][h, ] = apply(gMat2, 2, sum)
  }
  return(g)
}

############# Generate two GenoMatrix ##############
GenerateGenoMatrix <- function(n, freq, flg_scale = T)
{
  X = GenerateGeno(freq, n)
  if(flg_scale){
    X = scale(X)
  }
  return(X)
}

############# lab function ###############

lab_relations_true <- function(n1, n2, n_cp)
{
  lab <- array(9, dim = n1*n2)
  n1_0 = n1-sum(n_cp)
  n2_0 = n2-sum(n_cp)
  
  set_0 = c()
  set_1 = c()
  set_2 = c()
  set_3 = c()
  
  l0 = 1
  l1 = (n_cp[1])*n2+n_cp[1]+1
  l2 = (n_cp[1]+n_cp[2])*n2+n_cp[1]+n_cp[2]+1
  l3 = (n_cp[1]+n_cp[2]+n_cp[3])*n2+n_cp[1]+n_cp[2]+n_cp[3]+1
  
  r0 = (n_cp[1])*n2
  r1 = (n_cp[1]+n_cp[2])*n2
  r2 = (n_cp[1]+n_cp[2]+n_cp[3])*n2  
  r3 = (n_cp[1]+n_cp[2]+n_cp[3]+n_cp[4])*n2 
  
  # detected relations set
  if(l0<r0) set_0 <- seq(l0, r0, n2+1)
  if(l1<r1) set_1 <- seq(l1, r1, n2+1)
  if(l2<r2) set_2 <- seq(l2, r2, n2+1)
  if(l3<r3) set_3 <- seq(l3, r3, n2+1)
  
  lab[set_0] <- 0
  lab[set_1] <- 1
  lab[set_2] <- 2
  lab[set_3] <- 3
  return(lab)
}

#### random matrix ####
RandomMatrixEncryption <- function(Mat1, Mat2, M, K, ScaleRow = T)
{
  W <- matrix(rnorm(M*K,sd=sqrt(1/M)), M, K)
  RandomMat1 <- tcrossprod(Mat1, t(W))
  RandomMat2 <- tcrossprod(Mat2, t(W))
  if(ScaleRow)
  {
    # scale by row
    RandomMat1 <- t(apply(RandomMat1, 1, scale))
    RandomMat2 <- t(apply(RandomMat2, 1, scale))
  } else {
    # scale by col
    RandomMat1 <- scale(RandomMat1)
    RandomMat2 <- scale(RandomMat2) 
  }
  return(list(RandomMat1,RandomMat2))
}

RandomVecEncryption <- function(Mat, Vec, M, K, ScaleRow = T)
{
  W <- matrix(rnorm(M*K,sd=sqrt(1/M)), M, K)
  RandomMat <- tcrossprod(Mat, t(W))
  RandomVec <- tcrossprod(Vec, t(W))
  if(ScaleRow)
  {
    # scale by row
    RandomMat <- t(apply(RandomMat, 1, scale))
  } else {
    # scale by col
    RandomMat <- scale(RandomMat)
  }
  return(list(RandomMat,RandomVec))
}

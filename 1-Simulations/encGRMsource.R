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
lab_relations <- function(vec, n1, n2)
{
  lab <- array(9, dim = n1*n2)
  # detected relations set
  set_1 <- which(vec>=0.45 & vec<0.95)
  lab[set_1] <- 1
  return(lab)
}
lab_relations_true <- function(n1, n2, n_cp)
{
  lab <- array(0, dim = n1*n2)
  n1_0 = n1-sum(n_cp)
  n2_0 = n2-sum(n_cp)
  
  set_0 = c()
  if(n_cp[2]>0) set_1 = c()
  if(n_cp[3]>0) set_2 = c()
  if(n_cp[4]>0) set_3 = c()
  
  l0 = 1
  if(n_cp[2]>0) l1 = (n_cp[1])*n2+n_cp[1]+1
  if(n_cp[3]>0) l2 = (n_cp[1]+n_cp[2])*n2+n_cp[1]+n_cp[2]+1
  if(n_cp[4]>0) l3 = (n_cp[1]+n_cp[2]+n_cp[3])*n2+n_cp[1]+n_cp[2]+n_cp[3]+1
  
  r0 = (n_cp[1])*n2
  if(n_cp[2]>0) r1 = (n_cp[1]+n_cp[2])*n2
  if(n_cp[3]>0) r2 = (n_cp[1]+n_cp[2]+n_cp[3])*n2  
  if(n_cp[4]>0) r3 = (n_cp[1]+n_cp[2]+n_cp[3]+n_cp[4])*n2 
  
  # detected relations set
  if(l0<r0) set_0 <- seq(l0, r0, n2+1)
  if(n_cp[2]>0) { if(l1<r1) set_1 <- seq(l1, r1, n2+1) }
  if(n_cp[3]>0) { if(l2<r2) set_2 <- seq(l2, r2, n2+1) }
  if(n_cp[4]>0) { if(l3<r3) set_3 <- seq(l3, r3, n2+1) }
  
  lab[set_0] <- 1
  if(n_cp[2]>0) lab[set_1] <- 2
  if(n_cp[3]>0) lab[set_2] <- 3
  if(n_cp[4]>0) lab[set_3] <- 4
  return(lab)
}

########## random matrix ############
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

randsvd <- function(A, k)
{
  n = ncol(A)
  
  # Construct an appropriately random orthogonal matrix Omega.
  Omega = matrix(rnorm(k*n),k,n)
  Omega = qr.Q(qr(Omega))
  
  # Apply the appropriately random orthogonal matrix Omega
  # to every column of A^*, obtaining B.
  B = tcrossprod(Omega,A)
  
  # Construct a QR decomposition of B^*.
  QR = qr(t(B))
  Q = qr.Q(QR)
  R = qr.R(QR)
  
  # Calculate the SVD  U_R S V_R = R
  SVD = svd(R)
  U_R = SVD$u
  V_R = SVD$v
  sigma = SVD$d
  
  # Form U_A = Q U_R.
  U_A = tcrossprod(Q,t(U_R))
  # Form V_A = Omega^* V_R.
  V_A = crossprod(Omega,V_R)
  
  return(list(U = U_A, V = V_A, sigma = sigma))
}

#### Genomic Prediction function ####

GBLUP <- function(X1, X2, y1)
{
  if(ncol(X1)==ncol(X2)) M = ncol(X1)
  n1=nrow(X1)
  K11 = A.mat(X1)
  # Model:y=Zu+e
  dt = data.frame(y1, id=rownames(y1))
  colnames(dt)[1]="y1"
  mod = mmer(y1~1,
             random = ~vs(id, Gu=K11),
             rcov = ~units, data = dt, verbose = FALSE, date.warning = FALSE)
  sigmahat = mod$sigmaVector
  # V
  V = K11+diag(sigmahat[2]/sigmahat[1], n1, n1)
  solV = solve(V)
  # K21
  K21 = tcrossprod(X2,X1) / M
  y2hat = tcrossprod(K21,t(tcrossprod(solV, t(y1))))
  return(y2hat)
}

encGBLUP <- function(X1, X2, y1, K, Me)
{
  if(ncol(X1)==ncol(X2)) M = ncol(X1)
  n1 = nrow(X1)
  K11 = A.mat(X1)
  # Model:y=Zu+e
  dt = data.frame(y1, id=rownames(y1))
  colnames(dt)[1] = "y1"
  mod = mmer(y1~1,
             random = ~vs(id, Gu=K11),
             rcov = ~units, data = dt, verbose = FALSE, date.warning = FALSE)  
  sigmahat = mod$sigmaVector
  V = K11+diag(sigmahat[2]/sigmahat[1], n1, n1)
  solV = solve(V)
  
  # Encryption
  # the number of K
  A = RandomMatrixEncryption(X1, X2, M, K, FALSE) # scale by SNP
  K21hat = tcrossprod(A[[2]], A[[1]]) / K
  
  y2hat = tcrossprod(K21hat,t(tcrossprod(solV, t(y1))))
  return(y2hat)
}

gwasBLUP <- function(X1, X2, y1) # return y2hat
{
  bhat = apply(X1, 2, function(x) cov(x, y1)/var(x))
  y2hat = X2 %*% bhat
  return(y2hat)
}

encgwasBLUP <- function(X1, X2, y1, enc, Me)
{
  bhat = apply(X1, 2, function(x) cov(x, y1)/var(x))
  # Encryption
  # the number of K
  K = ceiling( Me / (1/enc-1) )
  A = RandomVecEncryption(X2, bhat, M, K, TRUE)
  y2hat = A[[1]] %*% t(A[[2]])
  return(y2hat)
}

GBLUP2 <- function(X1, X1_m2, X2, y1)
{
  if(ncol(X1)==ncol(X2)) M = ncol(X1)
  n1=nrow(X1)
  G11 = tcrossprod(X1) / M
  # Model:y=Zu+e
  dt = data.frame(y1, id=rownames(y1))
  colnames(dt)[1]="y1"
  mod = mmer(y1~1,
             random = ~vs(id, Gu=G11),
             rcov = ~units, data = dt, verbose = FALSE, date.warning = FALSE)
  sigmahat = mod$sigmaVector
  # V
  V = G11+diag(sigmahat[2]/sigmahat[1], n1, n1)
  solV = solve(V)
  # K21
  K21 = X2 %*% t(X1_m2) / M
  y2hat = K21 %*% solV %*% as.matrix(y1)
  return(y2hat)
}

#### GS ####

GBLUP_AddEnv <- function(dt, cv, A, E, EA)
{
  dt0=dt # store the original data
  env=length(levels(dt$env))
  n=length(levels(dt$var))
  colnames(dt)[1] = "X1"
  dt[cv, 1] = NA
  # y = miu + Zu_g + Zu_e + Zu_ge + e
  mod1 = mmer(X1~1,
              random = ~vs(var, Gu=A) + vs(env, Gu=E) + vs(env:var, Gu=EA),
              rcov = ~units,
              data = dt, 
              verbose = FALSE, date.warning = FALSE)
  # Breeding value
  BV = data.frame(rep(mod1$Beta$Estimate, n*env) + 
                    rep(mod1$U$`u:var`$X1, env) +
                    rep(mod1$U$`u:env`$X1, each = n) + 
                    mod1$U$`u:env:var`$X1,
                  dt0[order(dt0$env,dt0$var),1],
                  row.names = dt0[order(dt0$env,dt0$var),2])
  return(cor(BV[cv,1],BV[cv,2])^2)
}

GBLUP_AddDomEpi <- function(dt, cv, A, D, AA)
{
  dt0=dt # store the original data
  n=length(levels(dt$var))
  colnames(dt)[1] = "X1"
  dt[cv, 1] = NA
  
  ## GBLUP
  mod1 = mmer(X1~1,
              random = ~vs(var, Gu=A) + vs(var.1, Gu=D) + vs(var.2, Gu=AA) ,
              rcov = ~units,
              data = dt, 
              na.method.Y = "include",                    #impute y with the median value
              tolparinv = 1e-2,
              verbose = FALSE, date.warning = FALSE)
  # Breeding value
  BV = data.frame(rep(mod1$Beta$Estimate, n) + 
                    rep(mod1$U$`u:var`$X1) +
                    rep(mod1$U$`u:var.1`$X1) +
                    rep(mod1$U$`u:var.2`$X1),
                  dt0[order(dt0$var),1], 
                  row.names = dt0[order(dt0$var),2])
  return(cor(BV[cv,1],BV[cv,2])^2)
}

MMIBLUP_AddEnv <- function(dt, cv, fmla, A, E, EA)
{
  dt0=dt # store the original data
  env=length(levels(dt$env))
  n=length(levels(dt$var))
  colnames(dt)[1] = "X1"
  dt[cv, 1] = NA
  
  mod2 = mmer(fmla,
              random = ~vs(var, Gu=A) + vs(env, Gu=E) + vs(env:var, Gu=EA) ,
              rcov = ~units,
              data = dt,
              na.method.Y = "include",                           #impute y with the median value so that $fitted can be used
              verbose = FALSE, date.warning = FALSE)
  # fix = data.frame(as.matrix(cbind(rep(1,n*env), dt[,-c(1:4)])) %*% mod2$Beta$Estimate,
  #                  row.names = rownames(dt))
  fix = data.frame(mod2$fitted,row.names = rownames(dt))
  BV = data.frame(fix[order(row.names(fix)),]+
                    rep(mod2$U$`u:var`$X1, env) +
                    rep(mod2$U$`u:env`$X1, each = n) +
                    mod2$U$`u:env:var`$X1 ,
                  dt0[order(dt0$env,dt0$var),1],
                  row.names = dt0[order(dt0$env,dt0$var),2])
  return(cor(BV[cv,1],BV[cv,2])^2)
}

MMIBLUP_AddDomEpi <- function(dt, cv, fmla, A, D, AA)
{
  dt0=dt # store the original data
  n=length(levels(dt$var))
  colnames(dt)[1] = "X1"
  dt[cv, 1] = NA
  
  mod2 = mmer(fmla,
              random = ~vs(var, Gu=A) + vs(var.1, Gu=D) + vs(var.2, Gu=AA),
              rcov = ~units,
              data = dt,
              na.method.Y = "include",                      #impute y with the median value so that $fitted can be used
              tolparinv = 1e-2,
              verbose = FALSE, date.warning = FALSE)
  
  # Breeding value
  fix = data.frame(mod2$fitted,row.names = rownames(dt))
  BV = data.frame(fix[order(row.names(fix)),] +
                    mod2$U$`u:var`$X1 +
                    mod2$U$`u:var.1`$X1 +
                    mod2$U$`u:var.2`$X1 ,
                  dt0[order(dt0$var),1],
                  row.names = dt0[order(dt0$var),2])
  return(cor(BV[cv,1],BV[cv,2])^2)
}

GWAS_AddEnv <- function(dt, cv, Ga, A, E)
{
  dt0=dt # store the original data
  env=length(levels(dt$env))
  n=length(levels(dt$var))
  m=ncol(Ga)
  colnames(dt)[1] = "X1"
  dt[cv, 1] = NA
  
  # GWAS
  mod0 <- GWAS(X1~1,
               random = ~vs(var, Gu=A) + vs(env, Gu=E) + vs(ds(env), var, Gu=A),
               rcov = ~units,
               data = dt,
               M = Ga,
               gTerm = "u:var",
               verbose = FALSE, date.warning = FALSE)
  site_qtl <- which(mod0$scores>-log10(0.05/m) & mod0$scores != Inf) # bonferroni
  return(site_qtl)
}

GWAS_AddDomEpi <- function(dt, cv, Ga, A, D, AA)
{
  dt0=dt # store the original data
  n=length(levels(dt$var))
  m=ncol(Ga)
  colnames(dt)[1] = "X1"
  dt[cv, 1] = NA
  
  # GWAS
  mod0 = GWAS(X1~1,
              random = ~vs(var, Gu=A) + vs(var.1, Gu=D) + vs(var.2, Gu=AA),
              rcov = ~units,
              data = dt,
              M = Ga,
              gTerm = "u:var",
              tolparinv = 1e-2,
              verbose = FALSE, date.warning = FALSE)
  site_qtl <- which(mod0$scores>-log10(0.05/m) & mod0$scores != Inf) # bonferroni
  return(site_qtl)
}

GWASgxe_AddEnv <- function(dt, cv, Ga, A, E)
{  
  dt0=dt # store the original data
  env=length(levels(dt$env))
  n=length(levels(dt$var))
  m=ncol(Ga)
  colnames(dt)[1] = "X1"
  dt[cv, 1] = NA
  # GWASgxe
  site_env_qtl <- c()
  for(e in 1:3)
  {
    gterm = paste0(unique(y$env)[e],":var")
    mod0 <- GWAS(X1~1,
                 random = ~ vs(var, Gu=A) + vs(env, Gu=E) + vs(ds(env), var, Gu=A),
                 rcov = ~ units,
                 data = dt,
                 M = Ga,
                 gTerm = gterm,
                 verbose = FALSE, date.warning = FALSE)
    site_env_qtl <- c(site_env_qtl,which(mod0$scores>-log10(0.05/m) & mod0$scores != Inf)) # bonferroni
  }
  site_env_qtl <- unique(site_env_qtl)
  return(site_env_qtl)
}

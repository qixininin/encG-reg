library(ggplot2)
library(ggpointdensity)
library(egg)
library(viridis)

arg = commandArgs(T)
user1 = arg[1]
user2 = arg[2]
ref = "1KG-EUR"

# input
n1 = as.numeric(read.table(paste0(user1, ".n")))
n2 = as.numeric(read.table(paste0(user2, ".n")))
X0 = read.table(paste0(ref  , ".frq"), header = T)
X1 = read.table(paste0(user1, ".frq"), header = T)
X2 = read.table(paste0(user2, ".frq"), header = T)
rownames(X0) = X0$SNP
rownames(X1) = X1$SNP
rownames(X2) = X2$SNP

# take intersection
inter = intersect(X1$SNP, X2$SNP)
inter = intersect(inter, X0$SNP)
length(inter)

# line up with intersected SNPs
X0 = X0[inter,]
X1 = X1[inter,]
X2 = X2[inter,]

# remove flip
rs=cbind(as.character(X0$SNP),as.character(X1$SNP))
length(which(rs[,1]!=rs[,2]))
twoA1A2 = cbind(X0[,3:4],X1[,3:4])
rm1=which(twoA1A2[,1]!=twoA1A2[,3] & twoA1A2[,1]!=twoA1A2[,4])
rm2=which(twoA1A2[,2]!=twoA1A2[,3] & twoA1A2[,2]!=twoA1A2[,4])

rs=cbind(as.character(X0$SNP),as.character(X2$SNP))
length(which(rs[,1]!=rs[,2]))
twoA1A2 = cbind(X0[,3:4],X2[,3:4])
rm3=which(twoA1A2[,1]!=twoA1A2[,3] & twoA1A2[,1]!=twoA1A2[,4])
rm4=which(twoA1A2[,2]!=twoA1A2[,3] & twoA1A2[,2]!=twoA1A2[,4])
if(length(c(rm1,rm2,rm3,rm4))>0){
  rm = unique(rm1,rm2,rm3,rm4)
}else{
  rm = c()
}
if(length(rm)>1) inter = inter[-rm]
if(length(rm)>1) X1 = X1[-rm,]
if(length(rm)>1) X2 = X2[-rm,]

# generate data
dt = cbind(X1,X2)
M=nrow(dt)
# Change major and minor
changeA1 = which(dt[,3]!=dt[,9])
length(changeA1)
dtnew = dt[,c(5,11)]
dtnew[changeA1,2] = 1-dtnew[changeA1,2]
colnames(dtnew) = c("x","y")

# Density plot
png(paste0(user1,"-",user2,".maf.png"), width = 500, height = 1000)
ggplot(dtnew, aes(x=x,y=y))+
  geom_pointdensity(adjust=0.1) +
  scale_color_viridis()+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,0.5))+
  labs(x=paste0("Frequency in ",user1), y= paste0("Frequency in ",user2), color="density")+
  theme_bw()+
  theme(legend.position = "none", 
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank())
dev.off()

# M
alpha = 0.05 # type I error
beta = 0.1           # type II error two tail test
seta = 0.5*0.9       # 1-degree relatives
minM = ( ( qnorm(1-beta)*sqrt(1+seta^2) + qnorm(1-alpha/n1/n2) ) / seta )^2
M = ceiling(1.2 * minM)

# SNP list
interM = sample(inter, M, replace = F)
write.table(interM, file = paste0(user1,"-",user2,".snp"), quote = F, col.names = F, row.names = F)
write.table(X0[interM,c(2,3)], file = paste0(user1,"-",user2,".snpA1"), quote = F, col.names = F, row.names = F, sep = "\t")

# Me
plink = system("which plink",intern=T) # path to plink
system(paste0(plink," --bfile ",ref," --extract ",user1,"-",user2,".snp --make-bed --out ",ref,".", user1,"-",user2))
system(paste0(plink," --bfile ",ref,".", user1,"-",user2, " --make-grm-gz --out ",ref,".", user1,"-",user2))
gz = gzfile(paste0(ref,".",user1,"-",user2,".grm.gz"))
grm = read.table(gz, as.is = T)
offDiag = grm[grm[,1]!=grm[,2], 4]
Me = 1/var(offDiag, na.rm = TRUE)
Me = ceiling(Me)
write.table(Me, paste0(user1,"-",user2,".me"), quote = F, col.names = F, row.names = F)

# K
minK = ceiling(1/(1/(( ( qnorm(1-beta)*sqrt(1-seta^2) + qnorm(1-alpha/(n1*n2)) ) / seta )^2)-1/Me))
K = ceiling(minK)
write.table(K, paste0(user1,"-",user2,".k"), quote = F, col.names = F, row.names = F)

# seed
seed = n1 * n2
write.table(seed, paste0(user1,"-",user2,".seed"), quote = F, col.names = F, row.names = F)

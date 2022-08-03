# This R file is used to identify 1-degree relatives and identical samples based on .sscore results
# example input:
# Rscript encG-reg-step4.R ${user1} ${user2} 

# input parameter
arg = commandArgs(T)
user1=arg[1]
user2=arg[2]
coflag=paste0(user1,"-",user2)
ctl = read.table("1KG-EUR.txt")
Me = as.numeric(read.table(paste0(coflag,".me")))
encX1 = read.table(paste0(user1, ".", coflag, ".merged.sscore"))
encX2 = read.table(paste0(user2, ".", coflag, ".merged.sscore"))
n1 = nrow(encX1)
n2 = nrow(encX2)
name1 = as.character(encX1[,2])
name2 = as.character(encX2[,2])
count1 = encX1[,3]
count2 = encX2[,3]

# K0 and K1
alpha = 0.05
beta = 0.1
seta = 1*0.9
K0 = 1/(1/ ( ( ( qnorm(1-beta)*sqrt(1-seta^2) + qnorm(1-alpha/(n1-206)/(n2-206)) ) / seta )^2) -1/Me)
K0 = ceiling(K0) 
seta = 0.5*0.9
K1 = 1/(1/ ( ( ( qnorm(1-beta)*sqrt(1-seta^2) + qnorm(1-alpha/(n1-206)/(n2-206)) ) / seta )^2) -1/Me)
K1 = ceiling(K1) 


# Identical samples
encX1data = encX1[,c(5:(K0+4))]  ## the first K column
encX2data = encX2[,c(5:(K0+4))]
encX1data = t(apply(encX1data, 1, scale))
encX2data = t(apply(encX2data, 1, scale))
rownames(encX1data) = name1
rownames(encX2data) = name2
encG12 = tcrossprod(encX1data, encX2data) / K0

# determine threshold
thrd0 = qnorm(1-alpha/(n1-206)/(n2-206))*sqrt(1/Me+1/K0)
# locate individual id
sign = which(as.vector(encG12)>thrd0)
value = round(as.vector(encG12)[sign], 4)
row = sign %% n1 ; col = sign %/% n1 + 1
col[which(row==0)]=sign[which(row==0)]/n1 ; row[which(row==0)]=n1
# output individual id
id0 = matrix(NA, length(sign), 4)
id0[,1] = name1[row]
id0[,2] = count1[row]
id0[,3] = name2[col]
id0[,4] = count2[col]
colnames(id0) = c(user1,"Count1",user2,"Count2")
id0 = cbind(id0, value) 

####################################################
# first degree relatives
encX1data = encX1[,-c(1:4)]
encX2data = encX2[,-c(1:4)]
encX1data = t(apply(encX1data, 1, scale))
encX2data = t(apply(encX2data, 1, scale))
rownames(encX1data) = name1
rownames(encX2data) = name2
encG12 = tcrossprod(encX1data, encX2data) / K1

# determine threshold
thrd1 = qnorm(1-alpha/(n1-206)/(n2-206))*sqrt(1/Me+1/K1)
# locate individual id
sign = which(as.vector(encG12)>thrd1)
value = round(as.vector(encG12)[sign], 4)
row = sign %% n1 ; col = sign %/% n1 + 1
col[which(row==0)]=sign[which(row==0)]/n1 ; row[which(row==0)]=n1
# output individual id
id01 = matrix(NA, length(sign), 4)
id01[,1] = name1[row]
id01[,2] = count1[row]
id01[,3] = name2[col]
id01[,4] = count2[col]
colnames(id01) = c(user1,"Count1",user2,"Count2")
id01 = cbind(id01, value)

# remove identity pairs in 1-degree
df = data.frame(rbind(id0,id01))
colnames(df) = c(user1,"Count1", user2,"Count2", "value")
id1 = df[!(duplicated(df[c(user1, user2)]) | duplicated(df[c(user1, user2)], fromLast = TRUE)), ]
if( nrow(id0)+nrow(id1) != nrow(id01))  print("WARNING: there are identity pairs not identified by thresold for 1-degree, please check!")


# output list
rst = data.frame(rbind(id0,id1), flag = c(rep("Identital", nrow(id0)), rep("1st-degree", nrow(id1))))
# remove control
if(length(which(rst[,1] %in% ctl[,1]))!=206) print("WARNING: controls are not detected!")
rst = rst[-which(rst[,1] %in% ctl[,1]),] ; dim(rst)
write.table(rst, file = paste0(coflag, ".identified") , quote = F, row.names = F, col.names = T, sep = "\t")

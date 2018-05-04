library(bnlearn)
library(pcalg)
library(gRain)

setwd("C:\\Users\\Minh Tu Tran\\Documents\\test_git\\trunk")
data <- read.csv("BRCA_RNASeqv2_top50.csv")

dim(data)
sum(is.na(data))
str(data)
#remove class data
data_remove_class <- subset(data, select = -c(class) )

#1. Use a causal structure learning algorithm to find 
#the skeleton of the gene regulatory
#network using the gene expression data.
V <- colnames(data_remove_class) # labels aka node genes

n <- nrow    (data_remove_class)

## estimate Skeleton
skel.fit <- pcalg::skeleton(suffStat = list(C = cor(data_remove_class), n = n),
                     indepTest = gaussCItest,p=ncol(data_remove_class),
                     alpha = 0.01)
if (require(Rgraphviz)) {
  ## show estimated Skeleton
  par(mfrow=c(1,1))
  plot(skel.fit, main = "Estimated Skeleton")
}


#names

## estimate Skeleton
skel.fit <- pcalg::skeleton(suffStat = list(C = cor(data_remove_class), n = n),
                            indepTest = gaussCItest,labels = V,
                            alpha = 0.01)
if (require(Rgraphviz)) {
  ## show estimated Skeleton
  par(mfrow=c(1,1))
  plot(skel.fit, main = "Estimated Skeleton")
}
#2. Find the top 10 other genes that have strong causal effects on ABCA9 using a
#causal inference algorithm.
#find the index of ABCA9
grep("ABCA9", colnames(data_remove_class))

suffStat <- list(C = cor(data_remove_class), n = nrow(data_remove_class))
#get cpdag
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01,labels = V)
plot(pc.fit@graph)

#genes names
genes <- rownames(cov(data_remove_class))
genes=genes[-22]

#get effects of other nodes on the node 22
effects <- vector()
for (index in 1:50){
  if(index !=22){
   effect = ida(index, 22, cov(data_remove_class), pc.fit@graph)
   
   #the effect is the minimum of the absolute possible effects
   effects <- c(effects, min(abs(effect)))
  }
}

gene_effect = data.frame(genes,effects)

#sort effect from max to min
gene_effect_sort <- gene_effect[order(-effects),]
gene_effect_sort

#get top 10
head(gene_effect_sort,10)

#Use a local causal structure learning algorithm to find genes in the Markov blanket of
#ABCA9 from data.

MB.ABCA9=learn.mb(data_remove_class, 'ABCA9', method='iamb', alpha=0.01)
MB.ABCA9

#convert
data_new = data_remove_class
theMean = mean(as.matrix(data_remove_class))
for(i in 1:ncol(data_remove_class)){
  data_new[[i]] <- ifelse(data_remove_class[[i]] > theMean, 1, 0)
}
head(data_new)
data_new$class = data$class
str(data_new)
write.csv(data_new,"newData.csv")

#Use PC-simple algorithm (pcSelect) to find the parent and children set of the class
#variable.
#get class index
grep("class", colnames(data_new))
#convert to numeric to run pcSelect
data_new$class <- as.numeric(data_new$class)
#pcSelect
pcS <- pcSelect(data_new[,51], data_new[,-51], alpha=0.01)
pcS

#revert back to NC
data_new$class <- data$class



#Construct the conditional probability tables for the Bayesian network based on
#data.

str(data_new)

hl <- c("high","low")
nc <- c("normal","cancer")

mean(data_new$BTNL9==1)
mean(data_new$BTNL9==0)

sum(data_new$BTNL9==1)
sum(data_new$BTNL9==0)
BTNL <- cptable(~BTNL9, values=c(217,995),levels=hl)

mean(data_new$CD300LG[data_new$BTNL9==1]==1)
mean(data_new$CD300LG[data_new$BTNL9==0]==1)
mean(data_new$CD300LG[data_new$BTNL9==1]==0)
mean(data_new$CD300LG[data_new$BTNL9==0]==0)

sum(data_new$CD300LG[data_new$BTNL9==1]==1)
sum(data_new$CD300LG[data_new$BTNL9==0]==1)
sum(data_new$CD300LG[data_new$BTNL9==1]==0)
sum(data_new$CD300LG[data_new$BTNL9==0]==0)
BTNL_CD <- cptable(~CD300LG|BTNL9, values=c(139,78,8,987),levels=hl)

mean(data_new$class[data_new$CD300LG==1]=='C')
mean(data_new$class[data_new$CD300LG==0]=='C')
mean(data_new$class[data_new$CD300LG==1]=='N')
mean(data_new$class[data_new$CD300LG==0]=='N')

sum(data_new$class[data_new$CD300LG==1]=='C')
sum(data_new$class[data_new$CD300LG==0]=='C')
sum(data_new$class[data_new$CD300LG==1]=='N')
sum(data_new$class[data_new$CD300LG==0]=='N')

CD_CLASS <- cptable(~class|CD300LG, values=c(109,38,3,1062),levels=nc)

mean(data_new$IGSF10[data_new$class=='C']==1)
mean(data_new$IGSF10[data_new$class=='N']==1)
mean(data_new$IGSF10[data_new$class=='C']==0)
mean(data_new$IGSF10[data_new$class=='N']==0)

sum(data_new$IGSF10[data_new$class=='C']==1)
sum(data_new$IGSF10[data_new$class=='N']==1)
sum(data_new$IGSF10[data_new$class=='C']==0)
sum(data_new$IGSF10[data_new$class=='N']==0)
CLASS_IG <- cptable(~IGSF10|class, values=c(98,14,71,1029),levels=hl)

mean(data_new$ABCA9[data_new$BTNL9==1 & data_new$IGSF10==1]==1)
mean(data_new$ABCA9[data_new$BTNL9==1 & data_new$IGSF10==0]==1)
mean(data_new$ABCA9[data_new$BTNL9==0 & data_new$IGSF10==1]==1)
mean(data_new$ABCA9[data_new$BTNL9==0 & data_new$IGSF10==0]==1)
mean(data_new$ABCA9[data_new$BTNL9==1 & data_new$IGSF10==1]==0)
mean(data_new$ABCA9[data_new$BTNL9==1 & data_new$IGSF10==0]==0)
mean(data_new$ABCA9[data_new$BTNL9==0 & data_new$IGSF10==1]==0)
mean(data_new$ABCA9[data_new$BTNL9==0 & data_new$IGSF10==0]==0)

sum(data_new$ABCA9[data_new$BTNL9==1 & data_new$IGSF10==1]==1)
sum(data_new$ABCA9[data_new$BTNL9==1 & data_new$IGSF10==0]==1)
sum(data_new$ABCA9[data_new$BTNL9==0 & data_new$IGSF10==1]==1)
sum(data_new$ABCA9[data_new$BTNL9==0 & data_new$IGSF10==0]==1)
sum(data_new$ABCA9[data_new$BTNL9==1 & data_new$IGSF10==1]==0)
sum(data_new$ABCA9[data_new$BTNL9==1 & data_new$IGSF10==0]==0)
sum(data_new$ABCA9[data_new$BTNL9==0 & data_new$IGSF10==1]==0)
sum(data_new$ABCA9[data_new$BTNL9==0 & data_new$IGSF10==0]==0)

ABCA <- cptable(~ABCA9|BTNL9:IGSF10,values=c(114,4,11,40,
                                             26,73,16,928),levels=hl)

table(data_new$ABCA9,data_new$IGSF10)
plist <- compileCPT(list(BTNL,BTNL_CD,CD_CLASS,CLASS_IG,ABCA))
plist

plist$BTNL9
plist$CD300LG
plist$class
plist$IGSF10
plist$ABCA9

net1=grain(plist)
plot(net1)

#Estimate the probability of the four genes in the network having high expression
#levels.

querygrain(net1, nodes=c("BTNL9","CD300LG","IGSF10","ABCA9"), type="join")

querygrain(net1, nodes=c("class","CD300LG","BTNL9"), type="conditional")

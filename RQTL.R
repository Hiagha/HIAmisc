rm(list = ls())


##Data to run simulation on - IBM

Data <- read.csv("C:/Users/hia5y7/Desktop/R/EMP_genos.csv",stringsAsFactors=F,na.strings=c("-","NA"))
##Data1 <- read.csv("C:/Users/hia5y7/Desktop/R/EMP_Genotypes.csv",stringsAsFactors=F,na.strings=c("-","NA"))
#Data <- cbind(Data, Data1[,3:ncol(Data1)])
#Data[,1] <- paste("snp", 1:7333)

Data2 <- as.matrix(Data)
Data3 <- Data2[,1:5]
 
Data2[which(Data2 == "A")] <- 0
Data2[which(Data2 == "B")] <- 1
Data2[which(is.na(Data2) == T)] <- 0.5

Data4 <- cbind(Data3,Data2[,6:ncol(Data2)])
Data5 <- as.data.frame(Data4,stringsAsFactors=F)
Data5[,6:ncol(Data5)] <- apply(Data5[,6:ncol(Data5)],2,as.numeric)

### Run Data5 through simulations - obtain simulation data

setwd("C:/Users/hia5y7/Desktop/R")



### read matrix of all simulated data
allSims <- read.table(file = "Simulations/42_Add_QTN_4_Epi_QTN_h.2_0.8_add.eff_0.6_epis.eff_0.6_reps_100/Simulated.Data.100.Reps.Herit.0.8.txt",header=T,stringsAsFactors=F)
                                
### read genotypes
genos <- read.csv("IBM_Data/IBM_sample_final1.csv",header=T,stringsAsFactors=F,na.strings=c("-","NA"))


### Change directory to file in which you want data

setwd("C:/Users/hia5y7/Desktop/R/Simulations/42_Add_QTN_4_Epi_QTN_h.2_0.8_add.eff_0.6_epis.eff_0.6_reps_100/")

### Loop over analysis

library(qtl)
sigs <- read.table(file = "Genotypic.information.for.4.Additive.QTN.txt",stringsAsFactors=F,header=T)[,1]
twosigs <- read.table(file = "Genotypic.information.for.4.Epistatic.QTN.txt",stringsAsFactors=F,sep = "\t", header = T)[,1]
results <- matrix(NA,nrow=100,ncol=5)
colnames(results) <- c("threshold", sigs)
tworesults <- matrix(NA, nrow = 9, ncol = 100)
rownames(tworesults) <- c("threshold",twosigs)
for(i in 2:ncol(allSims)){
  genos.tmp <- genos
  genos.tmp[,1] <- c("","",allSims[,i])
  write.csv(genos.tmp,file=paste(i-1,".genos.tmp.csv"),row.names=F,quote=F)
  
  
  sim <- read.cross(format = "csv", file=paste(i-1,".genos.tmp.csv"), na.strings=c("-","NA"), genotypes = c("A","B"))
  sim2 <- jittermap(sim, 1e-6)
  
  sim3 <- convert2riself(sim2)
  
  sim4 <- calc.genoprob(sim3)
  
  CIM <- cim(sim4, n.marcov = 6, window = 10, method = "hk")
  
  CIMperms <- cim(sim4, n.marcov = 6, window= 10, method = "hk", n.perm = 100)
  
  thresh <- summary(CIMperms)[1,1]
  
  lod <- as.matrix(CIM$lod)
  rownames(lod) <- colnames(genos)[2:1340]
  first<- c(thresh,lod[rownames = sigs,])
  results[i-1,1:length(first)] <- first
}
  ###Scantwo search
  
  scantwo <- scantwo(sim4, chr = c(1,2,3,4,5,6,7,8,9,10), pheno.col = 1, model = "normal", method = "hk", addcovar = NULL, intcovar = NULL, weights = NULL, use = "all.obs")
  
  print(i) 
  
  scantwoperms <- scantwo(sim4, chr = c(1,2,3,4,5,6,7,8,9,10), pheno.col = 1, model = "normal", method = "hk", addcovar = NULL, intcovar = NULL, weights = NULL, use = "all.obs", n.perm = 100)
  
  scantwosum <- summary(scantwo, perm = scantwoperms, alpha = 0.05)

  twothresh <- summary(scantwoperms)[[1]][1]
  
  foo <- matrix(scantwo$lod, nrow = 1339, ncol = 1339)
  rownames(foo) <- rownames(CIM)
  colnames(foo) <- rownames(CIM)

  int1x2 <- foo[rownames = twosigs[1],colnames = twosigs[2]]
  int2x1 <- foo[rownames = twosigs[2],colnames = twosigs[1]]
  int3x4 <- foo[rownames = twosigs[3],colnames = twosigs[4]]
  int4x3 <- foo[rownames = twosigs[4],colnames = twosigs[3]]
  int5x6 <- foo[rownames = twosigs[5],colnames = twosigs[6]]
  int6x5 <- foo[rownames = twosigs[6],colnames = twosigs[5]]
  int7x8 <- foo[rownames = twosigs[7],colnames = twosigs[8]]
  int8x7 <- foo[rownames = twosigs[8],colnames = twosigs[7]]
  twocols <- rbind(twothresh, int1x2, int2x1, int3x4, int4x3, int5x6, int6x5, int7x8, int8x7) 
  
  tworesults[,i-1] <- twocols
  

  
  write.csv(tworesults, file = paste("results.epistatic.0.9.csv"))

  write.csv(results, file = paste("results.additive.0.9.csv"))


epistatic <- read.table("Simulated.Data.100.Reps.Herit.0.9.txt")
plot(epistatic[2:nrow(epistatic),1:2])
epistatic2 <- as.matrix(epistatic[2:211,2:101])
plot(epistatic2[1:ncol(QTN2),1], ylab = "Phenotype", xlab = "Genotype", main = "Simulated Phenotypes", col = col,pch=19, cex=2,ylim=c(0,0.75))


QTN <- read.csv(file = "../epistatic.qtn.genos.csv", stringsAsFactors = F)
QTN2 <- QTN[,6:ncol(QTN)]



for(i in 1:100){
  plot(epistatic2[,i])
  par(new = T)
}




col <- matrix("black",nrow = ncol(QTN2), ncol = 1)
  col[which((QTN2)[1,1:ncol(QTN2)] > 0 & (QTN2)[2,1:ncol(QTN2)] > 0)] <- "purple"

pch <- matrix(1,nrow = ncol(QTN2), ncol = 1)
  pch[which((QTN2)[1,1:ncol(QTN2)] > 0 & (QTN2)[2,1:ncol(QTN2)] > 0)] <- 0


cex <- matrix(1,nrow = ncol(QTN2), ncol = 1)
  cex[which((QTN2)[5,1:ncol(QTN2)] > 0 & (QTN2)[6,1:ncol(QTN2)] > 0)] <- 2
  
points <- matrix(NA,nrow = ncol(QTN2), ncol = 1)
for(i in 1:ncol(QTN2)){
  if((QTN2)[7,i] > 0 & (QTN2)[8,i] > 0){
    points[i] <- print(epistatic2[i,1])
  }
}
  
plot(epistatic2[1:ncol(QTN2),1], ylab = "Phenotype", xlab = "Genotype", 
     main = "Simulated Phenotypes", col = col, pch = pch, cex = cex,ylim=c(0,0.75))
  points(y = points, x = 1:ncol(QTN2), pch = 4, cex = 2)

legend("topright","c(x,y)",c("inter1","inter2","inter3","inter4"),col=c("purple","black","black","black"),pt.cex=c(1,1,2,1),pch=c(1,0,1,13), cex = 1.25)







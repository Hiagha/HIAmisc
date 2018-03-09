setwd("C:/Users/hia5y7/Desktop/R")

for(k in 1:9){

### read matrix of all simulated data
allSims <- read.table(file = paste("simulation_analyses_her0.",k,"/4_Add_QTN4_Epi_QTN_h.2_0.",k,"_add.eff_0.6_epis.eff_0.6_reps_100/Simulated.Data.100.Reps.Herit.0.",k,".txt"),header=T,stringsAsFactors=F)

### read genotypes
genos <- read.csv("IBM_Data/IBM_sample_final1.csv",header=T,stringsAsFactors=F,na.strings=c("-","NA"))


### Change directory to file in which you want data
setwd("C:/Users/hia5y7/Desktop/R/simulation_analyses_her0.9")

### Loop over analysis
library(qtl)

sigs <- read.table(file = "../simulation_analyses_her0.9/4_Add_QTN4_Epi_QTN_h.2_0.9_add.eff_0.6_epis.eff_0.6_reps_100/Genotypic.information.for.4.Additive.QTN.txt",stringsAsFactors=F,header=T)[,1:5]
sigs <- sigs[,1]
twosigs <- read.table(file = "significant.markers.epistatic.csv",stringsAsFactors=F,sep = ",", header = T)[,1]
results <- matrix(NA,nrow=100,ncol=4)
tworesults <- matrix(NA, nrow = 9, ncol = 100)
rownames(tworesults) <- c("threshold",twosigs)
for(i in 2:ncol(allSims)){
  genos.tmp <- genos
  genos.tmp[,1] <- c("","",allSims[,i])
  write.csv(genos.tmp,file=paste(i,".genos.tmp.csv"),row.names=F,quote=F)
  
  
  sim <- read.cross(format = "csv", ".", 
                    file=paste(i,".genos.tmp.csv"), na.strings=c("-","NA"), genotypes = c("A","B"))
  sim2 <- jittermap(sim, 1e-6)
  
  sim3 <- convert2riself(sim2)
  
  sim4 <- calc.genoprob(sim3)
  
  CIM <- cim(sim4, n.marcov = 6, window = 10, method = "hk")
  
  CIMperms <- cim(sim4, n.marcov = 6, window= 10, method = "hk", n.perm = 100)
  
  thresh <- summary(CIMperms)[1,1]
  
  first<- sigs[which(sigs %in% rownames(CIM[which(CIM[,3] >= thresh),]))]
  
  results[i-1,1:length(first)] <- first
  
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
  
}

write.csv(tworesults, file = paste("results.epistatic.0.9.csv"))

write.csv(results, file = paste("results.additive.0.9.csv"))

}
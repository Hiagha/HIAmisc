##setwd() to folder in which simulated data is located
setwd("../EMPandAdditive4_Add_QTN_4_Epi_QTN_h.2_0.9_add.eff_0.5_epis.eff_0.6_reps_100/")

##Read in and create a matrix of the simulated phenotypic data 
epistatic <- read.table("Simulated.Data.100.Reps.Herit.0.9.txt")
epistatic2 <- as.matrix(epistatic[2:nrow(epistatic),2:ncol(epistatic)])

##Read in genotypic data for simulated QTN
QTN <- read.table("Genotypic.information.for.4.Epistatic.QTN.txt", stringsAsFactors = F, header = T, row.names = NULL)
QTN2 <- QTN[,7:ncol(QTN)]

##Create matrix that identifies genotypes that carry first interacting QTN - Color
col <- matrix("black",nrow = ncol(QTN2), ncol = 1)
col[which((QTN2)[1,1:ncol(QTN2)] > 0 & (QTN2)[2,1:ncol(QTN2)])] <- "purple"

##Create matrix that identifies genotypes that carry second interacting QTN - Shape
pch <- matrix(1,nrow = ncol(QTN2), ncol = 1)
pch[which((QTN2)[1,1:ncol(QTN2)] > 0 & (QTN2)[2,1:ncol(QTN2)])] <- 0

##Create matrix that identifies genotypes that which carry third interacting QTN - Size
cex <- matrix(1,nrow = ncol(QTN2), ncol = 1)
cex[which((QTN2)[5,1:ncol(QTN2)] > 0 & (QTN2)[6,1:ncol(QTN2)])] <- 2

##Create matrix that identifies genotypes that carry fourth interacting QTN - Adds X over point
points <- matrix(NA,nrow = ncol(QTN2), ncol = 1)
for(i in 1:ncol(QTN2)){
  if(QTN2[7,i] > 0 & QTN2[8,i] > 0){
    points[i] <- print(epistatic2[i,1])}}

##Plot all phenotype points for genotypes - points with interactions present distinguished
plot(epistatic2[1:nrow(epistatic2),1], ylab = "Phenotype", xlab = "Genotype", 
     main = "Simulated Phenotypes", col = col, pch = pch, cex = cex)
points(y = points, x = 1:ncol(QTN2), pch = 4, cex = 2)

##Legend - if more than one interaction for a single genotype, plot attributes stack
legend("topright","c(x,y)",c("inter1","inter2","inter3","inter4"),
       col=c("purple","black","black","black"),pt.cex=c(1,1,2,1),pch=c(1,0,1,13), cex = 1.25)

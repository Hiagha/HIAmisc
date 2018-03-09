nils <- read.csv("C:/Users/hia5y7/Desktop/R/NILS.csv", stringsAsFactors = F, na.strings = c("N"))
nils <- nils[-which(nils[,4] == "C"),]
emp_genos <- matrix(NA, ncol = 2, nrow = nrow(nils))
emp_genos <- nils[,1:2]
nils <- as.matrix(nils)[,3:22]
nils[which(nils == "B")] <- 0
nils[which(nils == "M")] <- 1
nils[which(nils == "H")] <- 1
nils2 <- matrix(as.numeric(nils),ncol = 20)
nils2<-data.frame(nils2)
names(nils2) <- colnames(nils)

for(i in 1:20){
  for(j in 1:7313){
    if((is.na(nils2[j,i])) == T){
      for(k in 1:200){
        if(is.na(nils2[j+k,i]) == F){break}}
      nils2[j,i] <- ((nils2[j-1,i] + nils2[j+k,i])/2)
      for(l in 1:200){
        if(is.na(nils2[j+l,i]) == T){
          nils2[j+l, i] <- (nils2[j,i])
        }else{break}}
    }
  }  
  print(i)  }

nils2[(nils2 == 0.5)] <- 1

###This is an idealized EMP, all heterogenous locations on the DNA have been replaced with introgressions (MO17 DNA)
###Only to be used to run through R/qtl, not a believable genotype



for(i in 1:19){
  for(j in {i+1}:20){
    cols <- (nils2[,i] + nils2[,j])/2
    emp_genos <- cbind(emp_genos, cols)
    colnames(emp_genos)[ncol(emp_genos)] <- paste(colnames(nils2)[i], colnames(nils2)[j], sep = "x")
  }
}

for(i in 1:20){
  cols <- (nils2[,i])/2
  emp_genos <- cbind(emp_genos, cols)
  colnames(emp_genos)[ncol(emp_genos)] <- paste("B73",colnames(nils2)[i], sep = "x")
}
emp_genos2 <- t(emp_genos)
rownames(emp_genos2)[1:2] <- ""
colnames(emp_genos2) <- paste("snp", 1:7313, sep = " ")



write.csv(emp_genos2, file = "Idealized_EMP_genos.csv")

emp_genos1 <- emp_genos[,3:ncol(emp_genos)]

for(k in 1:210){
for(i in 1:7313){
both <- colnames(emp_genos1[which(emp_genos1[i,k] == 1)])
}
}
both <- matrix(NA, ncol = 2, nrow = 210)
for(i in 1:210){
  both[i,] <- range(emp_genos1[,i])
}
rownames(both) <- colnames(emp_genos1)
inters <- rownames(both[which(both[,2] != 1),])

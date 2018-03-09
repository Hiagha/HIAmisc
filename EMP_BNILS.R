nils <- read.csv("C:/Users/hia5y7/Desktop/R/NILS.csv", stringsAsFactors = F, na.strings = c("N", "C"))
emp_genos <- matrix(NA, ncol = 2, nrow = nrow(nils))
emp_genos <- nils[,1:2]
nils <- as.matrix(nils)[,3:22]
nils[which(nils == "B")] <- 0
nils[which(nils == "M")] <- 1
nils[which(nils == "H")] <- .5
nils2 <- matrix(as.numeric(nils),ncol = 20)
nils2<-data.frame(nils2)
names(nils2) <- colnames(nils)



for(i in 1:20){
    cols <- (nils2[,i])/2
    emp_genos <- cbind(emp_genos, cols)
    colnames(emp_genos)[ncol(emp_genos)] <- paste("B73",colnames(nils2)[i], sep = "x")
  }

write.csv(emp_genos, file = "C:/Users/hia5y7/Desktop/R/EMP_BNILS.csv")


dates <- read.csv("../Downloads/2017-08-08-01-02-54_mockEMP_table.csv", header = T, na.strings = "NA")
dates2 <- dates[,c(2,4:5)]
dates2 <- cbind(dates2, matrix(NA,nrow = nrow(dates), ncol = 1))
colnames(dates2) <- c("Genotype","Anthesis","Silk","Planting Date")
dates2[1:440,4]<- "5/17/2017"
dates2[441:880,4] <- "5/16/2017"
Anthesis <- as.Date(dates2[,2])
Silk <- as.Date(dates2[,3])
Planting <- as.Date(dates2[,4], format = "%m/%d/%Y" )
DaystoAnthesis <- Anthesis - Planting
DaystoSilk <- Silk - Planting
final <- matrix(c(dates$genotype,DaystoAnthesis,DaystoSilk), ncol = 3, nrow = nrow(dates2))
colnames(final) <- c("Days to Anthesis", "Days to Silk")
#write.csv(final, "../Desktop/R/maturity_data_emp.csv")

rm(list = ls())


setwd("/Users/hiagh/OneDrive/Documents/GitHub/sim_analyses/results/")

herit01 <- read.csv("results.epistatic.0.1.csv", header = T, row.names = 1)
foo01 <- (herit01[1,])

for(i in 2:800){
  if(!i %% 2){
  foo01 <- rbind(foo01, herit01[i,])
  
    }
}
rownames <- matrix(NA, nrow = 400, ncol = 1)
for(i in 2:800){
  if(!i %% 2){
  rownames[i/2,] <- paste(rownames(herit01[i,]), rownames(herit01[i+1,]), sep = "x")
}}
rownames(foo01) <- c("threshold", rownames)
results <- matrix(NA, nrow = 401, ncol = 100)
for(i in 1:100){
  for(j in 2:401){
    if(foo01[j,i]>foo01[1,i]){results[j,i] <- print(foo01[j,i])}
  }}

full02 <- read.csv("results.epistatic.0.2.csv")[c(1,2,4,6,8), 2:101]
fullmat02 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(full02[j,i]>full02[1,i]){fullmat02[j,i] <- print(full02[j,i])}
  }}


full03 <- read.csv("results.epistatic.0.3.csv")[c(1,2,4,6,8), 2:101]
fullmat03 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(full03[j,i]>full03[1,i]){fullmat03[j,i] <- print(full03[j,i])}
  }}

full04 <- read.csv("results.epistatic.0.4.csv")[c(1,2,4,6,8), 2:101]
fullmat04 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(full04[j,i]>full04[1,i]){fullmat04[j,i] <- print(full04[j,i])}
  }}

full05 <- read.csv("results.epistatic.0.5.csv")[c(1,2,4,6,8), 2:101]
fullmat05 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(full05[j,i]>full05[1,i]){fullmat05[j,i] <- print(full05[j,i])}
  }}

full06 <- read.csv("results.epistatic.0.6.csv")[c(1,2,4,6,8), 2:101]
fullmat06 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(full06[j,i]>full06[1,i]){fullmat06[j,i] <- print(full06[j,i])}
  }}

full07 <- read.csv("results.epistatic.0.7.csv")[c(1,2,4,6,8), 2:101]
fullmat07 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(full07[j,i]>full07[1,i]){fullmat07[j,i] <- print(full07[j,i])}
  }}

full08 <- read.csv("results.epistatic.0.8.csv")[c(1,2,4,6,8), 2:101]
fullmat08 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(full08[j,i]>full08[1,i]){fullmat08[j,i] <- print(full08[j,i])}
  }}
full09 <- read.csv("results.epistatic.0.9.csv")[c(1,2,4,6,8), 2:101]
fullmat09 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(full09[j,i]>full09[1,i]){fullmat09[j,i] <- print(full09[j,i])}
  }}

f01 <- length(which(fullmat01[5,] > 1)) + length(which(fullmat01[4,] > 1)) + length(which(fullmat01[3,] > 1)) + length(which(fullmat01[2,] > 1))
f02 <- length(which(fullmat02[5,] > 1)) + length(which(fullmat02[4,] > 1)) + length(which(fullmat02[3,] > 1)) + length(which(fullmat02[2,] > 1))
f03 <- length(which(fullmat03[5,] > 1)) + length(which(fullmat03[4,] > 1)) + length(which(fullmat03[3,] > 1)) + length(which(fullmat03[2,] > 1))
f04 <- length(which(fullmat04[5,] > 1)) + length(which(fullmat04[4,] > 1)) + length(which(fullmat04[3,] > 1)) + length(which(fullmat04[2,] > 1))
f05 <- length(which(fullmat05[5,] > 1)) + length(which(fullmat05[4,] > 1)) + length(which(fullmat05[3,] > 1)) + length(which(fullmat05[2,] > 1))
f06 <- length(which(fullmat06[5,] > 1)) + length(which(fullmat06[4,] > 1)) + length(which(fullmat06[3,] > 1)) + length(which(fullmat06[2,] > 1))
f07 <- length(which(fullmat07[5,] > 1)) + length(which(fullmat07[4,] > 1)) + length(which(fullmat07[3,] > 1)) + length(which(fullmat07[2,] > 1))
f08 <- length(which(fullmat08[5,] > 1)) + length(which(fullmat08[4,] > 1)) + length(which(fullmat08[3,] > 1)) + length(which(fullmat08[2,] > 1))
f09 <- length(which(fullmat09[5,] > 1)) + length(which(fullmat09[4,] > 1)) + length(which(fullmat09[3,] > 1)) + length(which(fullmat09[2,] > 1))

plot(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),c(f01/400,f02/400,f03/400,f04/400,f05/400,f06/400,f07/400,f08/400,f09/400),
     type = "o", main = "Epistatic Full Model", xlab = "Heritability", ylab = "Discovery Rate", ylim = 0:1)


add01 <- read.csv("results.epistatic.0.1.csv")[c(1,2,4,6,8), 2:101]
addmat01 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add01[j,i]>add01[1,i]){addmat01[j,i] <- print(add01[j,i])}
  }}

add02 <- read.csv("results.epistatic.0.2.csv")[c(1,2,4,6,8), 2:101]
addmat02 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add02[j,i]>add02[1,i]){addmat02[j,i] <- print(add02[j,i])}
  }}


add03 <- read.csv("results.epistatic.0.3.csv")[c(1,2,4,6,8), 2:101]
addmat03 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add03[j,i]>add03[1,i]){addmat03[j,i] <- print(add03[j,i])}
  }}

add04 <- read.csv("results.epistatic.0.4.csv")[c(1,2,4,6,8), 2:101]
addmat04 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add04[j,i]>add04[1,i]){addmat04[j,i] <- print(add04[j,i])}
  }}

add05 <- read.csv("results.epistatic.0.5.csv")[c(1,2,4,6,8), 2:101]
addmat05 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add05[j,i]>add05[1,i]){addmat05[j,i] <- print(add05[j,i])}
  }}

add06 <- read.csv("results.epistatic.0.6.csv")[c(1,2,4,6,8), 2:101]
addmat06 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add06[j,i]>add06[1,i]){addmat06[j,i] <- print(add06[j,i])}
  }}

add07 <- read.csv("results.epistatic.0.7.csv")[c(1,2,4,6,8), 2:101]
addmat07 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add07[j,i]>add07[1,i]){addmat07[j,i] <- print(add07[j,i])}
  }}

add08 <- read.csv("results.epistatic.0.8.csv")[c(1,2,4,6,8), 2:101]
addmat08 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add08[j,i]>add08[1,i]){addmat08[j,i] <- print(add08[j,i])}
  }}
add09 <- read.csv("results.epistatic.0.9.csv")[c(1,2,4,6,8), 2:101]
addmat09 <- matrix(NA, nrow = 5, ncol = 100)

for(i in 1:100){
  for(j in 2:5){
    if(add09[j,i]>add09[1,i]){addmat09[j,i] <- print(add09[j,i])}
  }}

a01 <- length(which(addmat01[5,] > 1)) + length(which(addmat01[4,] > 1)) + length(which(addmat01[3,] > 1)) + length(which(addmat01[2,] > 1))
a02 <- length(which(addmat02[5,] > 1)) + length(which(addmat02[4,] > 1)) + length(which(addmat02[3,] > 1)) + length(which(addmat02[2,] > 1))
a03 <- length(which(addmat03[5,] > 1)) + length(which(addmat03[4,] > 1)) + length(which(addmat03[3,] > 1)) + length(which(addmat03[2,] > 1))
a04 <- length(which(addmat04[5,] > 1)) + length(which(addmat04[4,] > 1)) + length(which(addmat04[3,] > 1)) + length(which(addmat04[2,] > 1))
a05 <- length(which(addmat05[5,] > 1)) + length(which(addmat05[4,] > 1)) + length(which(addmat05[3,] > 1)) + length(which(addmat05[2,] > 1))
a06 <- length(which(addmat06[5,] > 1)) + length(which(addmat06[4,] > 1)) + length(which(addmat06[3,] > 1)) + length(which(addmat06[2,] > 1))
a07 <- length(which(addmat07[5,] > 1)) + length(which(addmat07[4,] > 1)) + length(which(addmat07[3,] > 1)) + length(which(addmat07[2,] > 1))
a08 <- length(which(addmat08[5,] > 1)) + length(which(addmat08[4,] > 1)) + length(which(addmat08[3,] > 1)) + length(which(addmat08[2,] > 1))
a09 <- length(which(addmat09[5,] > 1)) + length(which(addmat09[4,] > 1)) + length(which(addmat09[3,] > 1)) + length(which(addmat09[2,] > 1))

plot(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),c(a01/400,a02/400,a03/400,a04/400,a05/400,a06/400,a07/400,a08/400,a09/400), 
     type = "o", main = "Epistatic Additive Model", xlab = "Heritability", ylab = "Discovery Rate", ylim = 0:1)






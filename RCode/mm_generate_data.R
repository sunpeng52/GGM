rm(list = ls())
library(dplyr)
library(spcov)
library(MASS)

library(parallel)
library(foreach)
library(doParallel)

source("functions.R")
source("functions_2.R")
options(dplyr.summarise.inform = FALSE)

numCores <- detectCores()
numCores

registerDoParallel(numCores)

########## generate data sets
generate.data2 <- function(R,G,n,N){
  n_all <- c(rep(n,99),N-99*n)
  da_list = list()
  p = nrow(R)
  for(j in 1:100){
    Yi_all = data.frame()
    m <-  100
    for(i in 1:m){
      ni = n_all[i]
      bi = matrix(mvrnorm(1,rep(0,p),G),ni,p,byrow = T)
      Yi = bi+mvrnorm(ni,rep(0,p),R)
      Yi_df = data.frame(i=rep(i,ni),j=1:ni,Yi)
      Yi_all = rbind(Yi_all,Yi_df)
    }
    da_list[[j]] = Yi_all
  }
  return(da_list)
}

set.seed(1234)
######### model 1, p =100
p <-  100
R.true <- banded.new(p,10,-1)
G.true <- banded.new(p,10)

da_m1_100_03 <- generate.data2(R.true,G.true,3,1000)
da_m1_100_04 <- generate.data2(R.true,G.true,4,1000)
da_m1_100_05 <- generate.data2(R.true,G.true,5,1000)
da_m1_100_06 <- generate.data2(R.true,G.true,6,1000)
da_m1_100_07 <- generate.data2(R.true,G.true,7,1000)
da_m1_100_08 <- generate.data2(R.true,G.true,8,1000)
da_m1_100_09 <- generate.data2(R.true,G.true,9,1000)
da_m1_100_10 <- generate.data2(R.true,G.true,10,1000)

######### model 1, p =200
p <-  200
R.true <- banded.new(p,10,-1)
G.true <- banded.new(p,10)

da_m1_200_03 <- generate.data2(R.true,G.true,3,1000)
da_m1_200_04 <- generate.data2(R.true,G.true,4,1000)
da_m1_200_05 <- generate.data2(R.true,G.true,5,1000)
da_m1_200_06 <- generate.data2(R.true,G.true,6,1000)
da_m1_200_07 <- generate.data2(R.true,G.true,7,1000)
da_m1_200_08 <- generate.data2(R.true,G.true,8,1000)
da_m1_200_09 <- generate.data2(R.true,G.true,9,1000)
da_m1_200_10 <- generate.data2(R.true,G.true,10,1000)

######### model 2, p =100
p <-  100
R.true <- ar.cov(p,-0.6)
G.true <- ar.cov(p,0.6)
min(eigen(R.true)$values)

da_m2_100_03 <- generate.data2(R.true,G.true,3,1000)
da_m2_100_04 <- generate.data2(R.true,G.true,4,1000)
da_m2_100_05 <- generate.data2(R.true,G.true,5,1000)
da_m2_100_06 <- generate.data2(R.true,G.true,6,1000)
da_m2_100_07 <- generate.data2(R.true,G.true,7,1000)
da_m2_100_08 <- generate.data2(R.true,G.true,8,1000)
da_m2_100_09 <- generate.data2(R.true,G.true,9,1000)
da_m2_100_10 <- generate.data2(R.true,G.true,10,1000)

######### model 2, p =200
p <-  200
R.true <- ar.cov(p,-0.6)
G.true <- ar.cov(p,0.6)
min(eigen(R.true)$values)

da_m2_200_03 <- generate.data2(R.true,G.true,3,1000)
da_m2_200_04 <- generate.data2(R.true,G.true,4,1000)
da_m2_200_05 <- generate.data2(R.true,G.true,5,1000)
da_m2_200_06 <- generate.data2(R.true,G.true,6,1000)
da_m2_200_07 <- generate.data2(R.true,G.true,7,1000)
da_m2_200_08 <- generate.data2(R.true,G.true,8,1000)
da_m2_200_09 <- generate.data2(R.true,G.true,9,1000)
da_m2_200_10 <- generate.data2(R.true,G.true,10,1000)

save(
  da_m1_100_03,da_m1_100_04,
  da_m1_100_05,da_m1_100_06,da_m1_100_07,
  da_m1_100_08,da_m1_100_09,da_m1_100_10,
  da_m2_100_03,da_m2_100_04,
  da_m2_100_05,da_m2_100_06,da_m2_100_07,
  da_m2_100_08,da_m2_100_09,da_m2_100_10,
  file = "data_mm_100.RData"
)

save(
  da_m1_200_03,da_m1_200_04,
  da_m1_200_05,da_m1_200_06,da_m1_200_07,
  da_m1_200_08,da_m1_200_09,da_m1_200_10,
  da_m2_200_03,da_m2_200_04,
  da_m2_200_05,da_m2_200_06,da_m2_200_07,
  da_m2_200_08,da_m2_200_09,da_m2_200_10,
  file = "data_mm_200.RData"
)

print(1)




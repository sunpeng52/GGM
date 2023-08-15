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

#numCores <- detectCores()
#numCores

#registerDoParallel(numCores)

########## generate data sets
generate.data <- function(R,G,n,m=50){
  #n_all <- c(rep(n,99),N-99*n)
  da_list = list()
  p = nrow(R)
  for(j in 1:100){
    Yi_all = data.frame()
    for(i in 1:m){
      ni = n[i]
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

#### scenario 1: opposite drection, 50 groups, n_i = 5
p <- 50
R.base <- banded.new(p,10,-1)
min(eigen(R.base)$values)
norm(R.base,"F")
G.base <- banded.new(p,10)
norm(G.base,"F")
min(eigen(G.base)$values)

n <- rep(5,50)

###
da_m3_50_n_00 <- generate.data(diag(1,p),G.base,n)
da_m3_50_n_01 <- generate.data(R.base,G.base,n)
da_m3_50_n_02 <- generate.data(2*R.base,G.base,n)
da_m3_50_n_03 <- generate.data(3*R.base,G.base,n)
da_m3_50_n_04 <- generate.data(4*R.base,G.base,n)
da_m3_50_n_05 <- generate.data(5*R.base,G.base,n)
da_m3_50_n_06 <- generate.data(6*R.base,G.base,n)
da_m3_50_n_07 <- generate.data(7*R.base,G.base,n)
da_m3_50_n_08 <- generate.data(8*R.base,G.base,n)
da_m3_50_n_09 <- generate.data(9*R.base,G.base,n)
da_m3_50_n_10 <- generate.data(10*R.base,G.base,n)


#### scenario 2: same drection, 50 groups, n_i = 5
p <- 50
R.base2 <- banded.new(p,10)
min(eigen(R.base2)$values)
norm(R.base2,"F")
G.base <- banded.new(p,10)
norm(G.base,"F")
min(eigen(G.base)$values)

n <- rep(5,50)
### a = 0
da_m3_50_p_01 <- generate.data(R.base2,G.base,n)
da_m3_50_p_02 <- generate.data(2*R.base2,G.base,n)
da_m3_50_p_03 <- generate.data(3*R.base2,G.base,n)
da_m3_50_p_04 <- generate.data(4*R.base2,G.base,n)
da_m3_50_p_05 <- generate.data(5*R.base2,G.base,n)
da_m3_50_p_06 <- generate.data(6*R.base2,G.base,n)
da_m3_50_p_07 <- generate.data(7*R.base2,G.base,n)
da_m3_50_p_08 <- generate.data(8*R.base2,G.base,n)
da_m3_50_p_09 <- generate.data(9*R.base2,G.base,n)
da_m3_50_p_10 <- generate.data(10*R.base2,G.base,n)


save(
  da_m3_50_n_00, 
  da_m3_50_n_01, da_m3_50_n_02, da_m3_50_n_03, da_m3_50_n_04, da_m3_50_n_05, 
  da_m3_50_n_06, da_m3_50_n_07, da_m3_50_n_08, da_m3_50_n_09, da_m3_50_n_10,
  da_m3_50_p_01, da_m3_50_p_02, da_m3_50_p_03, da_m3_50_p_04, da_m3_50_p_05, 
  da_m3_50_p_06, da_m3_50_p_07, da_m3_50_p_08, da_m3_50_p_09, da_m3_50_p_10,
  file = "data_mm_comparision.RData"
)




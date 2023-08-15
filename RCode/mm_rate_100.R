#setwd("/home/grad/sunpeng_duan/Group_Graphical_Model")
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

load("data_mm_100.RData")

p <-  100
R.true <- banded.new(p,10,-1)
G.true <- banded.new(p,10)

#fpr = function(res,true) sum(which(res != 0) %in% which(true == 0))/length(which(true == 0))
#tpr = function(res,true) sum(which(res != 0) %in% which(true != 0))/length(which(true != 0))


fpr2 = function(res,true) sum((res != 0)*(true == 0))/sum(true == 0)
tpr2 = function(res,true) sum((res != 0)*(true != 0))/sum(true != 0)

del = 10^(-5)
P <- matrix(1,p,p)
diag(P) <- 0

lam_vec = 10^seq(-2,0.5,0.05)
#maxit=500

rate_m1_100_10 <- foreach(j = 1:20, .combine = 'rbind') %dopar%{
  da <- da_m1_100_10[[j]]
  da_group = da %>% group_by(i)
  group_size = da_group %>% dplyr::summarise(N = n())
  group_size = group_size$N
  group_var = da_group %>% do(Var = var(.[,-c(1,2)]))
  group_var = group_var$Var
  group_mean = da_group %>% do(Mean = colMeans(.[,-c(1,2)]))
  group_mean = group_mean$Mean
  grand_mean = colMeans(da[,-c(1,2)])
  
  #### find moment estimation of R and G
  total_obs = Reduce("+",group_size)
  group_num = length(group_size)
  sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                          variance=group_var,size=group_size,SIMPLIFY = FALSE)
  )
  mse = sse/(total_obs - group_num)
  R_hat_n = mse
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_a <- A_hat_n - mean(1/group_size)*R_hat_n
  
  #### nested loop: compute fpr,tpr for a data simulation
  temp <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    
    R.admm = admm.cov(R_hat_n,del,lam,P)$Theta
    G.admm = admm.cov(G_hat_n,del,lam,P)$Theta
    Ga.admm = admm.cov(G_hat_a,del,lam,P)$Theta
    A.admm = admm.cov(A_hat_n,del,lam,P)$Theta
    
    c(j,
      fpr2(R.admm,R.true),tpr2(R.admm,R.true),
      fpr2(A.admm,R.true),tpr2(A.admm,R.true),
      fpr2(G.admm,G.true),tpr2(G.admm,G.true),
      fpr2(Ga.admm,G.true),tpr2(Ga.admm,G.true),
      fpr2(A.admm,G.true),tpr2(A.admm,G.true))
  }
  temp
}

rate_m1_100_09 <- foreach(j = 1:20, .combine = 'rbind') %dopar%{
  da <- da_m1_100_09[[j]]
  da_group = da %>% group_by(i)
  group_size = da_group %>% dplyr::summarise(N = n())
  group_size = group_size$N
  group_var = da_group %>% do(Var = var(.[,-c(1,2)]))
  group_var = group_var$Var
  group_mean = da_group %>% do(Mean = colMeans(.[,-c(1,2)]))
  group_mean = group_mean$Mean
  grand_mean = colMeans(da[,-c(1,2)])
  
  #### find moment estimation of R and G
  total_obs = Reduce("+",group_size)
  group_num = length(group_size)
  sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                          variance=group_var,size=group_size,SIMPLIFY = FALSE)
  )
  mse = sse/(total_obs - group_num)
  R_hat_n = mse
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_a <- A_hat_n - mean(1/group_size)*R_hat_n
  
  #### nested loop: compute fpr,tpr for a data simulation
  temp <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    
    R.admm = admm.cov(R_hat_n,del,lam,P)$Theta
    G.admm = admm.cov(G_hat_n,del,lam,P)$Theta
    Ga.admm = admm.cov(G_hat_a,del,lam,P)$Theta
    A.admm = admm.cov(A_hat_n,del,lam,P)$Theta
    
    c(j,
      fpr2(R.admm,R.true),tpr2(R.admm,R.true),
      fpr2(A.admm,R.true),tpr2(A.admm,R.true),
      fpr2(G.admm,G.true),tpr2(G.admm,G.true),
      fpr2(Ga.admm,G.true),tpr2(Ga.admm,G.true),
      fpr2(A.admm,G.true),tpr2(A.admm,G.true))
  }
  temp
}

rate_m1_100_08 <- foreach(j = 1:20, .combine = 'rbind') %dopar%{
  da <- da_m1_100_08[[j]]
  da_group = da %>% group_by(i)
  group_size = da_group %>% dplyr::summarise(N = n())
  group_size = group_size$N
  group_var = da_group %>% do(Var = var(.[,-c(1,2)]))
  group_var = group_var$Var
  group_mean = da_group %>% do(Mean = colMeans(.[,-c(1,2)]))
  group_mean = group_mean$Mean
  grand_mean = colMeans(da[,-c(1,2)])
  
  #### find moment estimation of R and G
  total_obs = Reduce("+",group_size)
  group_num = length(group_size)
  sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                          variance=group_var,size=group_size,SIMPLIFY = FALSE)
  )
  mse = sse/(total_obs - group_num)
  R_hat_n = mse
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_a <- A_hat_n - mean(1/group_size)*R_hat_n
  
  #### nested loop: compute fpr,tpr for a data simulation
  temp <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    
    R.admm = admm.cov(R_hat_n,del,lam,P)$Theta
    G.admm = admm.cov(G_hat_n,del,lam,P)$Theta
    Ga.admm = admm.cov(G_hat_a,del,lam,P)$Theta
    A.admm = admm.cov(A_hat_n,del,lam,P)$Theta
    
    c(j,
      fpr2(R.admm,R.true),tpr2(R.admm,R.true),
      fpr2(A.admm,R.true),tpr2(A.admm,R.true),
      fpr2(G.admm,G.true),tpr2(G.admm,G.true),
      fpr2(Ga.admm,G.true),tpr2(Ga.admm,G.true),
      fpr2(A.admm,G.true),tpr2(A.admm,G.true))
  }
  temp
}

rate_m1_100_07 <- foreach(j = 1:20, .combine = 'rbind') %dopar%{
  da <- da_m1_100_07[[j]]
  da_group = da %>% group_by(i)
  group_size = da_group %>% dplyr::summarise(N = n())
  group_size = group_size$N
  group_var = da_group %>% do(Var = var(.[,-c(1,2)]))
  group_var = group_var$Var
  group_mean = da_group %>% do(Mean = colMeans(.[,-c(1,2)]))
  group_mean = group_mean$Mean
  grand_mean = colMeans(da[,-c(1,2)])
  
  #### find moment estimation of R and G
  total_obs = Reduce("+",group_size)
  group_num = length(group_size)
  sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                          variance=group_var,size=group_size,SIMPLIFY = FALSE)
  )
  mse = sse/(total_obs - group_num)
  R_hat_n = mse
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_a <- A_hat_n - mean(1/group_size)*R_hat_n
  
  #### nested loop: compute fpr,tpr for a data simulation
  temp <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    
    R.admm = admm.cov(R_hat_n,del,lam,P)$Theta
    G.admm = admm.cov(G_hat_n,del,lam,P)$Theta
    Ga.admm = admm.cov(G_hat_a,del,lam,P)$Theta
    A.admm = admm.cov(A_hat_n,del,lam,P)$Theta
    
    c(j,
      fpr2(R.admm,R.true),tpr2(R.admm,R.true),
      fpr2(A.admm,R.true),tpr2(A.admm,R.true),
      fpr2(G.admm,G.true),tpr2(G.admm,G.true),
      fpr2(Ga.admm,G.true),tpr2(Ga.admm,G.true),
      fpr2(A.admm,G.true),tpr2(A.admm,G.true))
  }
  temp
}

rate_m1_100_06 <- foreach(j = 1:20, .combine = 'rbind') %dopar%{
  da <- da_m1_100_06[[j]]
  da_group = da %>% group_by(i)
  group_size = da_group %>% dplyr::summarise(N = n())
  group_size = group_size$N
  group_var = da_group %>% do(Var = var(.[,-c(1,2)]))
  group_var = group_var$Var
  group_mean = da_group %>% do(Mean = colMeans(.[,-c(1,2)]))
  group_mean = group_mean$Mean
  grand_mean = colMeans(da[,-c(1,2)])
  
  #### find moment estimation of R and G
  total_obs = Reduce("+",group_size)
  group_num = length(group_size)
  sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                          variance=group_var,size=group_size,SIMPLIFY = FALSE)
  )
  mse = sse/(total_obs - group_num)
  R_hat_n = mse
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_a <- A_hat_n - mean(1/group_size)*R_hat_n
  
  #### nested loop: compute fpr,tpr for a data simulation
  temp <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    
    R.admm = admm.cov(R_hat_n,del,lam,P)$Theta
    G.admm = admm.cov(G_hat_n,del,lam,P)$Theta
    Ga.admm = admm.cov(G_hat_a,del,lam,P)$Theta
    A.admm = admm.cov(A_hat_n,del,lam,P)$Theta
    
    c(j,
      fpr2(R.admm,R.true),tpr2(R.admm,R.true),
      fpr2(A.admm,R.true),tpr2(A.admm,R.true),
      fpr2(G.admm,G.true),tpr2(G.admm,G.true),
      fpr2(Ga.admm,G.true),tpr2(Ga.admm,G.true),
      fpr2(A.admm,G.true),tpr2(A.admm,G.true))
  }
  temp
}

rate_m1_100_05 <- foreach(j = 1:20, .combine = 'rbind') %dopar%{
  da <- da_m1_100_05[[j]]
  da_group = da %>% group_by(i)
  group_size = da_group %>% dplyr::summarise(N = n())
  group_size = group_size$N
  group_var = da_group %>% do(Var = var(.[,-c(1,2)]))
  group_var = group_var$Var
  group_mean = da_group %>% do(Mean = colMeans(.[,-c(1,2)]))
  group_mean = group_mean$Mean
  grand_mean = colMeans(da[,-c(1,2)])
  
  #### find moment estimation of R and G
  total_obs = Reduce("+",group_size)
  group_num = length(group_size)
  sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                          variance=group_var,size=group_size,SIMPLIFY = FALSE)
  )
  mse = sse/(total_obs - group_num)
  R_hat_n = mse
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_a <- A_hat_n - mean(1/group_size)*R_hat_n
  
  #### nested loop: compute fpr,tpr for a data simulation
  temp <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    
    R.admm = admm.cov(R_hat_n,del,lam,P)$Theta
    G.admm = admm.cov(G_hat_n,del,lam,P)$Theta
    Ga.admm = admm.cov(G_hat_a,del,lam,P)$Theta
    A.admm = admm.cov(A_hat_n,del,lam,P)$Theta
    
    c(j,
      fpr2(R.admm,R.true),tpr2(R.admm,R.true),
      fpr2(A.admm,R.true),tpr2(A.admm,R.true),
      fpr2(G.admm,G.true),tpr2(G.admm,G.true),
      fpr2(Ga.admm,G.true),tpr2(Ga.admm,G.true),
      fpr2(A.admm,G.true),tpr2(A.admm,G.true))
  }
  temp
}

rate_m1_100_04 <- foreach(j = 1:20, .combine = 'rbind') %dopar%{
  da <- da_m1_100_04[[j]]
  da_group = da %>% group_by(i)
  group_size = da_group %>% dplyr::summarise(N = n())
  group_size = group_size$N
  group_var = da_group %>% do(Var = var(.[,-c(1,2)]))
  group_var = group_var$Var
  group_mean = da_group %>% do(Mean = colMeans(.[,-c(1,2)]))
  group_mean = group_mean$Mean
  grand_mean = colMeans(da[,-c(1,2)])
  
  #### find moment estimation of R and G
  total_obs = Reduce("+",group_size)
  group_num = length(group_size)
  sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                          variance=group_var,size=group_size,SIMPLIFY = FALSE)
  )
  mse = sse/(total_obs - group_num)
  R_hat_n = mse
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_a <- A_hat_n - mean(1/group_size)*R_hat_n
  
  #### nested loop: compute fpr,tpr for a data simulation
  temp <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    
    R.admm = admm.cov(R_hat_n,del,lam,P)$Theta
    G.admm = admm.cov(G_hat_n,del,lam,P)$Theta
    Ga.admm = admm.cov(G_hat_a,del,lam,P)$Theta
    A.admm = admm.cov(A_hat_n,del,lam,P)$Theta
    
    c(j,
      fpr2(R.admm,R.true),tpr2(R.admm,R.true),
      fpr2(A.admm,R.true),tpr2(A.admm,R.true),
      fpr2(G.admm,G.true),tpr2(G.admm,G.true),
      fpr2(Ga.admm,G.true),tpr2(Ga.admm,G.true),
      fpr2(A.admm,G.true),tpr2(A.admm,G.true))
  }
  temp
}

rate_m1_100_03 <- foreach(j = 1:20, .combine = 'rbind') %dopar%{
  da <- da_m1_100_03[[j]]
  da_group = da %>% group_by(i)
  group_size = da_group %>% dplyr::summarise(N = n())
  group_size = group_size$N
  group_var = da_group %>% do(Var = var(.[,-c(1,2)]))
  group_var = group_var$Var
  group_mean = da_group %>% do(Mean = colMeans(.[,-c(1,2)]))
  group_mean = group_mean$Mean
  grand_mean = colMeans(da[,-c(1,2)])
  
  #### find moment estimation of R and G
  total_obs = Reduce("+",group_size)
  group_num = length(group_size)
  sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                          variance=group_var,size=group_size,SIMPLIFY = FALSE)
  )
  mse = sse/(total_obs - group_num)
  R_hat_n = mse
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_a <- A_hat_n - mean(1/group_size)*R_hat_n
  
  #### nested loop: compute fpr,tpr for a data simulation
  temp <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    
    R.admm = admm.cov(R_hat_n,del,lam,P)$Theta
    G.admm = admm.cov(G_hat_n,del,lam,P)$Theta
    Ga.admm = admm.cov(G_hat_a,del,lam,P)$Theta
    A.admm = admm.cov(A_hat_n,del,lam,P)$Theta
    
    c(j,
      fpr2(R.admm,R.true),tpr2(R.admm,R.true),
      fpr2(A.admm,R.true),tpr2(A.admm,R.true),
      fpr2(G.admm,G.true),tpr2(G.admm,G.true),
      fpr2(Ga.admm,G.true),tpr2(Ga.admm,G.true),
      fpr2(A.admm,G.true),tpr2(A.admm,G.true))
  }
  temp
}

save(rate_m1_100_10,rate_m1_100_09,rate_m1_100_08,rate_m1_100_07,
     rate_m1_100_06,rate_m1_100_05,rate_m1_100_04,rate_m1_100_03,
     file = "rate_m1_100.RData"
)



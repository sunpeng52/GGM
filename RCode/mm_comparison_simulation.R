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

load("data_mm_comparison.RData")

#### define moment estimation:
anova.est2 <- function(da){
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
  
  return(list(R=R_hat_n,G=G_hat_n,A=A_hat_n,Ga=G_hat_a))
}


p <- 50
G.true <- banded.m(p,10)

set.seed(1234)
########
A_m3_50_n_00 = list()
G_m3_50_n_00 = list()
lam_m3_50_n_00 = list()
for(j in 1:20){
  da = da_m3_50_n_00[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_00[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_00[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_00[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_01 = list()
G_m3_50_n_01 = list()
lam_m3_50_n_01 = list()
for(j in 1:20){
  da = da_m3_50_n_01[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_01[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_01[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_01[[j]] = lam.res
  print(paste0("DataSet",j))
}

mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m3_50_n_00))
mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m3_50_n_00))

mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m3_50_n_01))
mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m3_50_n_01))

########
A_m3_50_n_02 = list()
G_m3_50_n_02 = list()
lam_m3_50_n_02 = list()
for(j in 1:20){
  da = da_m3_50_n_02[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_02[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_02[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_02[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_03 = list()
G_m3_50_n_03 = list()
lam_m3_50_n_03 = list()
for(j in 1:20){
  da = da_m3_50_n_03[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_03[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_03[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_03[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_04 = list()
G_m3_50_n_04 = list()
lam_m3_50_n_04 = list()
for(j in 1:20){
  da = da_m3_50_n_04[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_04[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_04[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_04[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_05 = list()
G_m3_50_n_05 = list()
lam_m3_50_n_05 = list()
for(j in 1:20){
  da = da_m3_50_n_05[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_05[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_05[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_05[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_06 = list()
G_m3_50_n_06 = list()
lam_m3_50_n_06 = list()
for(j in 1:20){
  da = da_m3_50_n_06[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_06[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_06[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_06[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_07 = list()
G_m3_50_n_07 = list()
lam_m3_50_n_07 = list()
for(j in 1:20){
  da = da_m3_50_n_07[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_07[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_07[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_07[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_08 = list()
G_m3_50_n_08 = list()
lam_m3_50_n_08 = list()
for(j in 1:20){
  da = da_m3_50_n_08[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_08[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_08[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_08[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_09 = list()
G_m3_50_n_09 = list()
lam_m3_50_n_09 = list()
for(j in 1:20){
  da = da_m3_50_n_09[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_09[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_09[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_09[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_n_10 = list()
G_m3_50_n_10 = list()
lam_m3_50_n_10 = list()
for(j in 1:20){
  da = da_m3_50_n_10[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_n_10[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_n_10[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_n_10[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_01 = list()
G_m3_50_p_01 = list()
lam_m3_50_p_01 = list()
for(j in 1:20){
  da = da_m3_50_p_01[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_01[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_01[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_01[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_02 = list()
G_m3_50_p_02 = list()
lam_m3_50_p_02 = list()
for(j in 1:20){
  da = da_m3_50_p_02[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_02[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_02[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_02[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_03 = list()
G_m3_50_p_03 = list()
lam_m3_50_p_03 = list()
for(j in 1:20){
  da = da_m3_50_p_03[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_03[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_03[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_03[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_04 = list()
G_m3_50_p_04 = list()
lam_m3_50_p_04 = list()
for(j in 1:20){
  da = da_m3_50_p_04[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_04[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_04[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_04[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_05 = list()
G_m3_50_p_05 = list()
lam_m3_50_p_05 = list()
for(j in 1:20){
  da = da_m3_50_p_05[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_05[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_05[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_05[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_06 = list()
G_m3_50_p_06 = list()
lam_m3_50_p_06 = list()
for(j in 1:20){
  da = da_m3_50_p_06[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_06[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_06[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_06[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_07 = list()
G_m3_50_p_07 = list()
lam_m3_50_p_07 = list()
for(j in 1:20){
  da = da_m3_50_p_07[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_07[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_07[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_07[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_08 = list()
G_m3_50_p_08 = list()
lam_m3_50_p_08 = list()
for(j in 1:20){
  da = da_m3_50_p_08[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_08[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_08[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_08[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_09 = list()
G_m3_50_p_09 = list()
lam_m3_50_p_09 = list()
for(j in 1:20){
  da = da_m3_50_p_09[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_09[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_09[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_09[[j]] = lam.res
  print(paste0("DataSet",j))
}

########
A_m3_50_p_10 = list()
G_m3_50_p_10 = list()
lam_m3_50_p_10 = list()
for(j in 1:20){
  da = da_m3_50_p_10[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk2(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.a = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  
  mom = anova.est2(da)
  A_hat_n = mom$A
  G_hat_n = mom$G
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  A_m3_50_p_10[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  G_m3_50_p_10[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  lam_m3_50_p_10[[j]] = lam.res
  print(paste0("DataSet",j))
}



save(
  A_m3_50_n_00,G_m3_50_n_00,lam_m3_50_n_00,
  
  A_m3_50_n_01,G_m3_50_n_01,lam_m3_50_n_01,A_m3_50_n_02,G_m3_50_n_02,lam_m3_50_n_02,
  A_m3_50_n_03,G_m3_50_n_03,lam_m3_50_n_03,A_m3_50_n_04,G_m3_50_n_04,lam_m3_50_n_04,
  A_m3_50_n_05,G_m3_50_n_05,lam_m3_50_n_05,A_m3_50_n_06,G_m3_50_n_06,lam_m3_50_n_06,
  A_m3_50_n_07,G_m3_50_n_07,lam_m3_50_n_07,A_m3_50_n_08,G_m3_50_n_08,lam_m3_50_n_08,
  A_m3_50_n_09,G_m3_50_n_09,lam_m3_50_n_09,A_m3_50_n_10,G_m3_50_n_10,lam_m3_50_n_10,
  
  A_m3_50_p_01,G_m3_50_p_01,lam_m3_50_p_01,A_m3_50_p_02,G_m3_50_p_02,lam_m3_50_p_02,
  A_m3_50_p_03,G_m3_50_p_03,lam_m3_50_p_03,A_m3_50_p_04,G_m3_50_p_04,lam_m3_50_p_04,
  A_m3_50_p_05,G_m3_50_p_05,lam_m3_50_p_05,A_m3_50_p_06,G_m3_50_p_06,lam_m3_50_p_06,
  A_m3_50_p_07,G_m3_50_p_07,lam_m3_50_p_07,A_m3_50_p_08,G_m3_50_p_08,lam_m3_50_p_08,
  A_m3_50_p_09,G_m3_50_p_09,lam_m3_50_p_09,A_m3_50_p_10,G_m3_50_p_10,lam_m3_50_p_10,
  file = "moment_method_comparison_res.RData")






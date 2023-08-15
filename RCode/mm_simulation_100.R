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

set.seed(1234)
######### model 1, p =100
##### 03
R_m1_100_03 = list()
G_m1_100_03 = list()
Ga_m1_100_03 = list()
A_m1_100_03 = list()
lam_m1_100_03 = list()
for(j in 1:100){
  
  da = da_m1_100_03[[j]]
  
  group_num = length(unique(da[,1]))
  sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  lam_vec = 10^seq(-2,0.5,0.05)
  lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
    temp = matrix(nrow = 5,ncol = 4)
    for(chunkid in 1:5){
      temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
    }
    colMeans(temp)
  }
  lam.r = lam_vec[which.min(lam.res[,1])]
  lam.g = lam_vec[which.min(lam.res[,2])]
  lam.ga = lam_vec[which.min(lam.res[,3])]
  lam.a = lam_vec[which.min(lam.res[,4])]
  
  mom = anova.est2(da)
  R_hat_n = mom$R
  G_hat_n = mom$G
  G_hat_a = mom$Ga
  A_hat_n = mom$A
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  R_m1_100_03[[j]] = admm.cov(R_hat_n,del,lam.r,P)
  G_m1_100_03[[j]] = admm.cov(G_hat_n,del,lam.g,P)
  Ga_m1_100_03[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
  A_m1_100_03[[j]] = admm.cov(A_hat_n,del,lam.a,P)
  lam_m1_100_03[[j]] = lam.res
  print(paste0("DataSet",j))
}

##### 04
R_m1_100_04 = list()
G_m1_100_04 = list()
Ga_m1_100_04 = list()
A_m1_100_04 = list()
lam_m1_100_04 = list()
for(j in 1:100){
        
        da = da_m1_100_04[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m1_100_04[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m1_100_04[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m1_100_04[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m1_100_04[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m1_100_04[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 05
R_m1_100_05 = list()
G_m1_100_05 = list()
Ga_m1_100_05 = list()
A_m1_100_05 = list()
lam_m1_100_05 = list()
for(j in 1:100){
        
        da = da_m1_100_05[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m1_100_05[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m1_100_05[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m1_100_05[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m1_100_05[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m1_100_05[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 06
R_m1_100_06 = list()
G_m1_100_06 = list()
Ga_m1_100_06 = list()
A_m1_100_06 = list()
lam_m1_100_06 = list()
for(j in 1:100){
        
        da = da_m1_100_06[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m1_100_06[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m1_100_06[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m1_100_06[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m1_100_06[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m1_100_06[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 07
R_m1_100_07 = list()
G_m1_100_07 = list()
Ga_m1_100_07 = list()
A_m1_100_07 = list()
lam_m1_100_07 = list()
for(j in 1:100){
        
        da = da_m1_100_07[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m1_100_07[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m1_100_07[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m1_100_07[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m1_100_07[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m1_100_07[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 08
R_m1_100_08 = list()
G_m1_100_08 = list()
Ga_m1_100_08 = list()
A_m1_100_08 = list()
lam_m1_100_08 = list()
for(j in 1:100){
        
        da = da_m1_100_08[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m1_100_08[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m1_100_08[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m1_100_08[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m1_100_08[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m1_100_08[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 09
R_m1_100_09 = list()
G_m1_100_09 = list()
Ga_m1_100_09 = list()
A_m1_100_09 = list()
lam_m1_100_09 = list()
for(j in 1:100){
        
        da = da_m1_100_09[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m1_100_09[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m1_100_09[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m1_100_09[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m1_100_09[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m1_100_09[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 10
R_m1_100_10 = list()
G_m1_100_10 = list()
Ga_m1_100_10 = list()
A_m1_100_10 = list()
lam_m1_100_10 = list()
for(j in 1:100){
        
        da = da_m1_100_10[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m1_100_10[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m1_100_10[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m1_100_10[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m1_100_10[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m1_100_10[[j]] = lam.res
        print(paste0("DataSet",j))
}


######### model 2, p =100
##### 03
R_m2_100_03 = list()
G_m2_100_03 = list()
Ga_m2_100_03 = list()
A_m2_100_03 = list()
lam_m2_100_03 = list()
for(j in 1:100){
        
        da = da_m2_100_03[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m2_100_03[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m2_100_03[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m2_100_03[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m2_100_03[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m2_100_03[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 04
R_m2_100_04 = list()
G_m2_100_04 = list()
Ga_m2_100_04 = list()
A_m2_100_04 = list()
lam_m2_100_04 = list()
for(j in 1:100){
        
        da = da_m2_100_04[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m2_100_04[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m2_100_04[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m2_100_04[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m2_100_04[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m2_100_04[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 05
R_m2_100_05 = list()
G_m2_100_05 = list()
Ga_m2_100_05 = list()
A_m2_100_05 = list()
lam_m2_100_05 = list()
for(j in 1:100){
        
        da = da_m2_100_05[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m2_100_05[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m2_100_05[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m2_100_05[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m2_100_05[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m2_100_05[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 06
R_m2_100_06 = list()
G_m2_100_06 = list()
Ga_m2_100_06 = list()
A_m2_100_06 = list()
lam_m2_100_06 = list()
for(j in 1:100){
        
        da = da_m2_100_06[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m2_100_06[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m2_100_06[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m2_100_06[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m2_100_06[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m2_100_06[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 07
R_m2_100_07 = list()
G_m2_100_07 = list()
Ga_m2_100_07 = list()
A_m2_100_07 = list()
lam_m2_100_07 = list()
for(j in 1:100){
        
        da = da_m2_100_07[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m2_100_07[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m2_100_07[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m2_100_07[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m2_100_07[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m2_100_07[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 08
R_m2_100_08 = list()
G_m2_100_08 = list()
Ga_m2_100_08 = list()
A_m2_100_08 = list()
lam_m2_100_08 = list()
for(j in 1:100){
        
        da = da_m2_100_08[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m2_100_08[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m2_100_08[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m2_100_08[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m2_100_08[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m2_100_08[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 09
R_m2_100_09 = list()
G_m2_100_09 = list()
Ga_m2_100_09 = list()
A_m2_100_09 = list()
lam_m2_100_09 = list()
for(j in 1:100){
        
        da = da_m2_100_09[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m2_100_09[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m2_100_09[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m2_100_09[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m2_100_09[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m2_100_09[[j]] = lam.res
        print(paste0("DataSet",j))
}

##### 10
R_m2_100_10 = list()
G_m2_100_10 = list()
Ga_m2_100_10 = list()
A_m2_100_10 = list()
lam_m2_100_10 = list()
for(j in 1:100){
        
        da = da_m2_100_10[[j]]
        
        group_num = length(unique(da[,1]))
        sample.group.id = sample(cut(1:group_num, 5, labels=F))
        
        lam_vec = 10^seq(-2,0.5,0.05)
        lam.res <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
                temp = matrix(nrow = 5,ncol = 4)
                for(chunkid in 1:5){
                        temp[chunkid,] = do.chunk(chunkid, sample.group.id,da,lam)
                }
                colMeans(temp)
        }
        lam.r = lam_vec[which.min(lam.res[,1])]
        lam.g = lam_vec[which.min(lam.res[,2])]
        lam.ga = lam_vec[which.min(lam.res[,3])]
        lam.a = lam_vec[which.min(lam.res[,4])]
        
        mom = anova.est2(da)
        R_hat_n = mom$R
        G_hat_n = mom$G
        G_hat_a = mom$Ga
        A_hat_n = mom$A
        
        del = 10^(-5)
        p = dim(G_hat_n)[1]
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        R_m2_100_10[[j]] = admm.cov(R_hat_n,del,lam.r,P)
        G_m2_100_10[[j]] = admm.cov(G_hat_n,del,lam.g,P)
        Ga_m2_100_10[[j]] = admm.cov(G_hat_a,del,lam.ga,P)
        A_m2_100_10[[j]] = admm.cov(A_hat_n,del,lam.a,P)
        lam_m2_100_10[[j]] = lam.res
        print(paste0("DataSet",j))
}

save(
  R_m1_100_03, G_m1_100_03, A_m1_100_03, lam_m1_100_03, R_m1_100_04, G_m1_100_04, A_m1_100_04, lam_m1_100_04,
  R_m1_100_05, G_m1_100_05, A_m1_100_05, lam_m1_100_05, R_m1_100_06, G_m1_100_06, A_m1_100_06, lam_m1_100_06,
  R_m1_100_07, G_m1_100_07, A_m1_100_07, lam_m1_100_07, R_m1_100_08, G_m1_100_08, A_m1_100_08, lam_m1_100_08,
  R_m1_100_09, G_m1_100_09, A_m1_100_09, lam_m1_100_09, R_m1_100_10, G_m1_100_10, A_m1_100_10, lam_m1_100_10,
  
  R_m2_100_03, G_m2_100_03, A_m2_100_03, lam_m2_100_03, R_m2_100_04, G_m2_100_04, A_m2_100_04, lam_m2_100_04,
  R_m2_100_05, G_m2_100_05, A_m2_100_05, lam_m2_100_05, R_m2_100_06, G_m2_100_06, A_m2_100_06, lam_m2_100_06,
  R_m2_100_07, G_m2_100_07, A_m2_100_07, lam_m2_100_07, R_m2_100_08, G_m2_100_08, A_m2_100_08, lam_m2_100_08,
  R_m2_100_09, G_m2_100_09, A_m2_100_09, lam_m2_100_09, R_m2_100_10, G_m2_100_10, A_m2_100_10, lam_m2_100_10,
  
  Ga_m1_100_03, Ga_m1_100_04, Ga_m1_100_05, Ga_m1_100_06,
  Ga_m1_100_07, Ga_m1_100_08, Ga_m1_100_09, Ga_m1_100_10,
  Ga_m2_100_03, Ga_m2_100_04, Ga_m2_100_05, Ga_m2_100_06,
  Ga_m2_100_07, Ga_m2_100_08, Ga_m2_100_09, Ga_m2_100_10,
  
  file = "mm_simulation_100_res22.RData"
)

print(1)

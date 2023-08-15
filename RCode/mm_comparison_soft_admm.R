#### 
library(dplyr)
library(spcov)
library(MASS)

library(parallel)
library(foreach)
library(doParallel)

source("functions.R")
source("functions_2.R")
options(dplyr.summarise.inform = FALSE)
load("time_check_res.RData")

numCores <- detectCores()
numCores

registerDoParallel(numCores)

### define three functions for cross-validation

Estimate = function(da,lam,term){
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
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  if(term=="R"){
    return(admm.cov(R_hat_n,del,lam,P))
  }
  if(term=="G"){
    return(admm.cov(G_hat_n,del,lam,P))
  }
  
}

Validate = function(da,est,term){
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
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
  
  if(term=="R"){
    return(norm(est$Theta - R_hat_n,"F")^2)
  }
  if(term=="G"){
    return(norm(est$Theta - G_hat_n,"F")^2)
  }
  
}

Estimate.soft = function(da,lam,term){
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
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  if(term=="R"){
    return(SoftThreshold(R_hat_n, lam * P))
  }
  if(term=="G"){
    return(SoftThreshold(G_hat_n, lam * P))
  }
  
}

Validate.soft = function(da,est,term){
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
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
  
  if(term=="R"){
    return(norm(est - R_hat_n,"F")^2)
  }
  if(term=="G"){
    return(norm(est - G_hat_n,"F")^2)
  }
  
}

do.chunk = function(chunkid, chunkdef,da,lam){
  train_group = (chunkdef!=chunkid)
  test_group = !train_group
  
  group_num = length(unique(da[,1]))
  train = da[da$i %in% (1:group_num)[train_group],]
  test = da[da$i %in% (1:group_num)[test_group],]
  
  res.1 = c(Validate(test,Estimate(train,lam,"R"),"R"),Validate(test,Estimate(train,lam,"G"),"G"))
  res.2 = c(Validate.soft(test,Estimate.soft(train,lam,"R"),"R"),Validate.soft(test,Estimate.soft(train,lam,"G"),"G"))
  res = c(res.1,res.2)
  names(res) = c("r.error","g.error","r.error.soft","g.error.soft")
  return(res)
}



set.seed(12345)

##### model 1, p = 100
p <-  100
R.true <- banded.new(p,10,-1)
G.true <- banded.new(p,10)

R_m1_100 = list()
G_m1_100 = list()
R_m1_100.soft = list()
G_m1_100.soft = list()
lam_m1_100 = list()
for(j in 1:50){
  #group_num = length(unique(da[,1]))
  #sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  da = da_m1_100[[j]]
  
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
  lam.r.soft = lam_vec[which.min(lam.res[,3])]
  lam.g.soft = lam_vec[which.min(lam.res[,4])]
  
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
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  R_m1_100[[j]] = admm.cov(R_hat_n,del,lam.r)
  G_m1_100[[j]] = admm.cov(G_hat_n,del,lam.g)
  R_m1_100.soft[[j]] = SoftThreshold(R_hat_n, lam.r.soft * P)
  G_m1_100.soft[[j]] = SoftThreshold(G_hat_n, lam.g.soft * P)
  lam_m1_100[[j]] = c(lam.r,lam.g,lam.r.soft,lam.g.soft)
  
  print(paste0("DataSet",j))
}


##### model 1, p = 200
p <-  100
R.true <- banded.new(p,10,-1)
G.true <- banded.new(p,10)

R_m1_200 = list()
G_m1_200 = list()
R_m1_200.soft = list()
G_m1_200.soft = list()
lam_m1_200 = list()
for(j in 1:50){
  #group_num = length(unique(da[,1]))
  #sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  da = da_m1_200[[j]]
  
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
  lam.r.soft = lam_vec[which.min(lam.res[,3])]
  lam.g.soft = lam_vec[which.min(lam.res[,4])]
  
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
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  R_m1_200[[j]] = admm.cov(R_hat_n,del,lam.r)
  G_m1_200[[j]] = admm.cov(G_hat_n,del,lam.g)
  R_m1_200.soft[[j]] = SoftThreshold(R_hat_n, lam.r.soft * P)
  G_m1_200.soft[[j]] = SoftThreshold(G_hat_n, lam.g.soft * P)
  lam_m1_200[[j]] = c(lam.r,lam.g,lam.r.soft,lam.g.soft)
  
  print(paste0("DataSet",j))
}


##### model 2, p = 100
p <-  100
R.true <- ar.cov(p,-0.6)
G.true <- ar.cov(p,0.6)

R_m2_100 = list()
G_m2_100 = list()
R_m2_100.soft = list()
G_m2_100.soft = list()
lam_m2_100 = list()
for(j in 1:50){
  #group_num = length(unique(da[,1]))
  #sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  da = da_m2_100[[j]]
  
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
  lam.r.soft = lam_vec[which.min(lam.res[,3])]
  lam.g.soft = lam_vec[which.min(lam.res[,4])]
  
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
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  R_m2_100[[j]] = admm.cov(R_hat_n,del,lam.r)
  G_m2_100[[j]] = admm.cov(G_hat_n,del,lam.g)
  R_m2_100.soft[[j]] = SoftThreshold(R_hat_n, lam.r.soft * P)
  G_m2_100.soft[[j]] = SoftThreshold(G_hat_n, lam.g.soft * P)
  lam_m2_100[[j]] = c(lam.r,lam.g,lam.r.soft,lam.g.soft)
  
  print(paste0("DataSet",j))
}


##### model 2, p = 200
p <-  200
R.true <- ar.cov(p,-0.6)
G.true <- ar.cov(p,0.6)

R_m2_200 = list()
G_m2_200 = list()
R_m2_200.soft = list()
G_m2_200.soft = list()
lam_m2_200 = list()
for(j in 1:50){
  #group_num = length(unique(da[,1]))
  #sample.group.id = sample(cut(1:group_num, 5, labels=F))
  
  da = da_m2_200[[j]]
  
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
  lam.r.soft = lam_vec[which.min(lam.res[,3])]
  lam.g.soft = lam_vec[which.min(lam.res[,4])]
  
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
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
  
  del = 10^(-5)
  p = dim(G_hat_n)[1]
  P <- matrix(1,p,p)
  diag(P) <- 0
  
  R_m2_200[[j]] = admm.cov(R_hat_n,del,lam.r)
  G_m2_200[[j]] = admm.cov(G_hat_n,del,lam.g)
  R_m2_200.soft[[j]] = SoftThreshold(R_hat_n, lam.r.soft * P)
  G_m2_200.soft[[j]] = SoftThreshold(G_hat_n, lam.g.soft * P)
  lam_m2_200[[j]] = c(lam.r,lam.g,lam.r.soft,lam.g.soft)
  
  print(paste0("DataSet",j))
}


save(R_m1_100,
     G_m1_100,
     R_m1_100.soft,
     G_m1_100.soft,
     lam_m1_100,
     R_m1_200,
     G_m1_200,
     R_m1_200.soft,
     G_m1_200.soft,
     lam_m1_200,
     R_m2_100,
     G_m2_100,
     R_m2_100.soft,
     G_m2_100.soft,
     lam_m2_100,
     R_m2_200,
     G_m2_200,
     R_m2_200.soft,
     G_m2_200.soft,
     lam_m2_200,
     file = "time_check_res2.RData")



### model 1, p =100
p <-  100
R.true <- banded.new(p,10,-1)
G.true <- banded.new(p,10)

g_m1_100_f=mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100)
g_m1_100_f_soft=mapply(function(res) norm(res - G.true,"F"),res=G_m1_100.soft)
g_m1_100_l=mapply(function(res) norm(res$Theta - G.true,"2"),res=G_m1_100)
g_m1_100_l_soft=mapply(function(res) norm(res - G.true,"2"),res=G_m1_100.soft)

round(c(mean(g_m1_100_f_soft),mean(g_m1_100_f),
        mean(g_m1_100_l_soft),mean(g_m1_100_l)),4)

round(c(sd(g_m1_100_f_soft),sd(g_m1_100_f),
        sd(g_m1_100_l_soft),sd(g_m1_100_l))/sqrt(p),4)

sum(mapply(function(res) min(eigen(res$Theta)$values)>0,res=G_m1_100))/100
sum(mapply(function(res) min(eigen(res)$values)>0,res=G_m1_100.soft))/100

r_m1_100_f=mapply(function(res) norm(res$Theta - G.true,"F"),res=R_m1_100)
r_m1_100_f_soft=mapply(function(res) norm(res - G.true,"F"),res=R_m1_100.soft)
r_m1_100_l=mapply(function(res) norm(res$Theta - G.true,"2"),res=R_m1_100)
r_m1_100_l_soft=mapply(function(res) norm(res - G.true,"2"),res=R_m1_100.soft)

round(c(mean(r_m1_100_f_soft),mean(r_m1_100_f),
        mean(r_m1_100_l_soft),mean(r_m1_100_l)),4)

round(c(sd(r_m1_100_f_soft),sd(r_m1_100_f),
        sd(r_m1_100_l_soft),sd(r_m1_100_l))/sqrt(p),4)

sum(mapply(function(res) min(eigen(res$Theta)$values)>0,res=R_m1_100))/100
sum(mapply(function(res) min(eigen(res)$values)>0,res=R_m1_100.soft))/100

### model 1, p =200
p <-  200
R.true <- banded.new(p,10,-1)
G.true <- banded.new(p,10)

g_m1_200_f=mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200)
g_m1_200_f_soft=mapply(function(res) norm(res - G.true,"F"),res=G_m1_200.soft)
g_m1_200_l=mapply(function(res) norm(res$Theta - G.true,"2"),res=G_m1_200)
g_m1_200_l_soft=mapply(function(res) norm(res - G.true,"2"),res=G_m1_200.soft)

round(c(mean(g_m1_200_f_soft),mean(g_m1_200_f),
        mean(g_m1_200_l_soft),mean(g_m1_200_l)),4)

round(c(sd(g_m1_200_f_soft),sd(g_m1_200_f),
        sd(g_m1_200_l_soft),sd(g_m1_200_l))/sqrt(p),4)

sum(mapply(function(res) min(eigen(res$Theta)$values)>0,res=G_m1_200))/100
sum(mapply(function(res) min(eigen(res)$values)>0,res=G_m1_200.soft))/100

r_m1_200_f=mapply(function(res) norm(res$Theta - G.true,"F"),res=R_m1_200)
r_m1_200_f_soft=mapply(function(res) norm(res - G.true,"F"),res=R_m1_200.soft)
r_m1_200_l=mapply(function(res) norm(res$Theta - G.true,"2"),res=R_m1_200)
r_m1_200_l_soft=mapply(function(res) norm(res - G.true,"2"),res=R_m1_200.soft)

round(c(mean(r_m1_200_f_soft),mean(r_m1_200_f),
        mean(r_m1_200_l_soft),mean(r_m1_200_l)),4)

round(c(sd(r_m1_200_f_soft),sd(r_m1_200_f),
        sd(r_m1_200_l_soft),sd(r_m1_200_l))/sqrt(p),4)

sum(mapply(function(res) min(eigen(res$Theta)$values)>0,res=R_m1_200))/100
sum(mapply(function(res) min(eigen(res)$values)>0,res=R_m1_200.soft))/100


### model 2, p =100
p <-  100
R.true <- ar.cov(p,-0.6)
G.true <- ar.cov(p,0.6)

g_m2_100_f=mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100)
g_m2_100_f_soft=mapply(function(res) norm(res - G.true,"F"),res=G_m2_100.soft)
g_m2_100_l=mapply(function(res) norm(res$Theta - G.true,"2"),res=G_m2_100)
g_m2_100_l_soft=mapply(function(res) norm(res - G.true,"2"),res=G_m2_100.soft)

round(c(mean(g_m2_100_f_soft),mean(g_m2_100_f),
        mean(g_m2_100_l_soft),mean(g_m2_100_l)),4)

round(c(sd(g_m2_100_f_soft),sd(g_m2_100_f),
        sd(g_m2_100_l_soft),sd(g_m2_100_l))/sqrt(p),4)

sum(mapply(function(res) min(eigen(res$Theta)$values)>0,res=G_m2_100))/100
sum(mapply(function(res) min(eigen(res)$values)>0,res=G_m2_100.soft))/100

r_m2_100_f=mapply(function(res) norm(res$Theta - G.true,"F"),res=R_m2_100)
r_m2_100_f_soft=mapply(function(res) norm(res - G.true,"F"),res=R_m2_100.soft)
r_m2_100_l=mapply(function(res) norm(res$Theta - G.true,"2"),res=R_m2_100)
r_m2_100_l_soft=mapply(function(res) norm(res - G.true,"2"),res=R_m2_100.soft)

round(c(mean(r_m2_100_f_soft),mean(r_m2_100_f),
        mean(r_m2_100_l_soft),mean(r_m2_100_l)),4)

round(c(sd(r_m2_100_f_soft),sd(r_m2_100_f),
        sd(r_m2_100_l_soft),sd(r_m2_100_l))/sqrt(p),4)

sum(mapply(function(res) min(eigen(res$Theta)$values)>0,res=R_m2_100))/100
sum(mapply(function(res) min(eigen(res)$values)>0,res=R_m2_100.soft))/100



### model 2, p = 200
p <-  200
R.true <- ar.cov(p,-0.6)
G.true <- ar.cov(p,0.6)

g_m2_200_f=mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200)
g_m2_200_f_soft=mapply(function(res) norm(res - G.true,"F"),res=G_m2_200.soft)
g_m2_200_l=mapply(function(res) norm(res$Theta - G.true,"2"),res=G_m2_200)
g_m2_200_l_soft=mapply(function(res) norm(res - G.true,"2"),res=G_m2_200.soft)

round(c(mean(g_m2_200_f_soft),mean(g_m2_200_f),
        mean(g_m2_200_l_soft),mean(g_m2_200_l)),4)

round(c(sd(g_m2_200_f_soft),sd(g_m2_200_f),
        sd(g_m2_200_l_soft),sd(g_m2_200_l))/sqrt(p),4)

sum(mapply(function(res) min(eigen(res$Theta)$values)>0,res=G_m2_200))/100
sum(mapply(function(res) min(eigen(res)$values)>0,res=G_m2_200.soft))/100

r_m2_200_f=mapply(function(res) norm(res$Theta - G.true,"F"),res=R_m2_200)
r_m2_200_f_soft=mapply(function(res) norm(res - G.true,"F"),res=R_m2_200.soft)
r_m2_200_l=mapply(function(res) norm(res$Theta - G.true,"2"),res=R_m2_200)
r_m2_200_l_soft=mapply(function(res) norm(res - G.true,"2"),res=R_m2_200.soft)

round(c(mean(r_m2_200_f_soft),mean(r_m2_200_f),
        mean(r_m2_200_l_soft),mean(r_m2_200_l)),4)

round(c(sd(r_m2_200_f_soft),sd(r_m2_200_f),
        sd(r_m2_200_l_soft),sd(r_m2_200_l))/sqrt(p),4)

sum(mapply(function(res) min(eigen(res$Theta)$values)>0,res=R_m2_200))/100
sum(mapply(function(res) min(eigen(res)$values)>0,res=R_m2_200.soft))/100






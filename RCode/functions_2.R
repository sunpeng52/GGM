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
  if(term=="A"){
    return(admm.cov(A_hat_n,del,lam,P))
  }
  if(term=="Ga"){
    return(admm.cov(G_hat_a,del,lam,P))
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
  
  if(term=="R"){
    return(norm(est$Theta - R_hat_n,"F")^2)
  }
  if(term=="G"){
    return(norm(est$Theta - G_hat_n,"F")^2)
  }
  if(term=="A"){
    return(norm(est$Theta - A_hat_n,"F")^2)
  }
  if(term=="Ga"){
    return(norm(est$Theta - G_hat_a,"F")^2)
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
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
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
  if(term=="A"){
          return(SoftThreshold(A_hat_n, lam * P))
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
  
  ssa = Reduce("+",mapply(function(size,mean) size*matrix(mean-grand_mean,ncol = 1)%*%t(matrix(mean-grand_mean,ncol = 1)),
                          mean=group_mean,size=group_size,SIMPLIFY = FALSE)
  )
  msa = ssa/(group_num-1)
  
  ni_square = sum(mapply(function(x) x^2,group_size))
  n0 = (total_obs-ni_square/total_obs)/(group_num - 1)
  G_hat_n = (msa-mse)/n0
  
  da_avg <- Reduce("rbind",group_mean)
  A_hat_n <- cov(da_avg)
  
  
  if(term=="R"){
          return(norm(est - R_hat_n,"F")^2)
  }
  if(term=="G"){
          return(norm(est - G_hat_n,"F")^2)
  }
  if(term=="A"){
          return(norm(est - A_hat_n,"F")^2)
  }
  
}

do.chunk = function(chunkid, chunkdef,da,lam){
  train_group = (chunkdef!=chunkid)
  test_group = !train_group
  
  group_num = length(unique(da[,1]))
  train = da[da$i %in% (1:group_num)[train_group],]
  test = da[da$i %in% (1:group_num)[test_group],]
  
  res.1 = c(Validate(test,Estimate(train,lam,"R"),"R"),
            Validate(test,Estimate(train,lam,"G"),"G"),
            Validate(test,Estimate(train,lam,"Ga"),"Ga"),
            Validate(test,Estimate(train,lam,"A"),"A"))
  #res.2 = c(Validate.soft(test,Estimate.soft(train,lam,"R"),"R"),Validate.soft(test,Estimate.soft(train,lam,"G"),"G"))
  #res = c(res.1,res.2)
  #names(res) = c("r.error","g.error","r.error.soft","g.error.soft")
  names(res.1) = c("r.error","g.error","g.error.a","a.error")
  return(res.1)
}


do.chunk2 = function(chunkid, chunkdef,da,lam){
        train_group = (chunkdef!=chunkid)
        test_group = !train_group
        
        group_num = length(unique(da[,1]))
        train = da[da$i %in% (1:group_num)[train_group],]
        test = da[da$i %in% (1:group_num)[test_group],]
        
        res.1 = c(Validate(test,Estimate(train,lam,"A"),"A"),Validate(test,Estimate(train,lam,"G"),"G"))
        #res.2 = c(Validate.soft(test,Estimate.soft(train,lam,"R"),"R"),Validate.soft(test,Estimate.soft(train,lam,"G"),"G"))
        #res = c(res.1,res.2)
        #names(res) = c("r.error","g.error","r.error.soft","g.error.soft")
        names(res.1) = c("a.error","g.error")
        return(res.1)
}


do.chunksoft2 = function(chunkid, chunkdef,da,lam){
        train_group = (chunkdef!=chunkid)
        test_group = !train_group
        
        group_num = length(unique(da[,1]))
        train = da[da$i %in% (1:group_num)[train_group],]
        test = da[da$i %in% (1:group_num)[test_group],]
        
        #res.1 = c(Validate(test,Estimate(train,lam,"A"),"A"),Validate(test,Estimate(train,lam,"G"),"G"))
        res.2 = c(Validate.soft(test,Estimate.soft(train,lam,"A"),"A"),Validate.soft(test,Estimate.soft(train,lam,"G"),"G"))
        #res = c(res.1,res.2)
        #names(res) = c("r.error","g.error","r.error.soft","g.error.soft")
        names(res.2) = c("a.error","g.error")
        return(res.2)
}

banded.new <- function(n,b,same=1) {
        ## generate banded matrix 
        diff.ij <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                               (1:n - 1))
        (1-diff.ij/b)*(1 > (diff.ij/b))*(same)^diff.ij
}

ar.cov <- function(n,rho) {
        ## generate banded matrix 
        diff.ij <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                               (1:n - 1))
        (rho)^diff.ij
}

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
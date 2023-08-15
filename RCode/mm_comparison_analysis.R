library(dplyr)
library(spcov)
library(MASS)

#library(latex2exp)
#library(ggplot2)
#library(ggpubr)
#library(RColorBrewer)

source("functions.R")
source("functions_2.R")
options(dplyr.summarise.inform = FALSE)

load("data_mm_comparison.RData")
load("mm_comparison_res.RData")

p <- 50
G.true <- banded.m(p,10)

####
norm(diag(1,50),"F")/norm(G.true,"F")

R.base1 <- banded.new(50,10,-1)
R.base2 <- banded.new(50,10)

r.norm1 = c()
for(a in 1:10){
  r.norm1 = c(r.norm1,norm(a*R.base1,"F"))
}

r.norm2 = c()
for(a in 1:10){
  r.norm2 = c(r.norm2,norm(a*R.base2,"F"))
}

c(norm(diag(1,50),"F"),r.norm1)/norm(G.true,"F")
c(norm(diag(1,50),"F"),r.norm2)/norm(G.true,"F")


##### get moment estimation
f_norm_n_00 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_00[[i]])
  f_norm_n_00 = rbind(f_norm_n_00,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_01 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_01[[i]])
  f_norm_n_01 = rbind(f_norm_n_01,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_02 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_02[[i]])
  f_norm_n_02 = rbind(f_norm_n_02,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_03 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_03[[i]])
  f_norm_n_03 = rbind(f_norm_n_03,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_04 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_04[[i]])
  f_norm_n_04 = rbind(f_norm_n_04,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_05 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_05[[i]])
  f_norm_n_05 = rbind(f_norm_n_05,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_06 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_06[[i]])
  f_norm_n_06 = rbind(f_norm_n_06,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_07 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_07[[i]])
  f_norm_n_07 = rbind(f_norm_n_07,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_08 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_08[[i]])
  f_norm_n_08 = rbind(f_norm_n_08,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_09 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_09[[i]])
  f_norm_n_09 = rbind(f_norm_n_09,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n_10 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_n_10[[i]])
  f_norm_n_10 = rbind(f_norm_n_10,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_01 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_01[[i]])
  f_norm_p_01 = rbind(f_norm_p_01,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_02 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_02[[i]])
  f_norm_p_02 = rbind(f_norm_p_02,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_03 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_03[[i]])
  f_norm_p_03 = rbind(f_norm_p_03,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_04 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_04[[i]])
  f_norm_p_04 = rbind(f_norm_p_04,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_05 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_05[[i]])
  f_norm_p_05 = rbind(f_norm_p_05,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_06 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_06[[i]])
  f_norm_p_06 = rbind(f_norm_p_06,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_07 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_07[[i]])
  f_norm_p_07 = rbind(f_norm_p_07,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_08 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_08[[i]])
  f_norm_p_08 = rbind(f_norm_p_08,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_09 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_09[[i]])
  f_norm_p_09 = rbind(f_norm_p_09,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_p_10 = c()
for(i in 1:100){
  mom = anova.est2(da_m3_50_p_10[[i]])
  f_norm_p_10 = rbind(f_norm_p_10,c(norm(mom$G - G.true,"F"), norm(mom$A - G.true,"F")))
}

f_norm_n <- rbind(
  colMeans(f_norm_n_00[1:100,]),
  colMeans(f_norm_n_01[1:100,]),
  colMeans(f_norm_n_02[1:100,]),
  colMeans(f_norm_n_03[1:100,]),
  colMeans(f_norm_n_04[1:100,]),
  colMeans(f_norm_n_05[1:100,]),
  colMeans(f_norm_n_06[1:100,]),
  colMeans(f_norm_n_07[1:100,]),
  colMeans(f_norm_n_08[1:100,]),
  colMeans(f_norm_n_09[1:100,]),
  colMeans(f_norm_n_10[1:100,])
)

f_norm_p <- rbind(
  colMeans(f_norm_n_00[1:100,]),
  colMeans(f_norm_p_01[1:100,]),
  colMeans(f_norm_p_02[1:100,]),
  colMeans(f_norm_p_03[1:100,]),
  colMeans(f_norm_p_04[1:100,]),
  colMeans(f_norm_p_05[1:100,]),
  colMeans(f_norm_p_06[1:100,]),
  colMeans(f_norm_p_07[1:100,]),
  colMeans(f_norm_p_08[1:100,]),
  colMeans(f_norm_p_09[1:100,]),
  colMeans(f_norm_p_10[1:100,])
)



#P <- matrix(1,p,p)
#diag(P) = 0


g_err_n <- c(mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_00)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_01)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_02)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_03)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_04)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_05)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_06)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_07)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_08)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_09)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_10)))


a_err_n <- c(mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_00)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_01)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_02)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_03)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_04)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_05)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_06)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_07)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_08)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_09)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_10)))

g_err_p <- c(mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_n_00)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_01)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_02)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_03)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_04)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_05)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_06)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_07)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_08)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_09)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=G_m3_50_p_10)))


a_err_p <- c(mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_n_00)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_01)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_02)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_03)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_04)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_05)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_06)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_07)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_08)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_09)),
             mean(mapply(function(res) norm((res$Theta - G.true),"F"),res=A_m3_50_p_10)))

## model 4: negative; model 3: positive
da <- data.frame(
  a = rep(c(norm(diag(1,50),"F"),r.norm1)/norm(G.true,"F"),8),
  error = c(f_norm_p[,1],f_norm_p[,2],g_err_p,a_err_p,
            f_norm_n[,1],f_norm_n[,2],g_err_n,a_err_n),
  estimator = rep(c("Moment","Constrained"),each=22,times=2),
  scenario = rep(c("Model 3","Model 4"),each=44),
  #scenario = rep(c("Scenario 2","Scenario 1"),each=44),
  type = rep(c("Between Subject","Subject Aggregation"),each = 11,times=4)
)

write.csv(da,file="mm_comparison.csv",row.names = F)

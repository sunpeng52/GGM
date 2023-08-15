rm(list = ls())
library(latex2exp)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(cowplot)

library(dplyr)
library(spcov)
library(MASS)

source("functions.R")
source("functions_2.R")
options(dplyr.summarise.inform = FALSE)

balance = c()
for(x in 3:10){
        N = 1000
        n = c(rep(x,99),N-99*x)
        n_eff = (N-sum(n^2)/N)/(length(n)-1)
        balance = c(balance,max(n)/n_eff)
        
}

load("mm_simulation_100_res22.RData")

library(reshape2)

cv_m1_100_10_r <- sapply(lam_m1_100_10,function(res,col) res[,col],col=1)
cv_m1_100_10_r <- as.data.frame(cv_m1_100_10_r)
cv_m1_100_10_r$lam <- seq(-2,0.5,0.05)
cv_m1_100_10_r <- melt(cv_m1_100_10_r, id.vars = "lam")
cv_m1_100_10_r$balance <- rep("1.00",nrow(cv_m1_100_10_r))
cv_m1_100_10_r$type <- rep("SPD-MANOVA-W",nrow(cv_m1_100_10_r))


cv_m1_100_10_g <- sapply(lam_m1_100_10,function(res,col) res[,col],col=2)
cv_m1_100_10_g <- as.data.frame(cv_m1_100_10_g)
cv_m1_100_10_g$lam <- seq(-2,0.5,0.05)
cv_m1_100_10_g <- melt(cv_m1_100_10_g, id.vars = "lam")
cv_m1_100_10_g$balance <- rep("1.00",nrow(cv_m1_100_10_g))
cv_m1_100_10_g$type <- rep("SPD-MANOVA-B",nrow(cv_m1_100_10_g))

cv_m1_100_10_ga <- sapply(lam_m1_100_10,function(res,col) res[,col],col=3)
cv_m1_100_10_ga <- as.data.frame(cv_m1_100_10_ga)
cv_m1_100_10_ga$lam <- seq(-2,0.5,0.05)
cv_m1_100_10_ga <- melt(cv_m1_100_10_ga, id.vars = "lam")
cv_m1_100_10_ga$balance <- rep("1.00",nrow(cv_m1_100_10_ga))
cv_m1_100_10_ga$type <- rep("SPD-USS-B",nrow(cv_m1_100_10_ga))

cv_m1_100_10_a <- sapply(lam_m1_100_10,function(res,col) res[,col],col=4)
cv_m1_100_10_a <- as.data.frame(cv_m1_100_10_a)
cv_m1_100_10_a$lam <- seq(-2,0.5,0.05)
cv_m1_100_10_a <- melt(cv_m1_100_10_a, id.vars = "lam")
cv_m1_100_10_a$balance <- rep("1.00",nrow(cv_m1_100_10_a))
cv_m1_100_10_a$type <- rep("SPD-USS",nrow(cv_m1_100_10_a))

cv_m1_100_09_r <- sapply(lam_m1_100_09,function(res,col) res[,col],col=1)
cv_m1_100_09_r <- as.data.frame(cv_m1_100_09_r)
cv_m1_100_09_r$lam <- seq(-2,0.5,0.05)
cv_m1_100_09_r <- melt(cv_m1_100_09_r, id.vars = "lam")
cv_m1_100_09_r$balance <- rep("11.01",nrow(cv_m1_100_09_r))
cv_m1_100_09_r$type <- rep("SPD-MANOVA-W",nrow(cv_m1_100_09_r))


cv_m1_100_09_g <- sapply(lam_m1_100_09,function(res,col) res[,col],col=2)
cv_m1_100_09_g <- as.data.frame(cv_m1_100_09_g)
cv_m1_100_09_g$lam <- seq(-2,0.5,0.05)
cv_m1_100_09_g <- melt(cv_m1_100_09_g, id.vars = "lam")
cv_m1_100_09_g$balance <- rep("11.01",nrow(cv_m1_100_09_g))
cv_m1_100_09_g$type <- rep("SPD-MANOVA-B",nrow(cv_m1_100_09_g))

cv_m1_100_09_ga <- sapply(lam_m1_100_09,function(res,col) res[,col],col=3)
cv_m1_100_09_ga <- as.data.frame(cv_m1_100_09_ga)
cv_m1_100_09_ga$lam <- seq(-2,0.5,0.05)
cv_m1_100_09_ga <- melt(cv_m1_100_09_ga, id.vars = "lam")
cv_m1_100_09_ga$balance <- rep("11.01",nrow(cv_m1_100_09_ga))
cv_m1_100_09_ga$type <- rep("SPD-USS-B",nrow(cv_m1_100_09_ga))

cv_m1_100_09_a <- sapply(lam_m1_100_09,function(res,col) res[,col],col=4)
cv_m1_100_09_a <- as.data.frame(cv_m1_100_09_a)
cv_m1_100_09_a$lam <- seq(-2,0.5,0.05)
cv_m1_100_09_a <- melt(cv_m1_100_09_a, id.vars = "lam")
cv_m1_100_09_a$balance <- rep("11.01",nrow(cv_m1_100_09_a))
cv_m1_100_09_a$type <- rep("SPD-USS",nrow(cv_m1_100_09_a))

cv_m1_100_08_r <- sapply(lam_m1_100_08,function(res,col) res[,col],col=1)
cv_m1_100_08_r <- as.data.frame(cv_m1_100_08_r)
cv_m1_100_08_r$lam <- seq(-2,0.5,0.05)
cv_m1_100_08_r <- melt(cv_m1_100_08_r, id.vars = "lam")
cv_m1_100_08_r$balance <- rep("21.67",nrow(cv_m1_100_08_r))
cv_m1_100_08_r$type <- rep("SPD-MANOVA-W",nrow(cv_m1_100_08_r))


cv_m1_100_08_g <- sapply(lam_m1_100_08,function(res,col) res[,col],col=2)
cv_m1_100_08_g <- as.data.frame(cv_m1_100_08_g)
cv_m1_100_08_g$lam <- seq(-2,0.5,0.05)
cv_m1_100_08_g <- melt(cv_m1_100_08_g, id.vars = "lam")
cv_m1_100_08_g$balance <- rep("21.67",nrow(cv_m1_100_08_g))
cv_m1_100_08_g$type <- rep("SPD-MANOVA-B",nrow(cv_m1_100_08_g))

cv_m1_100_08_ga <- sapply(lam_m1_100_08,function(res,col) res[,col],col=3)
cv_m1_100_08_ga <- as.data.frame(cv_m1_100_08_ga)
cv_m1_100_08_ga$lam <- seq(-2,0.5,0.05)
cv_m1_100_08_ga <- melt(cv_m1_100_08_ga, id.vars = "lam")
cv_m1_100_08_ga$balance <- rep("21.67",nrow(cv_m1_100_08_ga))
cv_m1_100_08_ga$type <- rep("SPD-USS-B",nrow(cv_m1_100_08_ga))

cv_m1_100_08_a <- sapply(lam_m1_100_08,function(res,col) res[,col],col=4)
cv_m1_100_08_a <- as.data.frame(cv_m1_100_08_a)
cv_m1_100_08_a$lam <- seq(-2,0.5,0.05)
cv_m1_100_08_a <- melt(cv_m1_100_08_a, id.vars = "lam")
cv_m1_100_08_a$balance <- rep("21.67",nrow(cv_m1_100_08_a))
cv_m1_100_08_a$type <- rep("SPD-USS",nrow(cv_m1_100_08_a))

cv_m1_100_07_r <- sapply(lam_m1_100_07,function(res,col) res[,col],col=1)
cv_m1_100_07_r <- as.data.frame(cv_m1_100_07_r)
cv_m1_100_07_r$lam <- seq(-2,0.5,0.05)
cv_m1_100_07_r <- melt(cv_m1_100_07_r, id.vars = "lam")
cv_m1_100_07_r$balance <- rep("33.74",nrow(cv_m1_100_07_r))
cv_m1_100_07_r$type <- rep("SPD-MANOVA-W",nrow(cv_m1_100_07_r))


cv_m1_100_07_g <- sapply(lam_m1_100_07,function(res,col) res[,col],col=2)
cv_m1_100_07_g <- as.data.frame(cv_m1_100_07_g)
cv_m1_100_07_g$lam <- seq(-2,0.5,0.05)
cv_m1_100_07_g <- melt(cv_m1_100_07_g, id.vars = "lam")
cv_m1_100_07_g$balance <- rep("33.74",nrow(cv_m1_100_07_g))
cv_m1_100_07_g$type <- rep("SPD-MANOVA-B",nrow(cv_m1_100_07_g))

cv_m1_100_07_ga <- sapply(lam_m1_100_07,function(res,col) res[,col],col=3)
cv_m1_100_07_ga <- as.data.frame(cv_m1_100_07_ga)
cv_m1_100_07_ga$lam <- seq(-2,0.5,0.05)
cv_m1_100_07_ga <- melt(cv_m1_100_07_ga, id.vars = "lam")
cv_m1_100_07_ga$balance <- rep("33.74",nrow(cv_m1_100_07_ga))
cv_m1_100_07_ga$type <- rep("SPD-USS-B",nrow(cv_m1_100_07_ga))

cv_m1_100_07_a <- sapply(lam_m1_100_07,function(res,col) res[,col],col=4)
cv_m1_100_07_a <- as.data.frame(cv_m1_100_07_a)
cv_m1_100_07_a$lam <- seq(-2,0.5,0.05)
cv_m1_100_07_a <- melt(cv_m1_100_07_a, id.vars = "lam")
cv_m1_100_07_a$balance <- rep("33.74",nrow(cv_m1_100_07_a))
cv_m1_100_07_a$type <- rep("SPD-USS",nrow(cv_m1_100_07_a))

cv_m1_100_06_r <- sapply(lam_m1_100_06,function(res,col) res[,col],col=1)
cv_m1_100_06_r <- as.data.frame(cv_m1_100_06_r)
cv_m1_100_06_r$lam <- seq(-2,0.5,0.05)
cv_m1_100_06_r <- melt(cv_m1_100_06_r, id.vars = "lam")
cv_m1_100_06_r$balance <- rep("48.33",nrow(cv_m1_100_06_r))
cv_m1_100_06_r$type <- rep("SPD-MANOVA-W",nrow(cv_m1_100_06_r))


cv_m1_100_06_g <- sapply(lam_m1_100_06,function(res,col) res[,col],col=2)
cv_m1_100_06_g <- as.data.frame(cv_m1_100_06_g)
cv_m1_100_06_g$lam <- seq(-2,0.5,0.05)
cv_m1_100_06_g <- melt(cv_m1_100_06_g, id.vars = "lam")
cv_m1_100_06_g$balance <- rep("48.33",nrow(cv_m1_100_06_g))
cv_m1_100_06_g$type <- rep("SPD-MANOVA-B",nrow(cv_m1_100_06_g))

cv_m1_100_06_ga <- sapply(lam_m1_100_06,function(res,col) res[,col],col=3)
cv_m1_100_06_ga <- as.data.frame(cv_m1_100_06_ga)
cv_m1_100_06_ga$lam <- seq(-2,0.5,0.05)
cv_m1_100_06_ga <- melt(cv_m1_100_06_ga, id.vars = "lam")
cv_m1_100_06_ga$balance <- rep("48.33",nrow(cv_m1_100_06_ga))
cv_m1_100_06_ga$type <- rep("SPD-USS-B",nrow(cv_m1_100_06_ga))

cv_m1_100_06_a <- sapply(lam_m1_100_06,function(res,col) res[,col],col=4)
cv_m1_100_06_a <- as.data.frame(cv_m1_100_06_a)
cv_m1_100_06_a$lam <- seq(-2,0.5,0.05)
cv_m1_100_06_a <- melt(cv_m1_100_06_a, id.vars = "lam")
cv_m1_100_06_a$balance <- rep("48.33",nrow(cv_m1_100_06_a))
cv_m1_100_06_a$type <- rep("SPD-USS",nrow(cv_m1_100_06_a))

cv_m1_100_05_r <- sapply(lam_m1_100_05,function(res,col) res[,col],col=1)
cv_m1_100_05_r <- as.data.frame(cv_m1_100_05_r)
cv_m1_100_05_r$lam <- seq(-2,0.5,0.05)
cv_m1_100_05_r <- melt(cv_m1_100_05_r, id.vars = "lam")
cv_m1_100_05_r$balance <- rep("67.33",nrow(cv_m1_100_05_r))
cv_m1_100_05_r$type <- rep("SPD-MANOVA-W",nrow(cv_m1_100_05_r))


cv_m1_100_05_g <- sapply(lam_m1_100_05,function(res,col) res[,col],col=2)
cv_m1_100_05_g <- as.data.frame(cv_m1_100_05_g)
cv_m1_100_05_g$lam <- seq(-2,0.5,0.05)
cv_m1_100_05_g <- melt(cv_m1_100_05_g, id.vars = "lam")
cv_m1_100_05_g$balance <- rep("67.33",nrow(cv_m1_100_05_g))
cv_m1_100_05_g$type <- rep("SPD-MANOVA-B",nrow(cv_m1_100_05_g))

cv_m1_100_05_ga <- sapply(lam_m1_100_05,function(res,col) res[,col],col=3)
cv_m1_100_05_ga <- as.data.frame(cv_m1_100_05_ga)
cv_m1_100_05_ga$lam <- seq(-2,0.5,0.05)
cv_m1_100_05_ga <- melt(cv_m1_100_05_ga, id.vars = "lam")
cv_m1_100_05_ga$balance <- rep("67.33",nrow(cv_m1_100_05_ga))
cv_m1_100_05_ga$type <- rep("SPD-USS-B",nrow(cv_m1_100_05_ga))

cv_m1_100_05_a <- sapply(lam_m1_100_05,function(res,col) res[,col],col=4)
cv_m1_100_05_a <- as.data.frame(cv_m1_100_05_a)
cv_m1_100_05_a$lam <- seq(-2,0.5,0.05)
cv_m1_100_05_a <- melt(cv_m1_100_05_a, id.vars = "lam")
cv_m1_100_05_a$balance <- rep("67.33",nrow(cv_m1_100_05_a))
cv_m1_100_05_a$type <- rep("SPD-USS",nrow(cv_m1_100_05_a))

cv_m1_100_04_r <- sapply(lam_m1_100_04,function(res,col) res[,col],col=1)
cv_m1_100_04_r <- as.data.frame(cv_m1_100_04_r)
cv_m1_100_04_r$lam <- seq(-2,0.5,0.05)
cv_m1_100_04_r <- melt(cv_m1_100_04_r, id.vars = "lam")
cv_m1_100_04_r$balance <- rep("94.38",nrow(cv_m1_100_04_r))
cv_m1_100_04_r$type <- rep("SPD-MANOVA-W",nrow(cv_m1_100_04_r))


cv_m1_100_04_g <- sapply(lam_m1_100_04,function(res,col) res[,col],col=2)
cv_m1_100_04_g <- as.data.frame(cv_m1_100_04_g)
cv_m1_100_04_g$lam <- seq(-2,0.5,0.05)
cv_m1_100_04_g <- melt(cv_m1_100_04_g, id.vars = "lam")
cv_m1_100_04_g$balance <- rep("94.38",nrow(cv_m1_100_04_g))
cv_m1_100_04_g$type <- rep("SPD-MANOVA-B",nrow(cv_m1_100_04_g))

cv_m1_100_04_ga <- sapply(lam_m1_100_04,function(res,col) res[,col],col=3)
cv_m1_100_04_ga <- as.data.frame(cv_m1_100_04_ga)
cv_m1_100_04_ga$lam <- seq(-2,0.5,0.05)
cv_m1_100_04_ga <- melt(cv_m1_100_04_ga, id.vars = "lam")
cv_m1_100_04_ga$balance <- rep("94.38",nrow(cv_m1_100_04_ga))
cv_m1_100_04_ga$type <- rep("SPD-USS-B",nrow(cv_m1_100_04_ga))

cv_m1_100_04_a <- sapply(lam_m1_100_04,function(res,col) res[,col],col=4)
cv_m1_100_04_a <- as.data.frame(cv_m1_100_04_a)
cv_m1_100_04_a$lam <- seq(-2,0.5,0.05)
cv_m1_100_04_a <- melt(cv_m1_100_04_a, id.vars = "lam")
cv_m1_100_04_a$balance <- rep("94.38",nrow(cv_m1_100_04_a))
cv_m1_100_04_a$type <- rep("SPD-USS",nrow(cv_m1_100_04_a))

cv_m1_100_03_r <- sapply(lam_m1_100_03,function(res,col) res[,col],col=1)
cv_m1_100_03_r <- as.data.frame(cv_m1_100_03_r)
cv_m1_100_03_r$lam <- seq(-2,0.5,0.05)
cv_m1_100_03_r <- melt(cv_m1_100_03_r, id.vars = "lam")
cv_m1_100_03_r$balance <- rep("137.84",nrow(cv_m1_100_03_r))
cv_m1_100_03_r$type <- rep("SPD-MANOVA-W",nrow(cv_m1_100_03_r))


cv_m1_100_03_g <- sapply(lam_m1_100_03,function(res,col) res[,col],col=2)
cv_m1_100_03_g <- as.data.frame(cv_m1_100_03_g)
cv_m1_100_03_g$lam <- seq(-2,0.5,0.05)
cv_m1_100_03_g <- melt(cv_m1_100_03_g, id.vars = "lam")
cv_m1_100_03_g$balance <- rep("137.84",nrow(cv_m1_100_03_g))
cv_m1_100_03_g$type <- rep("SPD-MANOVA-B",nrow(cv_m1_100_03_g))

cv_m1_100_03_ga <- sapply(lam_m1_100_03,function(res,col) res[,col],col=3)
cv_m1_100_03_ga <- as.data.frame(cv_m1_100_03_ga)
cv_m1_100_03_ga$lam <- seq(-2,0.5,0.05)
cv_m1_100_03_ga <- melt(cv_m1_100_03_ga, id.vars = "lam")
cv_m1_100_03_ga$balance <- rep("137.84",nrow(cv_m1_100_03_ga))
cv_m1_100_03_ga$type <- rep("SPD-USS-B",nrow(cv_m1_100_03_ga))

cv_m1_100_03_a <- sapply(lam_m1_100_03,function(res,col) res[,col],col=4)
cv_m1_100_03_a <- as.data.frame(cv_m1_100_03_a)
cv_m1_100_03_a$lam <- seq(-2,0.5,0.05)
cv_m1_100_03_a <- melt(cv_m1_100_03_a, id.vars = "lam")
cv_m1_100_03_a$balance <- rep("137.84",nrow(cv_m1_100_03_a))
cv_m1_100_03_a$type <- rep("SPD-USS",nrow(cv_m1_100_03_a))



cv_m1_100_r <- rbind(cv_m1_100_10_r,cv_m1_100_07_r,cv_m1_100_04_r)
cv_m1_100_r$balance <- factor(cv_m1_100_r$balance,levels = c("1.00","33.74","94.38"),
                             labels = c(TeX("$max_i\\,n_i/n_0=1.00$"),
                                        TeX("$max_i\\,n_i/n_0=33.74$"),
                                        TeX("$max_i\\,n_i/n_0=94.38$"))
)
cv_m1_100_r$variable <- as.character(cv_m1_100_r$variable)

cv_m1_100_g <- rbind(cv_m1_100_10_g,cv_m1_100_07_g,cv_m1_100_04_g,
                     cv_m1_100_10_ga,cv_m1_100_07_ga,cv_m1_100_04_ga,
                     cv_m1_100_10_a,cv_m1_100_07_a,cv_m1_100_04_a)
cv_m1_100_g$balance <- factor(cv_m1_100_g$balance,levels = c("1.00","33.74","94.38"),
                              labels = c(TeX("$max_i\\,n_i/n_0=1.00$"),
                                         TeX("$max_i\\,n_i/n_0=33.74$"),
                                         TeX("$max_i\\,n_i/n_0=94.38$"))
)
cv_m1_100_g$variable <- as.character(cv_m1_100_g$variable)


cv_m1_100_r_min <- cv_m1_100_r %>% group_by(variable,balance,type) %>% slice(which.min(value))
cv_m1_100_r_min$variable <- as.character(cv_m1_100_r_min$variable)
cv_m1_100_g_min <- cv_m1_100_g %>% group_by(variable,balance,type) %>% slice(which.min(value))
cv_m1_100_g_min$variable <- as.character(cv_m1_100_g_min$variable)




da_roc <- read.csv("da_roc.csv")
da_roc$balance <- as.character(da_roc$balance)
da_roc$balance <- ifelse(da_roc$balance == "1","1.00",da_roc$balance)
da_cv <- read.csv("da_cv.csv")
da_cv$balance <- as.character(da_cv$balance)
da_cv$balance <- ifelse(da_cv$balance == "1","1.00",da_cv$balance)

newtype = c()
for(x in da_cv$type){
        if(x == 'Within Subject'){
                newtype = c(newtype,'SPD-MANOVA-W')
        }else if(x == 'Between Subject 1'){
                newtype = c(newtype,'SPD-MANOVA-B')
        }else if(x == 'Between Subject 2'){
                newtype = c(newtype,'SPD-USS-B')
        }else{
                newtype = c(newtype,'SPD-USS')
        }
}
da_cv$type = newtype

newtype = c()
for(x in da_roc$type){
        if(x == 'Within Subject'){
                newtype = c(newtype,'SPD-MANOVA-W')
        }else if(x == 'Between Subject 1'){
                newtype = c(newtype,'SPD-MANOVA-B')
        }else if(x == 'Between Subject 2'){
                newtype = c(newtype,'SPD-USS-B')
        }else{
                newtype = c(newtype,'SPD-USS')
        }
}
da_roc$type = newtype

###### roc curve
da_roc_sub4 <- subset(da_roc, group<=10 & balance %in% c("1.00","33.74","94.38") & Model == "Model 1" & p == "p=100")
da_cv_sub4 <- subset(da_cv,group<=10 & balance %in% c("1.00","33.74","94.38") & Model == "Model 1" & p == "p=100")

da_roc_sub4$balance <- factor(da_roc_sub4$balance,levels = c("1.00","33.74","94.38"),
                              labels = c(TeX("$max_i\\,n_i/n_0=1.00$"),
                                         TeX("$max_i\\,n_i/n_0=33.74$"),
                                         TeX("$max_i\\,n_i/n_0=94.38$"))
)

da_cv_sub4$balance <- factor(da_cv_sub4$balance,levels = c("1.00","33.74","94.38"),
                             labels = c(TeX("$max_i\\,n_i/n_0=1.00$"),
                                        TeX("$max_i\\,n_i/n_0=33.74$"),
                                        TeX("$max_i\\,n_i/n_0=94.38$"))
)





#brewer.pal(8,"Dark2")
#[1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
pr <- ggplot()+
        geom_line(data = filter(cv_m1_100_r,variable %in% paste('V',1:10,sep="")),
                  aes(x=lam,y=value,
                      group=variable),size=0.45,color="#E7298A")+#,alpha = 0.5)+
        facet_wrap(~balance,labeller=label_parsed,ncol=1) +
        
        geom_point(data = filter(cv_m1_100_r_min,variable %in% paste('V',1:10,sep="")),
                   mapping = aes(x=lam,y=value,
                                 group=variable),size=3.5,color="#E7298A",shape=18)+
        #xlim(0,1)+
        #ylim(0,1)+
        xlab(TeX("$log_{10}(\\lambda)$")) +
        #xlab(TeX("Log Lambda")) +
        ylab("Cross-Validation Criteria") +
        #facet_grid(cols=vars(p),rows=vars(Model)) +
        scale_color_brewer(palette="Dark2")+
        #guides(shape=guide_legend("Type of Correlations"),
        #       color=guide_legend("Type of Correlations"))+
        ggtitle(TeX("$\\Sigma_{\\epsilon}$: Cross-Validation Curves"))+
        theme_bw()+
        theme(plot.title = element_text(size=15,hjust = 0.5,family = "Times"),
              strip.text.x = element_text(size=15,family = "Times"),
              strip.text.y = element_text(size=15,family = "Times"),
              strip.background =element_rect(fill="white"),
              legend.position = "none",
              #legend.title=element_text(size=15),
              #legend.text=element_text(size=15),
              #legend.box = "horizontal",
              #legend.background = element_rect(fill = "white", color = "black"),
              axis.text=element_text(size=12,family = "Times"),
              axis.title=element_text(size=15,family = "Times"))
pr


pg <- ggplot()+
        geom_line(data = filter(cv_m1_100_g,variable %in% paste('V',1:10,sep="")),
                  aes(x=lam,y=value,
                      group=interaction(variable,type),
                      color = factor(type,levels = paste("SPD-",c("USS","USS-B","MANOVA-B"),sep = ""))),size=0.45)+#,alpha = 0.5)+
        facet_wrap(~balance,labeller=label_parsed,ncol=1,scales = "free_y") +
        
        geom_point(data = filter(cv_m1_100_g_min,variable %in% paste('V',1:10,sep="")),
                   mapping = aes(x=lam,y=value,
                                 group=interaction(variable,type),
                                 color = factor(type,levels = paste("SPD-",c("USS","USS-B","MANOVA-B"),sep = "")),
                                 shape=factor(type,levels = paste("SPD-",c("USS","USS-B","MANOVA-B"),sep = ""))),size=3.5)+
        #xlim(0,1)+
        #ylim(0,1)+
        xlab(TeX("$log_{10}(\\lambda)$")) +
        #xlab(TeX("Log Lambda")) +
        ylab("Cross-Validation Criteria") +
        #facet_grid(cols=vars(p),rows=vars(Model)) +
        scale_color_manual(values=c("#1B9E77","#D95F02","#7570B3"))+
        scale_shape_manual(values=c(15,16,17))+
        #guides(shape=guide_legend("Type of Correlations"),
        #       color=guide_legend("Type of Correlations"))+
        ggtitle(TeX("$\\Sigma_b$: Cross-Validation Curves"))+
        theme_bw()+
        theme(plot.title = element_text(size=15,hjust = 0.5,family = "Times"),
              strip.text.x = element_text(size=15,family = "Times"),
              strip.text.y = element_text(size=15,family = "Times"),
              strip.background =element_rect(fill="white"),
              legend.position = "none",
              #legend.title=element_text(size=15),
              #legend.text=element_text(size=15),
              #legend.box = "horizontal",
              #legend.background = element_rect(fill = "white", color = "black"),
              axis.text=element_text(size=12,family = "Times"),
              axis.title=element_text(size=15,family = "Times"))
pg

pc <- ggplot()+
        geom_line(data = da_roc_sub4,
                  aes(x=fpr,y=tpr,
                      group=interaction(group,type),
                      color = factor(type,levels = paste("SPD-",c("USS","USS-B","MANOVA-B","MANOVA-W"),sep = ""))),
                  size=0.45)+#,alpha = 0.5)+
        facet_wrap(~balance,labeller=label_parsed,ncol=1) +
        geom_abline(intercept = 0, slope = 1, linetype="dashed", color = "gray30",
                    size = 0.35) +
        geom_point(data = da_cv_sub4,
                   mapping = aes(x=fpr,y=tpr,
                                 group=interaction(group,type),
                                 color = factor(type,levels = paste("SPD-",c("USS","USS-B","MANOVA-B","MANOVA-W"),sep = "")),
                                 shape=factor(type,levels = paste("SPD-",c("USS","USS-B","MANOVA-B","MANOVA-W"),sep = ""))),size=3.5)+
        xlim(0,1)+
        ylim(0,1)+
        xlab("False Postive Rate") +
        ylab("True Positive Rate") +
        #facet_grid(cols=vars(p),rows=vars(Model)) +
        scale_color_brewer(palette="Dark2",
                           labels = c(TeX("$\\widehat{\\Sigma}^+$"),
                                      TeX("$\\widehat{\\Sigma}^+_b$"),
                                      TeX("$\\widetilde{\\Sigma}^+_b$"),
                                      TeX("$\\widehat{\\Sigma}^+_{\\epsilon}$")))+
        scale_shape_manual(values=c(15,16,17,18),
                           labels = c(TeX("$\\widehat{\\Sigma}^+$"),
                                      TeX("$\\widehat{\\Sigma}^+_b$"),
                                      TeX("$\\widetilde{\\Sigma}^+_b$"),
                                      TeX("$\\widehat{\\Sigma}^+_{\\epsilon}$")))+
        guides(shape=guide_legend("Type of Correlations"),
               color=guide_legend("Type of Correlations"))+
        ggtitle(TeX("ROC Curves for $\\Sigma_b$ and $\\Sigma_{\\epsilon}$"))+
        theme_bw()+
        theme(plot.title = element_text(size=15,hjust = 0.5,family = "Times"),
              strip.text.x = element_text(size=15,family = "Times"),
              strip.text.y = element_text(size=15,family = "Times"),
              strip.background =element_rect(fill="white"),
              legend.position = "none",#c(0.75, 0.06),
              #legend.title=element_text(size=15),
              #legend.text=element_text(size=15),
              #legend.box = "horizontal",
              #legend.background = element_rect(fill = "white", color = "black"),
              axis.text=element_text(size=12,family = "Times"),
              axis.title=element_text(size=15,family = "Times"))
pc





pall<-plot_grid(pr,pg,pc, align = "v",ncol=3)
pall

setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript("mm_roc2.eps", family="Times", width = 10, height = 9)
pall
dev.off()



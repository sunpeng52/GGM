rm(list = ls())

library(latex2exp)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

library(dplyr)
library(spcov)
library(MASS)

setwd("/Users/duanduan/Desktop/Research/Group Graphical Models - Covariance/04 GGM Coding")


source("functions.R")
source("functions_2.R")
options(dplyr.summarise.inform = FALSE)

data <- read.csv("mm_comparison.csv")

data <- data[data$a != min(data$a),]
data <- data[data$scenario == "Scenario 1",]
data$estimator <- ifelse(data$estimator != "Moment","Regularized Estimator","Moment Estimator")
data$type <- ifelse(data$type == "Between Subject", "Our Estimator", "USS")



#ggsave("mm_comparison.pdf", width = 12, height = 6)


p <-  100
R.true <- banded.new(p,10,-1)
G.true <- banded.new(p,10)

load("data_mm_100.RData")

i = 1
mom = anova.est2(da_m1_100_03[[i]])


##### get moment estimation
f_norm_03 = c()
for(i in 1:100){
        mom = anova.est2(da_m1_100_03[[i]])
        f_norm_03 = rbind(f_norm_03,c(norm(mom$G - G.true,"F"), norm(mom$Ga - G.true,"F")))
}

f_norm_04 = c()
for(i in 1:100){
        mom = anova.est2(da_m1_100_04[[i]])
        f_norm_04 = rbind(f_norm_04,c(norm(mom$G - G.true,"F"), norm(mom$Ga - G.true,"F")))
        print(i)
}

f_norm_05 = c()
for(i in 1:100){
        mom = anova.est2(da_m1_100_05[[i]])
        f_norm_05 = rbind(f_norm_05,c(norm(mom$G - G.true,"F"), norm(mom$Ga - G.true,"F")))
        print(i)
}

f_norm_06 = c()
for(i in 1:100){
        mom = anova.est2(da_m1_100_06[[i]])
        f_norm_06 = rbind(f_norm_06,c(norm(mom$G - G.true,"F"), norm(mom$Ga - G.true,"F")))
        print(i)
}

f_norm_07 = c()
for(i in 1:100){
        mom = anova.est2(da_m1_100_07[[i]])
        f_norm_07 = rbind(f_norm_07,c(norm(mom$G - G.true,"F"), norm(mom$Ga - G.true,"F")))
        print(i)
}

f_norm_08 = c()
for(i in 1:100){
        mom = anova.est2(da_m1_100_08[[i]])
        f_norm_08 = rbind(f_norm_08,c(norm(mom$G - G.true,"F"), norm(mom$Ga - G.true,"F")))
        print(i)
}

f_norm_09 = c()
for(i in 1:100){
        mom = anova.est2(da_m1_100_09[[i]])
        f_norm_09 = rbind(f_norm_09,c(norm(mom$G - G.true,"F"), norm(mom$Ga - G.true,"F")))
        print(i)
}

f_norm_10 = c()
for(i in 1:100){
        mom = anova.est2(da_m1_100_10[[i]])
        f_norm_10 = rbind(f_norm_10,c(norm(mom$G - G.true,"F"), norm(mom$Ga - G.true,"F")))
        print(i)
}


save(f_norm_03,f_norm_04,f_norm_05,
     f_norm_06,f_norm_07,f_norm_08,
     f_norm_09,f_norm_10, file="da_mm_100_f_norm.RData")
load("da_mm_100_f_norm.RData")

f_norm <- rbind(
        colMeans(f_norm_03),
        colMeans(f_norm_04),
        colMeans(f_norm_05),
        colMeans(f_norm_06),
        colMeans(f_norm_07),
        colMeans(f_norm_08),
        colMeans(f_norm_09),
        colMeans(f_norm_10)
)


load("mm_simulation_100_res22.RData")


g.g.error_m1 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100_03)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100_04)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100_05)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100_06)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100_07)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100_08)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100_09)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_100_10)))

g.ga.error_m1 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_100_03)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_100_04)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_100_05)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_100_06)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_100_07)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_100_08)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_100_09)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_100_10)))

balance = c()
for(x in 3:10){
        N = 1000
        n = c(rep(x,99),N-99*x)
        n_eff = (N-sum(n^2)/N)/(length(n)-1)
        balance = c(balance,max(n)/n_eff)
        
}



brewer.pal(4, "Dark2")
#"#1B9E77" "#D95F02" "#7570B3" "#E7298A"


error2 <- data.frame(
        im = rep(balance,4),
        error = c(f_norm[,1],f_norm[,2],g.g.error_m1,g.ga.error_m1),
        type = rep(c("MANOVA-B","Our Estimator"),each=8,time=2),
        estimator = rep(c("Moment Estimator","Regularized Estimator"), each=16)
)

f_norm[,2]

library(grid)
library(gridExtra)


p1 <- ggplot(data = data,
             aes(x=a,y=error,group=interaction(estimator,type),
                 color = factor(type,levels = c("USS","Our Estimator")),
                 shape = factor(type,levels = c("USS","Our Estimator")),
                 linetype=estimator))+
        geom_point(size=3)+
        geom_line(size=0.6)+
        #facet_wrap(~scenario)+
        #facet_grid(cols=vars(estimator)) +
        scale_color_manual(values=c("#1B9E77","#D95F02"))+
        scale_shape_manual(values=c(15,16))+
        guides(shape=guide_legend("Type of Covariances"),
               color=guide_legend("Type of Covariances"),
               linetype=guide_legend("Type of Estimators")
               )+
        ylim(0,45)+
        #ylab(TeX("$||\\Sigma^{est}-\\Sigma^{0}_b||_F$")) +
        #xlab(TeX("$|\\Sigma^{0}_{\\epsilon}|_{\\infty}/|\\Sigma^{0}_{b}|_{\\infty}$")) +
        ylab(TeX("Estimation Error")) +
        xlab(TeX("Inverse Signal-to-Noise Ratio")) +
        theme_bw()+
        theme(#plot.title = element_text(size=18),
                strip.text.x = element_text(size = 20),
                strip.background =element_rect(fill="white"),
                legend.position = c(0.37, 0.88),
                legend.title=element_text(size=16),
                legend.text=element_text(size=14),
                legend.box = "horizontal",
                legend.background = element_rect(fill = "white", color = "black"),
                axis.text=element_text(size=18),
                axis.title=element_text(size=20))
p1

p2 <- ggplot(data = error2,
             aes(x=im,y=error,group=interaction(type,estimator),
                 color = factor(type,levels = c("Our Estimator","MANOVA-B")),
                 shape = factor(type,levels = c("Our Estimator","MANOVA-B")),
                 linetype=estimator))+
        geom_point(size=3)+
        geom_line(size=0.6)+
        #facet_grid(cols=vars(estimator)) +
        scale_color_manual(values=c("#D95F02","#7570B3"),
                           labels = c(TeX("Our Estimator"),
                                      TeX("$MANOVA_b$")))+
        scale_shape_manual(values=c(16,17),
                           labels = c(TeX("Our Estimator"),
                                      TeX("$MANOVA_b$")))+
        guides(shape=guide_legend("Type of Covariances"),
               color=guide_legend("Type of Covariances"),
               linetype=guide_legend("Type of Estimators"))+
        ylim(0,45)+
        ylab(TeX("Estimation Error")) +
        xlab(TeX("Imbalance")) +
        theme_bw()+
        theme(#plot.title = element_text(size=18),
                strip.text.x = element_text(size = 20),
                strip.background =element_rect(fill="white"),
                legend.position = c(0.37, 0.88),
                legend.title=element_text(size=16),
                legend.text=element_text(size=14),
                legend.box = "horizontal",
                legend.background = element_rect(fill = "white", color = "black"),
                axis.text=element_text(size=18),
                axis.title=element_text(size=20))
p2


library(cowplot)
pp <- list(p1, p2)
plot_grid(plotlist=pp, ncol=2, align='h')

ggsave("mm_comparison_new.pdf", width = 15, height = 6)

###grid.arrange(gt1, gt2, ncol=2)


ggarrange(p1,p2,ncol = 2)


factor(c("MANOVA-B","MANOVA-W","USS-B","USS"),levels = c("USS","USS-B","MANOVA-B","MANOVA-W"))



############################
library(extrafont)

data <- read.csv("mm_comparison.csv")
data <- data[data$a != min(data$a),]
#data$scenario <- ifelse(data$scenario == "Scenario 2", "Model 3", "Model 4")
data$estimator <- ifelse(data$estimator != "Moment","Regularized Estimator","Moment Estimator")
data$type <- ifelse(data$type == "Between Subject", "Our Estimator", "USS")

p3 <- ggplot(data = data,
             aes(x=a,y=error,group=interaction(estimator,type),
                 color = factor(type,levels = c("USS","Our Estimator")),
                 shape = factor(type,levels = c("USS","Our Estimator")),
                 linetype=estimator))+
        geom_point(size=3)+
        geom_line(size=0.6)+
        facet_wrap(~scenario)+
        #facet_grid(cols=vars(estimator)) +
        scale_color_manual(values=c("#1B9E77","#D95F02"))+
        scale_shape_manual(values=c(15,16))+
        #guides(shape=guide_legend("Type of Covariances"),
        #       color=guide_legend("Type of Covariances"),
        #       linetype=guide_legend("Type of Estimators")
        #)+
        scale_x_continuous(breaks = c(1.0,2.5,5.0,7.5,10.0))+
        ylim(0,45)+
        ylab(TeX("$||\\Sigma^{est}-\\Sigma^{0}_b||_F$")) +
        xlab(TeX("$|\\Sigma^{0}_{\\epsilon}|_{\\infty}/|\\Sigma^{0}_{b}|_{\\infty}$")) +
        theme_bw()+
        theme(#plot.title = element_text(size=18),
                plot.margin=unit(rep(0.5,4), 'cm'),
                strip.text.x = element_text(size = 20,family = "Times"),
                strip.background =element_rect(fill="white"),
                legend.position = "none",
                #legend.position = c(0.18, 0.88),
                #legend.title=element_text(size=16),
                #legend.text=element_text(size=14),
                #legend.box = "horizontal",
                #legend.background = element_rect(fill = "white", color = "black"),
                axis.text=element_text(size=18,family = "Times"),
                axis.title=element_text(size=20,family = "Times"))
p3

setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript("mm_comparison.eps", family="Times", width = 9, height = 5)
p3
dev.off()

#ggsave("mm_comparison.eps", width = 15, height = 6)

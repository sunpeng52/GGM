rm(list =ls())

library(latex2exp)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

source("functions.R")
source("functions_2.R")
options(dplyr.summarise.inform = FALSE)

load("mm_simulation_100_res22.RData")

p <- 100
G.true <-  banded.m(p,10)
R.true <- banded.new(p,10,-1)


g.a.error_m1 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_100_03)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_100_04)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_100_05)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_100_06)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_100_07)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_100_08)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_100_09)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_100_10)))

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


r.a.error_m1 = c(mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_100_03)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_100_04)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_100_05)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_100_06)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_100_07)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_100_08)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_100_09)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_100_10)))

r.r.error_m1 = c(mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_100_03)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_100_04)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_100_05)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_100_06)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_100_07)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_100_08)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_100_09)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_100_10)))




p <-  100
R.true <- ar.cov(p,-0.6)
G.true <- ar.cov(p,0.6)


g.a.error_m2 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_100_03)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_100_04)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_100_05)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_100_06)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_100_07)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_100_08)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_100_09)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_100_10)))

g.g.error_m2 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100_03)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100_04)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100_05)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100_06)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100_07)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100_08)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100_09)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_100_10)))

g.ga.error_m2 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_100_03)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_100_04)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_100_05)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_100_06)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_100_07)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_100_08)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_100_09)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_100_10)))


r.a.error_m2 = c(mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_100_03)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_100_04)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_100_05)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_100_06)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_100_07)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_100_08)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_100_09)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_100_10)))

r.r.error_m2 = c(mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_100_03)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_100_04)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_100_05)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_100_06)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_100_07)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_100_08)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_100_09)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_100_10)))
balance = c()
for(x in 3:10){
        N = 1000
        n = c(rep(x,99),N-99*x)
        n_eff = (N-sum(n^2)/N)/(length(n)-1)
        balance = c(balance,max(n)/n_eff)
        
}

mm_res1 <- data.frame(
        im = rep(balance,10),
        error = c(g.a.error_m1,g.g.error_m1,g.ga.error_m1,r.r.error_m1,r.a.error_m1,
                  g.a.error_m2,g.g.error_m2,g.ga.error_m2,r.r.error_m2,r.a.error_m2),
        model = rep(c('Model 1', 'Model 2'),each = 5*8),
        type = rep(c("SPD-USS","SPD-MANOVA-B", "SPD-USS-B","SPD-MANOVA-W","SPD-USS"),each=8,time=2),
        comparison = rep(c(rep("Between Subject",3),rep("Within Subject",2)),each=8,time=2)
)
mm_res1$p = rep("p=100",nrow(mm_res1))

load("mm_simulation_200_res22.RData")

p <- 200
G.true <-  banded.m(p,10)
R.true <- banded.new(p,10,-1)


g.a.error_m1 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_200_03)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_200_04)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_200_05)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_200_06)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_200_07)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_200_08)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_200_09)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m1_200_10)))

g.g.error_m1 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200_03)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200_04)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200_05)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200_06)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200_07)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200_08)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200_09)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m1_200_10)))

g.ga.error_m1 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_200_03)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_200_04)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_200_05)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_200_06)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_200_07)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_200_08)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_200_09)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m1_200_10)))


r.a.error_m1 = c(mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_200_03)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_200_04)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_200_05)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_200_06)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_200_07)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_200_08)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_200_09)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m1_200_10)))

r.r.error_m1 = c(mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_200_03)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_200_04)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_200_05)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_200_06)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_200_07)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_200_08)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_200_09)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m1_200_10)))




p <-  200
R.true <- ar.cov(p,-0.6)
G.true <- ar.cov(p,0.6)


g.a.error_m2 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_200_03)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_200_04)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_200_05)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_200_06)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_200_07)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_200_08)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_200_09)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=A_m2_200_10)))

g.g.error_m2 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200_03)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200_04)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200_05)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200_06)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200_07)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200_08)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200_09)),
                 mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=G_m2_200_10)))

g.ga.error_m2 = c(mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_200_03)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_200_04)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_200_05)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_200_06)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_200_07)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_200_08)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_200_09)),
                  mean(mapply(function(res) norm(res$Theta - G.true,"F"),res=Ga_m2_200_10)))


r.a.error_m2 = c(mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_200_03)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_200_04)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_200_05)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_200_06)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_200_07)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_200_08)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_200_09)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=A_m2_200_10)))

r.r.error_m2 = c(mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_200_03)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_200_04)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_200_05)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_200_06)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_200_07)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_200_08)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_200_09)),
                 mean(mapply(function(res) norm(res$Theta - R.true,"F"),res=R_m2_200_10)))
balance = c()
for(x in 3:10){
        N = 1000
        n = c(rep(x,99),N-99*x)
        n_eff = (N-sum(n^2)/N)/(length(n)-1)
        balance = c(balance,max(n)/n_eff)
        
}



mm_res2 <- data.frame(
        im = rep(balance,10),
        error = c(g.a.error_m1,g.g.error_m1,g.ga.error_m1,r.r.error_m1,r.a.error_m1,
                  g.a.error_m2,g.g.error_m2,g.ga.error_m2,r.r.error_m2,r.a.error_m2),
        model = rep(c('Model 1', 'Model 2'),each = 5*8),
        type = rep(c("SPD-USS","SPD-MANOVA-B", "SPD-USS-B","SPD-MANOVA-W","SPD-USS"),each=8,time=2),
        comparison = rep(c(rep("Between Subject",3),rep("Within Subject",2)),each=8,time=2)
)

mm_res2$p = rep("p=200",nrow(mm_res2))


mm_res <- rbind(mm_res1,mm_res2)
#mm_res <- read.csv("mm_res.csv")

#mm_res2 <- filter(mm_res, !(type == "SPD-USS" & comparison == "Within Subject"))

#mm_res$model <- ifelse(mm_res$model == "Model 1", "Example 1", "Example 2")

library(facetscales)
scales_y <- list(
        `Model 1` = scale_y_continuous(limits = c(0, 46), breaks = seq(0, 46, 10)),
        `Model 2` = scale_y_continuous(limits = c(0, 22), breaks = seq(0, 22, 5))
)

error_f1 <- ggplot(data=mm_res, aes(x=im,y=error,group=interaction(model,type,comparison),
                                      color = factor(type,levels = paste("SPD-",c("USS","USS-B","MANOVA-B","MANOVA-W"),sep = "")),
                                      shape = factor(type,levels = paste("SPD-",c("USS","USS-B","MANOVA-B","MANOVA-W"),sep = "")),
                                      linetype = comparison))+
        geom_point(size = 4)+
        geom_line(size = 0.8)+
        facet_grid_sc(vars(model),vars(p),scales = list(y=scales_y))+
        #ylim(0,35)+
        scale_x_continuous(breaks=c(1,50, 100, 150),limits = c(1,150))+
        xlab(TeX("$max_i~n_i/n_0$"))+
        ylab(TeX("$||\\Sigma^{est}-\\Sigma^{0}||_F$"))+
        guides(shape=guide_legend("Type of Correlations"),
               color=guide_legend("Type of Correlations"),
               linetype=guide_legend("Comparison"))+
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
        theme_bw()+
        theme(#plot.title = element_text(size=18),
                plot.margin=unit(rep(0.5,4), 'cm'),
                strip.text.x = element_text(size = 20,family = "Times"),
                strip.text.y = element_text(size = 20,family = "Times"),
                strip.background =element_rect(fill="white"),
                legend.position = 'none',#c(0.21, 0.92),
                #legend.title=element_text(size=16),
                #legend.text=element_text(size=14),
                #legend.box = "horizontal",
                #legend.background = element_rect(fill = "white", color = "black"),
                axis.text=element_text(size=18,family = "Times"),
                axis.title=element_text(size=20,family = "Times"))
error_f1

setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")
postscript("mm_error.eps", family="Times", width = 9, height = 9)
error_f1
dev.off()

balance = c()
n.star= c()
for(x in 3:10){
        N = 1000
        n = c(rep(x,99),N-99*x)
        n_eff = (N-sum(n^2)/N)/(length(n)-1)
        balance = c(balance,max(n)/n_eff)
        n.star = c(n.star,mean(1/n))
        
}
(2-n.star)/n.star
balance

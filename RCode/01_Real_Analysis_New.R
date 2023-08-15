rm(list = ls())

# Load packages: 
library(dplyr)
library(tidyr)
library(data.table)
library(igraph)
library(latex2exp)
library(spcov)

library(ggplot2)
library(ggpubr)
library(gridExtra)

getwd()
setwd("/Users/duanduan/Desktop/Research/Group Graphical Models - Covariance/04 GGM Coding/Real Data2")

source("/Users/duanduan/Desktop/Research/Group Graphical Models - Covariance/04 GGM Coding/functions.R")
source("/Users/duanduan/Desktop/Research/Group Graphical Models - Covariance/04 GGM Coding/functions_2.R")
options(dplyr.summarise.inform = FALSE)

## read the data:
treatment <- fread("NIH_submitted_data/rad_016_067-01_DATA_Treatment_v001.csv")
patient <- fread("NIH_submitted_data/rad_016_067-01_DATA_PatientInfo_v001.csv")
labs <- fread("NIH_submitted_data/rad_016_067-01_DATA_Labs_v001.csv")

## join 3 tables with the random variable of interest
treat_sub <- treatment %>%
        dplyr::select(
                patient_id, FDD_to_txt_dt,
                pre_weight_lbs,post_weight_lbs, idwg,
                min_sbp, max_sbp, min_dbp, max_dbp, min_pulse, max_pulse
        ) %>%
        mutate(pre_weight = pre_weight_lbs*0.4536,
               post_weight = post_weight_lbs*0.4536) %>%
        dplyr::select(-pre_weight_lbs,-post_weight_lbs)

patient_sub <- patient %>%
        dplyr::select(patient_id,AgeAtFDD,gender,race,ethnicity,vintage,diabetes)

labs_sub <- labs %>% 
        dplyr::select(patient_id, FDD_to_draw_dt, covid_pcr) %>%
        mutate(FDD_to_txt_dt = FDD_to_draw_dt)

da <- patient_sub %>%
        right_join(treat_sub, by = "patient_id") %>%
        left_join(labs_sub, by = c("patient_id","FDD_to_txt_dt")) %>%
        mutate(ufv = pre_weight - post_weight)

## remove covid positive patients, data from the first year
tmp <- filter(da, covid_pcr == "Positive")
selid <- unique(tmp$patient_id)
da1 <- filter(da, !(patient_id %in% selid)) %>%
        filter(FDD_to_txt_dt>365)

## select male, white, non-hispanic, non-diabetic
da2 <- filter(da1, (gender==1)&(race==5)&(ethnicity==0)&(diabetes==0))

## select several variables, remove duplicates, remove missing values
da3 <- da2 %>%
        dplyr::select(
                patient_id,FDD_to_txt_dt,
                pre_weight, post_weight, ufv, idwg,
                min_sbp, max_sbp, min_dbp, max_dbp, min_pulse, max_pulse
        ) %>%
        unique() %>%
        na.omit()

summary(da3)

### comments: weight_loss has many outliers
### check weight
tmp <- da3 %>% 
        group_by(patient_id) %>%
        dplyr::summarise(change=(max(pre_weight,na.rm = T)-min(pre_weight,na.rm = T))/mean(pre_weight,na.rm = T))
summary(tmp$change)
hist(tmp$change)

mean(tmp$change > 0.1)
mean(tmp$change > 0.2)
mean(tmp$change > 0.5)

#patient05 = tmp[tmp$change > 0.5,]
#plist_05 = list()
#for(i in 1:nrow(patient05)){
#        patient <- filter(da3, patient_id == as.integer(patient05[i,1]))
#        plist_05[[i]] <- ggplot(patient, aes(x=FDD_to_txt_dt, y=pre_weight)) +
#                geom_line()+
#                xlab("FDD to Txt Date Time")+
#                ylab("Predialysis Weight") +
#                ggtitle(paste("Patient ID: ", patient05[i,1], 
#                              ", Change Ratio: ", round(patient05[i,2],2), sep=""))+
#                theme(text = element_text(size = 15))
#}
#grid.arrange(grobs = plist_05, nrow = 3) ## display plot
## patient05: 2
#ggsave(file = "pre_weight_ts.pdf", arrangeGrob(grobs = plist_05, ncol = 2),
#       width = 12, height = 6)

## recompute the time
first <- da3 %>% 
        group_by(patient_id) %>%
        dplyr::summarise(First = min(FDD_to_txt_dt))
da4 <- first %>%
        right_join(da3,by = "patient_id") %>%
        mutate(time = FDD_to_txt_dt-First) %>%
        arrange(patient_id,FDD_to_txt_dt) %>%
        group_by(patient_id) %>%
        dplyr::mutate(day_diff = ifelse(time==0,0,time-lag(time))) %>%
        filter(!(time!=0 & day_diff==0))


### patient: have at least 1 treatment every 30 days
tmp <- da4 %>% group_by(patient_id) %>% dplyr::summarise(max = max(day_diff,na.rm = T))
selid <- tmp$patient_id[tmp$max <= 30]
da5 <- filter(da4, patient_id %in% selid)

## at least three comlete observations
tmp <- da4 %>% group_by(patient_id) %>% dplyr::summarise(N=n())
selid <- tmp$patient_id[tmp$N>3]
da5 <- filter(da4, patient_id %in% selid)

length(unique(da5$patient_id))

## remove outliers
da5_sub <- da5 %>%
        dplyr::mutate(
                pre_weight_c = pre_weight-(lag(pre_weight)+lead(pre_weight))/2,
                post_weight_c = post_weight-(lag(post_weight)+lead(post_weight))/2,
                ufv_c = ufv-(lag(ufv)+lead(ufv))/2,
                min_sbp_c = min_sbp-(lag(min_sbp)+lead(min_sbp))/2,
                max_sbp_c = max_sbp-(lag(max_sbp)+lead(max_sbp))/2,
                min_dbp_c = min_dbp-(lag(min_dbp)+lead(min_dbp))/2,
                max_dbp_c = max_dbp-(lag(max_dbp)+lead(max_dbp))/2,
                min_pulse_c = min_pulse-(lag(min_pulse)+lead(min_pulse))/2,
                max_pulse_c = max_pulse-(lag(max_pulse)+lead(max_pulse))/2,
                idwg_c = idwg-(lag(idwg)+lead(idwg))/2,
        ) %>%
        dplyr::select(patient_id,time,pre_weight_c,post_weight_c,ufv_c,min_sbp_c,max_sbp_c,
               min_dbp_c,max_dbp_c,min_pulse_c,max_pulse_c,idwg_c)
summary(da5_sub)


outlier <- function(x){
        q2 = quantile(x,c(0.25,0.75),na.rm = T)
        iqr = IQR(x,na.rm = T)
        lower = q2[1]-1.5*iqr
        upper = q2[2]+1.5*iqr
        ifelse(is.na(x),T,x>=lower & x<=upper)
}


da5_sub2 <- da5_sub %>%
        dplyr::mutate(
                pre_weight_o = outlier(pre_weight_c),
                post_weight_o = outlier(post_weight_c),
                ufv_o =  outlier(ufv_c),
                min_sbp_o =  outlier(min_sbp_c),
                max_sbp_o =  outlier(max_sbp_c),
                min_dbp_o =  outlier(min_dbp_c),
                max_dbp_o =  outlier(max_dbp_c),
                min_pulse_o =  outlier(min_pulse_c),
                max_pulse_o =  outlier(max_pulse_c),
                idwg_o = outlier(idwg_c)
        ) %>%
        dplyr::select(patient_id,time,pre_weight_o,post_weight_o,min_sbp_o,max_sbp_o,
               min_dbp_o,max_dbp_o,min_pulse_o,max_pulse_o,idwg_o)

safe = as.logical(apply(da5_sub2[,-c(1,2)],1,prod))

da6 <- da5[safe,] %>%
        dplyr::select(patient_id,time,day_diff,
               ufv, idwg,
               min_sbp, max_sbp, min_dbp, max_dbp, min_pulse, max_pulse) %>%
        dplyr::mutate(day_diff = ifelse(time==0,0,time-lag(time))) %>%
        filter(!(time!=0 & day_diff==0))

names(da6)

#pdf(file="hist.pdf",height=7,width=9)
#for (i in 4:11) {
#        hist(da6[[i]],xlab=dimnames(da6)[[2]][i],col="gray", main="")
#        mtext(paste("Histogram of", dimnames(da6)[[2]][i]))
#}
#dev.off()

## ufv should not be very small
##

da6 <- da6[outlier(da6$ufv),]


### visualize the data
#pdf(file="hist.pdf",height=16,width=16)
#par(mfrow=c(4,4))
#for (i in 3:15) {
#        hist(da6[[i]],xlab=dimnames(da6)[[2]][i],col="gray", main="")
#        mtext(paste("Histogram of", dimnames(da4)[[2]][i]), cex=.5)
#}
#dev.off()

## at least three complete observations
tmp <- da6 %>%
        group_by(patient_id) %>%
        dplyr::summarise(diff_max = max(day_diff))
#tmp <- da6 %>% group_by(patient_id) %>% dplyr::summarise(N=n())
selid <- tmp$patient_id[tmp$diff_max<=30]
da7 <- filter(da6, patient_id %in% selid)

tmp <- da7 %>% group_by(patient_id) %>% dplyr::summarise(N=n())
min(tmp$N)
length(unique(da7$patient_id))

names(da7)

####### patient check:
#setwd("/Users/duanduan/Desktop/Research/Group Graphical Models - Covariance/04 GGM Coding/Real Data2/Patient_ts_plots")
#all_patient <- unique(da7$patient_id)
#for(i in 1:length(all_patient)){
#        patient <- filter(da7, patient_id == all_patient[i])
#        pp1 <- ggplot(patient, aes(x=time, y=ufv)) +
#                geom_line()+
#                xlab("")+
#                ylab("UFV")
        
#        pp2 <- ggplot(patient, aes(x=time, y=min_dbp)) +
#                geom_line()+
#                xlab("")+
#                ylab("Minimum DBP")
        
#        pp3 <- ggplot(patient, aes(x=time, y=max_dbp)) +
#                geom_line()+
#                xlab("")+
#                ylab("Maximum DBP")
        
#        pp4 <- ggplot(patient, aes(x=time, y=min_sbp)) +
#                geom_line()+
#                xlab("")+
#                ylab("Minimum SBP")
        
#        pp5 <- ggplot(patient, aes(x=time, y=max_sbp)) +
#                geom_line()+
#                xlab("")+
#                ylab("Maximum SBP")
        
#        pp6 <- ggplot(patient, aes(x=time, y=min_pulse)) +
#                geom_line()+
#                xlab("")+
#                ylab("Minimum Pulse")
        
#        pp7 <- ggplot(patient, aes(x=time, y=max_pulse)) +
#                geom_line()+
#                xlab("")+
#                ylab("Maximum Pulse")

#        pp8 <- ggplot(patient, aes(x=time, y=idwg)) +
#                geom_line()+
#                xlab("")+
#                ylab("IDWG")

#        ggpubr::ggarrange(
#                pp1,pp2,pp3,pp4,
#                pp5,pp6,pp7,pp8,
#                ncol=4,nrow = 2)
#        ggsave(paste0(i,"_",all_patient[i],".pdf"), width = 16*1.5, height = 8)
#}

patient_rm = c(120147,193989,208489,281677,
               284332,320146,343035,390551,
               405769,418917,453150,454144,
               486150,491489,528985,565838,
               617730,630490,656121,669104,
               682436,695804,735177,736002,
               793591,825979,856164,879071,
               928190)

setwd("/Users/duanduan/Desktop/Research/Group Graphical Models - Covariance/04 GGM Coding/Real Data2")
da8 <- da7 %>%
        filter(!(patient_id %in% patient_rm)) %>%
        dplyr::select(patient_id,ufv, idwg,
                      min_sbp, max_sbp, min_dbp, max_dbp, min_pulse, max_pulse)
summary(da8)


tmp <- da8 %>% group_by(patient_id) %>% dplyr::summarise(N=n())

ggplot(data=tmp, aes(x=N))+
        geom_histogram(binwidth=40,fill="white",color="black")+
        xlab(TeX("Count"))+
        ylab(TeX("Number of Group Observation"))+
        theme_bw()+
        theme(#plot.title = element_text(size=18),
                strip.text.x = element_text(size = 20),
                strip.text.y = element_text(size = 20),
                strip.background =element_rect(fill="white"),
                legend.position = c(0.21, 0.92),
                legend.title=element_text(size=16),
                legend.text=element_text(size=14),
                legend.box = "horizontal",
                legend.background = element_rect(fill = "white", color = "black"),
                axis.text=element_text(size=18),
                axis.title=element_text(size=20))

alln = tmp$N

nu=max(alln)
m=length(alln)
N = sum(alln)

n0 = (N-sum(alln^2)/N)/(m-1)

nu/n0

boxplot(da8$ufv)
#ggsave("histhist.pdf", width = 10, height = 10)

## how many patients? 276 patients
length(unique(da8$patient_id))

### cross validation
### define three functions for cross-validation
Estimate = function(da,lam,term){
        da_group = da %>% group_by(patient_id)
        group_size = da_group %>% dplyr::summarise(N = n())
        group_size = group_size$N
        group_var = da_group %>% do(Var = var(.[,-1]))
        group_var = group_var$Var
        group_mean = da_group %>% do(Mean = colMeans(.[,-1]))
        group_mean = group_mean$Mean
        
        grand_mean = colMeans(da[,-1])
        
        #### find moment estimation of R and G
        total_obs = Reduce("+",group_size)
        group_num = length(group_size)
        sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                                variance=group_var,size=group_size,SIMPLIFY = FALSE)
        )
        mse = sse/(total_obs - group_num)
        R_hat_n = mse
        R_diag = diag(diag(R_hat_n))
        R_hat_rho = solve(sqrt(R_diag))%*%R_hat_n%*%solve(sqrt(R_diag))
        
        da_avg <- Reduce("rbind",group_mean)
        A_hat_n <- cov(da_avg)
        G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
        G_diag = diag(diag(G_hat_n))
        G_hat_rho = solve(sqrt(G_diag))%*%G_hat_n%*%solve(sqrt(G_diag))
        
        A_hat_rho <- cor(da_avg)
        
        del = 10^(-5)
        p = length(grand_mean)
        P <- matrix(1,p,p)
        diag(P) <- 0
        
        if(term=="R"){
                return(admm.cov(R_hat_rho,del,lam,P))
        }
        if(term=="G"){
                return(admm.cov(G_hat_rho,del,lam,P))
        }
        if(term=="A"){
                return(admm.cov(A_hat_rho,del,lam,P))
        }
        
}

Validate = function(da,est,term){
        da_group = da %>% group_by(patient_id)
        group_size = da_group %>% dplyr::summarise(N = n())
        group_size = group_size$N
        group_var = da_group %>% do(Var = var(.[,-1]))
        group_var = group_var$Var
        group_mean = da_group %>% do(Mean = colMeans(.[,-1]))
        group_mean = group_mean$Mean
        
        grand_mean = colMeans(da[,-1])
        
        #### find moment estimation of R and G
        total_obs = Reduce("+",group_size)
        group_num = length(group_size)
        sse = Reduce("+",mapply(function(size,variance) variance*(size-1),
                                variance=group_var,size=group_size,SIMPLIFY = FALSE)
        )
        mse = sse/(total_obs - group_num)
        R_hat_n = mse
        R_diag = diag(diag(R_hat_n))
        R_hat_rho = solve(sqrt(R_diag))%*%R_hat_n%*%solve(sqrt(R_diag))
        
        da_avg <- Reduce("rbind",group_mean)
        A_hat_n <- cov(da_avg)
        G_hat_n <- A_hat_n - mean(1/group_size)*R_hat_n
        G_diag = diag(diag(G_hat_n))
        G_hat_rho = solve(sqrt(G_diag))%*%G_hat_n%*%solve(sqrt(G_diag))
        
        A_hat_rho <- cor(da_avg)
        
        if(term=="R"){
                return(norm(est$Theta - R_hat_rho,"F")^2)
        }
        if(term=="G"){
                return(norm(est$Theta - G_hat_rho,"F")^2)
        }
        if(term=="A"){
                return(norm(est$Theta - A_hat_rho,"F")^2)
        }
        
}

do.chunk = function(chunkid, chunkdef,da,lam){
        train_group = (chunkdef!=chunkid)
        test_group = !train_group
        
        #group_num = length(unique(da[,1]))
        train = da[da$patient_id %in% unique(da$patient_id)[train_group],]
        test = da[da$patient_id %in% unique(da$patient_id)[test_group],]
        
        res = c(Validate(test,Estimate(train,lam,"R"),"R"),
                Validate(test,Estimate(train,lam,"G"),"G"),
                Validate(test,Estimate(train,lam,"A"),"A"))
        names(res) = c("r.error","g.error","a.error")
        return(res)
}


library(parallel)
library(foreach)
library(doParallel)

numCores <- detectCores()
numCores

registerDoParallel(numCores)

### cross-validation
set.seed(12345)

group_num = length(unique(da8$patient_id))
sample.group.id = sample(cut(1:group_num, 10, labels=F))

lam_vec = 10^seq(-4,0,0.05)
lam.res.8 <- foreach(lam = lam_vec,.combine='rbind')%dopar%{
        temp = matrix(nrow = 10,ncol = 3)
        for(chunkid in 1:10){
                temp[chunkid,] = do.chunk(chunkid, sample.group.id,da8,lam)
        }
        c(colMeans(temp),apply(temp,2,sd))
}

lam_res.8 <- as.data.frame(lam.res.8)
names(lam_res.8) <- c('r.error','g.error','a.error','r.sd','g.sd','a.sd')

r.min = which.min(lam_res.8$r.error)
r.1se = lam_res.8$r.error[r.min]+lam_res.8$r.sd[r.min]/sqrt(10)
r.1se
which.max(which(lam_res.8$r.error<=r.1se))
lammda.r.min = seq(-4,0,0.05)[r.min]
lambda.r.1se = seq(-4,0,0.05)[which.max(which(lam_res.8$r.error<=r.1se))]

g.min = which.min(lam_res.8$g.error)
g.1se = lam_res.8$g.error[g.min]+lam_res.8$g.sd[g.min]/sqrt(10)
g.1se
which.max(which(lam_res.8$g.error<=g.1se))
lammda.g.min = seq(-4,0,0.05)[g.min]
lambda.g.1se = seq(-4,0,0.05)[which.max(which(lam_res.8$g.error<=g.1se))]

a.min = which.min(lam_res.8$a.error)
a.1se = lam_res.8$a.error[a.min]+lam_res.8$a.sd[a.min]/sqrt(10)
a.1se
which.max(which(lam_res.8$a.error<=a.1se))
lammda.a.min = seq(-4,0,0.05)[a.min]
lambda.a.1se = seq(-4,0,0.05)[which.max(which(lam_res.8$a.error<=a.1se))]


lammda.r.min
lammda.g.min
lammda.a.min

lambda.r.1se
lambda.g.1se
lambda.a.1se

### final estimation - cv

del = 10^(-5)
p = 8
P <- matrix(1,p,p)
diag(P) <- 0

R_final = Estimate(da8,10^lammda.r.min,"R")
G_final = Estimate(da8,10^lammda.g.min,"G")
A_final = Estimate(da8,10^lammda.a.min,"A")

RR = R_final$Theta
GG = G_final$Theta
AA = A_final$Theta

names(da8)
#"ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse"

colnames(RR) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")
row.names(RR) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")

colnames(GG) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")
row.names(GG) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")

colnames(AA) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")
row.names(AA) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")

melted_r <- reshape2::melt(RR)
melted_r$Var1 <- as.character(melted_r$Var1)
melted_r$Var2 <- as.character(melted_r$Var2)
melted_r = melted_r[melted_r$Var1!=melted_r$Var2 & melted_r$value!=0,]
melted_r$type = ifelse(sign(melted_r$value)==-1,1,2) ## 1: negative, 2: positive 
names(melted_r)[1:2] = c("from","to")

melted_g <- reshape2::melt(GG)
melted_g$Var1 <- as.character(melted_g$Var1)
melted_g$Var2 <- as.character(melted_g$Var2)
melted_g = melted_g[melted_g$Var1!=melted_g$Var2 & melted_g$value!=0,]
melted_g$type = ifelse(sign(melted_g$value)==-1,1,2) ## 1: negative, 2: positive 
names(melted_g)[1:2] = c("from","to")

melted_a <- reshape2::melt(AA)
melted_a$Var1 <- as.character(melted_a$Var1)
melted_a$Var2 <- as.character(melted_a$Var2)
melted_a = melted_a[melted_a$Var1!=melted_a$Var2 & melted_a$value!=0,]
melted_a$type = ifelse(sign(melted_a$value)==-1,1,2) ## 1: negative, 2: positive 
names(melted_a)[1:2] = c("from","to")

### final estimation - cv 1sd
R_final2 = Estimate(da8,10^lambda.r.1se,"R")
G_final2 = Estimate(da8,10^lambda.g.1se,"G")
A_final2 = Estimate(da8,10^lambda.a.1se,"A")

RR2 = R_final2$Theta
GG2 = G_final2$Theta
AA2 = A_final2$Theta

colnames(RR2) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")
row.names(RR2) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")

colnames(GG2) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")
row.names(GG2) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")

colnames(AA2) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")
row.names(AA2) = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse")

melted_r2 <- reshape2::melt(RR2)
melted_r2$Var1 <- as.character(melted_r2$Var1)
melted_r2$Var2 <- as.character(melted_r2$Var2)
melted_r2 = melted_r2[melted_r2$Var1!=melted_r2$Var2 & melted_r2$value!=0,]
melted_r2$type = ifelse(sign(melted_r2$value)==-1,1,2) ## 1: negative, 2: positive 
names(melted_r2)[1:2] = c("from","to")

melted_g2 <- reshape2::melt(GG2)
melted_g2$Var1 <- as.character(melted_g2$Var1)
melted_g2$Var2 <- as.character(melted_g2$Var2)
melted_g2 = melted_g2[melted_g2$Var1!=melted_g2$Var2 & melted_g2$value!=0,]
melted_g2$type = ifelse(sign(melted_g2$value)==-1,1,2) ## 1: negative, 2: positive 
names(melted_g2)[1:2] = c("from","to")

melted_a2 <- reshape2::melt(AA2)
melted_a2$Var1 <- as.character(melted_a2$Var1)
melted_a2$Var2 <- as.character(melted_a2$Var2)
melted_a2 = melted_a2[melted_a2$Var1!=melted_a2$Var2 & melted_a2$value!=0,]
melted_a2$type = ifelse(sign(melted_a2$value)==-1,1,2) ## 1: negative, 2: positive 
names(melted_a2)[1:2] = c("from","to")

nodes = data.frame(nodes = c("ufv","idwg","min_sbp","max_sbp","min_dbp","max_dbp","min_pulse","max_pulse"))

###
net_r <-  graph_from_data_frame(d=melted_r, vertices = nodes,directed=T)
net_r

net_g <-  graph_from_data_frame(d=melted_g, vertices = nodes,directed=T)
net_g

net_a <-  graph_from_data_frame(d=melted_a, vertices = nodes,directed=T)
net_a

###
net_r2 <-  graph_from_data_frame(d=melted_r2, vertices = nodes,directed=T)
net_r2

net_g2 <-  graph_from_data_frame(d=melted_g2, vertices = nodes,directed=T)
net_g2

net_a2 <-  graph_from_data_frame(d=melted_a2, vertices = nodes,directed=T)
net_a2

pdf(file="real_graph_all_pulse_3.pdf",height=12,width=18)
par(mfrow=c(2,3))
par(mar=c(4,4,4,4))
lr <- layout_in_circle(net_r)
plot(net_r, 
     vertex.color="white",
     #vertex.frame.color=V(net_r)$color,
     vertex.size = 30,
     vertex.label.color="gray10",
     vertex.label.cex=1.5,
     vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_r)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_r)$value)*4,
     layout=lr
)
title(main = TeX("(a) Denser Graph at Treatment Level"),#line=-3,
      cex.main = 2.5)

lg <- layout_in_circle(net_g)
plot(net_g, 
     vertex.color="white",
     #vertex.frame.color=V(net_g)$color,
     vertex.size = 30,
     vertex.label.color="gray10",
     vertex.label.cex=1.5,
     vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_g)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_g)$value)*4,
     layout=lg
)
title(main = TeX("(b) Denser Graph at Patient Level"),#line=-3,
      cex.main = 2.5)
#legend(0.5,-1,legend = c("Positive", "Negative"),
#       col = c("blue", "red"),lty = 1,lwd = 2,bty = "n",cex=1.5)

la <- layout_in_circle(net_a)
plot(net_a, 
     vertex.color="white",
     #vertex.frame.color=V(net_a)$color,
     vertex.size = 30,
     vertex.label.color="gray10",
     vertex.label.cex=1.5,
     vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_a)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_a)$value)*4,
     layout=la
)
title(main = TeX("(c) Subject Average, $\\lambda =10^{-1.55}$"),#line=-3,
      cex.main = 2.5)

lr2 <- layout_in_circle(net_r2)
plot(net_r2, 
     vertex.color="white",
     #vertex.frame.color=V(net_r)$color,
     vertex.size = 30,
     vertex.label.color="gray10",
     vertex.label.cex=1.5,
     vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_r2)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_r2)$value)*4,
     layout=lr2
)
title(main = TeX("(d) Sparser Graph at Treatment Level"),#line=-3,
      cex.main = 2.5)

lg2 <- layout_in_circle(net_g2)
plot(net_g2, 
     vertex.color="white",
     #vertex.frame.color=V(net_g)$color,
     vertex.size = 30,
     vertex.label.color="gray10",
     vertex.label.cex=1.5,
     vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_g2)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_g2)$value)*4,
     layout=lg2
)
title(main = TeX("(e) Sparser Graph at Patient Level$"),#line=-3,
      cex.main = 2.5)

la2 <- layout_in_circle(net_a2)
plot(net_a2, 
     vertex.color="white",
     #vertex.frame.color=V(net_a2)$color,
     vertex.size = 30,
     vertex.label.color="gray10",
     vertex.label.cex=1.5,
     vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_a2)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_a2)$value)*4,
     layout=la2
)
title(main = TeX("(f) Subject Average, $\\lambda =10^{-0.95}$"),#line=-3,
      cex.main = 3)
par(mfrow=c(1,1))
dev.off()



setEPS()
# naming the eps file
postscript(file="real22.eps",height=6,width=18,family="Times")
par(mfrow=c(1,3))
par(mar=c(2,2,2.5,2))
lr2 <- layout_in_circle(net_r2)
plot(net_r2, 
     vertex.color="white",
     #vertex.frame.color=V(net_r)$color,
     vertex.size = 50,
     vertex.label.color="gray10",
     vertex.label.cex=2.5,
     vertex.label.family = "Times",
     #vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_r2)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_r2)$value)*4,
     layout=lr2
)
title(main = TeX("Correlation Graph at Treatment Level ($\\widehat{\\Sigma}_{\\epsilon}^+$)"),line=-0.2,
      cex.main = 2.9,family="Times")

lg2 <- layout_in_circle(net_g2)
plot(net_g2, 
     vertex.color="white",
     #vertex.frame.color=V(net_g)$color,
     vertex.size = 50,
     vertex.label.color="gray10",
     vertex.label.cex=2.5,
     vertex.label.family = "Times",
     #vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_g2)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_g2)$value)*4,
     layout=lg2
)
title(main = TeX("Correlation Graph at Patient Level ($\\widehat{\\Sigma}_{b}^+$)"),line=-0.2,
      cex.main = 2.9,family="Times")

la2 <- layout_in_circle(net_a2)
plot(net_a2, 
     vertex.color="white",
     #vertex.frame.color=V(net_g)$color,
     vertex.size = 50,
     vertex.label.color="gray10",
     vertex.label.cex=2.5,
     vertex.label.family = "Times",
     #vertex.label.font=2,
     edge.color=c("blue","red")[(E(net_a2)$type==1)+1],
     edge.arrow.size=0.6,
     edge.width=abs(E(net_a2)$value)*4,
     layout=la2
)
title(main = TeX("Correlation Graph Using Aggregated Data ($\\bar{\\Sigma}^+$)"),line=-0.2,
      cex.main = 2.9,family="Times")
par(mfrow=c(1,1))
dev.off()


Prep = function(f1, X, Y, Vars){
 f1$Metric = rep(X,K)
 f1$Model = rep(Y,K)
 f1$Variable = Vars
 return(f1)
}

nchains = 3
iters = 5000

Vars = c("N_Males", "N_Females","Copulations","Season_Duration","Estrous_Duration","Estrous_Overlap", "Male_Dispersal", "Synchrony_index",
           "N_Males_Offset", "N_Females_Offset","Copulations_Offset","Season_Duration_Offset","Estrous_Duration_Offset","Estrous_Overlap_Offset")

############################################################################## Model M
M = standardize(kutsuskew_data$M_index)
N_m = sum(is.na(M))
Q_m = which(is.na(M))
M[is.na(M)] = -99
M_Raw = M

model_dat = list(
  N=N, K=K, J=J, M_Rough=M_Raw, S=S, Distance=Distance, Males=Males, Females_Rough=Females_Raw, Cops_Rough=Cops_Raw, Seas_Rough=Seas_Raw,EstDur_Rough=EstDur_Raw,
  EstOvO_Rough=EstOvO_Raw,  EstOvE_Rough=EstOvE_Raw,  Dispersal_Rough=Dispersal_Raw,  Sync_Rough=Sync_Raw,  N_m=N_m,  N_disp=N_disp,  N_sync=N_sync,  N_fems=N_fems,
  N_cops=N_cops,  N_seas=N_seas,  N_dur=N_dur,  N_ovo=N_ovo,  N_ove=N_ove,  Q_m=Q_m,  Q_fems=Q_fems,  Q_cops=Q_cops,  Q_seas=Q_seas,  Q_dur=Q_dur,  Q_sync=Q_sync,
  Q_disp=Q_disp,  Q_ovo=Q_ovo,  Q_ove=Q_ove, Pars=rep(1,K)) 

fit0 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)
fit0_M = fit0

model_dat$Pars = rep(0,K); model_dat$Pars[1] = 1; model_dat$Pars[9] = 1;
fit1 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[2] = 1; model_dat$Pars[10] = 1;
fit2 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[3] = 1; model_dat$Pars[11] = 1;
fit3 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[4] = 1; model_dat$Pars[12] = 1;
fit4 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[5] = 1; model_dat$Pars[13] = 1;
fit5 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[6] = 1; model_dat$Pars[14] = 1;
fit6 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[7] = 1;
fit7 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[8] = 1; 
fit8 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

f = vector("list",9)

f[[1]] = Prep(precis(fit1, pars="Beta", depth=2, prob=0.9)@output, "M", "Univariate", Vars)
f[[2]] = Prep(precis(fit2, pars="Beta", depth=2, prob=0.9)@output, "M", "Univariate", Vars)
f[[3]] = Prep(precis(fit3, pars="Beta", depth=2, prob=0.9)@output, "M", "Univariate", Vars)
f[[4]] = Prep(precis(fit4, pars="Beta", depth=2, prob=0.9)@output, "M", "Univariate", Vars)
f[[5]] = Prep(precis(fit5, pars="Beta", depth=2, prob=0.9)@output, "M", "Univariate", Vars)
f[[6]] = Prep(precis(fit6, pars="Beta", depth=2, prob=0.9)@output, "M", "Univariate", Vars)
f[[7]] = Prep(precis(fit7, pars="Beta", depth=2, prob=0.9)@output, "M", "Univariate", Vars)
f[[8]] = Prep(precis(fit8, pars="Beta", depth=2, prob=0.9)@output, "M", "Univariate", Vars)
f[[9]] = Prep(precis(fit0, pars="Beta", depth=2, prob=0.9)@output, "M", "Multivariate", Vars)

df = do.call('rbind', f)
colnames(df) = c("Mean", "StdDev", "L", "H", "n_eff", "Rhat", "Metric", "Model", "Variable")

df_M = df
print("M done")

############################################################################## Model B
M = standardize(kutsuskew_data$B_index)
N_m = sum(is.na(M))
Q_m = which(is.na(M))
M[is.na(M)] = -99
M_Raw = M

model_dat = list(
  N=N, K=K, J=J, M_Rough=M_Raw, S=S, Distance=Distance, Males=Males, Females_Rough=Females_Raw, Cops_Rough=Cops_Raw, Seas_Rough=Seas_Raw,EstDur_Rough=EstDur_Raw,
  EstOvO_Rough=EstOvO_Raw,  EstOvE_Rough=EstOvE_Raw,  Dispersal_Rough=Dispersal_Raw,  Sync_Rough=Sync_Raw,  N_m=N_m,  N_disp=N_disp,  N_sync=N_sync,  N_fems=N_fems,
  N_cops=N_cops,  N_seas=N_seas,  N_dur=N_dur,  N_ovo=N_ovo,  N_ove=N_ove,  Q_m=Q_m,  Q_fems=Q_fems,  Q_cops=Q_cops,  Q_seas=Q_seas,  Q_dur=Q_dur,  Q_sync=Q_sync,
  Q_disp=Q_disp,  Q_ovo=Q_ovo,  Q_ove=Q_ove, Pars=rep(1,K)) 

fit0 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)
fit0_B = fit0

model_dat$Pars = rep(0,K); model_dat$Pars[1] = 1; model_dat$Pars[9] = 1;
fit1 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[2] = 1; model_dat$Pars[10] = 1;
fit2 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[3] = 1; model_dat$Pars[11] = 1;
fit3 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[4] = 1; model_dat$Pars[12] = 1;
fit4 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[5] = 1; model_dat$Pars[13] = 1;
fit5 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[6] = 1; model_dat$Pars[14] = 1;
fit6 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[7] = 1;
fit7 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[8] = 1; 
fit8 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

f = vector("list",9)

f[[1]] = Prep(precis(fit1, pars="Beta", depth=2, prob=0.9)@output, "B", "Univariate", Vars)
f[[2]] = Prep(precis(fit2, pars="Beta", depth=2, prob=0.9)@output, "B", "Univariate", Vars)
f[[3]] = Prep(precis(fit3, pars="Beta", depth=2, prob=0.9)@output, "B", "Univariate", Vars)
f[[4]] = Prep(precis(fit4, pars="Beta", depth=2, prob=0.9)@output, "B", "Univariate", Vars)
f[[5]] = Prep(precis(fit5, pars="Beta", depth=2, prob=0.9)@output, "B", "Univariate", Vars)
f[[6]] = Prep(precis(fit6, pars="Beta", depth=2, prob=0.9)@output, "B", "Univariate", Vars)
f[[7]] = Prep(precis(fit7, pars="Beta", depth=2, prob=0.9)@output, "B", "Univariate", Vars)
f[[8]] = Prep(precis(fit8, pars="Beta", depth=2, prob=0.9)@output, "B", "Univariate", Vars)
f[[9]] = Prep(precis(fit0, pars="Beta", depth=2, prob=0.9)@output, "B", "Multivariate", Vars)

df = do.call('rbind', f)
colnames(df) = c("Mean", "StdDev", "L", "H", "n_eff", "Rhat", "Metric", "Model", "Variable")

df_B = df
print("B done")

############################################################################## Model Lambda
M = standardize(kutsuskew_data$lambda)
N_m = sum(is.na(M))
Q_m = which(is.na(M))
M[is.na(M)] = -99
M_Raw = M

model_dat = list(
  N=N, K=K, J=J, M_Rough=M_Raw, S=S, Distance=Distance, Males=Males, Females_Rough=Females_Raw, Cops_Rough=Cops_Raw, Seas_Rough=Seas_Raw,EstDur_Rough=EstDur_Raw,
  EstOvO_Rough=EstOvO_Raw,  EstOvE_Rough=EstOvE_Raw,  Dispersal_Rough=Dispersal_Raw,  Sync_Rough=Sync_Raw,  N_m=N_m,  N_disp=N_disp,  N_sync=N_sync,  N_fems=N_fems,
  N_cops=N_cops,  N_seas=N_seas,  N_dur=N_dur,  N_ovo=N_ovo,  N_ove=N_ove,  Q_m=Q_m,  Q_fems=Q_fems,  Q_cops=Q_cops,  Q_seas=Q_seas,  Q_dur=Q_dur,  Q_sync=Q_sync,
  Q_disp=Q_disp,  Q_ovo=Q_ovo,  Q_ove=Q_ove, Pars=rep(1,K)) 

fit0 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)
fit0_L = fit0
  
model_dat$Pars = rep(0,K); model_dat$Pars[1] = 1; model_dat$Pars[9] = 1;
fit1 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[2] = 1; model_dat$Pars[10] = 1;
fit2 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[3] = 1; model_dat$Pars[11] = 1;
fit3 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[4] = 1; model_dat$Pars[12] = 1;
fit4 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[5] = 1; model_dat$Pars[13] = 1;
fit5 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[6] = 1; model_dat$Pars[14] = 1;
fit6 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[7] = 1;
fit7 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[8] = 1; 
fit8 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

f = vector("list",9)

f[[1]] = Prep(precis(fit1, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Univariate", Vars)
f[[2]] = Prep(precis(fit2, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Univariate", Vars)
f[[3]] = Prep(precis(fit3, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Univariate", Vars)
f[[4]] = Prep(precis(fit4, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Univariate", Vars)
f[[5]] = Prep(precis(fit5, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Univariate", Vars)
f[[6]] = Prep(precis(fit6, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Univariate", Vars)
f[[7]] = Prep(precis(fit7, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Univariate", Vars)
f[[8]] = Prep(precis(fit8, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Univariate", Vars)
f[[9]] = Prep(precis(fit0, pars="Beta", depth=2, prob=0.9)@output, "Lambda", "Multivariate", Vars)

df = do.call('rbind', f)
colnames(df) = c("Mean", "StdDev", "L", "H", "n_eff", "Rhat", "Metric", "Model", "Variable")

df_L = df
print("L done")

############################################################################## Model MMP
M = standardize(kutsuskew_data$maximum_proportion_mating)
N_m = sum(is.na(M))
Q_m = which(is.na(M))
M[is.na(M)] = -99
M_Raw = M

model_dat = list(
  N=N, K=K, J=J, M_Rough=M_Raw, S=S, Distance=Distance, Males=Males, Females_Rough=Females_Raw, Cops_Rough=Cops_Raw, Seas_Rough=Seas_Raw,EstDur_Rough=EstDur_Raw,
  EstOvO_Rough=EstOvO_Raw,  EstOvE_Rough=EstOvE_Raw,  Dispersal_Rough=Dispersal_Raw,  Sync_Rough=Sync_Raw,  N_m=N_m,  N_disp=N_disp,  N_sync=N_sync,  N_fems=N_fems,
  N_cops=N_cops,  N_seas=N_seas,  N_dur=N_dur,  N_ovo=N_ovo,  N_ove=N_ove,  Q_m=Q_m,  Q_fems=Q_fems,  Q_cops=Q_cops,  Q_seas=Q_seas,  Q_dur=Q_dur,  Q_sync=Q_sync,
  Q_disp=Q_disp,  Q_ovo=Q_ovo,  Q_ove=Q_ove, Pars=rep(1,K)) 

fit0 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)
fit0_P = fit0
  
model_dat$Pars = rep(0,K); model_dat$Pars[1] = 1; model_dat$Pars[9] = 1;
fit1 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[2] = 1; model_dat$Pars[10] = 1;
fit2 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[3] = 1; model_dat$Pars[11] = 1;
fit3 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[4] = 1; model_dat$Pars[12] = 1;
fit4 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[5] = 1; model_dat$Pars[13] = 1;
fit5 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[6] = 1; model_dat$Pars[14] = 1;
fit6 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[7] = 1;
fit7 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

model_dat$Pars = rep(0,K); model_dat$Pars[8] = 1; 
fit8 = stan(file = "Code/Stan_Phylo.stan", data = model_dat, chains=nchains, iter=iters, control=list(adapt_delta=0.98), refresh=1)

f = vector("list",9)

f[[1]] = Prep(precis(fit1, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Univariate", Vars)
f[[2]] = Prep(precis(fit2, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Univariate", Vars)
f[[3]] = Prep(precis(fit3, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Univariate", Vars)
f[[4]] = Prep(precis(fit4, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Univariate", Vars)
f[[5]] = Prep(precis(fit5, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Univariate", Vars)
f[[6]] = Prep(precis(fit6, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Univariate", Vars)
f[[7]] = Prep(precis(fit7, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Univariate", Vars)
f[[8]] = Prep(precis(fit8, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Univariate", Vars)
f[[9]] = Prep(precis(fit0, pars="Beta", depth=2, prob=0.9)@output, "MMP", "Multivariate", Vars)

df = do.call('rbind', f)
colnames(df) = c("Mean", "StdDev", "L", "H", "n_eff", "Rhat", "Metric", "Model", "Variable")

df_P = df

print("P done")

############################################################## Plots
df = rbind(df_M,df_B,df_L,df_P)

df2 = df[which(df$H != 0 & df$L !=0),]

df2$Metric = factor(df2$Metric)
df2$Metric = factor(df2$Metric,levels(df2$Metric)[c(3,1,4,2)])

df3 = df2[which(df2$Variable %in% c("N_Males_Offset", "N_Females_Offset", "Copulations_Offset", "Season_Duration_Offset", 
                                     "Estrous_Duration_Offset", "Estrous_Overlap_Offset")),]

p1 = ggplot(data=df3,
    aes(x = Metric,y = Mean, ymin = L, ymax = H ))+
    geom_pointrange(aes(col=Metric))+
    geom_hline(yintercept = 0, linetype=2)+
    xlab(' ')+ ylab("Effect Size")+
    geom_errorbar(aes(ymin=L, ymax=H,col=Metric),width=0.25,cex=1.5)+ 
    facet_grid(Model~Variable,scales = "free") +
    theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold")) +
    coord_flip()+  scale_color_viridis(direction = -1,option="inferno",discrete=TRUE,begin=0.15, end=0.8)+   
    theme(strip.text.x = element_text(size=12,face="bold"), strip.text.y = element_text(size=12,face="bold")) +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=14,face="bold")) +
    theme(legend.title=element_text(size=12,face="bold"),legend.text=element_text(size=12))+ guides(color = guide_legend(reverse = TRUE))
 p1

p1b <- ggplot(df3,aes(x=Metric,y=Mean,ymin=L,ymax=H,color=Metric,linetype=Model))+ 
     geom_linerange(size=1,aes(color=Metric), position = position_dodge2(width=0.5,preserve = "total",padding = 1))+
     geom_point(size=2,aes(color=Metric), position = position_dodge2(width=0.5,preserve = "total",padding = 1))+
     facet_grid(~Variable,scales="free",space = "free_y")+geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
     labs(y="Regression parameters", x="") + theme(strip.text.x = element_text(size=14,face="bold"), 
     strip.text.y = element_text(size=14,face="bold"),axis.text=element_text(size=12),axis.title.y=element_text(size=14,
     face="bold"), axis.title.x=element_blank())+theme(strip.text.y = element_text(angle = 360))  + coord_flip() + theme(panel.spacing = unit(1, "lines")) +
     scale_color_viridis(direction = -1,option="inferno",discrete=TRUE,begin=0.15, end=0.8) +guides(color = FALSE)

 p1b  

ggsave("Phylo1.pdf",p1b,width=18,height=3)

df4 = df2[which(!df2$Variable %in% c("N_Males_Offset", "N_Females_Offset", "Copulations_Offset", "Season_Duration_Offset", 
                                     "Estrous_Duration_Offset", "Estrous_Overlap_Offset")),]

p2 = ggplot(data=df4,
    aes(x = Metric,y = Mean, ymin = L, ymax = H ))+
    geom_pointrange(aes(col=Metric))+
    geom_hline(yintercept = 0, linetype=2)+
    xlab(' ')+ ylab("Effect Size")+
    geom_errorbar(aes(ymin=L, ymax=H,col=Metric),width=0.25,cex=1.5)+ 
    facet_grid(Model~Variable,scales = "free") +
    theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold")) +
    coord_flip()+  scale_color_viridis(direction = -1,option="inferno",discrete=TRUE,begin=0.15, end=0.8)+   
    theme(strip.text.x = element_text(size=14,face="bold"), strip.text.y = element_text(size=14,face="bold")) +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=14,face="bold")) +
    theme(legend.title=element_text(size=14,face="bold"),legend.text=element_text(size=12))+ guides(color = guide_legend(reverse = TRUE))
 p2

p2b <- ggplot(df4,aes(x=Metric,y=Mean,ymin=L,ymax=H,color=Metric,linetype=Model))+ 
     geom_linerange(size=1,aes(color=Metric), position = position_dodge2(width=0.5,preserve = "total",padding = 1))+
     geom_point(size=2,aes(color=Metric), position = position_dodge2(width=0.5,preserve = "total",padding = 1))+
     facet_grid(~Variable,scales="free",space = "free_y")+geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
     labs(y="Regression parameters", x="") + theme(strip.text.x = element_text(size=14,face="bold"), 
     strip.text.y = element_text(size=14,face="bold"),axis.text=element_text(size=12),axis.title.y=element_text(size=14,
     face="bold"), axis.title.x=element_blank())+theme(strip.text.y = element_text(angle = 360))  + coord_flip() + theme(panel.spacing = unit(1, "lines")) +
     scale_color_viridis(direction = -1,option="inferno",discrete=TRUE,begin=0.15, end=0.8) +guides(color = FALSE)

 p2b  

ggsave("Phylo2.pdf",p2b,width=18,height=3)


ggsave("Trace_P_M.pdf", traceplot(fit0_M, pars=c("Beta")), height=4, width=12)   
ggsave("Trace_P_B.pdf", traceplot(fit0_B, pars=c("Beta")), height=4, width=12)   
ggsave("Trace_P_L.pdf", traceplot(fit0_L, pars=c("Beta")), height=4, width=12)   
ggsave("Trace_P_P.pdf", traceplot(fit0_P, pars=c("Beta")), height=4, width=12)   

print(fit0_M, pars=c("Beta"))
print(fit0_B, pars=c("Beta"))
print(fit0_L, pars=c("Beta"))
print(fit0_P, pars=c("Beta"))

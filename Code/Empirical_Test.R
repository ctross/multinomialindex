# Load library and attach data
data(SukumaMales) 
data(KipsigisMales) 
data(KipsigisFemales) 
data(ColombiaRS) 

d = ColombiaRS
d$age = d$age - 12 # Trim-off pre-reproductive life period, to get years of exposure to risk of RS
KipsigisMales$age = KipsigisMales$age - 12
KipsigisFemales$age = KipsigisFemales$age - 12

M_index_stan(d$rs[which(d$group=="AFROCOLOMBIAN" & d$sex=="M")], d$age[which(d$group=="AFROCOLOMBIAN" & d$sex=="M")], 
	          adapt_delta=0.995, samples = 8000, warmup = 4000, chains=4) 
M_post_A_male = extract(StanResults, pars="M_age")$M_age
M_point_A_male = M_index_age(model_dat$r, model_dat$t, Samples=500) 
StanResults_A_Male = StanResults

M_index_stan(d$rs[which(d$group=="AFROCOLOMBIAN" & d$sex=="F")], d$age[which(d$group=="AFROCOLOMBIAN" & d$sex=="F")],
	          adapt_delta=0.995, samples = 8000, warmup = 4000, chains=4)
M_post_A_female = extract(StanResults, pars="M_age")$M_age
M_point_A_female = M_index_age(model_dat$r, model_dat$t, Samples=500) 
StanResults_A_Female = StanResults

M_index_stan(d$rs[which(d$group=="EMBERA" & d$sex=="M")], d$age[which(d$group=="EMBERA" & d$sex=="M")],
	          adapt_delta=0.995, samples = 8000, warmup = 4000, chains=4)
M_post_E_male = extract(StanResults, pars="M_age")$M_age
M_point_E_male = M_index_age(model_dat$r, model_dat$t, Samples=500) 
StanResults_E_Male = StanResults

M_index_stan(d$rs[which(d$group=="EMBERA" & d$sex=="F")], d$age[which(d$group=="EMBERA" & d$sex=="F")],
	          adapt_delta=0.995, samples = 8000, warmup = 4000, chains=4)
M_post_E_female = extract(StanResults, pars="M_age")$M_age
M_point_E_female = M_index_age(model_dat$r, model_dat$t, Samples=500) 
StanResults_E_Female = StanResults

M_index_stan(KipsigisMales$rs, KipsigisMales$age, adapt_delta=0.995, samples = 8000, warmup = 4000, chains=4)
M_post_K_male = extract(StanResults, pars="M_age")$M_age
M_point_K_male = M_index_age(model_dat$r, model_dat$t, Samples=500) 
StanResults_K_Male = StanResults

M_index_stan(KipsigisFemales$rs, KipsigisFemales$age, adapt_delta=0.995, samples = 8000, warmup = 4000, chains=4)
M_post_K_female = extract(StanResults, pars="M_age")$M_age
M_point_K_female = M_index_age(model_dat$r, model_dat$t, Samples=500) 
StanResults_K_Female = StanResults

# Finally, plot the posterior estimates of M by group and sex
df1 = data.frame(M=M_post_A_male, Sex=rep("Male",length(M_post_A_male)), Group=rep("Afrocolombian",length(M_post_A_male)))
df2 = data.frame(M=M_post_A_female, Sex=rep("Female",length(M_post_A_female)), Group=rep("Afrocolombian",length(M_post_A_female)))
df3 = data.frame(M=M_post_E_male, Sex=rep("Male",length(M_post_E_male)), Group=rep("Embera",length(M_post_E_male)))
df4 = data.frame(M=M_post_E_female, Sex=rep("Female",length(M_post_E_female)), Group=rep("Embera",length(M_post_E_female)))
df6 = data.frame(M=M_post_K_male, Sex=rep("Male",length(M_post_K_male)), Group=rep("Kipsigis",length(M_post_K_male)))
df7 = data.frame(M=M_post_K_female, Sex=rep("Female",length(M_post_K_female)), Group=rep("Kipsigis",length(M_post_K_female)))
df = rbind(df1, df2, df3, df4, df6, df7)

dfb1 = data.frame(x0=M_point_A_male, Sex=rep("Male",1), Group=rep("Afrocolombian",1))
dfb2 = data.frame(x0=M_point_A_female, Sex=rep("Female",1), Group=rep("Afrocolombian",1))
dfb3 = data.frame(x0=M_point_E_male, Sex=rep("Male",1), Group=rep("Embera",1))
dfb4 = data.frame(x0=M_point_E_female, Sex=rep("Female",1), Group=rep("Embera",1))
dfb6 = data.frame(x0=M_point_K_male, Sex=rep("Male",1), Group=rep("Kipsigis",1))
dfb7 = data.frame(x0=M_point_K_female, Sex=rep("Female",1), Group=rep("Kipsigis",1))
dfb = rbind(dfb1, dfb2, dfb3, dfb4, dfb6, dfb7)

df$Group = factor(df$Group)
df$Group = factor(df$Group,levels(df$Group)[c(2,1,3)])
  
dfb$Group = factor(dfb$Group)
dfb$Group = factor(dfb$Group,levels(dfb$Group)[c(2,1,3)])

p4 = ggplot()  +
 stat_density_ridges(data=df, aes(x=M, y=Group, fill=0.5 - abs(0.5-..ecdf..)),
  geom = "density_ridges_gradient", calc_ecdf = TRUE, color="white") +
  scale_fill_viridis(name = "Tail probability", direction = -1, option="inferno")+ facet_grid(.~Sex) +   
  theme(strip.text.x = element_text(size=14, face="bold"), strip.text.y = element_text(size=14, face="bold")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  theme(legend.title=element_text(size=14),legend.text=element_text(size=12))+      
  geom_point(data=dfb, aes(x = x0, y = as.numeric(Group) + 0.25), color = "white",size=3,shape=18) +
  geom_segment(data=dfb, aes(x = x0, xend = x0, y = as.numeric(Group), yend = as.numeric(Group) + 0.25), color = "white",size=1.5) + 
  geom_point(data=dfb, aes(x = x0, y = as.numeric(Group) + 0.25), color = "royalblue3",size=2,shape=18) +
  geom_segment(data=dfb, aes(x = x0, xend = x0, y = as.numeric(Group), yend = as.numeric(Group) + 0.25), color = "royalblue3") + 
  theme_ridges(grid = TRUE, center = TRUE) +   geom_hline(yintercept=c(1,2,3), color="white") + ylab("")

p4

ggsave("EmpRes.pdf", p4, height=3.5, width=11.5)




ggsave("Trace_E_K_M.pdf", traceplot(StanResults_K_Male, 
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)    
ggsave("Trace_E_K_F.pdf", traceplot(StanResults_K_Female, 
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)       

ggsave("Trace_E_A_M.pdf", traceplot(StanResults_A_Male, 
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)       
ggsave("Trace_E_A_F.pdf", traceplot(StanResults_A_Female, 
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)      

ggsave("Trace_E_E_M.pdf", traceplot(StanResults_E_Male, 
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)   
ggsave("Trace_E_E_F.pdf", traceplot(StanResults_E_Female, 
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)   


print(StanResults_K_Male)
print(StanResults_K_Female)

print(StanResults_A_Male)
print(StanResults_A_Female)

print(StanResults_E_Male)
print(StanResults_E_Female)

print(mean(c(M_post_K_male-M_post_A_male)))
print(HPDI(c(M_post_K_male-M_post_A_male),0.9))

print(mean(c(M_post_K_male-M_post_E_male)))
print(HPDI(c(M_post_K_male-M_post_E_male),0.9))

print(mean(c(M_post_A_male-M_post_E_male)))
print(HPDI(c(M_post_A_male-M_post_E_male),0.9))

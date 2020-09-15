############################################################## Stan Sim
ResS = function(x){
           r = x$r
           t = x$t

           M_index_stan(r, t, chains=3, samples=6000, warmup=4000, adapt_delta=0.99) 
           M_post = extract(StanResults, pars="M_age")$M_age
           M_point = M_index_age(model_dat$r, model_dat$t, Samples=500)

           d = c(M_post,M_point)

           return(d)
           }

################################################################## RS Simulation
sim_rs = function(N=100, Rate=3.3, L=0, U=1, et=1, es=0.3){
           Time = 3*runif(N, L, 1)^U  # Age of each individual  

           B = 0.5              # Gamma scalar.
           Z = rgamma(N, 1*B, B)  # Random component of RS based on variation in quality.

           # Loop over individuals, and simulate RS given exposre time. Negative-binomial model, written here as a Gamma-Poisson mixture.
           r = rep(NA, N)       
           for(i in 1:N){
            mu = Rate*(Time[i]^et)*(Z[i]^es)  
            r[i] = rpois(1, mu) 
            }
           
           res = list(r=r, t=Time)
           
           return(res)
           }

#################################################################
 Q = 5   # Diferences
 S = 6001
 
 Pop = c(10,30,80,150,450) 
 
 Res1b = array(NA,c(Q,S))
 Res2b = array(NA,c(Q,S))
 Res3b = array(NA,c(Q,S))
 Res4b = array(NA,c(Q,S))
 Res5b = array(NA,c(Q,S))
 Res6b = array(NA,c(Q,S))
 Res7b = array(NA,c(Q,S))
 Res8b = array(NA,c(Q,S))
 Res9b = array(NA,c(Q,S))

 for(q in 1:Q){
 Res1b[q,] = ResS(sim_rs(N=Pop[q], Rate=1.0, L=0.2, U=1, et=1, es=0.01))
 Res2b[q,] = ResS(sim_rs(N=Pop[q], Rate=1.0, L=0.2, U=1, et=1, es=0.31))
 Res3b[q,] = ResS(sim_rs(N=Pop[q], Rate=1.0, L=0.2, U=1, et=1, es=0.61))

 Res4b[q,] = ResS(sim_rs(N=Pop[q], Rate=7.0, L=0.2, U=1, et=1, es=0.01))
 Res5b[q,] = ResS(sim_rs(N=Pop[q], Rate=7.0, L=0.2, U=1, et=1, es=0.31))
 Res6b[q,] = ResS(sim_rs(N=Pop[q], Rate=7.0, L=0.2, U=1, et=1, es=0.61))
 
 Res7b[q,] = ResS(sim_rs(N=Pop[q], Rate=20.0, L=0.2, U=1, et=1, es=0.01))
 Res8b[q,] = ResS(sim_rs(N=Pop[q], Rate=20.0, L=0.2, U=1, et=1, es=0.31))
 Res9b[q,] = ResS(sim_rs(N=Pop[q], Rate=20.0, L=0.2, U=1, et=1, es=0.61))
              }

 Parse2 = function(x,a,b,c){ 
  datF = vector("list",5)
   for(q in 1:5){
           dat = data.frame(M=t(x)[,q])
           dat$Lab1 = rep(a, dim(x)[2])
           dat$Lab2 = rep(b, dim(x)[2])
           dat$Lab3 = rep(c, dim(x)[2])
           dat$Lab4 = rep(Pop[q], dim(x)[2])
           dat$Lab5 = rep("Post", dim(x)[2])
           dat$Lab5[dim(x)[2]] = "Point"
           datF[[q]] = dat
           }        
  return(do.call('rbind', datF) )
  }

           
             
 dat1x = rbind(Parse2(Res1b, "Low Skew", "Low Mean", "Unequal Exposure"),
               Parse2(Res2b, "Mid Skew", "Low Mean", "Unequal Exposure"),
               Parse2(Res3b, "High Skew", "Low Mean", "Unequal Exposure"),

               Parse2(Res4b, "Low Skew", "Mid Mean", "Unequal Exposure"),
               Parse2(Res5b, "Mid Skew", "Mid Mean", "Unequal Exposure"),
               Parse2(Res6b, "High Skew", "Mid Mean", "Unequal Exposure"),

               Parse2(Res7b, "Low Skew", "High Mean", "Unequal Exposure"),
               Parse2(Res8b, "Mid Skew", "High Mean", "Unequal Exposure"),
               Parse2(Res9b, "High Skew", "High Mean", "Unequal Exposure")
  )
              
 dat1x$Lab2 = factor(dat1x$Lab2)
 dat1x$Lab2 = factor(dat1x$Lab2, levels(dat1x$Lab2)[c(1,3,2)])

 dat1x$Lab1 = factor(dat1x$Lab1)
 dat1x$Lab1 = factor(dat1x$Lab1, levels(dat1x$Lab1)[c(1,3,2)])

 dat2x = dat1x

 dat1x.p = dat2x[which(dat2x$Lab5=="Post"),]
 dat1x.po = dat2x[which(dat2x$Lab5=="Point"),]
 
 dat1x.p = dat1x.p[complete.cases(dat1x.p),]
 

 pB =  ggplot() +    
 stat_density_ridges(data= dat1x.p, aes(x=M, y=factor(Lab4), fill=0.5 - abs(0.5-..ecdf..)),
  geom = "density_ridges_gradient", calc_ecdf = TRUE, color="white") +
  scale_fill_viridis(name = "Tail probability", direction = -1,option="inferno")+   
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  theme(legend.title=element_text(size=14),legend.text=element_text(size=12))+      
  geom_point(data=dat1x.po, aes(x = M, y = as.numeric(factor(Lab4)) + 0.25), color = "white",size=3,shape=18) +
  geom_segment(data=dat1x.po, aes(x = M, xend = M, y = as.numeric(factor(Lab4)), yend = as.numeric(factor(Lab4)) + 0.25), color = "white",size=1.5) + 
  geom_point(data=dat1x.po, aes(x = M, y = as.numeric(factor(Lab4)) + 0.25), color = "royalblue3",size=2,shape=18) +
  geom_segment(data=dat1x.po, aes(x = M, xend = M, y = as.numeric(factor(Lab4)), yend = as.numeric(factor(Lab4)) + 0.25), color = "royalblue3") + 

  facet_wrap(Lab2~Lab1,scales="fixed") +  ylab("Sample size")  + coord_cartesian(xlim = c(-1.5, 2.7), expand = FALSE)
   
  pB
   
  ggsave("BayesCheck.pdf", pB, height=10, width=12)



x = sim_rs(N=Pop[1], Rate=7.0, L=0.2, U=1, et=1, es=0.31)
r = x$r
t = x$t
M_index_stan(r, t, chains=4, samples=6000, warmup=4000, adapt_delta=0.99) 
StanResults_Check1 = StanResults
ggsave("Trace_C_1.pdf", traceplot(StanResults_Check1,
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)   

x = sim_rs(N=Pop[3],Rate=7.0, L=0.2, U=1,et=1,es=0.31)
r = x$r
t = x$t
M_index_stan(r, t, chains=4, samples=6000, warmup=4000, adapt_delta=0.99) 
StanResults_Check2 = StanResults
ggsave("Trace_C_2.pdf", traceplot(StanResults_Check2,
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)   

x = sim_rs(N=Pop[5],Rate=7.0, L=0.2, U=1,et=1,es=0.31)
r = x$r
t = x$t
M_index_stan(r, t, chains=4, samples=6000, warmup=4000, adapt_delta=0.99) 
StanResults_Check3 = StanResults
ggsave("Trace_C_3.pdf", traceplot(StanResults_Check3, 
	pars=c("gamma", "Concentration", "M_raw", "M", "M_raw_age", "M_age", "alpha[1]", "alpha[5]", "alpha[10]")), height=4, width=12)   


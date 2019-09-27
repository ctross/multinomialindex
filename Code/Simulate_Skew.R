################################################################## RS Simulation
sim_rs <- function(N=100,Rate=3.3, L=0, U=1,et=1,es=0.3){
           Time <- 3*runif(N,L,1)^U  # Age of each individual  

           B<-.5
           Z <- rgamma(N,1*B,B)# Sexyness of each individual

           r <- rep(NA,N)        # RS
           for(i in 1:N)   {
           mu <- Rate*(Time[i]^et)*(Z[i]^es)  
           r[i] <- rpois(1, mu) }
           
           res <- list(r=r, t=Time)
           
           return(res)
           }

################################################################### Skew Metrics
B_index <- function(r,t) {   #r=rs,t=exposure
	         T <- sum(t)
	         Nbar <- T / max(t)
           R <- sum(r)
	         if(R>0){
           B = sum((r / R - t / T)^2) - (1 / R) * (1 - 1 / Nbar)
	         } else {
		       B = NA
	         }
	         return(B)
           }
           
I_index <- function(r,t) {
	         R <- sum(r)
           if(R>0){
           return(var(r)/(mean(r)^2))}
           else{
           return(NA)
           }
           }

Mraw_index <- function(r,t) {
	         R <- sum(r)
	         T <- sum(t)
           N <- length(t)
	         si <- ((r/R)-(t/T))^2
	         S <- sum(si)
	         C <- (N * S)
           if(R>0){
           return(C)}
           else{
           return(NA)
           }
           }

M_index <- function(r,t) {
           R <- sum(r)
           T <- sum(t)
           N <- length(t)
           
           alpha <-  t/T

           X <- rmultinom(1,sum(R),alpha)
         
           return(Mraw_index(r,t) - Mraw_index(X,t) )
     
           }

MMP_index <- function(r) { 
           R <- sum(r)
           if(R>0){
            X <- max(r)/R
            return(X)
            }else{
           return(NA)
           }
}

Lambda_index <- function(r) { # Cant really tell if this is correct, given the explaination in original paper
           R <- sum(r)
           N <- length(r)
           r <- rev(sort(r))
           y <- 1:N
           if(R>0){
    f = function(x){
            E_p_r <- (x*(1-x)^(y-1)) / (1 - (1-x)^N)
            O_p_r <- r/R
       return(sum((O_p_r-E_p_r)^2))
      }
   X<-optim( c(0.03), f,lower = 0, upper = 1 ,method="Brent")$par
            return(X)
            }else{
           return(NA)
           }
}

########################################################### Results Calculations
Res <- function(x){
           r <- x$r
           t <- x$t
           d<-c()
           d[1] <- I_index(r,t)
           d[2] <- I_index(r/t,t)
           d[3] <- B_index(r,t)
           d[4] <- Mraw_index(r,t)
           d[5] <- M_index(r,t)
           d[6] <- sum(r)
           d[7] <- sum(t)
           d[8] <- d[6]/d[7]
           d[9] <- length(t)
           d[10] <- MMP_index(r)
           d[11] <- Lambda_index(r)
           return(d)
           }
           
           
Parse <- function(Results,Lab1="A",Lab2="B",Lab3="C"){
           CV <- c(Results[,,1])
           rCV <- c(Results[,,2])
           B <- c(Results[,,3])
           Mraw <- c(Results[,,4])
           M <- c(Results[,,5])
           R <- c(Results[,,6])
           T <- c(Results[,,7])
           Rate <- c(Results[,,8])
           N <- c(Results[,,9])
           MMP <- c(Results[,,10])
           Lambda <- c(Results[,,11])
           
           CV_dat <- data.frame(Value=CV,N=N,R=R,T=T,Rate=Rate,Index=rep("I",length(CV)))  
           rCV_dat <- data.frame(Value=rCV,N=N,R=R,T=T,Rate=Rate,Index=rep("I_rate",length(CV)))   
           B_dat <- data.frame(Value=B,N=N,R=R,T=T,Rate=Rate,Index=rep("B",length(B))) 
           Mraw_dat <- data.frame(Value=Mraw,N=N,R=R,T=T,Rate=Rate,Index=rep("M_raw",length(Mraw))) 
           M_dat <- data.frame(Value=M,N=N,R=R,T=T,Rate=Rate,Index=rep("M",length(M))) 
           MMP_dat <- data.frame(Value=MMP,N=N,R=R,T=T,Rate=Rate,Index=rep("MMP",length(MMP))) 
           Lambda_dat <- data.frame(Value=Lambda,N=N,R=R,T=T,Rate=Rate,Index=rep("Lambda",length(Lambda))) 
  
           dat <- rbind(CV_dat, rCV_dat, B_dat, Mraw_dat, M_dat, MMP_dat, Lambda_dat)
           dat$Value <- ifelse(dat$Value==Inf | dat$Value== -Inf,NA,dat$Value)
           dat$Lab1 <- rep(Lab1,length(dat$Value))
           dat$Lab2 <- rep(Lab2,length(dat$Value))
           dat$Lab3 <- rep(Lab3,length(dat$Value))
           return(dat)
           }

############################################################## Run Simulation, N
 S <- 500  # Samples
 Q <- 30   # Diferences
 
 Pop <- round(exp(seq(2.2,7,length.out=Q)),0) 
 Res1 <- array(NA,c(S,Q,11))
 Res2 <- array(NA,c(S,Q,11))
 Res3 <- array(NA,c(S,Q,11))
 Res4 <- array(NA,c(S,Q,11))
 Res5 <- array(NA,c(S,Q,11))
 Res6 <- array(NA,c(S,Q,11))
 Res7 <- array(NA,c(S,Q,11))
 Res8 <- array(NA,c(S,Q,11))
 Res9 <- array(NA,c(S,Q,11))
 
 Res1b <- array(NA,c(S,Q,11))
 Res2b <- array(NA,c(S,Q,11))
 Res3b <- array(NA,c(S,Q,11))
 Res4b <- array(NA,c(S,Q,11))
 Res5b <- array(NA,c(S,Q,11))
 Res6b <- array(NA,c(S,Q,11))
 Res7b <- array(NA,c(S,Q,11))
 Res8b <- array(NA,c(S,Q,11))
 Res9b <- array(NA,c(S,Q,11))
 
for(s in 1:S){
 for(q in 1:Q){
 Res1[s,q,] <- Res(sim_rs(N=Pop[q],Rate=1.0, L=1, U=1,et=1,es=0.01))
 Res2[s,q,] <- Res(sim_rs(N=Pop[q],Rate=1.0, L=1, U=1,et=1,es=0.39))
 Res3[s,q,] <- Res(sim_rs(N=Pop[q],Rate=1.0, L=1, U=1,et=1,es=0.79))

 Res4[s,q,] <- Res(sim_rs(N=Pop[q],Rate=7.0, L=1, U=1,et=1,es=0.01))
 Res5[s,q,] <- Res(sim_rs(N=Pop[q],Rate=7.0, L=1, U=1,et=1,es=0.39))
 Res6[s,q,] <- Res(sim_rs(N=Pop[q],Rate=7.0, L=1, U=1,et=1,es=0.79))
 
 Res7[s,q,] <- Res(sim_rs(N=Pop[q],Rate=20.0, L=1, U=1,et=1,es=0.01))
 Res8[s,q,] <- Res(sim_rs(N=Pop[q],Rate=20.0, L=1, U=1,et=1,es=0.39))
 Res9[s,q,] <- Res(sim_rs(N=Pop[q],Rate=20.0, L=1, U=1,et=1,es=0.79))


 Res1b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=1.0, L=0.2, U=1,et=1,es=0.01))
 Res2b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=1.0, L=0.2, U=1,et=1,es=0.39))
 Res3b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=1.0, L=0.2, U=1,et=1,es=0.79))

 Res4b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=7.0, L=0.2, U=1,et=1,es=0.01))
 Res5b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=7.0, L=0.2, U=1,et=1,es=0.39))
 Res6b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=7.0, L=0.2, U=1,et=1,es=0.79))
 
 Res7b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=20.0, L=0.2, U=1,et=1,es=0.01))
 Res8b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=20.0, L=0.2, U=1,et=1,es=0.39))
 Res9b[s,q,] <- Res(sim_rs(N=Pop[q],Rate=20.0, L=0.2, U=1,et=1,es=0.79))
              }
                }
                
base_breaks <- function(n = 10){
    function(x) {
        axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
}

              
 dat1 <- rbind(Parse(Res1, "Low Skew", "Low Mean", "Equal Exposure"),
               Parse(Res2, "Mid Skew", "Low Mean", "Equal Exposure"),
               Parse(Res3, "High Skew", "Low Mean", "Equal Exposure"),

               Parse(Res4, "Low Skew", "Mid Mean", "Equal Exposure"),
               Parse(Res5, "Mid Skew", "Mid Mean", "Equal Exposure"),
               Parse(Res6, "High Skew", "Mid Mean", "Equal Exposure"),

               Parse(Res7, "Low Skew", "High Mean", "Equal Exposure"),
               Parse(Res8, "Mid Skew", "High Mean", "Equal Exposure"),
               Parse(Res9, "High Skew", "High Mean", "Equal Exposure")
  )
              
 dat2 <- rbind(Parse(Res1b, "Low Skew", "Low Mean", "Unequal Exposure"),
               Parse(Res2b, "Mid Skew", "Low Mean", "Unequal Exposure"),
               Parse(Res3b, "High Skew", "Low Mean", "Unequal Exposure"),

               Parse(Res4b, "Low Skew", "Mid Mean", "Unequal Exposure"),
               Parse(Res5b, "Mid Skew", "Mid Mean", "Unequal Exposure"),
               Parse(Res6b, "High Skew", "Mid Mean", "Unequal Exposure"),

               Parse(Res7b, "Low Skew", "High Mean", "Unequal Exposure"),
               Parse(Res8b, "Mid Skew", "High Mean", "Unequal Exposure"),
               Parse(Res9b, "High Skew", "High Mean", "Unequal Exposure")
  )
              
 dat<-rbind(dat1,dat2)

colnames(dat) <- c("Value", "N", "R", "T", "Rate", "Index", "Skew", "Mean", "Exposure")
  

dat$Mean <- factor(dat$Mean)
dat$Mean <- factor(dat$Mean,levels(dat$Mean)[c(1,3,2)])

dat$Skew <- factor(dat$Skew)
dat$Skew <- factor(dat$Skew,levels(dat$Skew)[c(1,3,2)])

dat$Index <- factor(dat$Index)
dat$Index <- factor(dat$Index,levels(dat$Index)[c(6,7,3,1,2,4,5)])

 (p1<- ggplot(dat, aes(N, Value, color=Mean, fill=Mean, lty=Exposure))+
  # geom_point()+
  stat_smooth()+
  facet_grid(Index~Skew,scales="free_y")  +  ylab("Index Value")+ xlab("Sample Size")
  + scale_x_continuous(trans='log10',breaks = base_breaks(n=9)) 
  +   theme(strip.text.x = element_text(size=14, face="bold"),
          strip.text.y = element_text(size=14, face="bold")) + theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+ 
          scale_color_viridis(discrete=TRUE,end=0.8, option="inferno",direction=-1)+ 
          scale_fill_viridis(discrete=TRUE,end=0.8, option="inferno",direction=-1) +
          theme(legend.position="bottom", legend.box = "vertical",legend.title = element_blank()) +
          theme(legend.text=element_text(size=12))+guides(fill = guide_legend(order=1),color = guide_legend(order=1),
         lty = guide_legend(order=2))+ guides(lty = guide_legend(override.aes = list(col = 'black')))
                                 )
  
ggsave("SimRes.pdf",p1,height=11.5,width=8.5)










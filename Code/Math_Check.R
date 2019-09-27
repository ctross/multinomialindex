




sim_rs <- function(N=100,Rate=3.3, L=0, U=1,et=1,es=0.3,j){
	set.seed(1)
           Time <- 3*runif(N,L,1)^U  # Age of each individual  
     set.seed(j)
           B<-.5
           Z <- rgamma(N,1*B,B)# Sexyness of each individual

           r <- rep(NA,N)        # RS
           for(i in 1:N)   {
           mu <- Rate*(Time[i]^et)*(Z[i]^es)  
           r[i] <- rpois(1, mu) }
           
           res <- list(r=r, t=Time)
           
           return(res)
           }
Reps <- 100000
res <- matrix(NA,nrow=Reps,ncol=8)

for( j in 1:Reps){
x <- sim_rs(N=1000,Rate=7.0, L=0.1, U=1,et=1,es=0.39,j=j+1)

r <- x$r
t <- x$t
N <- length(r)
R <- sum(r)
T <- sum(t)
t_hat <- t/T
r_hat <- r/R

res[j,1] <- ((1/N)*sum((r-(R*t_hat))^2))/(R^2/N^2)     # 3.2a-2
res[j,2] <- (((R^2)/N)*sum((r_hat-t_hat)^2))/(R^2/N^2) # 3.2a-1
res[j,3] <- N*sum((r_hat-t_hat)^2)                     # 3.2a
res[j,4] <- N*sum(((r/R)-(t/T))^2)                     # 3.2b
res[j,5] <- (N/R^2)*sum((r-R*t_hat)^2)                 # 3.2c and 3.2d

res[j,6] <- (N^2/R^2)*mean((r-R*t_hat)^2)              # 3.2e
res[j,7] <- (N^2/R^2)*(var(r)-var(R*t_hat))            # 3.2f
res[j,8] <- (N^2/R^2)*var(r)*(1-cor(r,R*t_hat)^2)      # 3.2g
}

print(colMeans(res))











 # MbF <- function (X, D, N, J){
 #  Q <- rep(NA,J)
 #  Res <- rep(NA,N)
    
 #     for(j in 1:J){ 
 #      ScrapV = 0;   
 #      ScrapS = 0;
   
 #      for( i in 1:N){
 #       ScrapV = ifelse(D[i]==j, X[i]+ScrapV, ScrapV)
 #       ScrapS = ifelse(D[i]==j, 1+ScrapS, ScrapS)   
 #      }
    
 #     Q[j] = ScrapV/ScrapS; 
 #    }
  
 #    for( i in 1:N){       
 #     Res[i] = Q[D[i]];
 #    }

 #    return(Q);    
 #    #return(Res);          
 # }


 # MbF(Males, S, N, J)
 # aggregate(Males~S, FUN=mean)
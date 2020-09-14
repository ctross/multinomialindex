sim_rs = function(N=100, Rate=3.3, L=0, U=1, et=1, es=0.3, j){
	set.seed(1) # Set seed for exposure time.
    Time = 3*runif(N,L,1)^U  # Simulate the age of each individual.  
	
    set.seed(j) # Reset seed for RS. Note that exposure time is held constant over replicates, while RS is allowed to vary.
           B = 0.5               # Gamma scalar.
           Z = rgamma(N,1*B,B)   # Random component of RS based on variation in quality.

           r = rep(NA,N)         # Initialize RS storage.

           # Loop over individuals, and simulate RS given exposre time. Negative-binomial model, written here as a Gamma-Poisson mixture.
           for(i in 1:N) {
            mu = Rate*(Time[i]^et)*(Z[i]^es)  
            r[i] = rpois(1, mu) 
           }
           
           res = list(r=r, t=Time) # Store results
           
           return(res)
           }

########################################################## Now check math in main paper
Reps = 10000 
res = matrix(NA,nrow=Reps,ncol=4)

for(j in 1:Reps){
 x = sim_rs(N=1000, Rate=7.0, L=0.1, U=1, et=1, es=0.39, j=j+1)

 r = x$r
 t = x$t
 N = length(r)
 R = sum(r)
 T = sum(t)
 t_hat = t/T
 r_hat = r/R

 res[j,1] = (N^2/R^2)*(1/N)*sum((r-R*t_hat)^2)         # 4.2a
 res[j,2] = N*sum((r_hat-t_hat)^2)                     # 4.2b
 res[j,3] = (N^2/R^2)*(var(r)-var(R*t_hat))            # 4.2c
 res[j,4] = (N^2/R^2)*var(r)*(1-cor(r,R*t_hat)^2)      # 4.2d
 }

print(colMeans(res))

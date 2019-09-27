functions{
 //# Function to calculate mean by factor
  vector MbF (vector X, int[] D, int N, int J){
  
   vector[J] Q;   //# Local storage, for each species in 1 to J
   vector[N] Res; //# Species mean for each data row i in 1 to N
  
   real ScrapV;   //# Local storage of values 
   real ScrapS;   //# Local storage of summands 
  
     for(j in 1:J){ //# For each species j
      ScrapV = 0;   //# Init counters at zero
      ScrapS = 0;
   
      for( i in 1:N){ //# Then loop through all data points and sum the values of the variable
       ScrapV = (D[i]==j)? X[i]+ScrapV: ScrapV; //# This is ifelse() function here
       ScrapS = (D[i]==j)? 1+ScrapS: ScrapS;   
      }
    
     Q[j] = ScrapV/ScrapS; //# For each species j, this is the mean value of X
    }
  
    for( i in 1:N){        //# Now, make each species mean back to the relevant data point
     Res[i] = Q[D[i]];
    }

    return(Res);           //# And return result
 }
}

data{ 
  int N;    //# Number of data points
  int J;    //# Number of species
  int K;    //# Number of regression parameters 
  
  int S[N]; //# Species ID of data point i

  vector[N] M_Rough; //# M (or B) values, with missing indicators
  
  matrix[J, J] Distance; //# Phylogenetic distance between species x and y
  
  vector[N] Males;           //# Count of males in group
  vector[N] Females_Rough;   //# Count of females in group, with missing indicators
  vector[N] Cops_Rough;      //# Count of copulations slash RS in group, with missing indicators
  vector[N] Seas_Rough;      //# Breeding season duration, with missing indicators
  vector[N] EstDur_Rough;    //# Estrous duration, with missing indicators
  vector[N] EstOvO_Rough;    //# Observed estrous overlap, with missing indicators
  vector[N] EstOvE_Rough;    //# Estimated estrous overlap, with missing indicators
  vector[N] Dispersal_Rough; //# Male dispersal, with missing indicators
  vector[N] Sync_Rough;      //# Synchrony index, with missing indicators
  
  int N_m;                  //# Number of missing values for each predictor
  int N_fems;               //#
  int N_cops;               //#
  int N_seas;               //#
  int N_dur;                //#
  int N_ovo;                //#
  int N_ove;                //#
  int N_disp;               //#
  int N_sync;               //#
  
  int Q_m [N_m];            //# Locations of missing values for each predictor
  int Q_fems [N_fems];      //#
  int Q_cops [N_cops];      //#
  int Q_seas [N_seas];      //#
  int Q_dur [N_dur];        //#
  int Q_ovo [N_ovo];        //#
  int Q_ove [N_ove];        //# 
  int Q_disp [N_disp];      //# 
  int Q_sync [N_sync];      //# 

  real Pars[K];             //# Which Pars in regression to include
 }
 
parameters{ 
  real Alpha;                 //# Intercept
  vector[K] Beta_Raw;         //# Slopes
  real<lower=0> Sigma;        //# Residual variance

  vector[N_m] M_I;            //# Missing data parameters
  vector[N_fems] Females_I;   //#
  vector[N_cops] Cops_I;      //#
  vector[N_seas] Seas_I;      //#
  vector[N_dur] EstDur_I;     //#
  vector[N_ovo] EstOvO_I;     //#
  vector[N_ove] EstOvE_I;     //#
  vector<lower=0,upper=1>[N_disp] Dispersal_I; //#
  vector[N_sync] Sync_I;      //#
  
  vector[J] Z_Rho;          //# Correlated random effects from phyogentic distance matrix, centered
  vector[J] Z_Gamma;        //# Species-level random effects, centered
  
  real<lower=0> SD_Rho;     //# Dispersion of random effects vectors
  real<lower=0> SD_Gamma;   //# 
   
  real<lower=0> Zeta;       //# Correlation decay rate
 } 

transformed parameters{ 
  vector[J] Gamma;          //# Correlated random effects from phyogentic distance matrix
  vector[J] Rho;            //# Species-level random effects

  vector[N] M;              //# Variables after missings have been replaced with parameters
  vector[N] Females;        //#
  vector[N] Cops;           //#
  vector[N] Seas;           //#
  vector[N] EstDur;         //#
  vector[N] EstOvO;         //#
  vector[N] EstOvE;         //#
  vector[N] Dispersal;      //#
  vector[N] Sync;           //#
  
  vector[N] Males_SM;       //# Species means for each predictor
  vector[N] Females_SM;     //#
  vector[N] Cops_SM;        //#
  vector[N] Seas_SM;        //#
  vector[N] EstDur_SM;      //#
  vector[N] EstOvO_SM;      //#
  vector[N] EstOvE_SM;      //#

  vector[K] Beta;           //# Slopes
  
  matrix[J,J] L;            //# Cholesky factor of phylogenetic covariance matrix
  
  for (i in 1:(J-1)){       //# Build L from Distance data and parameters
  for (j in (i+1):J){
                L[i,j] = exp( -Zeta * Distance[i,j]);     
                L[j,i] = L[i,j];                          
                       }}

 for (i in 1:J){
                L[i,i] = 1;                                             
                   }
  
  //############################# Now define the random effects vectors
  Gamma = SD_Gamma * Z_Gamma;
  Rho = SD_Rho * (cholesky_decompose(L) * Z_Rho);
  
  //############################# Here we map data and missing-data parameters into the final arrays
  M = M_Rough;
  for(i in 1:N_m) M[Q_m[i]] = M_I[i];
  
  Females = Females_Rough;
  for(i in 1:N_fems) Females[Q_fems[i]] = Females_I[i];
  
  Cops = Cops_Rough;
  for(i in 1:N_cops) Cops[Q_cops[i]] = Cops_I[i];
  
  Seas = Seas_Rough;
  for(i in 1:N_seas) Seas[Q_seas[i]] = Seas_I[i];
  
  EstDur = EstDur_Rough;
  for(i in 1:N_dur) EstDur[Q_dur[i]] = EstDur_I[i];
  
  EstOvO = EstOvO_Rough;
  for(i in 1:N_ovo) EstOvO[Q_ovo[i]] = EstOvO_I[i];

  EstOvE = EstOvE_Rough;
  for(i in 1:N_ove) EstOvE[Q_ove[i]] = EstOvE_I[i]; 

  Dispersal = Dispersal_Rough;
  for(i in 1:N_disp) Dispersal[Q_disp[i]] = Dispersal_I[i]; 

  Sync = Sync_Rough;
  for(i in 1:N_sync) Sync[Q_sync[i]] = Sync_I[i]; 
  
  //############################# Here we get species means accounting for missing data parameters
  Males_SM = MbF(Males, S, N, J);
  Females_SM = MbF(Females, S, N, J);
  Cops_SM = MbF(Cops, S, N, J);
  Seas_SM = MbF(Seas, S, N, J);
  EstDur_SM = MbF(EstDur, S, N, J);
  EstOvO_SM = MbF(EstOvO, S, N, J);
  EstOvE_SM = MbF(EstOvE, S, N, J);

//# Select pars to exclude from regresion
for(k in 1:K)
Beta[k] = Beta_Raw[k]*Pars[k] ;
 } 

model{ 
  vector[N] Mu;             //# Local storage                          

  Zeta ~ normal(0,1);       //# Decay rate prior is folded normal
  
  Z_Rho ~ normal(0,1);      //# Centered random effects are unit normals
  Z_Gamma ~ normal(0,1);
  
  SD_Rho ~ cauchy(0,1);     //# Random effects dispersion is folded Cauchy
  SD_Gamma ~ cauchy(0,1); 
  
  Alpha ~ normal(0,5);      //# Intercept prior
  Beta_Raw ~ normal(0,5);   //# Slope priors
  Sigma ~ cauchy(0,1);      //# Residual variance prior
  
  Females_I ~ normal(0,1);   //# Priors on normalized missing predictors
  Cops_I ~ normal(0,1);      //#
  Seas_I ~ normal(0,1);      //#
  EstDur_I ~ normal(0,1);    //#
  EstOvO_I ~ normal(0,1);    //#
  EstOvE_I ~ normal(0,1);    //#
  Dispersal_I ~ beta(1,1);   //#
  Sync_I ~ normal(0,1);      //#
  M_I ~ normal(0,1);         //#
  
//###################### Now just define a linear model
  for(i in 1:N){
  Mu[i] = Alpha + Gamma[S[i]] + Rho[S[i]] + Beta[1]*(Males_SM[i]) + Beta[2]*(Females_SM[i]) + 
                                            Beta[3]*(Cops_SM[i])  + Beta[4]*(Seas_SM[i]) +
                                            Beta[5]*(EstDur_SM[i])+ Beta[6]*(EstOvO_SM[i]) +
                                            Beta[7]*(Dispersal[i]) + Beta[8]*(Sync[i]) +
                                            Beta[9]*(Males[i]-Males_SM[i]) + Beta[10]*(Females[i]-Females_SM[i]) + 
                                            Beta[11]*(Cops[i]-Cops_SM[i])  + Beta[12]*(Seas[i]-Seas_SM[i]) +
                                            Beta[13]*(EstDur[i]-EstDur_SM[i])+ Beta[14]*(EstOvO[i]-EstOvO_SM[i]);
                                            }
  
//# And then model the outcomes
  M ~ normal(Mu, Sigma);
} 

generated quantities{
    vector[N] Resid;   

               for(i in 1:N){
  Resid[i] = M[i] -( Alpha + Gamma[S[i]] + Rho[S[i]] + Beta[1]*(Males_SM[i]) + Beta[2]*(Females_SM[i]) + 
                                            Beta[3]*(Cops_SM[i])  + Beta[4]*(Seas_SM[i]) +
                                            Beta[5]*(EstDur_SM[i])+ Beta[6]*(EstOvO_SM[i]) +
                                            Beta[7]*(Dispersal[i]) + Beta[8]*(Sync[i]) +
                                            Beta[9]*(Males[i]-Males_SM[i]) + Beta[10]*(Females[i]-Females_SM[i]) + 
                                            Beta[11]*(Cops[i]-Cops_SM[i])  + Beta[12]*(Seas[i]-Seas_SM[i]) +
                                            Beta[13]*(EstDur[i]-EstDur_SM[i])+ Beta[14]*(EstOvO[i]-EstOvO_SM[i]));
                                            }

}


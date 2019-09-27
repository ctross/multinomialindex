# load the data and phylogeny:
set.seed(1)
kutsuskew_data = read.csv("Data/KutsukakeNunn_SOMTable.csv")
kutsuskew_tree = read.nexus("Data/KutsukakeSkew_consensusTree_10kTrees_Primates_Version3.nex") 


for(i in 1:length(kutsuskew_data$B_index)){
if( !is.na(kutsuskew_data$N_copulations_observed[i]) & !is.na(kutsuskew_data$N_male[i]))
kutsuskew_data$M_index[i] <- M_index_from_B_index(kutsuskew_data$B_index[i],kutsuskew_data$N_copulations_observed[i],kutsuskew_data$N_male[i])
}

kutsuskew_tree_onechimp<- kutsuskew_tree
kutsuskew_tree_onechimp$tip.label[kutsuskew_tree_onechimp$tip.label=="Pan_troglodytes_verus"]<- "Pan_troglodytes"
kutsuskew_tree_onechimp<- drop.tip(kutsuskew_tree_onechimp, "Pan_troglodytes_schweinfurthii")

Distance <- cophenetic(kutsuskew_tree_onechimp)
Distance <- Distance/max(Distance)

S <- length(kutsuskew_data$Species)
for(i in 1:length(kutsuskew_data$Species))
S[i] <- which(kutsuskew_data$Species[i] == colnames(Distance))

standardize <- function(x){y = (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE);return(y)}

Females <- standardize(kutsuskew_data$N_female)
N_fems <- sum(is.na(Females))
Q_fems <- which(is.na(Females))
Females[is.na(Females)] <- -99
Females_Raw <- Females

Cops <- standardize(kutsuskew_data$N_copulations_observed)
N_cops <- sum(is.na(Cops))
Q_cops <- which(is.na(Cops))
Cops[is.na(Cops)] <- -99
Cops_Raw <- Cops

Seas <- standardize(kutsuskew_data$breeding_season_duration)
N_seas <- sum(is.na(Seas))
Q_seas <- which(is.na(Seas))
Seas[is.na(Seas)] <- -99
Seas_Raw <- Seas

EstDur <- standardize(kutsuskew_data$duration_estrous)
N_dur <- sum(is.na(EstDur))
Q_dur <- which(is.na(EstDur))
EstDur[is.na(EstDur)] <- -99
EstDur_Raw <- EstDur

EstOvO <- standardize(kutsuskew_data$observed_estrous_overlap)
N_ovo <- sum(is.na(EstOvO))
Q_ovo <- which(is.na(EstOvO))
EstOvO[is.na(EstOvO)] <- -99
EstOvO_Raw <- EstOvO

EstOvE <- standardize(kutsuskew_data$expected_estrus_overlap)
N_ove <- sum(is.na(EstOvE))
Q_ove <- which(is.na(EstOvE))
EstOvE[is.na(EstOvE)] <- -99
EstOvE_Raw <- EstOvE

Dispersal <- ifelse(kutsuskew_data$male_dispersal_pattern=="P",1,0)
N_disp <- sum(is.na(Dispersal))
Q_disp <- which(is.na(Dispersal))
Dispersal[is.na(Dispersal)] <- -99
Dispersal_Raw <- Dispersal

Sync <- standardize(kutsuskew_data$synchrony_index)
N_sync <- sum(is.na(Sync))
Q_sync <- which(is.na(Sync))
Sync[is.na(Sync)] <- -99
Sync_Raw <- Sync

Males <- kutsuskew_data$N_male

N <- length(Males)
K <- 14
J <- max(S)


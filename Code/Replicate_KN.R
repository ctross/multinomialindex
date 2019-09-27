# load the data and phylogeny:
set.seed(1)
kutsuskew_data = read.csv("Data/KutsukakeNunn_SOMTable.csv")

kutsuskew_tree = read.nexus("Data/KutsukakeSkew_consensusTree_10kTrees_Primates_Version3.nex") 
plot(kutsuskew_tree)

for(i in 1:length(kutsuskew_data$B_index)){
if( !is.na(kutsuskew_data$N_copulations_observed[i]) & !is.na(kutsuskew_data$N_male[i]))
kutsuskew_data$M_index[i] <- M_index_from_B_index(kutsuskew_data$B_index[i],kutsuskew_data$N_copulations_observed[i],kutsuskew_data$N_male[i])
}

## generate species averages
kutsuskew_data_means<- aggregate(kutsuskew_data, by=list(kutsuskew_data$Species), FUN=mean, na.rm=TRUE)

write.csv(cor(data.frame(kutsuskew_data_means$M_index,kutsuskew_data_means$B_index, kutsuskew_data_means$maximum_proportion_mating, kutsuskew_data_means$lambda),use="pairwise.complete"),"Cors.csv")

# merge male dispersal pattern back on
Disp <- ifelse(table(as.character(kutsuskew_data$Species),as.character(kutsuskew_data$male_dispersal_pattern))[,1]>0,"D","P")
kutsuskew_data_means <-  kutsuskew_data_means[,-which(names(kutsuskew_data_means) %in% c("Species","male_dispersal_pattern"))]
colnames(kutsuskew_data_means)[1]<- "Species"
Disp <- data.frame(Disp)
Disp$Species <- rownames(Disp)
colnames(Disp) <- c("male_dispersal_pattern","Species")

kutsuskew_data_means<- merge(kutsuskew_data_means, Disp, By=Species, all.x=TRUE)

## delete one chimp subspecies from tree, rename other
kutsuskew_tree_onechimp<- kutsuskew_tree
kutsuskew_tree_onechimp$tip.label[kutsuskew_tree_onechimp$tip.label=="Pan_troglodytes_verus"]<- "Pan_troglodytes"
kutsuskew_tree_onechimp<- drop.tip(kutsuskew_tree_onechimp, "Pan_troglodytes_schweinfurthii")

pdf("Tree_Full.pdf")
plot(kutsuskew_tree_onechimp)
dev.off()

## tree and species match?
setdiff(kutsuskew_tree_onechimp$tip.label,as.character(kutsuskew_data_means$Species))

## set all branch lengths to 1 and log transform all variables (see Kutsu paper p. 699)
kutsuskew_tree_onechimp_all1<- kutsuskew_tree_onechimp
kutsuskew_tree_onechimp_all1$edge.length<- rep(1, 59)

pdf("Tree_Ones.pdf")
plot(kutsuskew_tree_onechimp_all1)
dev.off()

kutsuskew_data_means$lognmale<- log(kutsuskew_data_means$N_male)
kutsuskew_data_means$lognfemale<- log(kutsuskew_data_means$N_female)
kutsuskew_data_means$loglambda<- log(kutsuskew_data_means$lambda)
kutsuskew_data_means$logmaxpop<- log(kutsuskew_data_means$maximum_proportion_mating)
kutsuskew_data_means$logrankmost<- log(kutsuskew_data_means$rank_most_successful_male)
kutsuskew_data_means$logncop<- log(kutsuskew_data_means$N_copulations_observed)
kutsuskew_data_means$logseas<- log(kutsuskew_data_means$breeding_season_duration)
kutsuskew_data_means$logestrdur<- log(kutsuskew_data_means$duration_estrous)
kutsuskew_data_means$logobsestrover<- log(kutsuskew_data_means$observed_estrous_overlap)
kutsuskew_data_means$logexpestrover<- log(kutsuskew_data_means$expected_estrus_overlap)
kutsuskew_data_means$M<- kutsuskew_data_means$M_index

#################################
### reproduce Kutsukake paper ###
#################################

## place tree and data in new R object as required by caper
kutsuskew = comparative.data(phy = kutsuskew_tree_onechimp_all1, data = kutsuskew_data_means, names.col = Species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
 

### Table 1: univariate ###
  
## using full datasets ##

## Lambda ##
lambda.male = crunch(loglambda ~ lognmale, data = kutsuskew, robust=3)
lambda.female = crunch(loglambda ~ lognfemale, data = kutsuskew, robust=3) 
lambda.cop = crunch(loglambda ~ logncop, data = kutsuskew, robust=3) 
lambda.seas = crunch(loglambda ~ logseas, data = kutsuskew, robust=3) 
lambda.estrdur = crunch(loglambda ~ logestrdur, data = kutsuskew, robust=3) 
lambda.estrexpover = crunch(loglambda ~ logexpestrover, data = kutsuskew, robust=3) 
lambda.estrobsover = crunch(loglambda ~ observed_estrous_overlap, data = kutsuskew, robust=3) 
lambda.synchrony = crunch(loglambda ~ synchrony_index, data = kutsuskew, robust=3) 
lambda.dispersal = brunch(loglambda ~ male_dispersal_pattern, data = kutsuskew, robust=3) 

res1 <- rbind(summary(lambda.male)$coefficients, summary(lambda.female)$coefficients, 
summary(lambda.cop)$coefficients, summary(lambda.seas)$coefficients, summary(lambda.estrdur)$coefficients,
summary(lambda.estrexpover)$coefficients, summary(lambda.estrobsover)$coefficients,
summary(lambda.synchrony)$coefficients, summary(lambda.dispersal)$coefficients)

res2 <- as.numeric(rbind(summary(lambda.male)$df[2], summary(lambda.female)$df[2], 
summary(lambda.cop)$df[2], summary(lambda.seas)$df[2], summary(lambda.estrdur)$df[2],
summary(lambda.estrexpover)$df[2], summary(lambda.estrobsover)$df[2],
summary(lambda.synchrony)$df[2], summary(lambda.dispersal)$df[2]))

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.", "Dur. estrous", 
                   "Exp. estr. overlap", "Obs. estr. overlap", "Synchrony", "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_lambda <- res


## logmaxpop ##
maxpop.male = crunch(logmaxpop ~ lognmale, data = kutsuskew, robust=3)
maxpop.female = crunch(logmaxpop ~ lognfemale, data = kutsuskew, robust=3) 
maxpop.cop = crunch(logmaxpop ~ logncop, data = kutsuskew, robust=3) 
maxpop.seas = crunch(logmaxpop ~ logseas, data = kutsuskew, robust=3) 
maxpop.estrdur = crunch(logmaxpop ~ logestrdur, data = kutsuskew, robust=3) 
maxpop.estrexpover = crunch(logmaxpop ~ logexpestrover, data = kutsuskew, robust=3) 
maxpop.estrobsover = crunch(logmaxpop ~ observed_estrous_overlap, data = kutsuskew, robust=3) 
maxpop.synchrony = crunch(logmaxpop ~ synchrony_index, data = kutsuskew, robust=3) 
maxpop.dispersal = brunch(logmaxpop ~ male_dispersal_pattern, data = kutsuskew, robust=3) 

res1 <- rbind(summary(maxpop.male)$coefficients, summary(maxpop.female)$coefficients, 
summary(maxpop.cop)$coefficients, summary(maxpop.seas)$coefficients, summary(maxpop.estrdur)$coefficients,
summary(maxpop.estrexpover)$coefficients, summary(maxpop.estrobsover)$coefficients,
summary(maxpop.synchrony)$coefficients, summary(maxpop.dispersal)$coefficients)

res2 <- as.numeric(rbind(summary(maxpop.male)$df[2], summary(maxpop.female)$df[2], 
summary(maxpop.cop)$df[2], summary(maxpop.seas)$df[2], summary(maxpop.estrdur)$df[2],
summary(maxpop.estrexpover)$df[2], summary(maxpop.estrobsover)$df[2],
summary(maxpop.synchrony)$df[2], summary(maxpop.dispersal)$df[2]))

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.", "Dur. estrous", 
                   "Exp. estr. overlap", "Obs. estr. overlap", "Synchrony", "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_mmp <- res



## B_index ##
B.male = crunch(B_index ~ lognmale, data = kutsuskew, robust=3)
B.female = crunch(B_index ~ lognfemale, data = kutsuskew, robust=3) 
B.cop = crunch(B_index ~ logncop, data = kutsuskew, robust=3) 
B.seas = crunch(B_index ~ logseas, data = kutsuskew, robust=3) 
B.estrdur = crunch(B_index ~ logestrdur, data = kutsuskew, robust=3) 
B.estrexpover = crunch(B_index ~ logexpestrover, data = kutsuskew, robust=3) 
B.estrobsover = crunch(B_index ~ observed_estrous_overlap, data = kutsuskew, robust=3) 
B.synchrony = crunch(B_index ~ synchrony_index, data = kutsuskew, robust=3) 
B.dispersal = brunch(B_index ~ male_dispersal_pattern, data = kutsuskew, robust=3) 

res1 <- rbind(summary(B.male)$coefficients, summary(B.female)$coefficients, 
summary(B.cop)$coefficients, summary(B.seas)$coefficients, summary(B.estrdur)$coefficients,
summary(B.estrexpover)$coefficients, summary(B.estrobsover)$coefficients,
summary(B.synchrony)$coefficients, summary(B.dispersal)$coefficients)

res2 <- as.numeric(rbind(summary(B.male)$df[2], summary(B.female)$df[2], 
summary(B.cop)$df[2], summary(B.seas)$df[2], summary(B.estrdur)$df[2],
summary(B.estrexpover)$df[2], summary(B.estrobsover)$df[2],
summary(B.synchrony)$df[2], summary(B.dispersal)$df[2]))

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.", "Dur. estrous", 
                   "Exp. estr. overlap", "Obs. estr. overlap", "Synchrony", "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_b <- res


## M_index ##
M.male = crunch(M_index ~ lognmale, data = kutsuskew, robust=3)
M.female = crunch(M_index ~ lognfemale, data = kutsuskew, robust=3) 
M.cop = crunch(M_index ~ logncop, data = kutsuskew, robust=3) 
M.seas = crunch(M_index ~ logseas, data = kutsuskew, robust=3) 
M.estrdur = crunch(M_index ~ logestrdur, data = kutsuskew, robust=3) 
M.estrexpover = crunch(M_index ~ logexpestrover, data = kutsuskew, robust=3) 
M.estrobsover = crunch(M_index ~ observed_estrous_overlap, data = kutsuskew, robust=3) 
M.synchrony = crunch(M_index ~ synchrony_index, data = kutsuskew, robust=3) 
M.dispersal = brunch(M_index ~ male_dispersal_pattern, data = kutsuskew, robust=3) 

res1 <- rbind(summary(M.male)$coefficients, summary(M.female)$coefficients, 
summary(M.cop)$coefficients, summary(M.seas)$coefficients, summary(M.estrdur)$coefficients,
summary(M.estrexpover)$coefficients, summary(M.estrobsover)$coefficients,
summary(M.synchrony)$coefficients, summary(M.dispersal)$coefficients)

res2 <- as.numeric(rbind(summary(M.male)$df[2], summary(M.female)$df[2], 
summary(M.cop)$df[2], summary(M.seas)$df[2], summary(M.estrdur)$df[2],
summary(M.estrexpover)$df[2], summary(M.estrobsover)$df[2],
summary(M.synchrony)$df[2], summary(M.dispersal)$df[2]))

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.", "Dur. estrous", 
                   "Exp. estr. overlap", "Obs. estr. overlap", "Synchrony", "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_m <- res


res <- cbind(
        c(rep("Lambda",length(rownames(res))),rep("MMP",length(rownames(res))),rep("B",length(rownames(res))),rep("M",length(rownames(res)))),
        round(rbind(res_lambda,res_mmp,res_b,res_m),2)
        )
        
res[,6] <- ifelse(res[,6]=="0","\\textbf{$<$0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.01","\\textbf{0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.02","\\textbf{0.02}",res[,6])
res[,6] <- ifelse(res[,6]=="0.03","\\textbf{0.03}",res[,6])
res[,6] <- ifelse(res[,6]=="0.04","\\textbf{0.04}",res[,6])
res[,6] <- ifelse(res[,6]=="0.05","\\textbf{0.05}",res[,6])

colnames(res)[1] <- "Metric"

write.csv(res, "Table1.csv", quote = FALSE)



## using only complete case dataset ##
kutsuskew_data_means_complete<- kutsuskew_data_means[!is.na(kutsuskew_data_means$M_index),]
kutsuskew.complete = comparative.data(phy = kutsuskew_tree_onechimp_all1, data = kutsuskew_data_means_complete, 
                                      names.col = Species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)


## Lambda ##
lambda.male = crunch(loglambda ~ lognmale, data = kutsuskew.complete, robust=3)
lambda.female = crunch(loglambda ~ lognfemale, data = kutsuskew.complete, robust=3) 
lambda.cop = crunch(loglambda ~ logncop, data = kutsuskew.complete, robust=3) 
lambda.seas = crunch(loglambda ~ logseas, data = kutsuskew.complete, robust=3) 
lambda.estrdur = crunch(loglambda ~ logestrdur, data = kutsuskew.complete, robust=3) 
lambda.estrexpover = crunch(loglambda ~ logexpestrover, data = kutsuskew.complete, robust=3) 
lambda.estrobsover = crunch(loglambda ~ observed_estrous_overlap, data = kutsuskew.complete, robust=3) 
lambda.synchrony = crunch(loglambda ~ synchrony_index, data = kutsuskew.complete, robust=3) 
lambda.dispersal = brunch(loglambda ~ male_dispersal_pattern, data = kutsuskew.complete, robust=3) 

res1 <- rbind(summary(lambda.male)$coefficients, summary(lambda.female)$coefficients, 
summary(lambda.cop)$coefficients, summary(lambda.seas)$coefficients, summary(lambda.estrdur)$coefficients,
summary(lambda.estrexpover)$coefficients, summary(lambda.estrobsover)$coefficients,
summary(lambda.synchrony)$coefficients, summary(lambda.dispersal)$coefficients)

res2 <- as.numeric(rbind(summary(lambda.male)$df[2], summary(lambda.female)$df[2], 
summary(lambda.cop)$df[2], summary(lambda.seas)$df[2], summary(lambda.estrdur)$df[2],
summary(lambda.estrexpover)$df[2], summary(lambda.estrobsover)$df[2],
summary(lambda.synchrony)$df[2], summary(lambda.dispersal)$df[2]))

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.", "Dur. estrous", 
                   "Exp. estr. overlap", "Obs. estr. overlap", "Synchrony", "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_lambda <- res


## logmaxpop ##
maxpop.male = crunch(logmaxpop ~ lognmale, data = kutsuskew.complete, robust=3)
maxpop.female = crunch(logmaxpop ~ lognfemale, data = kutsuskew.complete, robust=3) 
maxpop.cop = crunch(logmaxpop ~ logncop, data = kutsuskew.complete, robust=3) 
maxpop.seas = crunch(logmaxpop ~ logseas, data = kutsuskew.complete, robust=3) 
maxpop.estrdur = crunch(logmaxpop ~ logestrdur, data = kutsuskew.complete, robust=3) 
maxpop.estrexpover = crunch(logmaxpop ~ logexpestrover, data = kutsuskew.complete, robust=3) 
maxpop.estrobsover = crunch(logmaxpop ~ observed_estrous_overlap, data = kutsuskew.complete, robust=3) 
maxpop.synchrony = crunch(logmaxpop ~ synchrony_index, data = kutsuskew.complete, robust=3) 
maxpop.dispersal = brunch(logmaxpop ~ male_dispersal_pattern, data = kutsuskew.complete, robust=3) 

res1 <- rbind(summary(maxpop.male)$coefficients, summary(maxpop.female)$coefficients, 
summary(maxpop.cop)$coefficients, summary(maxpop.seas)$coefficients, summary(maxpop.estrdur)$coefficients,
summary(maxpop.estrexpover)$coefficients, summary(maxpop.estrobsover)$coefficients,
summary(maxpop.synchrony)$coefficients, summary(maxpop.dispersal)$coefficients)

res2 <- as.numeric(rbind(summary(maxpop.male)$df[2], summary(maxpop.female)$df[2], 
summary(maxpop.cop)$df[2], summary(maxpop.seas)$df[2], summary(maxpop.estrdur)$df[2],
summary(maxpop.estrexpover)$df[2], summary(maxpop.estrobsover)$df[2],
summary(maxpop.synchrony)$df[2], summary(maxpop.dispersal)$df[2]))

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.", "Dur. estrous", 
                   "Exp. estr. overlap", "Obs. estr. overlap", "Synchrony", "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_mmp <- res

## B_index ##
B.male = crunch(B_index ~ lognmale, data = kutsuskew.complete, robust=3)
B.female = crunch(B_index ~ lognfemale, data = kutsuskew.complete, robust=3) 
B.cop = crunch(B_index ~ logncop, data = kutsuskew.complete, robust=3) 
B.seas = crunch(B_index ~ logseas, data = kutsuskew.complete, robust=3) 
B.estrdur = crunch(B_index ~ logestrdur, data = kutsuskew.complete, robust=3) 
B.estrexpover = crunch(B_index ~ logexpestrover, data = kutsuskew.complete, robust=3) 
B.estrobsover = crunch(B_index ~ observed_estrous_overlap, data = kutsuskew.complete, robust=3) 
B.synchrony = crunch(B_index ~ synchrony_index, data = kutsuskew.complete, robust=3) 
B.dispersal = brunch(B_index ~ male_dispersal_pattern, data = kutsuskew.complete, robust=3) 

res1 <- rbind(summary(B.male)$coefficients, summary(B.female)$coefficients, 
summary(B.cop)$coefficients, summary(B.seas)$coefficients, summary(B.estrdur)$coefficients,
summary(B.estrexpover)$coefficients, summary(B.estrobsover)$coefficients,
summary(B.synchrony)$coefficients, summary(B.dispersal)$coefficients)

res2 <- as.numeric(rbind(summary(B.male)$df[2], summary(B.female)$df[2], 
summary(B.cop)$df[2], summary(B.seas)$df[2], summary(B.estrdur)$df[2],
summary(B.estrexpover)$df[2], summary(B.estrobsover)$df[2],
summary(B.synchrony)$df[2], summary(B.dispersal)$df[2]))

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.", "Dur. estrous", 
                   "Exp. estr. overlap", "Obs. estr. overlap", "Synchrony", "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_b <- res


## M_index ##
M.male = crunch(M_index ~ lognmale, data = kutsuskew.complete, robust=3)
M.female = crunch(M_index ~ lognfemale, data = kutsuskew.complete, robust=3) 
M.cop = crunch(M_index ~ logncop, data = kutsuskew.complete, robust=3) 
M.seas = crunch(M_index ~ logseas, data = kutsuskew.complete, robust=3) 
M.estrdur = crunch(M_index ~ logestrdur, data = kutsuskew.complete, robust=3) 
M.estrexpover = crunch(M_index ~ logexpestrover, data = kutsuskew.complete, robust=3) 
M.estrobsover = crunch(M_index ~ observed_estrous_overlap, data = kutsuskew.complete, robust=3) 
M.synchrony = crunch(M_index ~ synchrony_index, data = kutsuskew.complete, robust=3) 
M.dispersal = brunch(M_index ~ male_dispersal_pattern, data = kutsuskew.complete, robust=3) 

res1 <- rbind(summary(M.male)$coefficients, summary(M.female)$coefficients, 
summary(M.cop)$coefficients, summary(M.seas)$coefficients, summary(M.estrdur)$coefficients,
summary(M.estrexpover)$coefficients, summary(M.estrobsover)$coefficients,
summary(M.synchrony)$coefficients, summary(M.dispersal)$coefficients)

res2 <- as.numeric(rbind(summary(M.male)$df[2], summary(M.female)$df[2], 
summary(M.cop)$df[2], summary(M.seas)$df[2], summary(M.estrdur)$df[2],
summary(M.estrexpover)$df[2], summary(M.estrobsover)$df[2],
summary(M.synchrony)$df[2], summary(M.dispersal)$df[2]))


res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.", "Dur. estrous", 
                   "Exp. estr. overlap", "Obs. estr. overlap", "Synchrony", "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_m <- res

res <- cbind(
        c(rep("Lambda",length(rownames(res))),rep("MMP",length(rownames(res))),rep("B",length(rownames(res))),rep("M",length(rownames(res)))),
        round(rbind(res_lambda,res_mmp,res_b,res_m),2)
        )
        
res[,6] <- ifelse(res[,6]=="0","\\textbf{$<$0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.01","\\textbf{0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.02","\\textbf{0.02}",res[,6])
res[,6] <- ifelse(res[,6]=="0.03","\\textbf{0.03}",res[,6])
res[,6] <- ifelse(res[,6]=="0.04","\\textbf{0.04}",res[,6])
res[,6] <- ifelse(res[,6]=="0.05","\\textbf{0.05}",res[,6])

colnames(res)[1] <- "Metric"

write.csv(res, "Table2.csv", quote = FALSE)


### Table 2: multivariate ###

## using full data ##

## lambda
lambda.multi = crunch(loglambda ~ lognmale+lognfemale+logncop+logseas+logexpestrover+male_dispersal_pattern, data = kutsuskew, factor.action="allow", robust=3)
res1 <- summary(lambda.multi)$coefficients
res2 <- rep(summary(lambda.multi)$df[2], 6)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.",  
                   "Exp. estr. overlap",  "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_lambda <- res


## max prop
maxpop.multi = crunch(logmaxpop ~ lognmale+lognfemale+logncop+logseas+logexpestrover+male_dispersal_pattern, data = kutsuskew, factor.action="allow", robust=3)
res1 <- summary(maxpop.multi)$coefficients
res2 <- rep(summary(maxpop.multi)$df[2], 6)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.",  
                   "Exp. estr. overlap",  "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_mmp <- res


## B ##
B.multi = crunch(B_index ~ lognmale+lognfemale+logncop+logseas+logexpestrover+male_dispersal_pattern, data = kutsuskew, factor.action="allow", robust=3)
res1 <- summary(B.multi)$coefficients
res2 <- rep(summary(B.multi)$df[2], 6)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.",  
                   "Exp. estr. overlap",  "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_b <- res


## M ##
m.multi = crunch(M_index ~ lognmale+lognfemale+logncop+logseas+logexpestrover+male_dispersal_pattern, data = kutsuskew, factor.action="allow", robust=3)
res1 <- summary(m.multi)$coefficients
res2 <- rep(summary(m.multi)$df[2], 6)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.",  
                   "Exp. estr. overlap",  "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_m <- res


res <- cbind(
        c(rep("Lambda",length(rownames(res))),rep("MMP",length(rownames(res))),rep("B",length(rownames(res))),rep("M",length(rownames(res)))),
        round(rbind(res_lambda,res_mmp,res_b,res_m),2)
        )
        
res[,6] <- ifelse(res[,6]=="0","\\textbf{$<$0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.01","\\textbf{0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.02","\\textbf{0.02}",res[,6])
res[,6] <- ifelse(res[,6]=="0.03","\\textbf{0.03}",res[,6])
res[,6] <- ifelse(res[,6]=="0.04","\\textbf{0.04}",res[,6])
res[,6] <- ifelse(res[,6]=="0.05","\\textbf{0.05}",res[,6])

colnames(res)[1] <- "Metric"

write.csv(res, "Table3.csv", quote = FALSE)


## using only cases with all skew metrics ##

lambda.multi.complete = crunch(loglambda ~ lognmale+lognfemale+logncop+logseas+logexpestrover+male_dispersal_pattern, data = kutsuskew.complete, factor.action="allow", robust=3)
res1 <- summary(lambda.multi.complete)$coefficients
res2 <- rep(summary(lambda.multi.complete)$df[2], 6)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.",  
                   "Exp. estr. overlap",  "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_lambda <- res


maxpop.multi.complete = crunch(logmaxpop ~ lognmale+lognfemale+logncop+logseas+logexpestrover+male_dispersal_pattern, data = kutsuskew.complete, factor.action="allow", robust=3)
res1 <- summary(maxpop.multi.complete)$coefficients
res2 <- rep(summary(maxpop.multi.complete)$df[2], 6)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.",  
                   "Exp. estr. overlap",  "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_mmp <- res


B.multi.complete = crunch(B_index ~ lognmale+lognfemale+logncop+logseas+logexpestrover+male_dispersal_pattern, data = kutsuskew.complete, factor.action="allow", robust=3)
res1 <- summary(B.multi.complete)$coefficients
res2 <- rep(summary(B.multi.complete)$df[2], 6)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.",  
                   "Exp. estr. overlap",  "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_b <- res


M.multi.complete = crunch(M_index ~ lognmale+lognfemale+logncop+logseas+logexpestrover+male_dispersal_pattern, data = kutsuskew.complete, factor.action="allow", robust=3)
res1 <- summary(M.multi.complete)$coefficients
res2 <- rep(summary(M.multi.complete)$df[2], 6)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females", "N copulations", "Breed. Seas.",  
                   "Exp. estr. overlap",  "Male dispersal")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_m <- res


res <- cbind(
        c(rep("Lambda",length(rownames(res))),rep("MMP",length(rownames(res))),rep("B",length(rownames(res))),rep("M",length(rownames(res)))),
        round(rbind(res_lambda,res_mmp,res_b,res_m),2)
        )
        
res[,6] <- ifelse(res[,6]=="0","\\textbf{$<$0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.01","\\textbf{0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.02","\\textbf{0.02}",res[,6])
res[,6] <- ifelse(res[,6]=="0.03","\\textbf{0.03}",res[,6])
res[,6] <- ifelse(res[,6]=="0.04","\\textbf{0.04}",res[,6])
res[,6] <- ifelse(res[,6]=="0.05","\\textbf{0.05}",res[,6])

colnames(res)[1] <- "Metric"

write.csv(res, "Table4.csv", quote = FALSE)


### intra-specific analysis ###
kutsuskew_data_chimps <- kutsuskew_data[kutsuskew_data$Species=="Pan_troglodytes",]

lambda.chimps = lm(lambda ~ N_male+N_female, data = kutsuskew_data_chimps)
res1 <- summary(lambda.chimps)$coefficients[-1,]
res2 <- rep(summary(lambda.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_lambda1 <- res


maxpop.chimps = lm(maximum_proportion_mating ~ N_male+N_female, data = kutsuskew_data_chimps)
res1 <- summary(maxpop.chimps)$coefficients[-1,]
res2 <- rep(summary(maxpop.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_mmp1 <- res


B.chimps = lm(B_index ~ N_male+N_female, data = kutsuskew_data_chimps)
res1 <- summary(B.chimps)$coefficients[-1,]
res2 <- rep(summary(B.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_b1 <- res


M.chimps = lm(M_index ~ N_male+N_female, data = kutsuskew_data_chimps)
res1 <- summary(M.chimps)$coefficients[-1,]
res2 <- rep(summary(M.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_m1 <- res

res <- cbind(
        c(rep("Lambda",length(rownames(res))),rep("MMP",length(rownames(res))),rep("B",length(rownames(res))),rep("M",length(rownames(res)))),
        round(rbind(res_lambda1,res_mmp1,res_b1,res_m1),2)
        )
        
res[,6] <- ifelse(res[,6]=="0","\\textbf{$<$0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.01","\\textbf{0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.02","\\textbf{0.02}",res[,6])
res[,6] <- ifelse(res[,6]=="0.03","\\textbf{0.03}",res[,6])
res[,6] <- ifelse(res[,6]=="0.04","\\textbf{0.04}",res[,6])
res[,6] <- ifelse(res[,6]=="0.05","\\textbf{0.05}",res[,6])

resA <- res

colnames(resA)[1] <- "Metric"

lambda.chimps = lm(lambda ~ N_male+expected_estrus_overlap, data = kutsuskew_data_chimps)
res1 <- summary(lambda.chimps)$coefficients[-1,]
res2 <- rep(summary(lambda.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "Exp. estr. overlap")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_lambda2 <- res


maxpop.chimps = lm(maximum_proportion_mating ~ N_male+expected_estrus_overlap, data = kutsuskew_data_chimps)
res1 <- summary(maxpop.chimps)$coefficients[-1,]
res2 <- rep(summary(maxpop.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "Exp. estr. overlap")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_mmp2 <- res


B.chimps = lm(B_index ~ N_male+expected_estrus_overlap, data = kutsuskew_data_chimps)
res1 <- summary(B.chimps)$coefficients[-1,]
res2 <- rep(summary(B.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "Exp. estr. overlap")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_b2 <- res


M.chimps = lm(M_index ~ N_male+expected_estrus_overlap, data = kutsuskew_data_chimps)
res1 <- summary(M.chimps)$coefficients[-1,]
res2 <- rep(summary(M.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "Exp. estr. overlap")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_m2 <- res



res <- cbind(
        c(rep("Lambda",length(rownames(res))),rep("MMP",length(rownames(res))),rep("B",length(rownames(res))),rep("M",length(rownames(res)))),
        round(rbind(res_lambda2,res_mmp2,res_b2,res_m2),2)
        )
        
res[,6] <- ifelse(res[,6]=="0","\\textbf{$<$0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.01","\\textbf{0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.02","\\textbf{0.02}",res[,6])
res[,6] <- ifelse(res[,6]=="0.03","\\textbf{0.03}",res[,6])
res[,6] <- ifelse(res[,6]=="0.04","\\textbf{0.04}",res[,6])
res[,6] <- ifelse(res[,6]=="0.05","\\textbf{0.05}",res[,6])

resB <- res
colnames(resB)[1] <- "Metric"

write.csv(resA, "Table5a.csv", quote = FALSE)
write.csv(resB, "Table5b.csv", quote = FALSE)


## complete data only ##

kutsuskew_data_chimps_complete <- kutsuskew_data_chimps[!is.na(kutsuskew_data_chimps$M_index),]

lambda.chimps = lm(lambda ~ N_male+N_female, data = kutsuskew_data_chimps_complete)
res1 <- summary(lambda.chimps)$coefficients[-1,]
res2 <- rep(summary(lambda.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_lambda1 <- res


maxpop.chimps = lm(maximum_proportion_mating ~ N_male+N_female, data = kutsuskew_data_chimps_complete)
res1 <- summary(maxpop.chimps)$coefficients[-1,]
res2 <- rep(summary(maxpop.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_mmp1 <- res


B.chimps = lm(B_index ~ N_male+N_female, data = kutsuskew_data_chimps_complete)
res1 <- summary(B.chimps)$coefficients[-1,]
res2 <- rep(summary(B.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_b1 <- res


M.chimps = lm(M_index ~ N_male+N_female, data = kutsuskew_data_chimps_complete)
res1 <- summary(M.chimps)$coefficients[-1,]
res2 <- rep(summary(M.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "N females")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_m1 <- res

res <- cbind(
        c(rep("Lambda",length(rownames(res))),rep("MMP",length(rownames(res))),rep("B",length(rownames(res))),rep("M",length(rownames(res)))),
        round(rbind(res_lambda1,res_mmp1,res_b1,res_m1),2)
        )
        
res[,6] <- ifelse(res[,6]=="0","\\textbf{$<$0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.01","\\textbf{0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.02","\\textbf{0.02}",res[,6])
res[,6] <- ifelse(res[,6]=="0.03","\\textbf{0.03}",res[,6])
res[,6] <- ifelse(res[,6]=="0.04","\\textbf{0.04}",res[,6])
res[,6] <- ifelse(res[,6]=="0.05","\\textbf{0.05}",res[,6])

resA <- res

colnames(resA)[1] <- "Metric"

lambda.chimps = lm(lambda ~ N_male+expected_estrus_overlap, data = kutsuskew_data_chimps_complete)
res1 <- summary(lambda.chimps)$coefficients[-1,]
res2 <- rep(summary(lambda.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "Exp. estr. overlap")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_lambda2 <- res


maxpop.chimps = lm(maximum_proportion_mating ~ N_male+expected_estrus_overlap, data = kutsuskew_data_chimps_complete)
res1 <- summary(maxpop.chimps)$coefficients[-1,]
res2 <- rep(summary(maxpop.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "Exp. estr. overlap")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_mmp2 <- res


B.chimps = lm(B_index ~ N_male+expected_estrus_overlap, data = kutsuskew_data_chimps_complete)
res1 <- summary(B.chimps)$coefficients[-1,]
res2 <- rep(summary(B.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "Exp. estr. overlap")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_b2 <- res


M.chimps = lm(M_index ~ N_male+expected_estrus_overlap, data = kutsuskew_data_chimps_complete)
res1 <- summary(M.chimps)$coefficients[-1,]
res2 <- rep(summary(M.chimps)$df[2], 2)

res <- cbind(res1[,1:2],res2,res1[,3:4])

rownames(res) <- c("N males", "Exp. estr. overlap")
 
colnames(res) <- c("Slope", "SE", "DF", "T", "P")

res_m2 <- res



res <- cbind(
        c(rep("Lambda",length(rownames(res))),rep("MMP",length(rownames(res))),rep("B",length(rownames(res))),rep("M",length(rownames(res)))),
        round(rbind(res_lambda2,res_mmp2,res_b2,res_m2),2)
        )
        
res[,6] <- ifelse(res[,6]=="0","\\textbf{$<$0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.01","\\textbf{0.01}",res[,6])
res[,6] <- ifelse(res[,6]=="0.02","\\textbf{0.02}",res[,6])
res[,6] <- ifelse(res[,6]=="0.03","\\textbf{0.03}",res[,6])
res[,6] <- ifelse(res[,6]=="0.04","\\textbf{0.04}",res[,6])
res[,6] <- ifelse(res[,6]=="0.05","\\textbf{0.05}",res[,6])

resB <- res

colnames(resB)[1] <- "Metric"


write.csv(resA, "Table6a.csv", quote = FALSE)
write.csv(resB, "Table6b.csv", quote = FALSE)


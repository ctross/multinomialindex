# clean up workspace
rm(list=ls()) 

## load relevant packages
library(ape)
library(caper)
library(reshape2)
library(scales)
library(ggplot2)
library(SkewCalc)
library(ggridges)
library(viridis)


## set working directory
setwd("C:\\Users\\cody_ross\\Dropbox\\Completed and Published Projects\\1-papers\\Reproductive Skew\\Skew-Theory - V3\\Workflow")

# Check derivation
source("Code/Math_Check.R")

# Replicate K and N
source("Code/Replicate_KN.R")

# Simulate Figure 1
source("Code/Simulate_Skew.R")

# Simulate Figure 2
source("Code/Empirical_Test.R")

# Simulate Figure 3
source("Code/Bayes_Check.R")

# Last model
source("Code/Load_Phylo_Data.R")
library(rethinking)
source("Code/Run_Phylo.R")

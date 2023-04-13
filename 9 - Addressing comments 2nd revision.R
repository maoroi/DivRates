### Addressing reviewers comments from 15/3/23
library("readr")
library("ape")

tax <- read.csv("MDD_v1_6495species_JMamm.csv", header = TRUE)

alltrees <- read.nexus("MuHiSSE/trees/Data_S7_Mammalia_credibleTreeSets_tipDR/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_nexus.trees")

# 1. Assess advantage of the Cox 2021 dataset  ----------------------------

CoxSI2 <- read_csv("Cox2021_SupplementaryData2.csv")
# keep only AP columns
coxSI2 <- CoxSI2[,1:12]
# remove lines where activity data was imputed
cox <- coxSI2[which(coxSI2$Activity_source != "Imputed"),]
cox$Binomial <- paste(cox$Genus, cox$Species, sep='_')

# Retain unchanged binomials
unchanged <- which(cox$Binomial %in% tax$SciName)
cox$MDD <- NA
cox$MDD[unchanged] <- cox$Binomial[unchanged]

length(which(is.na(cox$MDD)))

# 2. Add data distribution figure -----------------------------------------


# 3. Add LTT plot (and more?) to explicitly account for time  -------------



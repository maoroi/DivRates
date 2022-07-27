setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")
# here::here()

library(hisse)  
#> packageVersion('hisse')
#[1] '1.9.13'
library(diversitree)
#> packageVersion('diversitree')
#[1] '0.9.11'
require(phytools)
#> packageVersion('phytools')
#[1] '0.6.99'

## reproducible evironemnt using {renv}: The packages used in this project are recorded into a lockfile, renv.lock.
## use renv::restore() to reinstall all of the packages as declared in the lockfile.
set.seed(88)

# 1 Estimating sampling fraction  -------------------------------------
# * 1.1 Proportion sampled out of all mammals --------------------------
## Proportion of sampled out of extant species by AP, assuming unbiased sampling. 

# total count in the entire dataset (inc species not in the phylogeny)
TotN <- length(which(data$AP == 'Nocturnal'))
TotC <- length(which(data$AP == 'Cathemeral'))
TotD <- length(which(data$AP == 'Diurnal'))
# proportions in data (sum up to 99%, crepuscular and others make the rest)
datN <- TotN / nrow(data)
datC <- TotC / nrow(data)
datD <- TotD / nrow(data)


# * 1.2 Assumptions ----------------------------------------------------
tax <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MDD_v1_6495species_JMamm.csv", header = TRUE)

# unsampled bats presumed nocturnal (n = 1096 in v1; 1122 in v1.1) 
bats_unsamp <- length(which(tax$Order == 'CHIROPTERA')) - length(which(data$Order == 'Chiroptera')) - 6  # unsampled bats minus 6 extinct spp

# unsampled squirrels (27 flying squirrels presumed noct; 87 other sciurids presumed diur)
SQRLS <- tax[which(tax$Family == 'SCIURIDAE'),]                         # all sciurids (all extant)
SQRLS_DATA <- data[which(data$Family == 'Sciuridae'),]                  # sampled sciurids 
flying <- length(which(SQRLS$Tribe == 'PTEROMYINI'))                    # 57 flying squirrels in MDD
noct_sqrls <- length(which(SQRLS_DATA$AP == "Nocturnal"))               # flying squirels are the only nocturnal sciurids
FLY_SQRLS_unsamp <- flying - noct_sqrls
Diur_sqrls_unsamp <- nrow(SQRLS) - nrow(SQRLS_DATA) - FLY_SQRLS_unsamp  # unsampled scuirids minus unsampled noct sciurids

# Simiiformes presumed diurnal (n = 148, excludes Aotus) 
PRIMS <- tax[which(tax$Order == "PRIMATES"),]
PRIMS <- PRIMS[which(PRIMS$extinct. == 0),]                             # only extant species
simians <- c("Atelidae", "Cebidae", "Cercopithecidae", "Pitheciidae", "Hominidae", "Hylobatidae")
simian_count <- length(which(PRIMS$Family %in% toupper(simians)))       # all simians recognised in MDD
SIMIAN_DATA <- data[which(data$Family %in% c(simians, "Aotidae")),]     # Aotidae demoted into Cebidae (Aotinae) in MDD but data$Family relates to MSW3   
# Out of 11 species of Aotus in MDD, 7 have AP info (no assumptions about remaining 4 spp)
Aotus_unsamp <- length(which(tax$Subfamily == "AOTINAE")) - length(which(data$Family == 'Aotidae'))
simi_diur_unsamp <- simian_count - nrow(SIMIAN_DATA) - Aotus_unsamp          # unsmapled (diur) simians minus unsampled (noct) Aotinae

# total unsampled (Noct=1123, Diur=235)
Nunsamp <- bats_unsamp + FLY_SQRLS_unsamp
Dunsamp <- Diur_sqrls_unsamp + simi_diur_unsamp
rm(bats_unsamp, FLY_SQRLS_unsamp, Diur_sqrls_unsamp, simi_diur_unsamp, PRIMS, simians, 
   simian_count, SIMIAN_DATA, Aotus_unsamp, noct_sqrls, flying, SQRLS_DATA, SQRLS)
##      ** end assumptions **       ##


# * 1.3 Estimate sampling fraction --------------------------------------
# species with no data or assumptions
all_extant <- length(unique(tax$SciName[which(tax$extinct. == 0)]))
unkn <- all_extant - (nrow(data) + Nunsamp + Dunsamp)                   # total number of extant species based on MDD v1
# expected AP distribution in unknown species if same proportion as in the data
expN <- datN * unkn
expC <- datC * unkn
expD <- datD * unkn
# proportion sampled out of likely total
Nprop <- TotN / (TotN + Nunsamp + expN)
Cprop <- TotC / (TotC + expC)
Dprop <- TotD / (TotD + Dunsamp + expD)

freq <- c(Nprop, Cprop, Dprop)

# Optional check: do nocturnal, cathemeral, and diurnal fractions account for 99% of mammals?
#AP_data <- length(which(data$AP %in% c("Nocturnal", "Cathemeral", "Diurnal"))) / nrow(data)
#AP_extrapol <- ((TotN + Nunsamp + expN) + (TotC + expC) + (TotD + Dunsamp + expD)) / all_extant
#round(AP_data, 2) == round(AP_extrapol, 2)

# 

# 2 Setting up the analysis (following MuHiSSE tutorial) --------------

for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
}

states <- data.frame(tree$tip.label, tree$tip.state, tree$tip.state)
states.trans <- states
for(i in 1:Ntip(tree)){
    if(states[i,2] == 1){
        states.trans[i,2] = 0
        states.trans[i,3] = 0
    }
    if(states[i,2] == 2){
        states.trans[i,2] = 0
        states.trans[i,3] = 1
    }
    if(states[i,2] == 3){
        states.trans[i,2] = 1
        states.trans[i,3] = 1
    }
}

## Setting up a MuSSE model using MuHiSSE
## As mentioned above, the number of free parameters in the model for both net turnover and extinction fraction
## are specified as index vectors provided to the function call. Each vector contains four entries that correspond
## to rates associated with the observed states (0 or 1) and the hidden states (A or B). They are always ordered
## as follows for a given hidden state, i: 00i, 01i, 10i, 11i. However, in this case we do not want any hidden
## states. But first let’s set up the “dull null” – i.e., turnover and extinction fraction are the same for all observed
## state combinations. 

# Note the “f” represents the sampling fraction for each observed state combination, which
# is a vector ordered in the same manner as for turnover and extinction fraction vectors:
f = c(freq[1], freq[2], 0, freq[3])

turnover <- c(1,1,0,1)
extinction.fraction <- c(1,1,0,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))
#Now, we can call MuHiSSE and estimate the parameters under this model using the default settings:
dull.null <- MuHiSSE(phy=tree, data=states.trans, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=FALSE,
                     trans.rate=trans.rate.mod, starting.vals=NULL)

# install code to allow using the console while models are running
remotes::install_github("lindeloev/job")

# MuSSE_1
turnover <- c(1,2,0,3)          
extinction.fraction <- c(1,1,0,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0, include.diagonals = FALSE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))   # disable transitions to/from fourth state

job::job(simple_MuSSE = {
    MuSSE_1 <- MuHiSSE(phy=tree, data=states.trans, f=f,
                       turnover=turnover,    
                       eps=extinction.fraction, 
                       hidden.states=FALSE, 
                       trans.rate=trans.rate.mod,
                       turnover.upper=1000, 
                       sann = TRUE)})
MuSSE_1 <- simple_MuSSE$MuSSE_1
save(MuSSE_1, file="MuSSE-1.Rdata")


# MuSSE_1A - extinction fraction free to vary
turnover <- c(1,2,0,3)          
extinction.fraction <-c(1,2,0,3) 
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0, include.diagonals = FALSE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))   # disable transitions to/from fourth state
MuSSE_eq <- MuHiSSE(phy=tree, data=states.trans, f=f, 
                    turnover=turnover,    
                    eps=extinction.fraction, 
                    hidden.states=FALSE, 
                    trans.rate=trans.rate.mod,
                    turnover.upper=1000, 
                    sann = TRUE)


# CID-8  (differential transitions among hidden states)
# *automating model specifications as function of hidden states* (the software can handle up to n=8 states)
for (n in 2:8) {
set_to_0 <- 4*(1:n)-1
#if (n <= 8) {
    drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)
#} else {stop(paste0("n too large"))}

# div rates for CID
for (z in 1:n) {reps <- rep(z, 4)
if (z == 1) {turnover <- reps} else {turnover <- c(turnover, reps)}}
extinction.fraction <- rep(1,4*n) 
turnover[set_to_0] <- extinction.fraction[set_to_0] <- 0

# trans rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = TRUE)
trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)])   
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA
}
# optional: identical transition rates across hidden states -- THINK CAREFULLY IF THIS MAKES SENSE!!
#EqTransAcrossHiddenStates <- function(rateMatrix, n_HiddenStates) {
#    lows <- 4*(1:(n_HiddenStates-1))+1
#    highs <- 4*(2:n_HiddenStates)
#    for (z in 1:(n_HiddenStates-1)) {
#        l <- lows[z]
#        h <- highs[z]
#        rateMatrix[l:h,l:h] <- rateMatrix[1:4,1:4]}
#    rateMatrix <- ParDrop(rateMatrix, c(0)) # update parameter indices (the indices are just names so this makes no difference to the model, just easier to read)
#}
#trans.rate.mod <- EqTransAcrossHiddenStates(trans.rate.mod, n) # identical transition rates across hidden states

# model 
CID_model <- MuHiSSE(phy=tree, data=states.trans, f=f,
                     turnover=turnover, 
                     eps=extinction.fraction, 
                     hidden.states=TRUE, 
                     trans.rate=trans.rate.mod, 
                     turnover.upper=1000, # excluding unreasonable parameter space to speed up (see new vignette)
                     sann = TRUE) # starting value optimiser on 
# MuHiSSE model
n <- 2
set_to_0 <- 4*(1:n)-1
TO_rates <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24)
drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)

# div rates for MuHiSSE
turnover <- TO_rates[1:(4*n)]
extinction.fraction <- rep(1,4*n) 
extinction.fraction[set_to_0] <- 0

# trans rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = TRUE)
trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)])   
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA


# 3 Extending rate matrix beyond 8 states ------------------------------------
# MuHiSSE model
n <- 2
set_to_0 <- 4*(1:n)-1
TO_rates <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24)
drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)

# div rates for MuHiSSE
turnover <- TO_rates[1:(4*n)]
extinction.fraction <- rep(1,4*n) 
extinction.fraction[set_to_0] <- 0

# trans rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = FALSE)
trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)])   
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA

## Adding onto the rate matrix
# columns
filler <- matrix(NA, ncol=4, nrow = nrow(trans.rate.mod))
trans.rate.app <- cbind(trans.rate.mod, filler)
colnames(trans.rate.app)[33:36] <- c("(00I)","(01I)","(10I)","(11I)")
# rows
filler <- matrix(NA, nrow=4, ncol = ncol(trans.rate.app))
trans.rate.app <- rbind(trans.rate.app, filler)
rownames(trans.rate.app)[33:36] <- c("(00I)","(01I)","(10I)","(11I)")

# claculating the correct indices is a pain if not starting from the source code 

## Additional analyses required by 2nd round of reviews NOT BASED ON COX 2021 DATASET

library(phytools)
library(readr)

setwd("C:/Users/rma10kg/Projects/CH2/AddMuHisse/DivRates")



# 1. Example of data distribution on a phylogeny --------------------------

tp <- read.nexus("tp5024.nex")
tp <- ladderize(tp, right=TRUE)
act <- read.table("APdata_2400spp.txt", header = TRUE)



# 2. Simmap reconstruction with LTT plot for timeline ---------------------



# 3. Run models of HiSSE with root.type="herr_als" ------------------------



# 4. Run a few examples of SecSSE to compare results  ---------------------



# 5. CorHMM reconstruction as well (if necessary) -------------------------



# 6. Plot evolutionary transition models (search gmail for code) ----------



# 7. Phylogenetic distances -----------------------------------------------

## if looking at the results of the best model from each tree does not provide good 
## enough a solution, phylogenetic distances can be compared based on Upham's backbone 
## tree (omitting tipwards "patches" that might hold most of the topological variation)






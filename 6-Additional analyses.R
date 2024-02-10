## Additional analyses required by 2nd round of reviews NOT BASED ON COX 2021 DATASET

library(phytools)
library(stringr)
library(readr)
library(dplyr)

setwd("C:/Users/rma10kg/Projects/CH2/AddMuHisse/DivRates")



# 1. Ensure data files are up to date -------------------------------------

tree <- read.nexus("tvar5024.nex")
tree <- ladderize(tree, right=TRUE)
# trim tip labels to binomial only
tree$tip.label <- str_extract(tree$tip.label,"[A-Z][a-z]*_[a-z]*")

act <- read.csv("ActivityData_MDD_v1.csv", header = TRUE)

## for consistency, the block below is copied from line 157 in file "2-Preparing for MuHiSSE.r"
# Matching binomials in phylogeny to taxonomy - Species names are updated in the column 'Phylo_name'
act$Phylo_name <- NA
act$Phylo_name_conforms_to <- NA

for (i in 1:nrow(act)){
   if (act$MDD[i] %in% tree$tip.label){ # if the MDD form is in the tree, it was retained
      act$Phylo_name[i] <- act$MDD[i]
      act$Phylo_name_conforms_to[i] <- "MDD_v1"
   } else if (act$MSW3[i] %in% tree$tip.label){ # if MDD isn't but the MSW3 form is - MSW3 was retained
      act$Phylo_name[i] <- act$MSW3[i]
      act$Phylo_name_conforms_to[i] <- "MSW3"
   } else {                # if neither taxonomies match, manually verified form was inserted (see below)
      act$Phylo_name[i] <- "-"
      act$Phylo_name_conforms_to[i] <- "-"
   }
}


# find the species for which both MSW3 and MDD are not in the tree
new_out <- which(!act$MDD %in% tree$tip.label)
old_out <- which(!act$MSW3 %in% tree$tip.label)


# correct tip labels in the tree, update the data table in the change
{act[which(act$MDD == "Bubalus_bubalis"), c(6,7)] <- c("not in tree", "_")
   tree$tip.label[which(tree$tip.label == "Nesotragus_moschatus")] <- "Neotragus_moschatus"    # typo
   act[which(act$MDD == "Neotragus_moschatus"), c(6,7)] <- c("Neotragus_moschatus", "typo")
   act[which(act$MDD == "Taurotragus_oryx"), c(6,7)] <- c("Tragelaphus_oryx", "MDD_v1.31")   # genus transfer in v1.31
   tree$tip.label[which(tree$tip.label ==  "Pseudalopex_culpaeus")] <- "Lycalopex_culpaeus"    # genus transfer v1
   act[which(act$MDD == "Lycalopex_culpaeus"), c(6,7)] <- c("Lycalopex_culpaeus", "not MDD_v1")
   tree$tip.label[which(tree$tip.label ==  "Pseudalopex_griseus")] <- "Lycalopex_griseus"      # genus transfer v1
   act[which(act$MDD == "Lycalopex_griseus"), c(6,7)] <- c("Lycalopex_griseus", "not MDD_v1")
   tree$tip.label[which(tree$tip.label ==  "Pseudalopex_gymnocercus")] <- "Lycalopex_gymnocercus"  # genus transfer v1
   act[which(act$MDD == "Lycalopex_gymnocercus"), c(6,7)] <- c("Lycalopex_gymnocercus", "not MDD_v1")
   tree$tip.label[which(tree$tip.label ==  "Pseudalopex_gymnocercus")] <- "Lycalopex_sechurae" # genus transfer v1
   act[which(act$MDD == "Lycalopex_sechurae"), c(6,7)] <- c("Lycalopex_sechurae", "not MDD_v1")
   tree$tip.label[which(tree$tip.label ==  "Pseudalopex_vetulus")] <- "Lycalopex_vetulus"      # genus transfer v1
   act[which(act$MDD == "Lycalopex_vetulus"), c(6,7)] <- c("Lycalopex_vetulus", "not MDD_v1")
   act[which(act$MDD == "Galerella_pulverulenta"), c(6,7)] <- c("Herpestes_pulverulentus", "MDD_v1.31")  # genus transfer in v1.31
   act[which(act$MDD == "Galerella_sanguinea"), c(6,7)] <- c("Herpestes_sanguineus", "MDD_v1.31")# genus transfer in v1.31
   tree$tip.label[which(tree$tip.label ==  "Lutra_maculicollis")] <- "Hydrictis_maculicollis"      # genus transfer v1
   act[which(act$MDD == "Hydrictis_maculicollis"), c(6,7)] <- c("Hydrictis_maculicollis", "not MDD_v1")
   tree$tip.label[which(tree$tip.label ==  "Monodelphis_unistriatus")] <- "Monodelphis_unistriata" # Latin grammar correction
   act[which(act$MDD == "Monodelphis_unistriata"), c(6,7)] <- c("Monodelphis_unistriata", "not MDD_v1")
   tree$tip.label[which(tree$tip.label ==  "Dactylonax_palpator")] <- "Dactylopsila_palpator"      # genus transfer v1
   act[which(act$MDD == "Dactylopsila_palpator"), c(6,7)] <- c("Dactylopsila_palpator", "not MDD_v1")
   tree$tip.label[which(tree$tip.label == "Elephantulus_revoilii")] <- "Elephantulus_revoili"      # spelling error
   act[which(act$MDD == "Elephantulus_revoili"), c(6,7)] <- c("Elephantulus_revoili", "typo")
   tree$tip.label[which(tree$tip.label == "Zaglossus_bruijnii")] <- "Zaglossus_bruijni"            # spelling error
   act[which(act$MDD == "Zaglossus_bruijni"), c(6,7)] <- c("Zaglossus_bruijni", "typo")
   tree$tip.label[which(tree$tip.label == "Equus_africanus")] <- "Equus_asinus"    # asinus is domestic E. africanus
   act[which(act$MDD == "Equus_asinus"), c(6,7)] <- c("Equus_asinus", "domestic form")
   act[which(act$MDD == "Cercopithecus_denti"), c(6,7)] <- c("not in tree", "split from C.pogonias")
   act[which(act$MDD == "Cercopithecus_wolfi"), c(6,7)] <- c("not in tree", "split from C.mitis")
   tree$tip.label[which(tree$tip.label ==  "Procolobus_badius")] <- "Piliocolobus_badius"      # genus transfer v1
   act[which(act$MDD == "Piliocolobus_badius"), c(6,7)] <- c("Piliocolobus_badius", " not MDD_v1")
   tree$tip.label[which(tree$tip.label ==  "Procolobus_kirkii")] <- "Piliocolobus_kirkii"      # genus transfer v1
   act[which(act$MDD == "Piliocolobus_kirkii"), c(6,7)] <- c("Piliocolobus_kirkii", "not MDD_v1")
   tree$tip.label[which(tree$tip.label == "Procolobus_pennantii")] <- "Piliocolobus_pennantii" # genus transfer v1
   act[which(act$MDD == "Piliocolobus_pennantii"), c(6,7)] <- c("Piliocolobus_pennantii", "not MDD_v1")
   tree$tip.label[which(tree$tip.label == "Procolobus_preussi")] <- "Piliocolobus_preussi"     # genus transfer v1
   act[which(act$MDD == "Piliocolobus_preussi"), c(6,7)] <- c("Piliocolobus_preussi", "not MDD_v1")
   tree$tip.label[which(tree$tip.label == "Procolobus_rufomitratus")] <- "Piliocolobus_rufomitratus"   # genus transfer v1
   act[which(act$MDD == "Piliocolobus_rufomitratus"), c(6,7)] <- c("Piliocolobus_rufomitratus", "not MDD_v1")
   act[which(act$MDD == "Loxodonta_cyclotis"), c(6,7)] <- c("not in tree", "_")
   tree$tip.label[which(tree$tip.label == "Pygeretmus_zhitkovi")] <- "Pygeretmus_shitkovi"     # correct spelling
   act[which(act$MDD == "Pygeretmus_shitkovi"), c(6,7)] <- c("Pygeretmus_shitkovi", "spelling")
   act[which(act$MDD == "Dephomys_eburneae"), c(6,7)] <- c("not in tree", "_")
   act[which(act$MDD == "Grammomys_poensis"), c(6,7)] <- c("not in tree", "_")
   tree$tip.label[which(tree$tip.label == "Aethomys_namaquensis")] <- "Micaelamys_namaquensis" # genus transfer v1
   act[which(act$MDD == "Micaelamys_namaquensis"), c(6,7)] <- c("Micaelamys_namaquensis", "not MDD_v1")
   tree$tip.label[which(tree$tip.label == "Tamiops_macclellandii")] <- "Tamiops_mcclellandii"  # typo
   act[which(act$MDD == "Tamiops_mcclellandii"), c(6,7)] <- c("Tamiops_mcclellandii", "typo")
}

# verify all names are matched
length(which(tree$tip.label %in% data$Phylo_name)) == 2424
length(which(data$Phylo_name %in% tree$tip.label)) == 2424


# 2. Example of data distribution on a phylogeny --------------------------



# 3. Simmap reconstruction with LTT plot for timeline ---------------------

####  REVELL CODE MODIFIED FROM FILES "9 - Addressing comments..."  #########

### making reproducible example for simmap LTT plot
packageVersion("phytools")
## [1] ‘1.5.1’

## I chose larger clades than in the blog example to see if individual (faint) LTT
## lines still appear or if thickness of the main line needs adjusting

# convert to named factor
AP <- factor(act$AP, levels = c("Nocturnal","Cathemeral","Diurnal"), exclude = "Crepuscular")
names(AP) <- act$MDD
AP <- AP[!is.na(AP)] # remove 24 crepuscular species
str(AP)



# correct tip names to match taxonomy - requires the data file loaded
CorTax <- function(tree){
   tree$tip.label[which(tree$tip.label == "Equus_africanus")]  <- "Equus_asinus"   # asinus is domestic E. africanus
   tree$tip.label[which(tree$tip.label == "Pygeretmus_zhitkovi")]  <- "Pygeretmus_shitkovi"    # correct spelling
   tree$tip.label[which(tree$tip.label == "Aethomys_namaquensis")] <- "Micaelamys_namaquensis" # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Tamiops_macclellandii")]<- "Tamiops_mcclellandii"   # typo
   tree$tip.label[which(tree$tip.label == "Nesotragus_moschatus")] <- "Neotragus_moschatus"    # typo
   tree$tip.label[which(tree$tip.label == "Pseudalopex_culpaeus")] <- "Lycalopex_culpaeus"     # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Pseudalopex_griseus")]  <- "Lycalopex_griseus"      # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Pseudalopex_gymnocercus")]  <- "Lycalopex_gymnocercus"  # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Pseudalopex_gymnocercus")]  <- "Lycalopex_sechurae" # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Pseudalopex_vetulus")]  <- "Lycalopex_vetulus"      # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Lutra_maculicollis")]   <- "Hydrictis_maculicollis" # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Monodelphis_unistriatus")]  <- "Monodelphis_unistriata" # Latin grammar correction
   tree$tip.label[which(tree$tip.label == "Dactylonax_palpator")]  <- "Dactylopsila_palpator"  # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Elephantulus_revoilii")]<- "Elephantulus_revoili"   # spelling error
   tree$tip.label[which(tree$tip.label == "Zaglossus_bruijnii")]   <- "Zaglossus_bruijni"      # spelling error
   tree$tip.label[which(tree$tip.label == "Procolobus_badius")]    <- "Piliocolobus_badius"    # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Procolobus_kirkii")]    <- "Piliocolobus_kirkii"    # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Procolobus_preussi")]   <- "Piliocolobus_preussi"   # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Procolobus_pennantii")] <- "Piliocolobus_pennantii" # genus transfer v1
   tree$tip.label[which(tree$tip.label == "Procolobus_rufomitratus")]  <- "Piliocolobus_rufomitratus"   # genus transfer v1
}

tvar <- CorTax(tvar)
tree <- drop.tip(tree, which(!tree$tip.label %in% data$Phylo_name))

# trim the phylogeny according to data
cld <- drop.tip(tvar, which(!tp$tip.label %in% names(AP)))

# specify transition model
qmat <- matrix(c(0, 1, 0, 2, 0, 3, 0, 4, 0), 3)

# run simmap
maps <- make.simmap(cld, AP, model=qmat, pi="fitzjohn", nsim=50) # this takes LONG
ltts <- ltt(maps)

options(scipen=10)

# plot
png("AP_ltt_tp5338.png",width=8,height=6,units="in",res=600)
par(mar=c(5.1,4.1,2.1,2.1))
cols <- setNames(c("dodgerblue2","#33CC33","gold"),levels(AP))
plot(ltts, colors=cols, show.total=FALSE, bty="n", cex.axis=0.8, las=1,
     xlab="millions of years (from the root)", axes=FALSE, log="y",
     ylim=c(1/length(maps),10000),alpha=0.1)
axis(1, at=round(seq(0,max(nodeHeights(cld)),length.out=9),1), cex.axis=0.8)
axis(2,las=1,cex.axis=0.8)
clip(0,max(nodeHeights(cld)),1/length(maps),1000)
grid()
dev.off()


####  END REVELL CODE  #########


# 4. Run models of HiSSE with root.type="herr_als" ------------------------



# 5. Run a few examples of SecSSE to compare results  ---------------------



# 6. CorHMM reconstruction as well (if necessary) -------------------------



# 7. Plot evolutionary transition models (search gmail for code) ----------



# 8. Phylogenetic distances -----------------------------------------------

## if looking at the results of the best model from each tree does not provide good
## enough a solution, phylogenetic distances can be compared based on Upham's backbone
## tree (omitting tipwards "patches" that might hold most of the topological variation)






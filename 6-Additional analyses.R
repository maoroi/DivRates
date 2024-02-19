## Additional analyses required by 2nd round of reviews NOT BASED ON COX 2021 DATASET

library(phytools)
library(stringr)
library(readr)
library(dplyr)

setwd("C:/Users/rma10kg/Projects/CH2/AddMuHisse/DivRates")



# 1. Ensure data files are up to date -------------------------------------

# function to correct tip names to match taxonomy and trim - requires the data file loaded
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

   tree <- drop.tip(tree, which(!tree$tip.label %in% act$Phylo_name))
}

## all variants that were used in the models
allvar <- c("0028","0612","0677","1030","1166","1774","1845","2966","3865","4496","5024","5221",
            "5333","5455","5684","5844","6259","6375","6404","6668","7198","7306","8981","9128",
            "0320","0369","0400","0647","0675","0772","0911","0917","1017","1233","1322","1350",
            "1564","1575","1666","1777","1803","1805","1816","1823","1904","2154","2316","2417",
            "2555","2953","2981","3019","3046","3120","3152","3168","3229","3349","3490","3599",
            "3734","3756","3766","3824","4012","4029","4224","4304","4387","4420","4423","4439",
            "4467","4568","4683","4831","4945","5039","5133","5242","5272","5338","5340","5568",
            "5594","5976","6255","6345","6747","6760","6774","6865","6914","6935","6947","6972",
            "7009","7051","7072","7074","7180","7444","7549","7627","7983","8004","8156","8307",
            "8377","8394","8429","8573","8705","8775","8974","9011","9023","9455","9572","9743",
            "9873","9923","9989","9997")


## loop over trees and reprocess to ensure the data and trees match perfectly
for(j in 1:length(allvar)){
   # start with a fresh copy of the data to ensure consistency
   act <- read.csv("ActivityData_MDD_v1.csv", header = TRUE)

   tree <- read.nexus(paste0("tvar",allvar[j],".nex"))
   tree <- ladderize(tree, right=TRUE)
   # trim tip labels to binomial only
   tree$tip.label <- str_extract(tree$tip.label,"[A-Z][a-z]*_[a-z]*")


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

   # correct tip labels in the tree, update the data table in the change
   {act[which(act$MDD == "Bubalus_bubalis"), c(6,7)] <- c("not in tree", "-")
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
      act[which(act$MDD == "Piliocolobus_badius"), c(6,7)] <- c("Piliocolobus_badius", "not MDD_v1")
      tree$tip.label[which(tree$tip.label ==  "Procolobus_kirkii")] <- "Piliocolobus_kirkii"      # genus transfer v1
      act[which(act$MDD == "Piliocolobus_kirkii"), c(6,7)] <- c("Piliocolobus_kirkii", "not MDD_v1")
      tree$tip.label[which(tree$tip.label == "Procolobus_pennantii")] <- "Piliocolobus_pennantii" # genus transfer v1
      act[which(act$MDD == "Piliocolobus_pennantii"), c(6,7)] <- c("Piliocolobus_pennantii", "not MDD_v1")
      tree$tip.label[which(tree$tip.label == "Procolobus_preussi")] <- "Piliocolobus_preussi"     # genus transfer v1
      act[which(act$MDD == "Piliocolobus_preussi"), c(6,7)] <- c("Piliocolobus_preussi", "not MDD_v1")
      tree$tip.label[which(tree$tip.label == "Procolobus_rufomitratus")] <- "Piliocolobus_rufomitratus"   # genus transfer v1
      act[which(act$MDD == "Piliocolobus_rufomitratus"), c(6,7)] <- c("Piliocolobus_rufomitratus", "not MDD_v1")
      act[which(act$MDD == "Loxodonta_cyclotis"), c(6,7)] <- c("not in tree", "-")
      tree$tip.label[which(tree$tip.label == "Pygeretmus_zhitkovi")] <- "Pygeretmus_shitkovi"     # correct spelling
      act[which(act$MDD == "Pygeretmus_shitkovi"), c(6,7)] <- c("Pygeretmus_shitkovi", "spelling")
      act[which(act$MDD == "Dephomys_eburneae"), c(6,7)] <- c("not in tree", "-")
      act[which(act$MDD == "Grammomys_poensis"), c(6,7)] <- c("not in tree", "-")
      tree$tip.label[which(tree$tip.label == "Aethomys_namaquensis")] <- "Micaelamys_namaquensis" # genus transfer v1
      act[which(act$MDD == "Micaelamys_namaquensis"), c(6,7)] <- c("Micaelamys_namaquensis", "not MDD_v1")
      tree$tip.label[which(tree$tip.label == "Tamiops_macclellandii")] <- "Tamiops_mcclellandii"  # typo
      act[which(act$MDD == "Tamiops_mcclellandii"), c(6,7)] <- c("Tamiops_mcclellandii", "typo")
   }

   tree <- CorTax(tree)

   # verify all names are matched
   if(length(which(tree$tip.label %in% act$Phylo_name)) == 2424){
      if(length(which(act$Phylo_name %in% tree$tip.label)) == 2424){
         write.nexus(tree, file=paste0("ntree",allvar[j],".nex"))
      }
   }
   # now that the tree is aligned with taxonomy, align the data to tree tips
   act <- act[which(act$Phylo_name %in% tree$tip.label),]
}


# 2. Example of data distribution on a phylogeny --------------------------


phy1 <- read.tree("UltMCC_2424spp_4AP_5kvars.nex")
tree <- ladderize(phy1)

label <- character(length(tree$tip.label))
names(label) <- act[which(act$Phylo_name %in% tree$tip.label),6]
label[act$AP=="Nocturnal"] <- "dodgerblue3"
label[act$AP=="Cathemeral"] <- "green3"
label[act$AP=="Diurnal"] <- "gold"
label[act$AP=="Crepuscular"] <- "magenta2"
label[which(label[] == '')] <- "white"

# ** change this to use original tip labels for taxonomic labeling
{artio <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Artiodactyla')])))
   chiro <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Chiroptera')])))
   carni <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Carnivora')])))
   prim <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Primates')])))
   roden <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Rodentia')])))
   eulip <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order %in% c('Erinaceomorpha','Soricomorpha'))])))
   xenar <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order %in% c('Cingulata','Pholidota','Pilosa'))])))
   afro <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order %in% c('Afrosoricida','Hyracoidea','Macroscelidea','Proboscidea','Sirenia','Tubulidentata'))])))
}
orders <- c(afro, xenar, eulip, chiro, roden, artio, prim, carni)
labels <- c("Afrotheria","Xenarthra","Eulipotyphla","Chiroptera","Rodentia","Artiodactyla","Primates","Carnivora")

pdf(file="DataMCC.pdf",width=40, height=40)
plot(tree, type = "fan", lwd = 0.4, label.offset = 0.4, cex = 0.3, show.tip.label = TRUE,
     tip.color = "white")#, edge.color = BemEDGE[,3])
ring(20, tree, style = "ring", offset = 5, col = label[tree$tip.label])
# adding arcs for orders
for(i in 1:4){#length(orders)) {
   arc.cladelabels(tree, labels[i], orders[i], cex=6,
                   MoreArgs = list(mark.node=FALSE, lab.offset=25, wing.length=0.5))
}
#ring(, tree, style = "ring", offset = -65.5, col = alpha("gray50", 0.25))
dev.off()

##OPTIONAL:
# colour main order brnaches
colours <- c("brown","cyan","green4","blue","darkorange2","goldenrod1","purple4","gray64")
for (k in 1:length(orders)) {
   AncNode <- orders[k]
   CladeSpp <- Descendants(tree, AncNode, type="all")
   edges <- BemEDGE[which(BemEDGE[,2] %in% CladeSpp),]
   BemEDGE[which(BemEDGE[,2] == AncNode),3] <- as.character(colours[k]) ## colouring the order root branch
   for (i in 1:length(edges[,2])) {
      edges[i,3] <- as.character(colours[k])
      BemEDGE[which(BemEDGE[,2] == edges[i,2]),3] <- edges[i,3]
   }
}

# 3. Simmap reconstruction with LTT plot for timeline ---------------------

####  REVELL CODE MODIFIED FROM FILES "9 - Addressing comments..."  #########

### making reproducible example for simmap LTT plot
packageVersion("phytools")
## [1] ‘1.5.1’

## I chose larger clades than in the blog example to see if individual (faint) LTT
## lines still appear or if thickness of the main line needs adjusting

# convert to named factor
AP <- factor(act$AP, levels = c("Nocturnal","Cathemeral","Diurnal"), exclude = "Crepuscular")
names(AP) <- act$Phylo_name
AP <- AP[!is.na(AP)] # remove 24 crepuscular species
str(AP)


# specify transition model
qmat <- matrix(c(0, 1, 0, 2, 0, 3, 0, 4, 0), 3)

for(j in 1:length(allvar)){
   tree <- read.nexus(paste0("SIMMAP/ntree",allvar[j],".nex"))
   tree <- drop.tip(tree, which(!tree$tip.label %in% names(AP)))

   # run SIMMAP and calculate LTT - these two lines take LONG
   maps <- make.simmap(tree, AP, model=qmat, pi="fitzjohn", nsim=100)
   ltts <- ltt(maps)

   saveRDS(maps, file=paste0("SIMMAP/simmap",allvar[j],".RDS"))
   saveRDS(ltts, file=paste0("SIMMAP/LTTs_",allvar[j],".RDS"))
}


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

library(hisse)
set.seed(88)

# start with a fresh copy of the data to ensure consistency
act <- read.csv("ActivityData_FEB24.csv", header = TRUE)
# remove crepusculars
length(which(act$AP == 'Crepuscular')) == 24
act3 <- act[-which(act$AP == 'Crepuscular'),]

# convert to numbers
act3$AP <- as.character(act3$AP)
act3$AP[act3$AP == 'Nocturnal'] <- 1
act3$AP[act3$AP == 'Cathemeral'] <- 2
act3$AP[act3$AP == 'Diurnal'] <- 3


## setup the unvarying part of the model
f <- read.csv("f.csv")
f <- f[,2]

# MuHiSSE model params
n <- 4
set_to_0 <- 4*(1:n)-1
TO_rates <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24)
drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)

# diversification rates for MuHiSSE
turnover <- TO_rates[1:(4*n)]
extinction.fraction <- rep(1,4*n)
extinction.fraction[set_to_0] <- 0

# transition matrix
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = TRUE)
trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)])
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA

## iterate over all phylogenetic trees
allvar <- c("0028","0612","0677","1030","1166","1774","1845","2966","3865","4496","5024","5221",
            "5333","5455","5684","5844","6259","6375","6404","6668","7198","7306","8981","9128",
            "0320","0369","0400","0647","0675","0772","0911","0917","1017","1233","1322","1350",
            "1564","1575","1666","1777","1803","1805","1816","1823","1904","2154","2316","2417",
            "2555","2953","2981","3019","3046","3120","3152","3168","3229","3349","3490","3599",
            "3734","3756","3766","3824","4012","4029","4224","4304","4387","4420","4423","4439",
            "4467","4568","4683","4831","4945","5039","5133","5242","5272","5338","5340","5568",
            "5594","5976","6255","6345","6747","6760","6774","6865","6914","6935","6947","6972",
            "7009","7051","7072","7074","7180","7444","7549","7627","7983","8004","8156","8307",
            "8377","8394","8429","8573","8705","8775","8974","9011","9023","9455","9572","9743",
            "9873","9923","9989","9997")

if(length(allvar) != 124) {allvar <- numeric()}

for(j in 1:length(allvar)){
   tree <- read.nexus(paste0("SIMMAP/ntree",allvar[j],".nex"))
   tree <- drop.tip(tree, which(!tree$tip.label %in% act3$Phylo_name))


   # transform to binary data
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

   # model
   model <- MuHiSSE(phy=tree, data=states.trans, f=f,
                    turnover=turnover,
                    eps=extinction.fraction,
                    hidden.states=TRUE,
                    trans.rate=trans.rate.mod,
                    root.type="herr_als",
                    turnover.upper=1000,
                    sann = TRUE)

   saveRDS(model, file=paste0("MuHiSSE",n,"_herr-als_tree", allvar[j],".RDS"))
}


# 5. Run models of HiSSE with ext_frac instead of turnover_rates ----------

## preparation steps are the same regardless of root.type

library(hisse)
set.seed(88)

# start with a fresh copy of the data to ensure consistency
act <- read.csv("ActivityData_FEB24.csv", header = TRUE)
# remove crepusculars
length(which(act$AP == 'Crepuscular')) == 24
act3 <- act[-which(act$AP == 'Crepuscular'),]

# convert to numbers
act3$AP <- as.character(act3$AP)
act3$AP[act3$AP == 'Nocturnal'] <- 1
act3$AP[act3$AP == 'Cathemeral'] <- 2
act3$AP[act3$AP == 'Diurnal'] <- 3

## sampled proportion by AP
f <- read.csv("f.csv")
f <- f[,2]


# MuHiSSE model params
n <- 4
set_to_0 <- 4*(1:n)-1
ext_frac <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24)
drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)

# diversification rates for MuHiSSE
turnover <- rep(1,4*n)
turnover[set_to_0] <- 0
extinction.fraction <- ext_frac[1:(4*n)]

# transition matrix
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = TRUE)
trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)])
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA


## 5.1 Estimating ext_frac while holding turnover rates constant ----------

for(j in 1:length(allvar)){
   tree <- read.nexus(paste0("SIMMAP/ntree",allvar[j],".nex"))
   tree <- drop.tip(tree, which(!tree$tip.label %in% act3$Phylo_name))


   # transform to binary data
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

   # model
   model <- MuHiSSE(phy=tree, data=states.trans, f=f,
                    turnover=turnover,
                    eps=extinction.fraction,
                    hidden.states=TRUE,
                    trans.rate=trans.rate.mod,
                    root.type="madfitz",
                    eps.upper=0.999,
                    sann = TRUE)

   saveRDS(model, file=paste0("MuHiSSE",n,"_madfitz_ExtFrac_tree", allvar[j],".RDS"))
}


## 5.2 Same as 5.1 but w root.type="herr_als" for comparison --------------

for(j in 1:length(allvar)){
   tree <- read.nexus(paste0("SIMMAP/ntree",allvar[j],".nex"))
   tree <- drop.tip(tree, which(!tree$tip.label %in% act3$Phylo_name))


   # transform to binary data
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

   # model
   model <- MuHiSSE(phy=tree, data=states.trans, f=f,
                    turnover=turnover,
                    eps=extinction.fraction,
                    hidden.states=TRUE,
                    trans.rate=trans.rate.mod,
                    root.type="herr_als",
                    eps.upper=0.999,
                    sann = TRUE)

   saveRDS(model, file=paste0("MuHiSSE",n,"_herrals_ExtFrac_tree", allvar[j],".RDS"))
}

# 6. CorHMM reconstruction as well (if necessary) -------------------------



# 7. Plot evolutionary transition models (search gmail for code) ----------



# 8. Phylogenetic distances -----------------------------------------------

## if looking at the results of the best model from each tree does not provide good
## enough a solution, phylogenetic distances can be compared based on Upham's backbone
## tree (omitting tipwards "patches" that might hold most of the topological variation)






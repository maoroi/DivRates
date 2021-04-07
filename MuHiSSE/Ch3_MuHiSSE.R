
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")
# here::here()

library(hisse)  
#> packageVersion('hisse')
#[1] ‘1.9.13’
library(diversitree)
#> packageVersion('diversitree')
#[1] ‘0.9.11’
require(phytools)
#> packageVersion('phytools')
#[1] ‘0.6.99’
library(phangorn)
#> packageVersion('phangorn')
#[1] ‘2.5.5’
require(scales)
#> packageVersion('scales')
#[1] ‘1.0.0’

## reproducible evironemnt using {renv}: The packages used in this project are recorded into a lockfile, renv.lock.
## use renv::restore() to reinstall all of the packages as declared in the lockfile.
set.seed(88)

# 1. Data preparation -----------------------------------------------------
# * 1.1 Functions to iron out phylogeny -----------------------------------
# the Upham 2019 trees aren't precisely ultrametric, so I used this function from Liam Revell: 
# http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html
force.ultrametric <- function(tree,method=c("nnls","extend")){
    method<-method[1]
    if(method=="nnls") tree <- nnls.tree(cophenetic(tree),tree, rooted=TRUE,trace=0)
    else if(method=="extend"){
        h<-diag(vcv(tree))
        d<-max(h)-h
        ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
                   y=tree$edge[,2])
        tree$edge.length[ii]<-tree$edge.length[ii]+d
    } else 
        cat("method not recognized: returning input tree\n\n")
    tree
}

# automatically correct tip names to match taxonomy - requires the data file loaded (based on section 1.4 below)
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
    tree <- drop.tip(tree, which(!tree$tip.label %in% act3$Phylo_name))
}

# * 1.2 Phylogenetic data (see also 1.4) ------------------------------------------
phy1 <- read.tree("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre")
tree <- ladderize(phy1)

# plot tree for inspection if needed
#pdf(file='tree.pdf', height=40, width=40)
#plot(phy1, type='fan', cex=0.15, label.offset = 0.05)
#dev.off()

## random sample of 250 trees out of 10k from Upham 2019
indices <- sprintf("%04d", sample(0:9999, 5000, replace=FALSE)) # sprintf makes all samples indices 4-digit long

for (i in 1:length(indices)){
    phyfile <- paste0("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree",indices[i],".tre")
    treevar <- read.tree(phyfile)
    treevar <- CorTax(treevar)
    #write.tree(treevar, file=paste0('treevar',indices[i],'.nex'))
    if (i == 1) {trees <- treevar} else {trees <- ape:::c.phylo(trees, treevar)} # group trees to one object
}
# MCC tree made of 5000 variants
clades <- prop.part(trees)
MCCtree <- maxCladeCred(trees, part = clades)


# * 1.3 Activity data  ----------------------------------------------------
read.csv("ActivityData_MDD_v1.csv") -> data   # see file "Uploading taxonomy.r" for how this file was made

# making sure names in data and in tree match
# checks below use the "%in%" syntax becasue setdiff() returns unique entries only, messing up the totals 
length(which(!data$MDD %in% tree$tip.label)) == 153     # checking names match - they don't
length(which(!data$MSW3 %in% tree$tip.label)) == 136    # somehow MSW3 has fewer mismatches


# * 1.4 Matching names in the data to phylogeny ---------------------------
# new columns for the names that are actually in the phylogeny, and where they're taken from
data$Phylo_name <- NA
data$Phylo_name_conforms_to <- NA

# * 1.4.1 Matching binomials in phylogeny to taxonomy ---------------------
# Species names were updated in the column 'Phylo_name'
for (i in 1:nrow(data)){ 
    if (as.character(data$MDD[i]) %in% tree$tip.label){ # if the MDD form is in the tree, it was retained
        data$Phylo_name[i] <- as.character(data$MDD[i])
        data$Phylo_name_conforms_to[i] <- "MDD_v1"
    } else if (as.character(data$MSW3[i]) %in% tree$tip.label){ # if MDD isn't but the MSW3 form is - MSW3 was retained
        data$Phylo_name[i] <- as.character(data$MSW3[i])
        data$Phylo_name_conforms_to[i] <- "MSW3"
    } else {                # if neither taxonomies match, manually verified form was inserted (see below)
        data$Phylo_name[i] <- "_"
        data$Phylo_name_conforms_to[i] <- "_"
    } 
}


# * 1.4.2 Correcting mismatches manually ----------------------------------
# (it's the safest solution)

# find the species for which both MSW3 and MDD are not in the tree
new_out <- which(!data$MDD %in% tree$tip.label)
old_out <- which(!data$MSW3 %in% tree$tip.label)

# list species for later checks and consistency
data$MDD[intersect(new_out, old_out)]
#[1] "Bubalus_bubalis"           "Neotragus_moschatus"       "Taurotragus_oryx"          "Lycalopex_culpaeus"       
#[5] "Lycalopex_griseus"         "Lycalopex_gymnocercus"     "Lycalopex_sechurae"        "Lycalopex_vetulus"        
#[9] "Galerella_pulverulenta"    "Galerella_sanguinea"       "Urva_edwardsii"            "Hydrictis_maculicollis"   
#[13] "Dermanura_azteca"          "Hypsugo_ariel"             "Marmosa_paraguayana"       "Marmosa_isthmica"         
#[17] "Marmosa_germana"           "Monodelphis_unistriata"    "Dactylopsila_palpator"     "Elephantulus_revoili"     
#[21] "Zaglossus_bruijni"         "Equus_asinus"              "Cercopithecus_denti"       "Cercopithecus_wolfi"      
#[25] "Piliocolobus_badius"       "Piliocolobus_kirkii"       "Piliocolobus_pennantii"    "Piliocolobus_preussi"     
#[29] "Piliocolobus_rufomitratus" "Paragalago_zanzibaricus"   "Loxodonta_cyclotis"        "Pygeretmus_shitkovi"      
#[33] "Dephomys_eburneae"         "Grammomys_poensis"         "Micaelamys_namaquensis"    "Tamiops_mcclellandii" 

# corrected names
{data[which(data$MDD == "Bubalus_bubalis"), c(6,7)] <- c("not in tree", "_")
    tree$tip.label[which(tree$tip.label == "Nesotragus_moschatus")] <- "Neotragus_moschatus"    # typo
    data[which(data$MDD == "Neotragus_moschatus"), c(6,7)] <- c("Neotragus_moschatus", "typo")   
    data[which(data$MDD == "Taurotragus_oryx"), c(6,7)] <- c("Tragelaphus_oryx", "MDD_v1.31")   # genus transfer in v1.31         
    tree$tip.label[which(tree$tip.label ==  "Pseudalopex_culpaeus")] <- "Lycalopex_culpaeus"    # genus transfer v1
    data[which(data$MDD == "Lycalopex_culpaeus"), c(6,7)] <- c("Lycalopex_culpaeus", "not MDD_v1") 
    tree$tip.label[which(tree$tip.label ==  "Pseudalopex_griseus")] <- "Lycalopex_griseus"      # genus transfer v1
    data[which(data$MDD == "Lycalopex_griseus"), c(6,7)] <- c("Lycalopex_griseus", "not MDD_v1")    
    tree$tip.label[which(tree$tip.label ==  "Pseudalopex_gymnocercus")] <- "Lycalopex_gymnocercus"  # genus transfer v1
    data[which(data$MDD == "Lycalopex_gymnocercus"), c(6,7)] <- c("Lycalopex_gymnocercus", "not MDD_v1")  
    tree$tip.label[which(tree$tip.label ==  "Pseudalopex_gymnocercus")] <- "Lycalopex_sechurae" # genus transfer v1
    data[which(data$MDD == "Lycalopex_sechurae"), c(6,7)] <- c("Lycalopex_sechurae", "not MDD_v1")   
    tree$tip.label[which(tree$tip.label ==  "Pseudalopex_vetulus")] <- "Lycalopex_vetulus"      # genus transfer v1
    data[which(data$MDD == "Lycalopex_vetulus"), c(6,7)] <- c("Lycalopex_vetulus", "not MDD_v1")   
    data[which(data$MDD == "Galerella_pulverulenta"), c(6,7)] <- c("Herpestes_pulverulentus", "MDD_v1.31")  # genus transfer in v1.31   
    data[which(data$MDD == "Galerella_sanguinea"), c(6,7)] <- c("Herpestes_sanguineus", "MDD_v1.31")# genus transfer in v1.31
    tree$tip.label[which(tree$tip.label ==  "Lutra_maculicollis")] <- "Hydrictis_maculicollis"      # genus transfer v1
    data[which(data$MDD == "Hydrictis_maculicollis"), c(6,7)] <- c("Hydrictis_maculicollis", "not MDD_v1")
    tree$tip.label[which(tree$tip.label ==  "Monodelphis_unistriatus")] <- "Monodelphis_unistriata" # Latin grammar correction
    data[which(data$MDD == "Monodelphis_unistriata"), c(6,7)] <- c("Monodelphis_unistriata", "not MDD_v1")   
    tree$tip.label[which(tree$tip.label ==  "Dactylonax_palpator")] <- "Dactylopsila_palpator"      # genus transfer v1
    data[which(data$MDD == "Dactylopsila_palpator"), c(6,7)] <- c("Dactylopsila_palpator", "not MDD_v1")   
    tree$tip.label[which(tree$tip.label == "Elephantulus_revoilii")] <- "Elephantulus_revoili"      # spelling error
    data[which(data$MDD == "Elephantulus_revoili"), c(6,7)] <- c("Elephantulus_revoili", "typo")   
    tree$tip.label[which(tree$tip.label == "Zaglossus_bruijnii")] <- "Zaglossus_bruijni"            # spelling error
    data[which(data$MDD == "Zaglossus_bruijni"), c(6,7)] <- c("Zaglossus_bruijni", "typo")     
    tree$tip.label[which(tree$tip.label == "Equus_africanus")] <- "Equus_asinus"    # asinus is domestic E. africanus
    data[which(data$MDD == "Equus_asinus"), c(6,7)] <- c("Equus_asinus", "domestic form")     
    data[which(data$MDD == "Cercopithecus_denti"), c(6,7)] <- c("not in tree", "split from C.pogonias")                
    data[which(data$MDD == "Cercopithecus_wolfi"), c(6,7)] <- c("not in tree", "split from C.mitis")
    tree$tip.label[which(tree$tip.label ==  "Procolobus_badius")] <- "Piliocolobus_badius"      # genus transfer v1
    data[which(data$MDD == "Piliocolobus_badius"), c(6,7)] <- c("Piliocolobus_badius", " not MDD_v1")   
    tree$tip.label[which(tree$tip.label ==  "Procolobus_kirkii")] <- "Piliocolobus_kirkii"      # genus transfer v1
    data[which(data$MDD == "Piliocolobus_kirkii"), c(6,7)] <- c("Piliocolobus_kirkii", "not MDD_v1")   
    tree$tip.label[which(tree$tip.label == "Procolobus_pennantii")] <- "Piliocolobus_pennantii" # genus transfer v1
    data[which(data$MDD == "Piliocolobus_pennantii"), c(6,7)] <- c("Piliocolobus_pennantii", "not MDD_v1")   
    tree$tip.label[which(tree$tip.label == "Procolobus_preussi")] <- "Piliocolobus_preussi"     # genus transfer v1
    data[which(data$MDD == "Piliocolobus_preussi"), c(6,7)] <- c("Piliocolobus_preussi", "not MDD_v1")   
    tree$tip.label[which(tree$tip.label == "Procolobus_rufomitratus")] <- "Piliocolobus_rufomitratus"   # genus transfer v1
    data[which(data$MDD == "Piliocolobus_rufomitratus"), c(6,7)] <- c("Piliocolobus_rufomitratus", "not MDD_v1")   
    data[which(data$MDD == "Loxodonta_cyclotis"), c(6,7)] <- c("not in tree", "_")
    tree$tip.label[which(tree$tip.label == "Pygeretmus_zhitkovi")] <- "Pygeretmus_shitkovi"     # correct spelling
    data[which(data$MDD == "Pygeretmus_shitkovi"), c(6,7)] <- c("Pygeretmus_shitkovi", "spelling")   
    data[which(data$MDD == "Dephomys_eburneae"), c(6,7)] <- c("not in tree", "_")
    data[which(data$MDD == "Grammomys_poensis"), c(6,7)] <- c("not in tree", "_")
    tree$tip.label[which(tree$tip.label == "Aethomys_namaquensis")] <- "Micaelamys_namaquensis" # genus transfer v1
    data[which(data$MDD == "Micaelamys_namaquensis"), c(6,7)] <- c("Micaelamys_namaquensis", "not MDD_v1")   
    tree$tip.label[which(tree$tip.label == "Tamiops_macclellandii")] <- "Tamiops_mcclellandii"  # typo
    data[which(data$MDD == "Tamiops_mcclellandii"), c(6,7)] <- c("Tamiops_mcclellandii", "typo")}

# verify all names are matched
length(which(tree$tip.label %in% data$Phylo_name)) == 2424
length(which(data$Phylo_name %in% tree$tip.label)) == 2424

# write to file just in case (this wont help to skip the matching step, because some tip names have to be 
# change in the phylogeny, which would be repeated for evert tree variant uploaded) 
#write.csv(data, file="ActivityData_MDD_v1_match.csv", row.names = FALSE)


# * 1.5 Final data processing steps ---------------------------------------
# read full data table (to save repeating the previous steps)
data <- read.csv(file="ActivityData_MDD_v1_match.csv")
# remove crepusculars
act <- data
length(which(act$AP == 'Crepuscular')) == 24
act3 <- act[-which(act$AP == 'Crepuscular'),]

# convert to numbers 
act3$AP <- as.character(act3$AP)
act3$AP[act3$AP == 'Nocturnal'] <- 1
act3$AP[act3$AP == 'Cathemeral'] <- 2
act3$AP[act3$AP == 'Diurnal'] <- 3

# keep only species that have activity AND phylogentic data
act3 <- act3[which(act3$Phylo_name %in% tree$tip.label),]
tree <- drop.tip(tree, which(!tree$tip.label %in% act3$Phylo_name))

# force ultrametric - placing this step after trimming the tree to save processing time
is.ultrametric(tree)
tree <- force.ultrametric(tree, "nnls")
is.ultrametric(tree)
#write.tree(tree, file='UltMCC.nex')


# * 1.6 Transforming to binary data (optional) ----------------------------
### THIS COULD ALSO BE DONE BY SETTING PARAMETER IDENTITIES IN THE MODELS
#  * 1.6.1 Any daytime activity (N <-> CD dichotomy) ----------------------


#  * 1.6.2 Only daytime activity (NC <-> D dichotomy ----------------------


#dat$C <- rep.int(0, nrow(dat))
dat$ND <- rep.int(0, nrow(dat))
dat$C <- rep.int(0, nrow(dat))
colnames(dat) <- c('Binomial','AP','ND','Generalist')
# transform to binary data ()
for (i in 1:nrow(dat)){
    if (dat[i,2]==1){dat[i,c(3:4)] <- c(0,0)} 
    if (dat[i,2]==2){dat[i,c(3:4)] <- c(0,1)}
    if (dat[i,2]==3){dat[i,c(3:4)] <- c(1,0)}
}



# 2. Analysis -------------------------------------------------------------
# * 2.1 Estimating sampling fraction  -------------------------------------
#  * 2.1.1 Proportion sampled out of all mammals --------------------------
## Proportion of sampled out of extant species by AP, assuming unbiased sampling. 

# total count in the entire dataset (inc species not in the phylogeny)
TotN <- length(which(data$AP == 'Nocturnal'))
TotC <- length(which(data$AP == 'Cathemeral'))
TotD <- length(which(data$AP == 'Diurnal'))
# proportions in data (sum up to 99%, crepuscular and others make the rest)
datN <- TotN / nrow(data)
datC <- TotC / nrow(data)
datD <- TotD / nrow(data)


#  * 2.1.2 Assumptions ----------------------------------------------------
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


# * 2.1.3 Estimate sampling fraction --------------------------------------
# species with no data or assumptions
all_extant <- length(unique(tax$SciName[which(tax$extinct. == 0)]))
unkn <- all_extant - (nrow(data) + Nunsamp + Dunsamp)                   # total number of extant species based on MDD v1
# expected species out of unknown if same proportion as in the data
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


# * 2.2 Setting up the analysis (following MuHiSSE tutorial) --------------

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
# As mentioned above, the number of free parameters in the model for both net turnover and extinction fraction
# are specified as index vectors provided to the function call. Each vector contains four entries that correspond
# to rates associated with the observed states (0 or 1) and the hidden states (A or B). They are always ordered
# as follows for a given hidden state, i: 00i, 01i, 10i, 11i. However, in this case we do not want any hidden
# states. But first let’s set up the “dull null” – i.e., turnover and extinction fraction are the same for all observed
# state combinations. Note the “f” represents the sampling fraction for each observed state combination, which
# is a vector ordered in the same manner as for turnover and extinction fraction vectors:
f = c(freq[1], freq[2], 0, freq[3])

turnover <- c(1,1,0,1)
extinction.fraction <- c(1,1,0,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))
#Now, we can call MuHiSSE and estimate the parameters under this model using the default settings:
dull.null <- MuHiSSE(phy=tree, data=states.trans, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=FALSE,
                     trans.rate=trans.rate.mod)


# CID-8  (differential transitions among hidden states)
# *automating model specifications as function of hidden states*
n <- 8                  
set_to_0 <- 4*(1:n)-1
if (n < 9) {   # For models w >8 hidden states, transition matrix has to be coded manually to bypass TransMatMakerMuHiSSE() limitations
    drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)
} else {drops <- NA}

# div rates formulae
for (z in 1:n) {reps <- rep(z, 4)
    if (z == 1) {turnover <- reps} else {turnover <- c(turnover, reps)}}
extinction.fraction <- rep(1,4*n) 
turnover[set_to_0] <- extinction.fraction[set_to_0] <- 0
# trans rates formulae
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = TRUE)
trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)])   
trans.rate.mod[,set_to_0] = trans.rate.mod[set_to_0,] <- 0

# optional: identical transition rates across hidden states
EqTransAcrossHiddenStates <- function(rateMatrix, n_HiddenStates) {
    lows <- 4*(1:(n_HiddenStates-1))+1
    highs <- 4*(2:n_HiddenStates)
    for (z in 1:(n_HiddenStates-1)) {
        l <- lows[z]
        h <- highs[z]
        rateMatrix[l:h,l:h] <- rateMatrix[1:4,1:4]}
    rateMatrix <- ParDrop(rateMatrix, c(0)) # update parameter indices (the indices are just names so this makes no difference to the model, just easier to read)
}
trans.rate.modmod <- EqTransAcrossHiddenStates(trans.rate.mod, n) # identical transition rates across hidden states

# model 
MuHiSSE <- MuHiSSE(phy=tree, data=states.trans, f=f, turnover=turnover, eps=extinction.fraction, 
                   hidden.states=TRUE, trans.rate=trans.rate.mod)


# ------------------------------------------------------------------ #
# enough demonstrations, below is the part that was actually done:   #
# ------------------------------------------------------------------ #

# -load tree (from section 1.2 above)
# -load data (from section 1.3 above)
# -calculate sampling fraction (section 2.1 above)

states.trans <- states
for(i in 1:Ntip(phy)){
    if(states[i,1] == 1){
        states.trans[i,1] = 0
        states.trans[i,2] = 0
    }
    if(states[i,1] == 2){
        states.trans[i,1] = 0
        states.trans[i,2] = 1
    }
    if(states[i,1] == 3){
        states.trans[i,1] = 1
        states.trans[i,2] = 1
    }
}

## Setting up a three-state model using MuHiSSE. This may be beneficial with rather large trees.
## "we are going to assume we do not have a fourth state"
turnover <- c(1,2,3,0)
extinction.fraction <- c(1,1,1,0)
f = freq

## generate transition rate matrix (unordered trait model):
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
print(trans.rate)
## Oredered model: use the ParDrop() to remove transitions to and from the fourth state in the model:
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))
print(trans.rate.mod)
## From there the function call is all the same:
MuHiSSE <- MuHiSSE(phy=tree, data=states, f=f, turnover=turnover, eps=extinction.fraction, 
                   hidden.states=FALSE, trans.rate=trans.rate.mod)

## If you wanted to set up a three-state model, but include hidden states, you would do the following:
turnover <- c(1,2,3,0, 4,5,6,0)
extinction.fraction <- c(1,1,1,0, 1,1,1,0)
f = c(0.6742257, 0.3854165, 0.5727366, 0.0000000)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1)
trans.rate.mod <- ParDrop(trans.rate, c(4,6,7,8,12,14,15,16))
## run analysis
MuHiSSE1 <- MuHiSSE(phy=datatree, data=states, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=TRUE,
                    trans.rate=trans.rate.mod)

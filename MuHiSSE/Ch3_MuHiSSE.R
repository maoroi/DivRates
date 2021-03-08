
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
require(scales)
#> packageVersion('scales')
#[1] ‘1.0.0’

## reproducible evironemnt using {renv}: The packages used in this project are recorded into a lockfile, renv.lock.
## use renv::restore() to reinstall all of the packages as declared in the lockfile.
set.seed(88)

# 1. Data preparation -----------------------------------------------------
# * 1.1 Function to iron out phylogeny ------------------------------------
# the Upham 2019 trees aren't precisely ultrametric, so I used code from Liam Revell: 
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

# * 1.2 Phylogenetic data (see also 1.4) ------------------------------------------
phy1 <- read.tree("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre")
phy1 <- ladderize(phy1)

# force ultrametric
is.ultrametric(phy1)
tree <- force.ultrametric(phy1, "nnls")
is.ultrametric(phy1)
#write.tree(phy1, file='tree????.nex')


# plot tree for inspection if needed
#pdf(file='tree.pdf', height=40, width=40)
#plot(phy1, type='fan', cex=0.15, label.offset = 0.05)
#dev.off()

## random sample of 250 trees out of 10k from Upham 2019
indices <- sprintf("%04d", sample(0:9999, 250, replace=FALSE)) # sprintf makes all indices 4-digit long

#clades <- prop.part(tree)
#MCCtree <- maxCladeCred(tree, part = clades)


# * 1.2 Activity data  ----------------------------------------------------
read.csv("ActivityData_MDD_v1.csv") -> data   # see file "Uploading taxonomy.r" for how this file was made

# making sure names in data and in tree match
# checks below use the "%in%" syntax becasue setdiff() returns unique entries only, messing up the totals 
length(which(!data$MDD %in% tree$tip.label)) == 153     # checking names match - they don't
length(which(!data$MSW3 %in% tree$tip.label)) == 136    # somehow MSW3 has fewer mismatches


# * 1.3 Matching names in the data to phylogeny ---------------------------
# new columns for the names that are actually in the phylogeny, and where they're taken from
data$Phylo_name <- NA
data$Phylo_name_note <- NA

# * 1.3.1 Matching binomials in phylogeny to taxonomy ---------------------
# Species names were updated in the column 'Phylo_name'
for (i in 1:nrow(data)){ 
    if (as.character(data$MDD[i]) %in% tree$tip.label){ # if the MDD form is in the tree, it was retained
        data$Phylo_name[i] <- as.character(data$MDD[i])
        data$Phylo_name_note[i] <- "MDD_v1"
    } else if (as.character(data$MSW3[i]) %in% tree$tip.label){ # if MDD isn't but the MSW3 form is - MSW3 was retained
        data$Phylo_name[i] <- as.character(data$MSW3[i])
        data$Phylo_name_note[i] <- "MSW3"
    } else {                # if neither taxonomies match, manually verified form was inserted (see below)
        data$Phylo_name[i] <- "_"
        data$Phylo_name_note[i] <- "_"
    } 
}


# * 1.3.2 Correcting mismatches manually ----------------------------------
# it's the safest solution

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
write.csv(data, file="ActivityData_MDD_v1_match.csv", row.names = FALSE)


# * 1.4 Final data processing steps ---------------------------------------
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


# 2. Analysis -------------------------------------------------------------


# * 2.1 Estimating sampling fraction  -------------------------------------

#  * 2.1.1 Proportion sampled out of all mammals --------------------------
## Proportion of sampled out of extant species by AP, assuming unbiased sampling. 

# total count in the entire dataset (inc species not in Faurby 2015 phylogeny)
TotN <- length(which(data$AP == 'Nocturnal'))
TotC <- length(which(data$AP == 'Cathemeral'))
TotD <- length(which(data$AP == 'Diurnal'))
# proportions in data
datN <- TotN / nrow(data)
datC <- TotC / nrow(data)
datD <- TotD / nrow(data)


# * 2.1.2 Assumptions -----------------------------------------------------
# unsampled nocturnal bats (n = 1096)
bats_unsamp <- length(which(tax$Order == 'CHIROPTERA')) - length(which(data$Order == 'Chiroptera')) - 6  # unsampled bats minus 6 extinct spp

# unsampled squirrels (27 noct flying squirrels; 87 diur sciurids)
SQRLS <- tax[which(tax$Family == 'SCIURIDAE'),]                         # all sciurids (all extant)
SQRLS_DATA <- data[which(data$Family == 'Sciuridae'),]                  # sampled sciurids 
flying <- length(which(SQRLS$Tribe == 'PTEROMYINI'))                    # 57 flying squirrels in MDD
noct_sqrls <- length(which(SQRLS_DATA$AP == "Nocturnal"))               # flying squirels are the only nocturnal sciurids
FLY_SQRLS_unsamp <- flying - noct_sqrls
Diur_sqrls_unsamp <- nrow(SQRLS) - nrow(SQRLS_DATA) - FLY_SQRLS_unsamp  # unsampled scuirids minus unsampled noct sciurids

# diurnal Simiiformes (n = 148, excl. Aotus) 
PRIMS <- tax[which(tax$Order == "PRIMATES"),]
PRIMS <- PRIMS[which(PRIMS$extinct. == 0),]                             # only extant species
simians <- c("Atelidae", "Cebidae", "Cercopithecidae", "Pitheciidae", "Hominidae", "Hylobatidae")
simian_count <- length(which(PRIMS$Family %in% toupper(simians)))       # all simians recognised in MDD
SIMIAN_DATA <- data[which(data$Family %in% c(simians, "Aotidae")),]     # Aotidae demoted into Cebidae (Aotinae) in MDD but data$Family relates to MSW3   
length(which(tax$Subfamily == "AOTINAE")) == 11                         # 11 species in MDD (no assumption of AP of unknown spp)
length(which(data$Family == 'Aotidae')) == 7                            # night monkeys in data 
simi_unsamp <- simian_count - nrow(SIMIAN_DATA) - (11-7)                # unsmapled (diur) simians minus unsampled (noct) Aotinae

# total unsampled (N=1123, D=235)
Nunsamp <- bats_unsamp + FLY_SQRLS_unsamp
Dunsamp <- Diur_sqrls_unsamp + simi_unsamp
rm(bats_unsamp, FLY_SQRLS_unsamp, Diur_sqrls_unsamp, simi_unsamp, PRIMS, simians, 
   simian_count, SIMIAN_DATA, noct_sqrls, flying, SQRLS_DATA, SQRLS)
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

# check that proportions were retained so that nocturnal, cathemeral, and diurnal species account for ~99% of mammals
AP_data <- length(which(data$AP %in% c("Nocturnal", "Cathemeral", "Diurnal"))) / nrow(data)
AP_extrapol <- ((TotN + Nunsamp + expN) + (TotC + expC) + (TotD + Dunsamp + expD)) / all_extant
round(AP_data, 2) == round(AP_extrapol, 2)


#states <- as.numeric(act3$AP)
#names(states) <- act3$Phylo_name

# * 2.2 Setting up the analysis following steps from MuHiSSE tutor --------

for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
}

states <- data.frame(tree$tip.state, tree$tip.state, row.names=names(tree$tip.state))
states <- states[tree$tip.label,]
#
#
#
#
#
#
# * 999 Sampling fraction -----------------------------------------------------




## simulate parameters for three-state model
pars <- c(.1, .15, .2, # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0, # q12, q13
          .05, .05, # q21, q23
          0, .05) # q31, q32

phy <- tree.musse(pars, 30, x0=1)
states <- phy$tip.state
lik <- make.musse(phy, states, 3)
lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1, mu2 ~ mu1, mu3 ~ mu1,
                      q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0.03, q32 ~ q12)

##  Here is the likelihood for a maximally constrained model using diversitree:
diversitree.constrained = lik.base(c(.1, .03, .05))
print(diversitree.constrained)

## Now, let's compare the likelihood using a three-state model in GeoHiSSE:
states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
states <- states[phy$tip.label,]
states[states[,1]==3,] = 4
pars.hisse <- c(0.1, 0.1, 0.1, 0.03, 0.03, 0.03, 0.05, 0, 0.05, 0.05, 0, 0.05)
model.vec = rep(0,120)
model.vec[1:12] = pars.hisse
phy$node.label = NULL
cache <- hisse:::ParametersToPassMuSSE(phy, states[,1], model.vec, f=c(1,1,1), hidden.states="TEST1")
geosse.constrained <- hisse:::DownPassMusse(phy, cache, hidden.states=FALSE, root.type="madfitz", 
                                            condition.on.survival=TRUE, root.p=NULL)
comparison <- identical(round(geosse.constrained,4), round(diversitree.constrained,4))
print(comparison)

## We can show that you can obtain the same likelihood using MuHiSSE. Warning: as before it's a bit convoluted to 
## format the data properly, but it can be done:
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
    if(states[i,1] == 4){
        states.trans[i,1] = 1
        states.trans[i,2] = 0
    }
}
pars.hisse <- c(0.1+0.03,0.1+0.03,0.1+0.03,0,
                0.03/0.1,0.03/0.1,0.03/0.1,0,
                0.05,0,0, 0.05,0.05,0, 0.03,0.05,0, 0,0,0)
model.vec = rep(0,384)
model.vec[1:20] = pars.hisse
phy$node.label = NULL
cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=TRUE, nb.tip=Ntip(phy), 
                                         nb.node=Nnode(phy), bad.likelihood=exp(-500), ode.eps=0)
gen <- hisse:::FindGenerations(phy)
dat.tab <- hisse:::OrganizeData(states.trans, phy, f=c(1,1,1,0), hidden.states=TRUE)
muhisse.constrained <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache, root.type="madfitz", 
                                               condition.on.survival=TRUE, root.p=NULL)
comparison <- identical(round(muhisse.constrained,4), round(diversitree.constrained,4))
print(comparison)

# ------------------------------------------------------------------ #
# enough demonstrations, below is the part that was actually done:   #
# ------------------------------------------------------------------ #

## load data (taken as is from "MuSSE MCC tree ALL.R")
newdata <- read.csv('EvolRatesData.csv')
newdata <- newdata[,c(1,2,5,6)]
newdata[,3] <- gsub(' ', '_', newdata[,3])
# cleaning the data
newdata[which(newdata[,3] == 'Acomys_russatus'),4] <- 'Cathemeral'
newdata <- newdata[-which(newdata[,3] == 'Acomys_russatus')[1],]    ## remove one of two indentical Acomys_russtaus entries
newdata[which(newdata[,3] == 'Stenella_longirostris'),4] <- 'Cathemeral'
newdata <- newdata[-which(newdata[,3] == 'Stenella_longirostris')[1],]
# remove species with multiple entries
dups <- as.character(newdata[which(duplicated(newdata[,3]) == TRUE),3])
RTdata <- newdata[-which(newdata[,3] %in% dups),]
RTdata[,4] <- gsub('Nocturnal', '1', RTdata[,4])
RTdata[,4] <- gsub('Cathemeral', '2', RTdata[,4])
RTdata[,4] <- gsub('Diurnal', '3', RTdata[,4])
# unify activity descriptions
aberr <- unique(RTdata[,4])
RTdata[which(RTdata[,4] %in% aberr[c(4,12,13)]),4] <- 1
RTdata[which(RTdata[,4] %in% aberr[c(7,8,10,11)]),4] <- 2
RTdata[which(RTdata[,4] %in% aberr[5]),4] <- 3
# remove species with other activity descriptions
rawdat <- RTdata[which(RTdata[,4] %in% c(1,2,3)),]   ## remove 24 crepuscular species, two 'diurnal or cathemeral'
# make data input table
dat <- data.frame(rawdat[,c(3,4)])

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

# load tree
tree <- read.tree("MCCtree.tre")
nodata <- which(!tree$tip.label %in% dat$Binomial)
datatree <- drop.tip(tree,nodata)

## calculate expected state frequencies 
# proportion of sampled species corrected for KNOWN biases (836 bats 800 rodents presumed nocturnal, 71 primates and 
# 93 squirrels presumed diurnal, and assuming unbiased sampling for all other species
## the below section is written this way in case any assumptions about the AP of missing species are needed
TotC <- length(which(dat$AP == 2))
TotN <- length(which(dat$AP == 1))
TotD <- length(which(dat$AP == 3))
# species with no data nor assumptions
unkn <- 5416-(nrow(dat)+1636+164)   
# proportions in data
datC <- TotC / nrow(dat)
datN <- TotN / nrow(dat)
datD <- TotD / nrow(dat)
# expected species out of unknown
expC <- datC * unkn
expN <- datN * unkn
expD <- datD * unkn
# proportion sampled out of likely total
Cprop <- TotC/(TotC + expC)
Nprop <- TotN/(TotN + 1636 + expN)
Dprop <- TotD/(TotD + 164 + expD)
# expected proportions
freq = c(Cprop,Nprop,Dprop,0)
# remove from data species not in the tree (done after sampling proportion estimation)
dat <- dat[which(dat$Binomial %in% datatree$tip.label),]
states <- data.frame(dat$Binomial, dat[,c(3:4)], row.names=dat$Binomial)
#states <- states[which(states$dat.Binomial %in% datatree$tip.label),]

## Setting up a three-state model using MuHiSSE. This may be beneficial with rather large trees.
## "we are going to assume we do not have a fourth state"
turnover <- c(1,2,3,0)
extinction.fraction <- c(1,1,1,0)
f = freq

## generate, and then modify the transition rate matrix (unordered trait model):
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
print(trans.rate)
## Oredered model: use the ParDrop() to remove transitions to and from the fourth state in the model:
trans.rate.mod <- ParDrop(trans.rate, c(4,6,7,8))
print(trans.rate.mod)

## From there the function call is all the same:
MuHiSSE <- MuHiSSE(phy=datatree, data=states, f=f, turnover=turnover, eps=extinction.fraction, 
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

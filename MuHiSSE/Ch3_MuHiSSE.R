
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
set.seed(98)

# 1. Data preparation -----------------------------------------------------

## random tree sample from Upham et al. 2019
indices <- sprintf("%04d", sample(0:9999, 250, replace=FALSE)) # sprintf makes all indices 4-digit long

# * 1.1 Taxonomy ----------------------------------------------------------
read.csv("MamTaxMDD_v1.31_6513sp.csv") -> TAX   # taxonomy from Mammal Diversity Database [https://www.mammaldiversity.org/], last accessed 23/02/21
TAX[,c(1,8,13,16,17,21,26,27,35:39)] -> tax      # relevant columns only (also, some species are wrongly marked EXTINCT, e.g. Alcelaphus buselaphus)


# * 1.2 Phylogeny (see also 1.4) ------------------------------------------
Btree <- read.tree("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre")
tree <- ladderize(Btree)

# plot tree for inspection if needed
pdf(file='tree.pdf', height=40, width=40)
plot(tree, type='fan', cex=0.15, label.offset = 0.05)
dev.off()

#clades <- prop.part(tree)
#MCCtree <- maxCladeCred(tree, part = clades)


# * 1.3 Activity data -----------------------------------------------------
read.csv('Appendix_1-ActivityData.csv') -> DATA
DATA[,c(1:4)] -> DATA
DATA[,3] <- gsub(' ', '_', DATA[,3], ignore.case = FALSE)
DATA <- cbind(DATA[,c(1:3)], DATA[,c(3:4)])
colnames(DATA) <- c(colnames(DATA[1:3]), "MSW3", colnames(DATA[ncol(DATA)]))
DATA$Binomial <- NA
# P. magdalenae and P. poliopus dropped from Burgin 2018 without clue why/where, so I omitted from data to not trip the "manual" loop 
data <- DATA

# retain unchanged binomials
unchanged <- which(data$MSW3 %in% tax$sciName)
data$Binomial[unchanged] <- data$MSW3[unchanged]
# consistency test
length(which(is.na(data$Binomial))) == 301


# * 1.3.1 Converting taxonomy to Burgin et al. 2018 -----------------------
## There are several types of changes: genus transfer, split/elevation to full species, spelling errors, historical errors

# First, match species based on MDD authors specifications
for (i in 1:nrow(data)){
    if (is.na(data$Binomial[i])){
        try(data$Binomial[i] <- as.character(tax$sciName[which(tax$MSW3_sciName == data$MSW3[i])]), silent = TRUE)
    }
}
# consistency test
length(which(is.na(data$Binomial))) == 45


# this was done manually because automated text searches give false hits. 
{data$Binomial[which(data$MSW3 == "Capricornis_milneedwardsii")] <- "Capricornis_sumatraensis"
    data$Binomial[which(data$MSW3 == "Prionailurus_iriomotensis")] <- "Prionailurus bengalensis" # searched IUCN Red List data
    data$Binomial[which(data$MSW3 == "Bassaricyon_beddardi")] <- "Bassaricyon_alleni"              
    data$Binomial[which(data$MSW3 == "Pteropus_giganteus")] <- "Pteropus_medius" 
    data$Binomial[which(data$MSW3 == "Hypsugo_bodenheimeri")] <- "Hypsugo_ariel"
    data$Binomial[which(data$MSW3 == "Scotophilus_borbonicus")] <- "Scotophilus_trujilloi"
    data$Binomial[which(data$MSW3 == "Monodelphis_sorex")] <- "Monodelphis_dimidiata"
    data$Binomial[which(data$MSW3 == "Eulemur_albocollaris")] <- "Eulemur_cinereiceps" 
    data$Binomial[which(data$MSW3 == "Cratogeomys_neglectus")] <- "Cratogeomys_fumosus"
    data$Binomial[which(data$MSW3 == "Monodelphis_maraxina")] <- "Monodelphis_glirina"
    data$Binomial[which(data$MSW3 == "Monodelphis_rubida")] <- "Monodelphis_americana" 
    data$Binomial[which(data$MSW3 == "Akodon_molinae")] <- "Akodon_dolores"
    data$Binomial[which(data$MSW3 == "Orthogeomys_thaeleri")] <- "Heterogeomys_dariensis"
    data$Binomial[which(data$MSW3 == "Microtus_bavaricus")] <- "Microtus_liechtensteini"
    data$Binomial[which(data$MSW3 == "Orthogeomys_cuniculus")] <- "Orthogeomys_grandis"
    data$Binomial[which(data$MSW3 == "Sorex_bairdi")] <- "Sorex_bairdii"
    data$Binomial[which(data$MSW3 == "Proechimys_magdalenae")] <- "Proechimys_chrysaeolus"
    data$Binomial[which(data$MSW3 == "Nyctophilus_timoriensis")] <- "Nyctophilus_major" # taxon split
    data$Binomial[which(data$MSW3 == "Micoureus_regina")] <- "Marmosa_isthmica" # genus transfer and split
}   
    
## insert added species (splits) in Nyctophilus and Marmosa (the last two lines above)
incumbent <- which(data$Binomial == "Nyctophilus_major")
template <- data[incumbent,]
template$Binomial <- "Nyctophilus_corbeni" # taxon split
insert <- template
template$Binomial <- "Nyctophilus_sherrini" # taxon split
insert <- rbind(insert, template)
template$Binomial <- "Nyctophilus_shirleyae" # taxon split
insert <- rbind(insert, template)

# merge into data
data <- rbind(data[1:incumbent,], insert, data[(incumbent+1):nrow(data),])

## and again for Marmosa
incumbent <- which(data$Binomial == "Marmosa_isthmica")
template <- data[incumbent,]
template$Binomial <- "Marmosa_germana" # genus transfer and split

# merge into data
data <- rbind(data[1:incumbent,], template, data[(incumbent+1):nrow(data),])

## removing redundancies after name update (merging results in multiple entries under the same new name)
# start by removing all entries not updated to burgin 2018 - these are duplicates
data <- data[-which(is.na(data$Binomial)),]
# the remaining are two species that merge different APs and will be made cathemeral
dups <- data$Binomial[which(duplicated(data$Binomial))] # find all duplicated names 
whittle <- data[which(data$Binomial %in% dups),]
# first remove excess records
rm1 <- which(data$MSW3 == "Monodelphis_rubida")
rm2 <- which(data$MSW3 == "Orthogeomys_thaeleri")
data <- data[-c(rm1, rm2),]
# then change the AP of the correspondigng record to cathemeral
rm1 <- which(data$Binomial == "Monodelphis_americana")
rm2 <- which(data$Binomial == "Heterogeomys_dariensis")
data$Diel.activity.pattern[c(rm1,rm2)] <- "Cathemeral"
# save 
write.csv(data, file="ActivityData_2018Tax.csv", row.names = FALSE)









## no need actually update these species because they were all merged into existing entries w same AP
{#data$Binomial[which(data$MSW3 == "Alcelaphus_lichtensteinii")] <- "Alcelaphus_buselaphus"
   #data$Binomial[which(data$MSW3 == "Pseudois_schaeferi")] <- "Pseudois_nayaur"          
   #data$Binomial[which(data$MSW3 == "Galidictis_grandidieri")] <- "Galidictis_fasciata"
   #data$Binomial[which(data$MSW3 == "Conepatus_humboldtii")] <- "Conepatus_chinga"
   #data$Binomial[which(data$MSW3 == "Bassaricyon_lasius")] <- "Bassaricyon_gabbii"          
   #data$Binomial[which(data$MSW3 == "Bassaricyon_pauli")] <- "Bassaricyon_gabbii"
   #data$Binomial[which(data$MSW3 == "Chaetophractus_nationi")] <- "Chaetophractus_vellerosus"
   #data$Binomial[which(data$MSW3 == "Sminthopsis_aitkeni")] <- "Sminthopsis_fuliginosus"
   #data$Binomial[which(data$MSW3 == "Sminthopsis_griseoventer")] <- "Sminthopsis_fuliginosus"
   #data$Binomial[which(data$MSW3 == "Marmosops_cracens")] <- "Marmosops_fuscatus"
   #data$Binomial[which(data$MSW3 == "Marmosops_dorothea")] <- "Marmosops_noctivagus"
   #data$Binomial[which(data$MSW3 == "Monodelphis_theresa")] <- "Monodelphis_scalops"   # searched IUCN Red List data
   #data$Binomial[which(data$MSW3 == "Trichosurus_arnhemensis")] <- "Trichosurus_vulpecula" 
   #data$Binomial[which(data$MSW3 == "Aotus_hershkovitzi")] <- "Aotus_lemurinus"
   #data$Binomial[which(data$MSW3 == "Alouatta_coibensis")] <- "Alouatta_palliata"
   #data$Binomial[which(data$MSW3 == "Cercopithecus_albogularis")] <- "Cercopithecus_mitis"
   #data$Binomial[which(data$MSW3 == "Callicebus_dubius")] <- "Plecturocebus_caligatus"
   #data$Binomial[which(data$MSW3 == "Lagidium_peruanum")] <- "Lagidium_viscacia"           # searched IUCN Red List data
   #data$Binomial[which(data$MSW3 == "Microtus_breweri")] <- "Microtus_pennsylvanicus"
   #data$Binomial[which(data$MSW3 == "Neotoma_martinensis")] <- "Neotoma_bryanti"
   #data$Binomial[which(data$MSW3 == "Oligoryzomys_delticola")] <- "Oligoryzomys_nigripes"
   #data$Binomial[which(data$MSW3 == "Dasyprocta_cristata")] <- "Dasyprocta_leporina"
   #data$Binomial[which(data$MSW3 == "Proechimys_poliopus")] <- "Proechimys_guairae"
   #data$Binomial[which(data$MSW3 == "Sphiggurus_villosus")] <- "Coendou_spinosus" # genus transfer
   #data$Binomial[which(data$MSW3 == "Pseudomys_pilligaensis")] <- "Pseudomys_delicatulus"
   #data$Binomial[which(data$MSW3 == "Petinomys_sagitta")] <- "Hylopetes_sagitta"           # based on notes in Burgin 2018
} 
## THE BELOW REFERS TO WORKING WITH THE FILE "asm-species-2020-02-04.csv" BUT NOT WITH THE CURRENT "MDD" FILE
{# match binomials after genus transfer - this worked 
#transfer <- which(data$MSW3[which(is.na(data$Binomial))] %in% tax$IfTransfer_oldSciName) # this wrongly includes "Leptonycteris_yerbabuenae" which is in MSW3
#for (i in 1:length(transfer)){
#    data$Binomial[transfer[i]] <- as.character(tax$SciName[which(tax$IfTransfer_oldSciName == data$MSW3[transfer[i]])])
#}
    
# For some transfers the previous binomial is never mentioned in the new taxonomy table (e.g. Bison_bison). 
# SOLUTION: finding matches for sp. name in its family (that's what the loop below does). For families w multiple 
# matches I checked taxonomy and IUCN website and corrected manually. Manual corrections are placed before the 
# loop so they don't trip it    
   
#data$Binomial[which(data$MSW3 == "Vulpes_rueppellii")] <- "Vulpes_rueppelli" 
   #data$Binomial[which(data$MSW3 == "Leopardus_colocolo")] <- "Leopardus_colocola"
   #data$Binomial[which(data$MSW3 == "Herpestes_brachyurus")] <- "Urva_brachyura"
   #data$Binomial[which(data$MSW3 == "Herpestes_javanicus")] <- "Urva_javanica"
   #data$Binomial[which(data$MSW3 == "Herpestes_semitorquatus")] <- "Urva_semitorquata"
   #data$Binomial[which(data$MSW3 == "Herpestes_vitticollis")] <- "Urva_vitticolla"
   #data$Binomial[which(data$MSW3 == "Aonyx_cinerea")] <- "Aonyx_cinereus"
   #data$Binomial[which(data$MSW3 == "Triaenops_rufus")] <- "Triaenops_menamena"
   #data$Binomial[which(data$MSW3 == "Sminthopsis_fuliginosus")] <- "Sminthopsis_griseoventer"
   #data$Binomial[which(data$MSW3 == "Lepus_microtis")] <- "Lepus_victoriae"                # searched IUCN Red List data
   #data$Binomial[which(data$MSW3 == "Heterocephalus_glaber")] <- "Heterocephalus_glaber"   # wrongly named glader in Burgin 2018
   #data$Binomial[which(data$MSW3 == "Abrothrix_olivaceus")] <- "Abrothrix_olivacea"
   #data$Binomial[which(data$MSW3 == "Phyllomys_blainvillii")] <- "Phyllomys_blainvillii"   # missing one 'l' in Burgin 2018  
   #data$Binomial[which(data$MSW3 == "Hylopetes_lepidus")] <- "Hylopetes_sagitta"           # based on notes in Burgin 2018
   #data$Binomial[which(data$MSW3 == "Callicebus_dubius")] <- "Plecturocebus_caligatus"
   #data$Binomial[which(data$MSW3 == "Vampyressa_bidens")] <- "Vampyriscus_bidens"          # the loop fails due to >1 candidates in the family
   #data$Binomial[which(data$MSW3 == "Cercopithecus_preussi")] <- "Allochrocebus_preussi"   # the loop fails due to >1 candidates in the family
}

# finding species-name matches in the family (for when the genus transfer notes fall short)
manual <- data$MSW3[which(is.na(data$Binomial))]
manual <- manual[which(manual %in% tree$tip.label)] 
for (i in 1:length(manual)){
    spp <- strsplit(manual[i], '_', fixed = TRUE)[[1]][2]           # get species name
    fam <- as.character(TAX$SciName[which(TAX$Family == toupper(data$Family[which(data$MSW3 == manual[i])]))]) # all sp. in the family
    # get all species names in the family
    confams <- as.character()
    for (j in 1:length(fam)){
        confams <- c(confams, strsplit(fam[j], '_', fixed = TRUE)[[1]][2])
    }
    if (length(which(confams == spp)) == 1){ # if exactly one match
        data$Binomial[which(data$MSW3 == manual[i])] <- fam[which(confams == spp)] # this should be it
    } else {print[i]} # flag to screen for manual verification
    print(c(data$MSW3[which(data$MSW3 == manual[i])], fam[which(confams == spp)])) # eyeballing results, looks good-ish
}
length(which(is.na(data$Binomial)))


# * 1.3.2 Remove species not in phylogeny or AP not in evol model ------------

# remove 24 crepuscular species
length(which(data$Diel.activity.pattern == 'Crepuscular')) == 24
dat <- data[-which(data$Diel.activity.pattern == 'Crepuscular'),]
# remove 402 species not in phylogeny
length(which(!dat$Binomial %in% tree$tip.label)) == 208 # making sure all trees using this process have the same sample
length(which(!dat$MSW3 %in% tree$tip.label)) == 141     # .. and double checking
act <- dat[-which(!dat$MSW3 %in% tree$tip.label),]

act$Diel.activity.pattern <- as.character(act$Diel.activity.pattern)
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Nocturnal'] <- 1
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Cathemeral'] <- 2
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Diurnal'] <- 3


# * 1.4 Prepare MCC tree ------------------------------------------------
# the MCC tree is not precisely ultrametric, so I use code from Liam Revell's phytools blog: 
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

is.ultrametric(MCCtree)
MCCphy <- force.ultrametric(MCCtree, "nnls")
is.ultrametric(MCCphy)
#write.tree(MCCphy, file='MCCphy.nex')

MCCphy <- drop.tip(MCCphy, MCCphy$tip.label[which(!MCCphy$tip.label %in% act$MSW3)])



# * 999 Sampling fraction -----------------------------------------------------
## calculate proportions of AP in sampled species accounting for *KNOWN* biases (incl. species not 
## in the phylogeny) given the updated taxonomy and assuming unbiased sampling for all other species. 
## The below section is written this way in case I need to make assumptions about AP of missing species.

# total count in the entire dataset (inc species not in Faurby 2015 phylogeny)
TotN <- length(which(data$Diel.activity.pattern == 'Nocturnal'))
TotC <- length(which(data$Diel.activity.pattern == 'Cathemeral'))
TotD <- length(which(data$Diel.activity.pattern == 'Diurnal'))
# proportions in data
datN <- TotN / (TotN+TotC+TotD)
datC <- TotC / (TotN+TotC+TotD)
datD <- TotD / (TotN+TotC+TotD)

## ** assumptions **
# nocturnal bats (n=1162)
Nunsamp <- length(which(tax$Order == 'CHIROPTERA')) - length(which(data$Order == 'Chiroptera')) - 6  # unsampled bats minus 6 extinct spp

# diurnal squirrels (n=113)
SQRLS <- length(which(tax$Family == 'SCIURIDAE'))                       # all sciurids (all extant)
SQRLS_DATA <- data[which(data$Family == 'Sciuridae'),]                  # sampled sciurids 
length(which(TAX$Tribe == 'PTEROMYINI')) == 57                          # all flying squirrels (nocturnal)
length(which(SQRLS_DATA$Diel.activity.pattern == "Nocturnal")) == 31    # flying squirels in data
FLY_SQRLS_unsamp <- 57-31

Dunsamp <- SQRLS - nrow(SQRLS_DATA) - FLY_SQRLS_unsamp                  # unsampled scuirids minus unsampled noct sciurids
Nunsamp <- Nunsamp + FLY_SQRLS_unsamp

# diurnal Haplorhini (n=146) 
Dprims <- c("Atelidae", "Cebidae", "Cercopithecidae", "Pitheciidae", "Hominidae", "Hylobatidae") # Simian families
TAX[which(TAX$Family %in% toupper(Dprims)),] -> primD                   # all Simiiformes INCL. AOTIDAE
length(which(primD$extinct. == 0)) == 359                               # only extant species
SIMIAN_DATA <- data[which(data$Family %in% c(Dprims, 'Aotidae')),]      # Aotidae demoted into Cebidae (Aotinae) in Burgin 2018 but data$Family relates to MSW3   
length(which(TAX$Genus == 'Aotus')) == 11                               # 11 night monkies (not all nocturnal so no assumption made)
length(which(data$Family == 'Aotidae')) == 8                            # night monkeys in data
# unsampled nocturnal primates (all extant)
length(which(TAX$Family == "LORISIDAE")) - length(which(data$Family == "Lorisidae")) == 8
length(which(TAX$Family == "TARSIIDAE")) - length(which(data$Family == "Tarsiidae")) == 8
length(which(TAX$Family == "GALAGIDAE")) - length(which(data$Family == "Galagidae")) == 9
length(which(TAX$Family == "CHEIROGALEIDAE")) - length(which(data$Family == "Cheirogaleidae")) == 26
Nprims <- 8+8+9+26

Dunsamp <- Dunsamp + 359 - nrow(SIMIAN_DATA) - (11-8)
Nunsamp <- Nunsamp + Nprims
##      ** end assumptions **       ##

# species with no data or assumptions
unkn <- 6399-(nrow(data)+Nunsamp+Dunsamp)  # total number of extant species taken from Burgin 2018 abstract
# expected species out of unknown if same proportion as in the data
expN <- datN * unkn
expC <- datC * unkn
expD <- datD * unkn
# proportion sampled out of likely total
Nprop <- TotN / (TotN + Nunsamp + expN)
Cprop <- TotC / (TotC + expC)
Dprop <- TotD / (TotD + Dunsamp + expD)

freq <- c(Nprop, Cprop, Dprop)



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

##################################################################################################################
### enough demonstrations, below is the part that was actually done:
##################################################################################################################

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

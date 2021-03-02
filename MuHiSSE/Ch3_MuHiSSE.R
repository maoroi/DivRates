
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
# 1.1 Phylogeny (see also 1.4) ------------------------------------------
Btree <- read.tree("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre")
tree <- ladderize(Btree)

tree$tip.label[which(tree$tip.label == "Nesotragus_moschatus")] <- "Neotragus_moschatus"    # typo
# plot tree for inspection if needed
pdf(file='tree.pdf', height=40, width=40)
plot(tree, type='fan', cex=0.15, label.offset = 0.05)
dev.off()

## random sample of 250 trees out of 10k from Upham 2019
indices <- sprintf("%04d", sample(0:9999, 250, replace=FALSE)) # sprintf makes all indices 4-digit long

#clades <- prop.part(tree)
#MCCtree <- maxCladeCred(tree, part = clades)


# * 1.2 Activity data  ----------------------------------------------------
read.csv("ActivityData_MDD_v1.1.csv") -> data   #species names updated to MDD v1.1

# remove crepusculars
length(which(data$AP == 'Crepuscular')) == 24
dat3 <- data[-which(data$AP == 'Crepuscular'),]

# remove 402 species not in phylogeny
length(setdiff(dat3$MDD, tree$tip.label)) == 152
length(which(!dat3$MSW3 %in% tree$tip.label)) == 141     # .. and double checking
act <- dat3[-which(!dat3$MSW3 %in% tree$tip.label),]

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

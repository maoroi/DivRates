# Updating Chapter 3 to include Faurby & Svenning (2015) phylogeny and Burgin et al. (2018) taxonomy
require("ape")
require("diversitree")
require("phangorn")
require("phytools")

### 1. Data preparation
##  1.1 Taxonomy
read.csv("MamTax2018.csv") -> TAX
TAX[,c(3,4,8,21,34,35)] -> tax      # relevant columns only (also, some species are wrongly marked EXTINCT, e.g. Alcelaphus buselaphus)

##  1.2 Phylogeny (see 1.4)
FStree <- read.tree("Small_phylogeny_4125_species.nex")
clades <- prop.part(FStree)
MCCtree <- maxCladeCred(FStree, part = clades)

##  1.3 Activity data
read.csv('Appendix_1-ActivityData.csv') -> DATA
DATA[,c(1:4)] -> dat
dat[,3] <- gsub(' ', '_', dat[,3], ignore.case = FALSE)
dat <- cbind(dat[,c(1:3)], dat[,c(3:4)])
colnames(dat) <- c(colnames(dat[1:3]), "MSW3", colnames(dat[ncol(dat)]))
dat$Binomial <- "NA"
data <- dat

# data verification
length(which(data$MSW3 %in% tax$SciName)) == 2217
length(which(MCCtree$tip.label %in% tax$SciName)) == 3613
length(which(data$MSW3 %in% MCCtree$tip.label)) == 2060

# retain unchanged binomials
unchanged <- which(data$MSW3 %in% tax$SciName)
data$Binomial[unchanged] <- data$MSW3[unchanged]

# match updated binomials (there are a few classes of changes)
updated <- which(data$MSW3 %in% tax$IfTransfer_oldSciName) # this wrongly includes "Leptonycteris_yerbabuenae" which is in MSW3
for (i in 1:length(updated)){
    data$Binomial[updated[i]] <- as.character(tax$SciName[which(tax$IfTransfer_oldSciName == data$MSW3[updated[i]])])
}
# this is a real bitch because the 'genus transfer' and 'taxonomy notes' columns are inconsistent. For some 
# transfers the previous binomial is never mentioned in the table, e.g. Bison_bison). 
# SOLUTION: finding matches for sp. name in its family (that's what the loop below does). For families w >1 match 
# I checked taxonomy and IUCN website and corrected manually)
{data$Binomial[which(data$MSW3 == "Sus_salvanius")] <- "Porcula_salvania"
data$Binomial[which(data$MSW3 == "Vulpes_rueppellii")] <- "Vulpes_rueppelli"
data$Binomial[which(data$MSW3 == "Leopardus_colocolo")] <- "Leopardus_colocola"
data$Binomial[which(data$MSW3 == "Leopardus_jacobitus")] <- "Leopardus_jacobita"
data$Binomial[which(data$MSW3 == "Prionailurus_iriomotensis")] <- "Prionailurus bengalensis" # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Herpestes_brachyurus")] <- "Urva_brachyura"
data$Binomial[which(data$MSW3 == "Herpestes_edwardsi")] <- "Urva_edwardsii"
data$Binomial[which(data$MSW3 == "Herpestes_javanicus")] <- "Urva_javanica"
data$Binomial[which(data$MSW3 == "Herpestes_naso")] <- "Xenogale_naso"
data$Binomial[which(data$MSW3 == "Herpestes_urva")] <- "Urva_urva"
data$Binomial[which(data$MSW3 == "Herpestes_semitorquatus")] <- "Urva_semitorquata"
data$Binomial[which(data$MSW3 == "Herpestes_smithii")] <- "Urva_smithii"
data$Binomial[which(data$MSW3 == "Herpestes_vitticollis")] <- "Urva_vitticolla"
data$Binomial[which(data$MSW3 == "Conepatus_humboldtii")] <- "Conepatus_chinga"
data$Binomial[which(data$MSW3 == "Aonyx_cinerea")] <- "Aonyx_cinereus"
data$Binomial[which(data$MSW3 == "Bassaricyon_beddardi")] <- "Bassaricyon_alleni"   # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Bassaricyon_lasius")] <- "Bassaricyon_gabbii"     # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Bassaricyon_pauli")] <- "Bassaricyon_gabbii"      # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Physeter_catodon")] <- "Physeter_macrocephalus"   # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Triaenops_rufus")] <- "Triaenops_menamena"
data$Binomial[which(data$MSW3 == "Neoromicia_nanus")] <- "Neoromicia_nana"
data$Binomial[which(data$MSW3 == "Chaetophractus_nationi")] <- "Chaetophractus_vellerosus"
data$Binomial[which(data$MSW3 == "Neophascogale_lorentzi")] <- "Neophascogale_lorentzii"
data$Binomial[which(data$MSW3 == "Ningaui_yvonnae")] <- "Ningaui_yvonneae"
data$Binomial[which(data$MSW3 == "Sminthopsis_fuliginosus")] <- "Sminthopsis_griseoventer"
data$Binomial[which(data$MSW3 == "Galeopterus_variegates")] <- "Galeopterus_variegatus"
data$Binomial[which(data$MSW3 == "Marmosops_dorothea")] <- "Marmosops_noctivagus"   # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Micoureus_paraguayanus")] <- "Marmosa_paraguayana"# searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Monodelphis_sorex")] <- "Monodelphis_dimidiata"   # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Monodelphis_theresa")] <- "Monodelphis_scalops"   # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Lepus_microtis")] <- "Lepus_victoriae"            # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Heterocephalus_glaber")] <- "Heterocephalus_glaber"   # wrongly named glader in Burgin 2018
data$Binomial[which(data$MSW3 == "Lagidium_peruanum")] <- "Lagidium_viscacia"       # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Abrothrix_olivaceus")] <- "Abrothrix_olivacea"
data$Binomial[which(data$MSW3 == "Akodon_molinae")] <- "Akodon_dolores"
data$Binomial[which(data$MSW3 == "Microtus_bavaricus")] <- "Microtus_liechtensteini"
data$Binomial[which(data$MSW3 == "Phyllomys_blainvillii")] <- "Phyllomys_blainvillii"   # missing one 'l' in Burgin 2018  
data$Binomial[which(data$MSW3 == "Sphiggurus_villosus")] <- "Coendou_spinosus"      # searched IUCN Red List data
data$Binomial[which(data$MSW3 == "Orthogeomys_cuniculus")] <- "Orthogeomys_grandis"
data$Binomial[which(data$MSW3 == "Pseudomys_pilligaensis")] <- "Pseudomys_pilligaensis" # not in Burgin but valid in IUCN
data$Binomial[which(data$MSW3 == "Hylopetes_lepidus")] <- "Hylopetes_sagitta"
data$Binomial[which(data$MSW3 == "Sorex_bairdi")] <- "Sorex_bairdii"
data$Binomial[which(data$MSW3 == "Abrothrix_andinus")] <- "Abrothrix_andina"
data$Binomial[which(data$MSW3 == "Kunia_fronto")] <- "Gyldenstolpia_fronto"
data$Binomial[which(data$MSW3 == "Cercopithecus_preussi")] <- "Allochrocebus_preussi"
data$Binomial[which(data$MSW3 == "Vampyressa_bidens")] <- "Vampyriscus_bidens"
data$Binomial[which(data$MSW3 == "Naemorhedus_caudatus")] <- "Naemorhedus_caudatus" # 'Nemorhaedus' in Burgin but I kept this name to comply with phylogeny
data$Binomial[which(data$MSW3 == "Naemorhedus_goral")] <- "Naemorhedus_goral"}      # 'Nemorhaedus' in Burgin but I kept this name to comply with phylogeny

manual <- data$MSW3[which(data$Binomial == "NA")]
manual <- manual[which(manual %in% MCCtree$tip.label)]    # this is really shit so only updating species included in the phylogeny
for (i in 1:length(manual)){
    spp <- strsplit(manual[i], '_', fixed = TRUE)[[1]][2]           # get species name
    fam <- as.character(TAX$SciName[which(TAX$Family == toupper(data$Family[which(data$MSW3 == manual[i])]))]) # all sp. in the family
    # get all species names in the family
    confams <- as.character()
    for (j in 1:length(fam)){
        confams <- c(confams, strsplit(fam[j], '_', fixed = TRUE)[[1]][2])
    }
    data$Binomial[which(data$MSW3 == manual[i])] <- fam[which(confams == spp)] # this should be it
}

# remove 24 crepuscular species
length(which(data$Diel.activity.pattern == 'Crepuscular')) == 24
data <- data[-which(data$Diel.activity.pattern == 'Crepuscular'),]
# remove 402 species not in phylogeny
length(which(!data$Binomial %in% MCCtree$tip.label)) == 539 # Faurby 2015 phylogeny follows MSW3 binomials
length(which(!data$MSW3 %in% MCCtree$tip.label)) == 402
act <- data[-which(!data$MSW3 %in% MCCtree$tip.label),]

act$Diel.activity.pattern <- as.character(act$Diel.activity.pattern)
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Nocturnal'] <- 1
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Cathemeral'] <- 2
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Diurnal'] <- 3

##  1.4 iron out Phylogeny
# the MCC tree is not precisely ultrametric, so I use code from Liam Revell's phytools blog: http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html
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

MCCphy <- drop.tip(MCCphy, MCCphy$tip.label[which(!MCCphy$tip.label %in% act$MSW3)])

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

MCCphy <- drop.tip(MCCphy, MCCphy$tip.label[which(!MCCphy$tip.label %in% act$MSW3)])

### 2. Preparing MuSSE arguments
states <- as.numeric(act$Diel.activity.pattern)
names(states) <- act$MSW3

## calculate proportion of sampled species correcting for *KNOWN* biases given updated taxonomy and assuming 
## unbiased sampling for all other species. 
## The below section is written this way in case any assumptions about the AP of missing species are needed

# total count in the entire dataset (inc species not in Faurby 2015 phylogeny)
TotN <- length(which(dat$Diel.activity.pattern == 'Nocturnal'))
TotC <- length(which(dat$Diel.activity.pattern == 'Cathemeral'))
TotD <- length(which(dat$Diel.activity.pattern == 'Diurnal'))
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
# diurnal Haplorhini (n=348)
Dprims <- c("ATELIDAE","CEBIDAE","CERCOPITHECIDAE","PITHECIIDAE","HOMINIDAE","HYLOBATIDAE") # Simian families
TAX[which(TAX$Family %in% Dprims),] -> primD    # all Simiiformes
length(which(primD$extinct. == 0)) == 359       # only extant species
length(which(TAX$Genus == 'Aotus')) == 11       # 11 night monkies (not all nocturnal so no assumption made)
# unsampled nocturnal primates
length(which(TAX$Family == "LORISIDAE")) - length(which(data$Family == "Lorisidae")) == 8
length(which(TAX$Family == "TARSIIDAE")) - length(which(data$Family == "Tarsiidae")) == 8
length(which(TAX$Family == "GALAGIDAE")) - length(which(data$Family == "Galagidae")) == 9
length(which(TAX$Family == "CHEIROGALEIDAE")) - length(which(data$Family == "Cheirogaleidae")) == 26

Dunsamp <- Dunsamp + (359-11)
##      ** end assumptions **       ##

# species with no data or assumptions
unkn <- 6399-(nrow(dat)+Nunsamp+Dunsamp)  # total number of extant species taken from Burgin 2018 abstract
# expected species out of unknown if same proportion as in the data
expN <- datN * unkn
expC <- datC * unkn
expD <- datD * unkn
# proportion sampled out of likely total
Nprop <- TotN / (TotN + Nunsamp + expN)
Cprop <- TotC / (TotC + expC)
Dprop <- TotD / (TotD + Dunsamp + expD)

freq <- c(Nprop, Cprop, Dprop)


### 3. MuSSE analyses (3 state dataset)

## TODO BAMM analysis

##  3.2 Mammalia wide models
lik <- make.musse(MCCphy, states, 3, sampling.f=freq, strict=TRUE, control=list())
p <- starting.point.musse(MCCphy, 3)
prior <- make.prior.exponential(1/(2 * (p[1] - p[4])))

##  3.2.1 Null model: diversification equal for all states, transition rates follow ordered model
lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1, mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q31 ~ 0)
fit.base <- find.mle(lik.base, p[argnames(lik.base)])
# running an MCMC instead of ML
samples.b <- mcmc(lik.base, coef(fit.base), nstep = 1000, w = 1, prior = prior, print.every = 50)
# saving to file
write.table(samples.b, file=paste('MCC_ALL_MuSSE_transitions_only.csv', sep=''), row.names = FALSE, 
            col.names = TRUE, sep=',', append = FALSE, quote = TRUE, eol = "\n", na = "NA", dec = ".", 
            qmethod = c("escape", "double"), fileEncoding = "")
pdf(file=paste('MCC_ALL_MuSSE_transitions_only.pdf', sep=''), height=6, width=8)
profiles.plot(samples.b[2:7], lwd = 1, col.line = c('grey10','grey60','firebrick','cornflowerblue','cyan','purple'), 
              col.fill = alpha(c('grey10','grey60','firebrick','cornflowerblue','cyan','purple'), alpha=0.8), 
              opacity = 0.2, n.br = 120)
dev.off()

##  3.2.2 Diversification unconstrained, ordered model (AP distrib patterns due to diversification rates alone, not transition rates)
lik.div <- constrain(lik, q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
fit.div <- find.mle(lik.div, p[argnames(lik.div)])
samples.d <- mcmc(lik.div, coef(fit.div), nstep = 1000, w = 1, prior = prior, print.every = 50)
write.table(samples.d, file=paste('MCC_ALL_MuSSE_diversification_only.csv', sep=''), row.names = FALSE, 
            col.names = TRUE, sep=',', append = FALSE, quote = TRUE, eol = "\n", na = "NA", dec = ".", 
            qmethod = c("escape", "double"), fileEncoding = "")
pdf(file=paste('MCC_ALL_MuSSE_diversification_only.pdf', sep=''), height=6, width=8)
profiles.plot(samples.d[2:4], lwd = 1, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), 
              col.fill = alpha(c('dodgerblue3','#22dd11','darkgoldenrod2'), alpha=0.8), opacity = 0.2, n.br = 80)
dev.off()

##  3.2.3 Diversification unconstrained, ordered character-change model
lik.free <- constrain(lik, q13 ~ 0, q31 ~ 0)        
fit.free <- find.mle(lik.free, p[argnames(lik.free)])
samples.f <- mcmc(lik.free, coef(fit.free), nstep = 1000, w = 1, prior = prior, print.every = 50)
write.table(samples.f, file=paste('MCC_ALL_MuSSE_transition_diversification.csv', sep=''), row.names = FALSE,
            col.names = TRUE, sep=',', append = FALSE, quote = TRUE, eol = "\n", na = "NA", dec = ".", 
            qmethod = c("escape", "double"), fileEncoding = "")
pdf(file=paste('MCC_ALL_MuSSE_transition_diversification.pdf', sep=''), height=6, width=8)
profiles.plot(samples.f[2:4], lwd = 1, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), 
              col.fill = alpha(c('dodgerblue3','#22dd11','darkgoldenrod2'), alpha=0.8), opacity = 0.2, n.br = 100)
dev.off()

# TODO: add legends, plot in a tiled figure.

#################################################################################################################
## this loop reads all data files for a taxon and plots their likelihood distributions to find the most likely  #
                                                                                                                #
clade <- c("ALL","Primates","Carnivora","Rodentia","Artiodactyla")                                              #
type <- c('transitions_only','diversification_only','transition_diversification')                               #
for (k in 1:5) {                                                                                                #
    for (i in 1:3) {                                                                                            #
        sam <- read.csv(paste('MCC_',clade[k],'_MuSSE_',type[i],'.csv',sep=''))                                 #
        curve <- density(sam[,ncol(sam)])                                                                       #
        if (i == 1) {                                                                                           #
            plot(curve, col=i, lwd=2, xlim=c(min(curve$x)-100,max(curve$x)+100))                                #
        } else {                                                                                                #
            lines(curve, col=i, lwd=2)                                                                          #
        }                                                                                                       #
    }                                                                                                           #
}                                                                                                               #
#################################################################################################################

##  3.3  Separate analyses for primates, Carnivora, Rodentia, and Artiodactlya
#   making a list of ancestral nodes for make.musse.split
breaks <- getMRCA(MCCphy, c('Hapalemur_griseus','Pan_paniscus'))                # primates 
breaks <- c(breaks, getMRCA(MCCphy, c('Nandinia_binotata','Lutra_lutra')))      # Carnivora 
breaks <- c(breaks, getMRCA(MCCphy, c('Vicugna_vicugna','Ovis_ammon')))         # Cetartiodactyla
    #breaks <- c(breaks, getMRCA(datatree, c('Tragulus_napu','Cervus_elaphus')))        # Ruminantia
    #breaks <- c(breaks, getMRCA(datatree, c('Phocoena_phocoena','Eschrichtius_robustus'))) # Cetacea
breaks <- c(breaks, getMRCA(MCCphy, c('Felovia_vae','Akodon_dayi')))                  # Rodentia 
    #breaks <- c(breaks, getMRCA(datatree, c('Cratogeomys_neglectus','Castor_fiber')))        # Castorimorpha
    #breaks <- c(breaks, getMRCA(datatree, c('Felovia_vae','Dasyprocta_mexicana')))         # Hystricomorpha
    #breaks <- c(breaks, getMRCA(datatree, c('Platacanthomys_lasiurus','Jaculus_jaculus'))) # Myomorpha
#breaks <- c(breaks, getMRCA(MCCphy, c('Eidolon_helvum','Lavia_frons')))              # Chiroptera
names(breaks) <- c("Primates", "Carnivora","Artiodactyla","Rodentia")

for (k in 1:length(breaks)) {
    Partdat <- dat[which(dat$Order == breaks[k]),]          # data for the order
    
    OrdTax <- TAX[which(TAX$Order == toupper(breaks[k])),]  # total extant species in new taxonomy
    totspp <- length(which(OrdTax$extinct. == 0))
    
    # total count in the entire dataset (inc species not in Faurby 2015 phylogeny)
    TotN <- length(which(Partdat$Diel.activity.pattern == 'Nocturnal'))
    TotC <- length(which(Partdat$Diel.activity.pattern == 'Cathemeral'))
    TotD <- length(which(Partdat$Diel.activity.pattern == 'Diurnal'))
    # proportions in data
    datN <- TotN / (TotN+TotC+TotD)
    datC <- TotC / (TotN+TotC+TotD)
    datD <- TotD / (TotN+TotC+TotD)
    
    ## ** assumptions **
    Nunsamp = 0
    Dunsamp = 0
    if (k == 1) {   # primates
        # diurnal simiiformes (n=348) calculated in the all-mammalia section
        Dunsamp <- 348
        length(which(TAX$Family == 'TARSIIDAE')) == 13  # 13 nocturnal tarsiers
        Nunsamp <- 13
        
    }
    if (k == 4) {   # rodentia
        # diurnal squirrels (n=144)
        Dunsamp <- length(which(tax$Family == 'SCIURIDAE')) - length(which(data$Family == 'Sciuridae'))      # unsampled sciurids (all extant)
    }
    # species with no data or assumptions
    unkn <- totspp-(nrow(Partdat)+Nunsamp+Dunsamp)  # total number of extant species taken from Burgin 2018 abstract
    # expected species out of unknown if same proportion as in the data
    expN <- datN * unkn
    expC <- datC * unkn
    expD <- datD * unkn
    # proportion sampled out of likely total
    Nprop <- TotN / (TotN + Nunsamp + expN)
    Cprop <- TotC / (TotC + expC)
    Dprop <- TotD / (TotD + Dunsamp + expD)
    
    freq <- c(Nprop, Cprop, Dprop)
    
}    


##  3.4  Partitioned analyses based on the rates shifts found in the BAMM analysis

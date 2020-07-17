# Updating Chapter 3 to include Faurby & Svenning (2015) phylogeny and Burgin et al. (2018) taxonomy
require("ape")
require("diversitree")
require("phangorn")
require("phytools")
require("scales")

# 1. Data preparation -----------------------------------------------------

# * 1.1 Taxonomy ----------------------------------------------------------
read.csv("MamTax2018.csv") -> TAX
TAX[,c(3,4,8,21,34,35)] -> tax      # relevant columns only (also, some species are wrongly marked EXTINCT, e.g. Alcelaphus buselaphus)

# * 1.2 Phylogeny (see also 1.4) ------------------------------------------
FStree <- read.tree("Faurby_Small_4125_species.nex")
clades <- prop.part(FStree)
MCCtree <- maxCladeCred(FStree, part = clades)

# * 1.3 Activity data -----------------------------------------------------
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
dat <- data[-which(data$Diel.activity.pattern == 'Crepuscular'),]
# remove 402 species not in phylogeny
length(which(!dat$Binomial %in% MCCtree$tip.label)) == 539 # Faurby 2015 phylogeny follows MSW3 binomials
length(which(!dat$MSW3 %in% MCCtree$tip.label)) == 402
act <- dat[-which(!dat$MSW3 %in% MCCtree$tip.label),]

act$Diel.activity.pattern <- as.character(act$Diel.activity.pattern)
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Nocturnal'] <- 1
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Cathemeral'] <- 2
act$Diel.activity.pattern[act$Diel.activity.pattern == 'Diurnal'] <- 3

# * 1.4 Iron out phylogeny ------------------------------------------------
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


# 2. MuSSE arguments ------------------------------------------------------
# * 2.1 STATES argument ---------------------------------------------------
states <- as.numeric(act$Diel.activity.pattern)
names(states) <- act$MSW3


# * 2.2 FREQ argument -----------------------------------------------------
## calculate proportions of APs in sampled species accounting for *KNOWN* biases (incl. species not in the 
## phylogeny) given the updated taxonomy and assuming unbiased sampling for all other species. 
## The below section is written this way in case any assumptions about the AP of missing species are needed

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


# 3. MuSSE analyses (3 state dataset)--------------------------------------

# * 3.1 BAMM analysis (extrenal code) -------------------------------------
# TODO insert BAMM input and output file names here
#### SCRAPPED for now - shifts inserted arbitrarily in roots of main mammal orders/clades ####


# * 3.2 Mammalia-wide models ----------------------------------------------
lik <- make.musse(MCCphy, states, 3, sampling.f=freq, strict=TRUE, control=list())
p <- starting.point.musse(MCCphy, 3)
prior <- make.prior.exponential(1/(2 * (p[1] - p[4])))


# ** 3.2.1 Null model: equal diversification for all APs -----------------
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
legend("topright", legend=c('Speciation','Extinction','Noct -> Cath','Cath -> Noct','Cath -> Diur','Diur  -> Cath'),
       border = NA, fill=alpha(c('grey10','grey60','firebrick','cornflowerblue','cyan','purple'), alpha=0.8), 
       cex = 1.2, bty="n")
dev.off()


# ** 3.2.2 Free diversification, all transition rates equal  -------------
#   (AP distrib patterns due to diversification alone, not due to transitions)
lik.div <- constrain(lik, q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
fit.div <- find.mle(lik.div, p[argnames(lik.div)])
samples.d <- mcmc(lik.div, coef(fit.div), nstep = 1000, w = 1, prior = prior, print.every = 50)
write.table(samples.d, file=paste('MCC_ALL_MuSSE_diversification_only.csv', sep=''), row.names = FALSE, 
            col.names = TRUE, sep=',', append = FALSE, quote = TRUE, eol = "\n", na = "NA", dec = ".", 
            qmethod = c("escape", "double"), fileEncoding = "")
pdf(file=paste('MCC_ALL_MuSSE_diversification_only.pdf', sep=''), height=6, width=8)
profiles.plot(samples.d[2:4], lwd = 1, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), 
              col.fill = alpha(c('dodgerblue4','#22cc11','darkgoldenrod2'), alpha=0.8), opacity = 0.2, n.br = 40)
dev.off()


# ** 3.2.3 Free diversification, free ordered transitions ----------------
lik.free <- constrain(lik, q13 ~ 0, q31 ~ 0)        
fit.free <- find.mle(lik.free, p[argnames(lik.free)])
samples.f <- mcmc(lik.free, coef(fit.free), nstep = 1000, w = 1, prior = prior, print.every = 50)
write.table(samples.f, file=paste('MCC_ALL_MuSSE_transition_diversification.csv', sep=''), row.names = FALSE,
            col.names = TRUE, sep=',', append = FALSE, quote = TRUE, eol = "\n", na = "NA", dec = ".", 
            qmethod = c("escape", "double"), fileEncoding = "")
pdf(file=paste('MCC_ALL_MuSSE_transition_diversification.pdf', sep=''), height=6, width=8)
profiles.plot(samples.f[2:4], lwd = 1, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), 
              col.fill = alpha(c('dodgerblue4','#22bb11','darkgoldenrod2'), alpha=0.8), opacity = 0.2, n.br = 40)
dev.off()

# TODO: add legends, plot in a tiled figure.


# * 3.3 Separate analysis for main orders ---------------------------------
# (Primates, Carnivora, Rodentia, Artiodactlya)

#   get ancestral nodes
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
    Partdat <- data[which(data$Order == names(breaks)[k]),]         # data for the order
    OrdTax <- TAX[which(TAX$Order == toupper(names(breaks)[k])),]   # total extant species in new taxonomy
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
    Dunsamp = 0
    Nunsamp = 0
    
    if (k == 1) {   # primates
        # diurnal simiiformes (n=146)
        Dunsamp <- 359 - nrow(SIMIAN_DATA) - (11-8)
        # nocturnal tarsiers, galagos & lorises (n=51)
        Nunsamp <- Nprims
    }
    if (k == 4) {   # rodentia
        # diurnal squirrels (n=87)
        Dunsamp <- SQRLS - nrow(SQRLS_DATA) - FLY_SQRLS_unsamp      # unsampled sciurids (all extant)
        # nocturnal flying squirrels (n=26) 
        Nunsamp <- FLY_SQRLS_unsamp
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
    
    
    ##  MuSSE models
    phy <- extract.clade(MCCphy, breaks[k])
    lik <- make.musse(phy, states, 3, sampling.f=freq, strict=TRUE, control=list())
    p <- starting.point.musse(phy, 3)
    prior <- make.prior.exponential(1/(2 * (p[1] - p[4])))
    
    lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1, mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q31 ~ 0)
    fit.base <- find.mle(lik.base, p[argnames(lik.base)])
    samples.b <- mcmc(lik.base, coef(fit.base), nstep = 1000, w = 1, prior = prior, print.every = 50)
    write.table(samples.b, file=paste('MCC_', names(breaks[k]), '_MuSSE_transitions_only.csv', sep=''), 
                row.names = FALSE, col.names = TRUE, sep=',', append = FALSE, quote = TRUE, eol = "\n", 
                na = "NA", dec = ".", qmethod = c("escape", "double"), fileEncoding = "")
    pdf(file=paste('MCC_', names(breaks[k]), '_MuSSE_transitions_only.pdf', sep=''), height=6, width=8)
    profiles.plot(samples.b[2:7], lwd = 1.5, col.line = c('grey10','grey40','firebrick','cornflowerblue','cyan','purple'), 
                  col.fill = alpha(c('grey10','grey60','firebrick','cornflowerblue','cyan','purple'), alpha=0.8), 
                  opacity = 0.2, n.br = 40)
    #if (k<4) {
    #    if (k==2) {
            legend("topright", legend=c('Speciation','Extinction','Noct -> Cath','Cath -> Noct','Cath -> Diur','Diur  -> Cath'),
                   border = NA, fill=alpha(c('grey10','grey60','firebrick','cornflowerblue','cyan','purple'), alpha=0.8), 
                   cex = 1.2, bty="n")
     #   } else {title(ylab = 'Probability density', line=2)}
    #}  
    dev.off()
    
    ##  Diversification unconstrained, ordered model (AP distrib patterns due to diversification rates alone, not transition rates)
    lik.div <- constrain(lik, q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
    fit.div <- find.mle(lik.div, p[argnames(lik.div)])
    samples.d <- mcmc(lik.div, coef(fit.div), nstep = 1000, w = 1, prior = prior, print.every = 50)
    write.table(samples.d, file=paste('MCC_', names(breaks[k]), '_MuSSE_diversification_only.csv', sep=''), 
                row.names = FALSE, col.names = TRUE, sep=',', append = FALSE, quote = TRUE, eol = "\n", 
                na = "NA", dec = ".", qmethod = c("escape", "double"), fileEncoding = "")
    pdf(file=paste('MCC_', names(breaks[k]), '_MuSSE_diversification_only.pdf', sep=''), height=6, width=8)
    profiles.plot(samples.d[2:4], lwd = 1.5, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), n.br = 40, 
                  col.fill = alpha(c('dodgerblue4','#22cc11','darkgoldenrod2'), alpha=0.8), opacity = 0.2)
    dev.off()
    
    ##  Diversification unconstrained, ordered character-change model
    lik.free <- constrain(lik, q13 ~ 0, q31 ~ 0)        
    fit.free <- find.mle(lik.free, p[argnames(lik.free)])
    samples.f <- mcmc(lik.free, coef(fit.free), nstep = 1000, w = 1, prior = prior, print.every = 50)
    write.table(samples.f, file=paste('MCC_', names(breaks[k]), '_MuSSE_transition_diversification.csv', sep=''), 
                row.names = FALSE, col.names = TRUE, sep=',', append = FALSE, quote = TRUE, eol = "\n", na = "NA", 
                dec = ".", qmethod = c("escape", "double"), fileEncoding = "")
    pdf(file=paste('MCC_', names(breaks[k]), '_MuSSE_transition_diversification.pdf', sep=''), height=6, width=8)
    profiles.plot(samples.f[2:4], lwd = 1.5, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), n.br = 40, 
                  col.fill = alpha(c('dodgerblue4','#22cc11','darkgoldenrod2'), alpha=0.8), opacity = 0.2)
    dev.off()
    
    ## this loop reads all data files for a taxon and plots their likelihood distributions to find the most likely  
    type <- c('transitions_only','diversification_only','transition_diversification')                               
    for (i in 1:3) {                                                                                            
        sam <- read.csv(paste('MCC_',names(breaks[k]),'_MuSSE_',type[i],'.csv',sep=''))                                 
        curve <- density(sam[,ncol(sam)])                                                                       
        if (i == 1) {                                                                                           
            plot(curve, col=i, lwd=2, xlim=c(min(curve$x)-20,max(curve$x)+50))                                
        } else {                                                                                                
            lines(curve, col=i, lwd=2)                                                                          
        }                                                                                                       
    }                                                                                                           
}    

# * 3.4 Likelihood comparisons --------------------------------------------
# this loop reads all data files for a taxon and plots their likelihood 
# distributions to find the most likely

clade <- c("ALL","Primates","Carnivora","Rodentia","Artiodactyla")              
type <- c('transitions_only','diversification_only','transition_diversification')  
for (k in 1:5) {                                                                
    for (i in 1:3) {                                                            
        sam <- read.csv(paste('MCC_',clade[k],'_MuSSE_',type[i],'.csv',sep='')) 
        curve <- density(sam[,ncol(sam)])                                       
        if (i == 1) {plot(curve, col=i, lwd=2, xlim=c(min(curve$x)-100,max(curve$x)+100))
        } else {lines(curve, col=i, lwd=2)}                                                                       
    }                                                                           
}                                                                               

# * 3.5 Partitioned analyses ----------------------------------------------

# ** 3.5.1 Define partitions ----------------------------------------------

## making a list of ancestral nodes for make.musse.split (consider adding shifts one branch above and 
## below each one of these to reflect placement uncertainty)
breaks <- getMRCA(MCCphy, c('Galago_moholi','Hylobates_lar'))                       # primates
#breaks <- c(breaks, getMRCA(datatree, c('Tarsius_bancanus','Pan_paniscus')))           # Haplorrhini
#getMRCA(datatree, c('Papio_anubis','Pan_paniscus'))     # catarrhini - Old world
#getMRCA(datatree, c('Alouatta_pigra','Callithrix_kuhlii'))  # Platyrrhini - New world
breaks <- c(breaks, getMRCA(MCCphy, c('Elephas_maximus','Neamblysomus_gunningi')))  # Afrotheria
#breaks <- c(breaks, getMRCA(datatree, c('Elephantulus_fuscus','Setifer_setosus')))     # Macroscelidea+Afrosoricida 
breaks <- c(breaks, getMRCA(MCCphy, c('Camelus_dromedarius','Ovis_ammon')))         # Cetartiodactyla
breaks <- c(breaks, getMRCA(MCCphy, c('Solenodon_cubanus','Sorex_hoyi')))           # Eulypotyphla
breaks <- c(breaks, getMRCA(MCCphy, c('Tragulus_napu','Cervus_elaphus')))           # Ruminantia
breaks <- c(breaks, getMRCA(MCCphy, c('Phocoena_phocoena','Eschrichtius_robustus')))# Cetacea
breaks <- c(breaks, getMRCA(MCCphy, c('Nyctalus_azoreum','Hipposideros_ater')))     # Chiroptera
breaks <- c(breaks, getMRCA(MCCphy, c('Panthera_leo','Lutra_lutra')))               # Carnivora
breaks <- c(breaks, getMRCA(MCCphy, c('Glis_glis','Castor_fiber')))                 # Rodentia 
breaks <- c(breaks, getMRCA(MCCphy, c('Microtus_arvalis','Castor_fiber')))          # castorimorpha
#breaks <- c(breaks, getMRCA(MCCphy, c('Platacanthomys_lasiurus','Jaculus_jaculus')))   # Myomorpha
breaks <- c(breaks, getMRCA(MCCphy, c('Glis_glis','Felovia_vae')))                  # Hystricomorpha
#breaks <- c(breaks, getMRCA(datatree, c('Glis_glis','Aplodontia_rufa')))               # Sciuromorpha
#breaks <- c(breaks, getMRCA(datatree, c('Sciurus_alleni','Cynomys_leucurus'))          # Sciuridae + Aplodontidae                    # Sciuridae


## make likelihood function 
lik.s <- make.musse.split(MCCphy, states, 3, nodes=breaks, split.t=Inf, sampling.f=freq, strict=TRUE)
p <- starting.point.musse(MCCphy, 3)
p.s <- rep(p,(length(breaks)+1))
names(p.s) <- argnames(lik.s)


# ** 3.5.2 Equal diversification among states  ----------------------------
# Div. rates independent in each tree partition, transition rates follow ordered model (this means AP distribution 
# is entirely due to character transitions and not speciation or extinction) 

# automate model formulation (this will massively pay off with usage)    
formula.base <- character()
qs.base <- character()

for (i in 1:(length(breaks)+1)) {
    for (param in c("lambda","mu")) {
        formula.base <- paste0(formula.base, paste0(', ',param,'1.',i,' ~ ',param,'2.',i,', ',param,'3.',i,' ~ ',param,'2.',i)) # the formula for equal diversification and extinction among all states in partition i 
    }
    qs.base <- paste0(qs.base, paste0(', q13.',i,'~0, q31.',i,'~0')) # disable nocturnal<->diurnal direct transitions
}
formula.base <- paste0('lik.s', formula.base, qs.base)
## Haven't found a quicker solution than a stupid copy-paste to remove the quotes in the formula
formula.base
# "lik.s, lambda1.1 ~ lambda2.1, lambda3.1 ~ lambda2.1, mu1.1 ~ mu2.1, mu3.1 ~ mu2.1, lambda1.2 ~ lambda2.2, lambda3.2 ~ lambda2.2, mu1.2 ~ mu2.2, mu3.2 ~ mu2.2, lambda1.3 ~ lambda2.3, lambda3.3 ~ lambda2.3, mu1.3 ~ mu2.3, mu3.3 ~ mu2.3, lambda1.4 ~ lambda2.4, lambda3.4 ~ lambda2.4, mu1.4 ~ mu2.4, mu3.4 ~ mu2.4, lambda1.5 ~ lambda2.5, lambda3.5 ~ lambda2.5, mu1.5 ~ mu2.5, mu3.5 ~ mu2.5, lambda1.6 ~ lambda2.6, lambda3.6 ~ lambda2.6, mu1.6 ~ mu2.6, mu3.6 ~ mu2.6, lambda1.7 ~ lambda2.7, lambda3.7 ~ lambda2.7, mu1.7 ~ mu2.7, mu3.7 ~ mu2.7, lambda1.8 ~ lambda2.8, lambda3.8 ~ lambda2.8, mu1.8 ~ mu2.8, mu3.8 ~ mu2.8, lambda1.9 ~ lambda2.9, lambda3.9 ~ lambda2.9, mu1.9 ~ mu2.9, mu3.9 ~ mu2.9, lambda1.10 ~ lambda2.10, lambda3.10 ~ lambda2.10, mu1.10 ~ mu2.10, mu3.10 ~ mu2.10, lambda1.11 ~ lambda2.11, lambda3.11 ~ lambda2.11, mu1.11 ~ mu2.11, mu3.11 ~ mu2.11, lambda1.12 ~ lambda2.12, lambda3.12 ~ lambda2.12, mu1.12 ~ mu2.12, mu3.12 ~ mu2.12, q13.1~0, q31.1~0, q13.2~0, q31.2~0, q13.3~0, q31.3~0, q13.4~0, q31.4~0, q13.5~0, q31.5~0, q13.6~0, q31.6~0, q13.7~0, q31.7~0, q13.8~0, q31.8~0, q13.9~0, q31.9~0, q13.10~0, q31.10~0, q13.11~0, q31.11~0, q13.12~0, q31.12~0"

# compute model likelihood
lik.s.base <- constrain(lik.s, 
                        lambda1.1 ~ lambda2.1, lambda3.1 ~ lambda2.1, mu1.1 ~ mu2.1, mu3.1 ~ mu2.1, 
                        lambda1.2 ~ lambda2.2, lambda3.2 ~ lambda2.2, mu1.2 ~ mu2.2, mu3.2 ~ mu2.2, 
                        lambda1.3 ~ lambda2.3, lambda3.3 ~ lambda2.3, mu1.3 ~ mu2.3, mu3.3 ~ mu2.3, 
                        lambda1.4 ~ lambda2.4, lambda3.4 ~ lambda2.4, mu1.4 ~ mu2.4, mu3.4 ~ mu2.4, 
                        lambda1.5 ~ lambda2.5, lambda3.5 ~ lambda2.5, mu1.5 ~ mu2.5, mu3.5 ~ mu2.5, 
                        lambda1.6 ~ lambda2.6, lambda3.6 ~ lambda2.6, mu1.6 ~ mu2.6, mu3.6 ~ mu2.6, 
                        lambda1.7 ~ lambda2.7, lambda3.7 ~ lambda2.7, mu1.7 ~ mu2.7, mu3.7 ~ mu2.7, 
                        lambda1.8 ~ lambda2.8, lambda3.8 ~ lambda2.8, mu1.8 ~ mu2.8, mu3.8 ~ mu2.8, 
                        lambda1.9 ~ lambda2.9, lambda3.9 ~ lambda2.9, mu1.9 ~ mu2.9, mu3.9 ~ mu2.9, 
                        lambda1.10 ~ lambda2.10, lambda3.10 ~ lambda2.10, mu1.10 ~ mu2.10, mu3.10 ~ mu2.10, 
                        lambda1.11 ~ lambda2.11, lambda3.11 ~ lambda2.11, mu1.11 ~ mu2.11, mu3.11 ~ mu2.11, 
                        lambda1.12 ~ lambda2.12, lambda3.12 ~ lambda2.12, mu1.12 ~ mu2.12, mu3.12 ~ mu2.12, 
                        q13.1~0, q31.1~0, q13.2~0, q31.2~0, q13.3~0, q31.3~0, q13.4~0, q31.4~0, 
                        q13.5~0, q31.5~0, q13.6~0, q31.6~0, q13.7~0, q31.7~0, q13.8~0, q31.8~0, 
                        q13.9~0, q31.9~0, q13.10~0, q31.10~0, q13.11~0, q31.11~0, q13.12~0, q31.12~0)                   
fit.s.base <- find.mle(lik.s.base, p.s, control=list(maxit=20000))

## Diversification unconstrained, single character transition rate
# (AP distribution due to diversification rates alone, not transition rates)
lik.s.div <- constrain(lik.s, q12.1 ~ q21.1, q23.1 ~ q21.1, q32.1 ~ q21.1, q12.2 ~ q21.2, q23.2 ~ q21.2, q32.2 ~ q21.2, 
                       q12.3 ~ q21.3, q23.3 ~ q21.3, q32.3 ~ q21.3, q12.4 ~ q21.4, q23.4 ~ q21.4, q32.4 ~ q21.4, 
                       q12.5 ~ q21.5, q23.5 ~ q21.5, q32.5 ~ q21.5, q12.6 ~ q21.6, q23.6 ~ q21.6, q32.6 ~ q21.6, 
                       q12.7 ~ q21.7, q23.7 ~ q21.7, q32.7 ~ q21.7, 
                       q13.1~0, q31.1~0, q13.2~0, q31.2~0, q13.3~0, q31.3~0, q13.4~0, q31.4~0, q13.5~0, q31.5~0, q13.6~0, q31.6~0, q13.7~0, q31.7~0)
fit.s.div <- find.mle(lik.s.div, p.s, control=list(maxit=20000))

## diversification unconstrained, ordered character transition model 
# (AP distrib patterns due to both diversification and transitions)
lik.s.free <- constrain(lik.s, q13.1~0, q31.1~0, q13.2~0, q31.2~0, q13.3~0, q31.3~0, q13.4~0, q31.4~0, q13.5~0, q31.5~0, 
                        q13.6~0, q31.6~0, q13.7~0, q31.7~0)
fit.s.free <- find.mle(lik.s.free, p.s, control=list(maxit=20000))



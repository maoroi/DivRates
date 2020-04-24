# Updating Chapter 3 to include Faurby & Svenning (2015) phylogeny and Burgin et al. (2018) taxonomy
require("ape")
require("phangorn")
require("phytools")

## preparing data
read.csv("MamTax2018.csv") -> TAX
TAX[,c(3,4,8,21,34,35)] -> tax      # relevant columns only (also, some species are wrongly marked EXTINCT, e.g. Alcelaphus buselaphus)

read.tree("Small_phylogeny_4125_species.nex") -> FStree
clades <- prop.part(FStree)
MCCtree <- maxCladeCred(FStree, part = clades)

read.csv('Appendix_1-ActivityData.csv') -> DATA
DATA[,c(1:4)] -> dat
dat[,3] <- gsub(' ', '_', dat[,3], ignore.case = FALSE)
dat <- cbind(dat[,c(1:3)], dat[,c(3:4)])
colnames(dat) <- c(colnames(dat[1:3]), "MSW3", colnames(dat[ncol(dat)]))
dat$Binomial <- "NA"
data <- dat

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
data$Binomial[which(data$MSW3 == "Sus_salvanius")] <- "Porcula_salvania"
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
data$Binomial[which(data$MSW3 == "Naemorhedus_goral")] <- "Naemorhedus_goral"       # 'Nemorhaedus' in Burgin but I kept this name to comply with phylogeny

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

## preparing MuSSE arguments
states <- as.numeric(act$Diel.activity.pattern)
names(states) <- act$MSW3

## calculate proportion of sampled species correcting for *KNOWN* biases given updated taxonomy (836 bats presumed 
## nocturnal, 94 primates and 93 squirrels presumed diurnal, and assuming unbiased sampling for all other species
## the below section is written this way in case any assumptions about the AP of missing species are needed

# total count in the entire dataset (inc species not in Faurby 2015 phylogeny)
TotN <- length(which(dat$Diel.activity.pattern == 'Nocturnal'))
TotC <- length(which(dat$Diel.activity.pattern == 'Cathemeral'))
TotD <- length(which(dat$Diel.activity.pattern == 'Diurnal'))
# proportions in data
datN <- TotN / (TotN+TotC+TotD)
datC <- TotC / (TotN+TotC+TotD)
datD <- TotD / (TotN+TotC+TotD)

## ** assumptions **
# nocturnal bats
Nunsamp <- length(which(tax$Order == 'CHIROPTERA')) - length(which(act$Order == 'Chiroptera')) - 6  # unsampled bats minus 6 extinct spp
# diurnal squirrels
Dunsamp <- length(which(tax$Family == 'SCIURIDAE')) - length(which(act$Family == 'Sciuridae'))      # unsampled sciurids (all extant)
# diurnal Haplorhini
Dprims <- c("ATELIDAE","CEBIDAE","CERCOPITHECIDAE","PITHECIIDAE","HOMINIDAE","HYLOBATIDAE")         # Haplorhine families
TAX[which(TAX$Family %in% Dprims),] -> primD    # all Haplorhines
length(which(primD$extinct. == 0)) == 359       # only extant species
length(which(TAX$Genus == 'Aotus')) == 11       # 11 night monkies to remove
Dunsamp <- Dunsamp + 348

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

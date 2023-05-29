### Addressing reviewers comments from 15/3/23
library("readr")
library("ape")
library("phytools")


# 1. Update phylogeny tips to taxonomy ------------------------------------

## get taxonomic and phylogenetic data
Tax <- read.csv("MDD_v1_6495species_JMamm.csv", header = TRUE)
Tax$phyl <- paste(Tax$SciName, Tax$Family, Tax$Order, sep="_")
tax <- Tax[which(Tax$extinct. == 0),]

# download phylogeny from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.tb03d03 (on 7th April?)
#alltrees <- read.nexus("MuHiSSE/trees/Data_S7_Mammalia_credibleTreeSets_tipDR/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_nexus.trees")

## all variants that were eventually used
allvar <- c("0028","0612","0677","1030","1166","1774","1845","2966","3865","4496","5024","8981","9128","0320","0369","0400","0647","0675","0772","0911","0917","1017","1233","1322","1350","1564","1575","1666","1777","1803","1805","1816","1823","1904","2154","2316","2417","2555","2953","2981","3019","3046","3120","3152","3168","3229","3349","3490",
            "3599","3734","3756","3766","3824","4012","4029","4224","4304","4387","4420","4423","4439","4467","4568","4683","4831","4945","5039","5133","5242","5272","5338","5340","5568","5594","5976","6255","6345","6747","6760","6774","6865","6914","6935","6947","6972","7009","7051","7072","7074","7180","7444","7549","7627","7983","8004","8156",
            "8307","8377","8394","8429","8573","8705","8775","8974","9011","9023","9455","9572","9743","9873","9923","9989","9997")

# write individual trees
#for(i in 111:1110){#000){
#   write.nexus(alltrees[i], file=paste0("raw",sprintf("%04d", i-1),".nex"))
#}

# write the relevant trees as independent files
for(i in allvar){
   tree <- alltrees[[as.numeric(i)]]
   tree <- drop.tip(tree, which(!tree$tip.label %in% tax$phyl))
   write.nexus(tree, file=paste0("tvar",i,".nex"))
}


# 2. Assess difference of the Cox 2021 dataset  ---------------------------

CoxSI2 <- read_csv("Cox2021_SupplementaryData2.csv")
# keep only AP columns
coxSI2 <- CoxSI2[,1:12]
# remove lines where activity data was imputed
cox <- coxSI2[which(coxSI2$Activity_source != "Imputed"),]
cox$Binomial <- paste(cox$Genus, cox$Species, sep='_')


# get binomials to check overlap with binomials in AP data
labs <- tree$tip.label
for(a in 1:length(labs)){
   tx <- strsplit(labs[a],"_")
   labs[a] <- paste(tx[[1]][1], tx[[1]][[2]], sep="_")
}
labs <- as.data.frame(labs)
# get species name for partial matching with data in subsection 2.1
labs$sp <- NA
for(i in 1:nrow(labs)){
   labs$sp[i] <- strsplit(labs$labs[i], "_")[[1]][2]
   }

cox$inphy <- NA
cox$inphy[which(cox$Binomial %in% labs$labs)] <- cox$Binomial[which(cox$Binomial %in% labs$labs)]
#length(which(cox$Binomial %in% labs$labs))
# [1] 4792  -- THIS IS THE NUMBER FOR "tvarxxxx" TREES. IF YOU LOAD A "tpxxxx" TREE (see section 2.2) THE NUMBER WILL BE 4108.

# Retain unchanged binomials
unchanged <- which(cox$Binomial %in% tax$SciName)
cox$MDD <- NA
cox$MDD[unchanged] <- cox$Binomial[unchanged]

#length(which(is.na(cox$MDD)))


## 2.1 Updating taxonomy -------------------------------------------------

# match names from data to (not up to date) names in the phylogeny
notinphy <- cox$Binomial[which(is.na(cox$inphy))]

# 3 useful code lines - lookup missing species among congeners in taxonomy and phylogeny
#tree$tip.label[which(str_detect(tree$tip.label, regex("Vampyressa")))]
#tax$SciName[which(str_detect(tax$IfTransfer_oldSciName, regex("juldaschi")))]
#tax[which(tax$Genus == "Vampyrodes"),]


# updating names in phylogeny
{
   tree$tip.label[which(tree$tip.label == "Dermanura_aztecus_PHYLLOSTOMIDAE_CHIROPTERA")] <- "Dermanura_azteca_PHYLLOSTOMIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Dermanura_cinereus_PHYLLOSTOMIDAE_CHIROPTERA")] <- "Dermanura_cinerea_PHYLLOSTOMIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Dermanura_glaucus_PHYLLOSTOMIDAE_CHIROPTERA")] <- "Dermanura_glauca_PHYLLOSTOMIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Dermanura_gnomus_PHYLLOSTOMIDAE_CHIROPTERA")] <- "Dermanura_gnoma_PHYLLOSTOMIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Dermanura_toltecus_PHYLLOSTOMIDAE_CHIROPTERA")] <- "Dermanura_tolteca_PHYLLOSTOMIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Tamias_sibiricus_SCIURIDAE_RODENTIA")] <- "Eutamias_sibiricus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Orthogeomys_cherriei_GEOMYIDAE_RODENTIA")] <- "Heterogeomys_cherriei_GEOMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Orthogeomys_dariensis_GEOMYIDAE_RODENTIA")] <- "Heterogeomys_dariensis_GEOMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Liomys_adspersus_HETEROMYIDAE_RODENTIA")] <- "Heteromys_adspersus_HETEROMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Liomys_irroratus_HETEROMYIDAE_RODENTIA")] <- "Heteromys_irroratus_HETEROMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Liomys_pictus_HETEROMYIDAE_RODENTIA")] <- "Heteromys_pictus_HETEROMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Liomys_salvini_HETEROMYIDAE_RODENTIA")] <- "Heteromys_salvini_HETEROMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Liomys_spectabilis_HETEROMYIDAE_RODENTIA")] <- "Heteromys_spectabilis_HETEROMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Lutra_maculicollis_MUSTELIDAE_CARNIVORA")] <- "Hydrictis_maculicollis_MUSTELIDAE_CARNIVORA"
   tree$tip.label[which(tree$tip.label == "Hypsugo_anthonyi_VESPERTILIONIDAE_CHIROPTERA")] <- "Pipistrellus_anthonyi_VESPERTILIONIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Hypsugo_joffrei_VESPERTILIONIDAE_CHIROPTERA")] <- "Pipistrellus_joffrei_VESPERTILIONIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Hypsugo_kitcheneri_VESPERTILIONIDAE_CHIROPTERA")] <- "Pipistrellus_kitcheneri_VESPERTILIONIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Hypsugo_lophurus_VESPERTILIONIDAE_CHIROPTERA")] <- "Pipistrellus_lophurus_VESPERTILIONIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Hypsugo_macrotis_VESPERTILIONIDAE_CHIROPTERA")] <- "Pipistrellus_macrotis_VESPERTILIONIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Hypsugo_savii_VESPERTILIONIDAE_CHIROPTERA")] <- "Pipistrellus_savii_VESPERTILIONIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Hypsugo_vordermanni_VESPERTILIONIDAE_CHIROPTERA")] <- "Pipistrellus_vordermanni_VESPERTILIONIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Lepus_microtis_LEPORIDAE_LAGOMORPHA")] <- "Lepus_victoriae_LEPORIDAE_LAGOMORPHA"
   tree$tip.label[which(tree$tip.label == "Scotonycteris_occidentalis_PTEROPODIDAE_CHIROPTERA")] <- "Lophostoma_occidentalis_PTEROPODIDAE_CHIROPTERA"
   tree$tip.label[which(tree$tip.label == "Pseudalopex_culpaeus_CANIDAE_CARNIVORA")] <- "Lycalopex_culpaeus_CANIDAE_CARNIVORA"
   tree$tip.label[which(tree$tip.label == "Pseudalopex_fulvipes_CANIDAE_CARNIVORA")] <- "Lycalopex_fulvipes_CANIDAE_CARNIVORA"
   tree$tip.label[which(tree$tip.label == "Pseudalopex_griseus_CANIDAE_CARNIVORA")] <- "Lycalopex_griseus_CANIDAE_CARNIVORA"
   tree$tip.label[which(tree$tip.label == "Pseudalopex_gymnocercus_CANIDAE_CARNIVORA")] <- "Lycalopex_gymnocercus_CANIDAE_CARNIVORA"
   tree$tip.label[which(tree$tip.label == "Pseudalopex_sechurae_CANIDAE_CARNIVORA")] <- "Lycalopex_sechurae_CANIDAE_CARNIVORA"
   tree$tip.label[which(tree$tip.label == "Pseudalopex_vetulus_CANIDAE_CARNIVORA")] <- "Lycalopex_vetulus_CANIDAE_CARNIVORA"
   tree$tip.label[which(tree$tip.label == "Callibella_humilis_CALLITRICHIDAE_PRIMATES")] <- "Mico_humilis_CALLITRICHIDAE_PRIMATES"
   tree$tip.label[which(tree$tip.label == "Monodelphis_unistriatus_DIDELPHIDAE_DIDELPHIMORPHIA")] <- "Monodelphis_unistriata_DIDELPHIDAE_DIDELPHIMORPHIA"
   tree$tip.label[which(tree$tip.label == "Mysateles_garridoi_CAPROMYIDAE_RODENTIA")] <- "Mysateles_garridoi_ECHIMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Mysateles_gundlachi_CAPROMYIDAE_RODENTIA")] <- "Mysateles_gundlachi_ECHIMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Mysateles_meridionalis_CAPROMYIDAE_RODENTIA")] <- "Mysateles_meridionalis_ECHIMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Mysateles_prehensilis_CAPROMYIDAE_RODENTIA")] <- "Mysateles_prehensilis_ECHIMYIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Akodon_bogotensis_CRICETIDAE_RODENTIA")] <- "Neomicroxus_bogotensis_CRICETIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Akodon_latebricola_CRICETIDAE_RODENTIA")] <- "Neomicroxus_latebricola_CRICETIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_alpinus_SCIURIDAE_RODENTIA")] <- "Neotamias_alpinus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_amoenus_SCIURIDAE_RODENTIA")] <- "Neotamias_amoenus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_bulleri_SCIURIDAE_RODENTIA")] <- "Neotamias_bulleri_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_canipes_SCIURIDAE_RODENTIA")] <- "Neotamias_canipes_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_cinereicollis_SCIURIDAE_RODENTIA")] <- "Neotamias_cinereicollis_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_dorsalis_SCIURIDAE_RODENTIA")] <- "Neotamias_dorsalis_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_durangae_SCIURIDAE_RODENTIA")] <- "Neotamias_durangae_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_merriami_SCIURIDAE_RODENTIA")] <- "Neotamias_merriami_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_minimus_SCIURIDAE_RODENTIA")] <- "Neotamias_minimus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_obscurus_SCIURIDAE_RODENTIA")] <- "Neotamias_obscurus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_ochrogenys_SCIURIDAE_RODENTIA")] <- "Neotamias_ochrogenys_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_palmeri_SCIURIDAE_RODENTIA")] <- "Neotamias_palmeri_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_panamintinus_SCIURIDAE_RODENTIA")] <- "Neotamias_panamintinus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_quadrimaculatus_SCIURIDAE_RODENTIA")] <- "Neotamias_quadrimaculatus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_quadrivittatus_SCIURIDAE_RODENTIA")] <- "Neotamias_quadrivittatus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_ruficaudus_SCIURIDAE_RODENTIA")] <- "Neotamias_ruficaudus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_rufus_SCIURIDAE_RODENTIA")] <- "Neotamias_rufus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_senex_SCIURIDAE_RODENTIA")] <- "Neotamias_senex_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_siskiyou_SCIURIDAE_RODENTIA")] <- "Neotamias_siskiyou_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_sonomae_SCIURIDAE_RODENTIA")] <- "Neotamias_sonomae_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_speciosus_SCIURIDAE_RODENTIA")] <- "Neotamias_speciosus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_townsendii_SCIURIDAE_RODENTIA")] <- "Neotamias_townsendii_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Tamias_umbrinus_SCIURIDAE_RODENTIA")] <- "Neotamias_umbrinus_SCIURIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Paralomys_gerbillus_CRICETIDAE_RODENTIA")] <- "Phyllotis_gerbillus_CRICETIDAE_RODENTIA"
   tree$tip.label[which(tree$tip.label == "Procolobus_badius_CERCOPITHECIDAE_PRIMATES")] <- "Piliocolobus_badius_CERCOPITHECIDAE_PRIMATES"
   tree$tip.label[which(tree$tip.label == "Procolobus_gordonorum_CERCOPITHECIDAE_PRIMATES")] <- "Piliocolobus_gordonorum_CERCOPITHECIDAE_PRIMATES"
   tree$tip.label[which(tree$tip.label == "Procolobus_kirkii_CERCOPITHECIDAE_PRIMATES")] <- "Piliocolobus_kirkii_CERCOPITHECIDAE_PRIMATES"
   tree$tip.label[which(tree$tip.label == "Procolobus_pennantii_CERCOPITHECIDAE_PRIMATES")] <- "Piliocolobus_pennantii_CERCOPITHECIDAE_PRIMATES"
   tree$tip.label[which(tree$tip.label == "Procolobus_preussi_CERCOPITHECIDAE_PRIMATES" )] <- "Piliocolobus_preussi_CERCOPITHECIDAE_PRIMATES"
   tree$tip.label[which(tree$tip.label == "Procolobus_rufomitratus_CERCOPITHECIDAE_PRIMATES")] <- "Piliocolobus_rufomitratus_CERCOPITHECIDAE_PRIMATES"
}

# updating names in Cox data
{
   cox$inphy[which(cox$Binomial == "Allochrocebus_lhoesti")] <- "Cercopithecus_lhoesti"
   cox$inphy[which(cox$Binomial == "Allochrocebus_preussi")] <- "Cercopithecus_preussi"
   cox$inphy[which(cox$Binomial == "Allochrocebus_solatus")] <- "Cercopithecus_solatus"
   cox$inphy[which(cox$Binomial == "Aonyx_cinereus")] <- "Aonyx_cinerea"
   cox$inphy[which(cox$Binomial ==  "Austronomus_australis")] <- "Tadarida_australis"
   cox$inphy[which(cox$Binomial ==  "Austronomus_kuboriensis")] <- "Tadarida_kuboriensis"
   cox$inphy[which(cox$Binomial ==  "Dermanura_rosenbergi")] <- "Dermanura_rosenbergii"
   cox$inphy[which(cox$Binomial ==  "Diclidurus_isabella")] <- "Diclidurus_isabellus"
   cox$inphy[which(cox$Binomial ==  "Eropeplus_1.5O")] <- "Eropeplus_canus"
   cox$inphy[which(cox$Binomial ==  "Galagoides_demidoff")] <- "Galagoides_demidovii"
   cox$inphy[which(cox$Binomial ==  "Ictonyx_libycus")] <- "Ictonyx_libyca"
   cox$inphy[which(cox$Binomial ==  "Lissonycteris_angolensis")] <- "Myonycteris_angolensis"
   cox$inphy[which(cox$Binomial ==  "Mico_chrysoleucos")] <- "Mico_chrysoleucus"
   cox$inphy[which(cox$Binomial ==  "Neodon_juldaschi")] <- "Blanfordimys_juldaschi" # Microtus in MDD
   cox$inphy[which(cox$Binomial ==  "Niviventer_langbianis")] <- "Chiromyscus_langbianis"
   cox$inphy[which(cox$Binomial ==  "Ochotona_pallasii")] <- "Ochotona_pallasi"
   cox$inphy[which(cox$Binomial ==  "Parahyaena_brunnea")] <- "Hyaena_brunnea"
   cox$inphy[which(cox$Binomial ==  "Paremballonura_atrata")] <- "Emballonura_atrata"
   cox$inphy[which(cox$Binomial ==  "Paremballonura_tiavato")] <- "Emballonura_tiavato"
   cox$inphy[which(cox$Binomial ==  "Perognathus_alticola")] <- "Perognathus_alticolus"
   cox$inphy[which(cox$Binomial ==  "Phaiomys_leucurus")] <- "Neodon_leucurus"
   cox$inphy[which(cox$Binomial ==  "Pipistrellus_anchietae")] <- "Hypsugo_anchietae"
   cox$inphy[which(cox$Binomial ==  "Pipistrellus_hesperus")] <- "Parastrellus_hesperus"
   cox$inphy[which(cox$Binomial ==  "Proechimys_trinitatis")] <- "Proechimys_trinitatus"
   cox$inphy[which(cox$Binomial ==  "Pteropus_leucopterus")] <- "Desmalopex_leucopterus"
   cox$inphy[which(cox$Binomial ==  "Saguinus_fuscicollis")] <- "Leontocebus_fuscicollis"
   cox$inphy[which(cox$Binomial ==  "Saguinus_melanoleucus")] <- "Leontocebus_fuscicollis"
   cox$inphy[which(cox$Binomial ==  "Saguinus_nigricollis")] <- "Leontocebus_nigricollis"
   cox$inphy[which(cox$Binomial ==  "Scotonycteris_ophiodon")] <- "Casinycteris_ophiodon"
   cox$inphy[which(cox$Binomial ==  "Vampyriscus_bidens")] <- "Vampyressa_bidens"
   cox$inphy[which(cox$Binomial ==  "Vampyriscus_brocki")] <- "Vampyressa_brocki"
   cox$inphy[which(cox$Binomial ==  "Vampyriscus_nymphaea")] <- "Vampyressa_nymphaea"
}

# species not in the phylogeny
{## "Apodemus_avicennicus"  # new species 2006
## "Arctonyx_albogularis"  # split from A.collaris
## "Arctonyx_hoevenii"     # split from A.collaris
## "Artibeus_schwartzi"    # split from A.jamaicensis
## "Cebus_cuscinus"        # split_from_C._albifrons
## "Cercocebus_lunulatus"  # split_from_atys
## "Cercopithecus_denti"   # split_from_pogonias
## "Cercopithecus_lowei"   # split_from_campbelli
## "Cercopithecus_roloway" # split_from_diana
## "Cerradomys_andersoni"  # no data
## "Cervus_canadensis"     # split_from_elaphus
## "Chaerephon_jobimena"   # not data
## "Cynomops_milleri"      # split_from_paranus
## "Epomophorus_minor"     # split_from_labiatus
## "Eupleres_major"        # no data
## "Galea_leucoblephara"   # newly described
## "Gazella_marica"        # no data
## "Herpestes_auropunctatus"  # Now in Urva, not on phylogeny
## "Jaculus_thaleri"       # no data
## "Lasiopodomys_fuscus"   # no data
## "Leopardus_guttulus"    # split_from_tigrina (tigrinus)
## "Miniopterus_africanus" # split_from_inflatus
## "Mormopterus_kalinowskii" # recently elevated
## "Myonycteris_leptodon"  # recently_recognized_as_a_distinct_species
## "Myotis_nyctor"         # recently_recognized_as_a_distinct_species
## "Mysateles_melanurus"   # no data
## "Ochotona_opaca"        # no data
## "Otomys_saundersiae"    # included in irroratus
## "Papio_kindae"          # split_from_cynocephalus
## "Pattonomys_carrikeri"  # no data
## "Pattonomys_flavidus"   # no data
## "Pattonomys_punctatus"  # no data
## "Perodicticus_edwardsi" # split_from_P._potto
## "Perodicticus_ibeanus"  # split_from_P._potto
## "Piliocolobus_bouvieri" # no data
## "Piliocolobus_epieni"   # no data
## "Piliocolobus_oustaleti"# no data
## "Platyrrhinus_aquilus"  # no data
## "Promops_davisoni"      # re-elevated_Species
## "Pteronotus_mesoamericanus" # elevated_subspecies
## "Pteronotus_rubiginosus"# elevated_subspecies
## "Saguinus_tripartitus"  # no data
## "Spermophilus_musicus"  # no data
## "Tympanoctomys_aureus"  # no data
## "Vampyrodes_major"      # new species
}


## 2.2 Trim tree and data to the shared set of species --------------------

# make labels like in phylogeny
cox$treelab <- NA
for(i in 1:nrow(cox)){
   if(is.na(cox$inphy[i])){
      cox$treelab[i] <- NA
   } else {
   cox$treelab[i] <- paste(cox$inphy[i], toupper(cox$Family[i]), toupper(cox$Order[i]), sep="_")
   }
}

# remove tips with no data
tips_to_drop <- which(!tree$tip.label %in% cox$treelab)
# write pruned trees as independent files
for(i in allvar){
   tree <- alltrees[[as.numeric(i)]]
   pruned <- drop.tip(tree, tips_to_drop)
   write.nexus(pruned, file=paste0("tp",i,".nex"))
}

# remove species not in the tree
cdat <- cox[which(cox$treelab %in% pruned$tip.label),]
write.csv(cdat,"Cox_data_trimmed.csv", row.names = FALSE)


# 3. Add data distribution figure -----------------------------------------

tree <- read.nexus("tp9997.nex")
cdat <- read_csv("Cox_data_trimmed.csv")

labs <- as.data.frame(tree$tip.label)
labs$tcol <- cdat$Activity_DD[match(tree$tip.label, cdat$treelab)]

table(labs$tcol)
# Cathemeral   Crepuscular    Diurnal     Nocturnal
# 350          56             654         3048
# 8.52%        1.36%          15.92%      74.19%

table(cdat$Activity_DD)/nrow(cdat)  # for reference
# Cathemeral   Crepuscular    Diurnal     Nocturnal
# 0.09579832   0.01911765     0.18172269  0.70336134


# 4. Add LTT plot (and more?) to explicitly account for time  -------------

####  REVELL CODE (from email on 10/8/22) TO ALTER  #########

### making reproducible example for simmap LTT plot
library(dplyr)
library(phytools)
packageVersion("phytools")
## [1] ‘1.5.1’

## I chose larger clades than in the blog example to see if individual (faint) LTT
## lines still appear or if thickness of the main line needs adjusting

act <- cdat %>% select(c("treelab", "Activity_DD"))
#act <- act[which(act$Activity_DD != "Crepuscular"),] # remove 56 crepuscular species
#act <- read.table("APprim.txt")  # ("APcar.txt")

AP <- factor(act$Activity_DD, levels = c("Nocturnal","Cathemeral","Diurnal"), exclude = "Crepuscular") # remove 56 crepuscular species
names(AP) <- act$treelab
str(AP)

# trim the phylogeny according to data
cld <- read.nexus(file="tp5338.nex") # file="Carnivora.nex")
cld <- drop.tip(cld, which(!cld$tip.label %in% names(AP)))

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

# 5. Summarise descriptive tree statistics --------------------------------

# generate some descriptive stats of tree sample
treesum <- as.data.frame(matrix(data = NA, nrow = length(allvar), ncol=10))
colnames(treesum) <- c("tree","Root age","Total BrLen","Mean BrLen","BrLen Var","Min BrLen","BrLrn Q1","BrLrn Q2","BrLrn Q3","Max BrLen")
for(i in 1:length(allvar)){
   phyfile <- paste0("tvar",allvar[i],".nex")
   tree <- read.nexus(phyfile)
   treesum$tree[i] <- allvar[i]
   treesum$`Root age`[i] <- max(nodeHeights(tree))
   treesum$`Total BrLen`[i] <- sum(tree$edge.length)
   treesum$`Mean BrLen`[i] <- mean(tree$edge.length)
   treesum$`BrLen Var`[i] <- var(tree$edge.length)
   treesum[i,c(5:9)] <- quantile(tree$edge.length)[1:5]
}
write.csv(treesum, file="treesum.csv", row.names = FALSE)


# 6. Repeat MuHiSSE analyses for Cox2021 dataset --------------------------

library(hisse)
#> packageVersion('hisse')
#[1] '2.1.11'
library(diversitree)
#> packageVersion('diversitree')
#[1] '0.9.18'
set.seed(88)


## 6.1 Estimating sampling fraction  --------------------------------------

# using the entire dataset (inc species not in the phylogeny)
valid <- cox[!is.na(cox$treelab),]

# total counts in data
TotN <- length(which(valid$Activity_DD == 'Nocturnal'))
TotC <- length(which(valid$Activity_DD == 'Cathemeral'))
TotD <- length(which(valid$Activity_DD == 'Diurnal'))

# proportions in data
datN <- TotN/nrow(valid)
datC <- TotC/nrow(valid)
datD <- TotD/nrow(valid)

# unsampled bats presumed nocturnal (n = 1096 in v1; 1122 in v1.1)
bats_unsamp <- length(which(tax$Order == 'CHIROPTERA')) - length(which(valid$Order == 'Chiroptera')) - 6  # unsampled bats minus 6 extinct spp

# unsampled squirrels (27 flying squirrels assumed noct; 87 other sciurids assumed diur)
SQRLS <- tax[which(tax$Family == 'SCIURIDAE'),]                         # all sciurids (all extant)
SQRLS_DATA <- valid[which(valid$Family == 'Sciuridae'),]                # sampled sciurids
flying <- length(which(SQRLS$Tribe == 'PTEROMYINI'))                    # 57 flying squirrels in MDD
noct_sqrls <- length(which(SQRLS_DATA$Activity_DD == "Nocturnal"))      # flying squirrels are the only nocturnal sciurids
FLY_SQRLS_unsamp <- flying - noct_sqrls
Diur_sqrls_unsamp <- nrow(SQRLS) - nrow(SQRLS_DATA) - FLY_SQRLS_unsamp  # unsampled scuirids minus unsampled noct sciurids


# Simiiformes presumed diurnal (n = 148, excludes Aotus)
PRIMS <- tax[which(tax$Order == "PRIMATES"),]
PRIMS <- PRIMS[which(PRIMS$extinct. == 0),]                             # only extant species
simians <- c("Atelidae", "Cebidae", "Cercopithecidae", "Pitheciidae", "Hominidae", "Hylobatidae")
simian_count <- length(which(PRIMS$Family %in% toupper(simians)))       # all simians recognised in MDD
SIMIAN_DATA <- valid[which(valid$Family %in% c(simians, "Aotidae")),]     # Aotidae demoted into Cebidae (Aotinae) in MDD but data$Family relates to MSW3
# Out of 11 species of Aotus in MDD, 7 have AP info (no assumptions about remaining 4 spp)
Aotus_unsamp <- length(which(tax$Subfamily == "AOTINAE")) - length(which(valid$Family == 'Aotidae'))
simi_diur_unsamp <- simian_count - nrow(SIMIAN_DATA) - Aotus_unsamp          # unsmapled (diur) simians minus unsampled (noct) Aotinae

# total unsampled (Noct=1123, Diur=235)
Nunsamp <- bats_unsamp + FLY_SQRLS_unsamp
Dunsamp <- Diur_sqrls_unsamp + simi_diur_unsamp
rm(bats_unsamp, FLY_SQRLS_unsamp, Diur_sqrls_unsamp, simi_diur_unsamp, PRIMS, simians,
   simian_count, SIMIAN_DATA, Aotus_unsamp, noct_sqrls, flying, SQRLS_DATA, SQRLS)
##      ** end assumptions **       ##


### Sampling fraction calculation
all_extant <- length(unique(tax$SciName[which(tax$extinct. == 0)]))

# species with no data or assumptions
unkn <- all_extant - (nrow(valid) + Nunsamp + Dunsamp)

# expected AP distribution in unknown species if same proportion as in the data
expN <- datN * unkn
expC <- datC * unkn
expD <- datD * unkn

# proportion sampled out of expected total
Nprop <- TotN / (TotN + Nunsamp + expN)
Cprop <- TotC / (TotC + expC)
Dprop <- TotD / (TotD + Dunsamp + expD)

freq <- c(Nprop, Cprop, Dprop)


## 6.2 Setting up MuHiSSE analysis ----------------------------------------

# removing crepuscular because SSE models cannot reliably deal with phenotypes rarer than ~10% of the data
act3 <- cdat[which(cdat$Activity_DD != "Crepuscular"),] # 91 spp removed

# trim the phylogeny to the data
pruned <- read.nexus(file="tp5338.nex")
tree <- drop.tip(pruned, which(!pruned$tip.label %in% act3$treelab))

# convert to numbers
act3$AP <- NA
act3$AP[act3$Activity_DD == 'Nocturnal'] <- 1
act3$AP[act3$Activity_DD == 'Cathemeral'] <- 2
act3$AP[act3$Activity_DD == 'Diurnal'] <- 3

# assign AP to tree tips
for (i in 1:length(tree$tip.label)){
   tree$tip.state[i] <- as.numeric(act3$AP[which(act3$treelab == tree$tip.label[i])])
}

# express AP in binary annotation: noct-00, cath-01, diur-11
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

# Note: “f” is the sampling fraction for each observed state, and is a vector
# ordered the same way as the turnover and extinction fraction vectors:
f = c(freq[1], freq[2], 0, freq[3])

turnover <- c(1,2,0,3)
extinction.fraction <- c(1,1,0,1)
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))  # disable transitions to/from fourth state

# run MuHiSSE to estimate parameters
mod_musse <- MuHiSSE(phy=tree, data=states.trans, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=FALSE,
                     trans.rate=trans.rate.mod, starting.vals=NULL)

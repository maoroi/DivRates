### Updating activity data to the Burgin et al. (2018, J Mamm) taxonomy. 
# This taxonomy is linked to the Mammal Diversity Database [https://www.mammaldiversity.org/] which is continually 
# updated, so newer versions of the taxonomy don't match the phylogeny of Upham et al. (2019, PLoS Biol)


#  1. Activity data -------------------------------------------------------
read.csv('Appendix_1-ActivityData.csv') -> DATA
DATA[,c(1:4)] -> DATA
colnames(DATA) <- c("Order", "Family", "MSW3", "AP")
DATA$MSW3 <- gsub(' ', '_', DATA$MSW3, ignore.case = FALSE)
data <- DATA


#  2. Updating to MDD -----------------------------------------------------
#tax_v1 <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MDD_v1_6495species_JMamm.csv", header = TRUE)
#tax_v1.1 <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MamTax2018_v1.1.csv", header = TRUE)
#tax_v1.31 <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MDD_v1.31_6513species.csv") # taxonomy from Mammal Diversity Database [https://www.mammaldiversity.org/], last accessed 23/02/21

## ** didn't have time to put the entire rocedure into this loop so focusing on MDD_v1 because that's the best
## fit to the phylogeny ** 
for (z in c("v1","v1.1","v1.31")){
    if (z == "v1"){ # MDD v1 data accessed 28/02/2021
        data <- DATA # reset data table every time this loop is run
        tax <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MDD_v1_6495species_JMamm.csv", header = TRUE)
        tax <- tax[,c(3,4,8:11,21,34,35)]
    } else if (z == "v1.1"){ #MDD v1.1 data accessed 04/02/2020)
        tax <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MamTax2018_v1.1.csv", header = TRUE)
        tax <- tax[,c(3,4,8:11,21,34,35)]
    } else if (z == "v1.31") { # MDD v1.31 data accessed 08/01/2021, 23/02/2021)
        tax <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MDD_v1.31_6513species.csv") # taxonomy from Mammal Diversity Database [https://www.mammaldiversity.org/], last accessed 23/02/21
        tax <- tax[,c(1,8,13,16,17,21,26,27,35:39)]
    }
    if ("sciName" %in% colnames(tax)) { # making column names consistent across files
        colnames(tax)[which(colnames(tax) == "sciName")] <- "SciName"
    }
    # Retain unchanged binomials
    unchanged <- which(data$MSW3 %in% tax$SciName)
    data$new_col[unchanged] <- data$MSW3[unchanged]
    # consistency test
    print(length(which(is.na(data$new_col)))) # 248, 245, 301 respectively
    
    #colnames(data)[which(colnames(data) == "new_col")] <- paste0("MDD_",z)
    #curr_col <- which(colnames(data) == paste0("MDD_",z))
}


# * 2.2 Retain unchanged binomials ----------------------------------------
unchanged <- which(data$MSW3 %in% tax$SciName)
data$MDD[unchanged] <- data$MSW3[unchanged]

# consistency test
length(which(is.na(data$MDD))) == 248 # 90% of binomials remain unchanged

##  NOTE: The change of Pteropus_giganteus to Pteropus_medius is misrepresented in the data, so the number of 
# changed binomials should be 246.

# * 2.3 Converting MSW3 taxonomy to Burgin et al. 2018 (MDD v1)------------
## There are several types of changes: genus transfer, split/elevation to full species, spelling errors, historical errors

# * 2.3.1 Match species based on MDD name-change records ------------------
updated <- which(data$MSW3 %in% tax$IfTransfer_oldSciName) # this wrongly includes "Leptonycteris_yerbabuenae" which is in MSW3
for (i in 1:length(updated)){
    data$MDD[updated[i]] <- as.character(tax$SciName[which(tax$IfTransfer_oldSciName == data$MSW3[updated[i]])])
}
# consistency test
length(which(is.na(data$MDD))) == 131


# * 2.3.2 Manual matching for exceptions ---------------------------------
# done manually because these entries trip the following loop of automatic text search
{data$MDD[which(data$MSW3 == "Monodelphis_maraxina")] <- "Monodelphis_glirina"
data$MDD[which(data$MSW3 == "Monodelphis_rubida")] <- "Monodelphis_americana" 
data$MDD[which(data$MSW3 == "Capricornis_milneedwardsii")] <- "Capricornis_sumatraensis"
data$MDD[which(data$MSW3 == "Prionailurus_iriomotensis")] <- "_" # removing from data to not alter AP of all P. bengalensis
data$MDD[which(data$MSW3 == "Bassaricyon_beddardi")] <- "Bassaricyon_alleni"              
data$MDD[which(data$MSW3 == "Pteropus_giganteus")] <- "Pteropus_medius" 
data$MDD[which(data$MSW3 == "Hypsugo_bodenheimeri")] <- "Hypsugo_ariel"
data$MDD[which(data$MSW3 == "Scotophilus_borbonicus")] <- "Scotophilus_trujilloi"
data$MDD[which(data$MSW3 == "Monodelphis_sorex")] <- "Monodelphis_dimidiata"
data$MDD[which(data$MSW3 == "Eulemur_albocollaris")] <- "Eulemur_cinereiceps" 
data$MDD[which(data$MSW3 == "Cratogeomys_neglectus")] <- "Cratogeomys_fumosus"
data$MDD[which(data$MSW3 == "Akodon_molinae")] <- "Akodon_dolores"
data$MDD[which(data$MSW3 == "Orthogeomys_thaeleri")] <- "Heterogeomys_dariensis"
data$MDD[which(data$MSW3 == "Microtus_bavaricus")] <- "Microtus_liechtensteini"
data$MDD[which(data$MSW3 == "Orthogeomys_cuniculus")] <- "Orthogeomys_grandis"
data$MDD[which(data$MSW3 == "Sorex_bairdi")] <- "Sorex_bairdii"
data$MDD[which(data$MSW3 == "Proechimys_magdalenae")] <- "Proechimys_chrysaeolus"
data$MDD[which(data$MSW3 == "Nyctophilus_timoriensis")] <- "Nyctophilus_major"  # taxon split
data$MDD[which(data$MSW3 == "Micoureus_regina")] <- "Marmosa_isthmica"          # genus transfer and split
## all species in the next block are name changes that are not flagged as such MDD_v1.1, but some are noted in v1.31     
data$MDD[which(data$MSW3 == "Sus_salvanius")] <- "Porcula_salvania"     # genus transfer and epithet change
data$MDD[which(data$MSW3 == "Leopardus_jacobitus")] <- "Leopardus_jacobita"
data$MDD[which(data$MSW3 == "Herpestes_edwardsi")] <- "Urva_edwardsii"  # new genus and epithet spelling
data$MDD[which(data$MSW3 == "Physeter_catodon")] <-"Physeter_macrocephalus"
data$MDD[which(data$MSW3 == "Neoromicia_nanus")] <-"Neoromicia_nana"    # noted as genus transfer from Pippistrellus
data$MDD[which(data$MSW3 == "Neophascogale_lorentzi")] <-"Neophascogale_lorentzii"
data$MDD[which(data$MSW3 == "Ningaui_yvonnae")] <-"Ningaui_yvonneae"
data$MDD[which(data$MSW3 == "Galeopterus_variegates")] <- "Galeopterus_variegatus"
data$MDD[which(data$MSW3 == "Micoureus_paraguayanus")] <- "Marmosa_paraguayana"
data$MDD[which(data$MSW3 == "Callithrix_argentata")] <- "Mico_argentatus"       # elevated to full genus 
data$MDD[which(data$MSW3 == "Callithrix_chrysoleuca")] <- "Mico_chrysoleucus"   # elevated to full genus 
data$MDD[which(data$MSW3 == "Callithrix_humeralifera")] <- "Mico_humeralifer"   # elevated to full genus 
data$MDD[which(data$MSW3 == "Callithrix_intermedia")] <- "Mico_intermedius"     # elevated to full genus 
data$MDD[which(data$MSW3 == "Callithrix_nigriceps")] <- "Mico_nigriceps"        # elevated to full genus 
data$MDD[which(data$MSW3 == "Galago_demidoff")] <- "Galagoides_demidovii"       # genus transfer 
data$MDD[which(data$MSW3 == "Oryzomys_macconnelli")] <- "Euryoryzomys_macconnelli"  # located by tribe (2 hits on epithet in the family)
data$MDD[which(data$MSW3 == "Oryzomys_melanotis")] <- "Handleyomys_melanotis"       # located by tribe (2 hits on epithet in the family)
data$MDD[which(data$MSW3 == "Coendou_rothschildi")] <- "Coendou_quichua"            # notes in MDD v1.31
# technically there's no need to update the below species because these were merged into existing entries w same 
# ActPat but it's easier to update and delete than to manually bypass the loop every time one of them trips it
data$MDD[which(data$MSW3 == "Alcelaphus_lichtensteinii")] <- "Alcelaphus_buselaphus"
data$MDD[which(data$MSW3 == "Pseudois_schaeferi")] <- "Pseudois_nayaur"          
data$MDD[which(data$MSW3 == "Galidictis_grandidieri")] <- "Galidictis_fasciata"
data$MDD[which(data$MSW3 == "Conepatus_humboldtii")] <- "Conepatus_chinga"
data$MDD[which(data$MSW3 == "Bassaricyon_lasius")] <- "Bassaricyon_gabbii"          
data$MDD[which(data$MSW3 == "Bassaricyon_pauli")] <- "Bassaricyon_gabbii"
data$MDD[which(data$MSW3 == "Chaetophractus_nationi")] <- "Chaetophractus_vellerosus"
data$MDD[which(data$MSW3 == "Sminthopsis_aitkeni")] <- "Sminthopsis_fuliginosus"
data$MDD[which(data$MSW3 == "Sminthopsis_griseoventer")] <- "Sminthopsis_fuliginosus"
data$MDD[which(data$MSW3 == "Marmosops_cracens")] <- "Marmosops_fuscatus"
data$MDD[which(data$MSW3 == "Marmosops_dorothea")] <- "Marmosops_noctivagus"
data$MDD[which(data$MSW3 == "Monodelphis_theresa")] <- "Monodelphis_scalops"# searched IUCN Red List data
data$MDD[which(data$MSW3 == "Trichosurus_arnhemensis")] <- "Trichosurus_vulpecula" 
data$MDD[which(data$MSW3 == "Aotus_hershkovitzi")] <- "Aotus_lemurinus"
data$MDD[which(data$MSW3 == "Alouatta_coibensis")] <- "Alouatta_palliata"
data$MDD[which(data$MSW3 == "Cercopithecus_albogularis")] <- "Cercopithecus_mitis"
data$MDD[which(data$MSW3 == "Callicebus_dubius")] <- "Plecturocebus_caligatus"
data$MDD[which(data$MSW3 == "Lagidium_peruanum")] <- "Lagidium_viscacia"    # searched IUCN Red List data
data$MDD[which(data$MSW3 == "Microtus_breweri")] <- "Microtus_pennsylvanicus"
data$MDD[which(data$MSW3 == "Neotoma_martinensis")] <- "Neotoma_bryanti"
data$MDD[which(data$MSW3 == "Oligoryzomys_delticola")] <- "Oligoryzomys_nigripes"
data$MDD[which(data$MSW3 == "Dasyprocta_cristata")] <- "Dasyprocta_leporina"
data$MDD[which(data$MSW3 == "Proechimys_poliopus")] <- "Proechimys_guairae"
data$MDD[which(data$MSW3 == "Sphiggurus_villosus")] <- "Coendou_spinosus"   # genus transfer
data$MDD[which(data$MSW3 == "Pseudomys_pilligaensis")] <- "Pseudomys_delicatulus"
data$MDD[which(data$MSW3 == "Petinomys_sagitta")] <- "Hylopetes_sagitta"    # based on notes in Burgin 2018
# For some transfers the previous binomial is never mentioned in the new taxonomy table (e.g. Bison_bison). 
# SOLUTION: finding matches for sp. name in its family (that's what the loop below does). For families w multiple 
# matches I checked taxonomy and IUCN website and corrected manually. Manual corrections are placed before the 
# loop so they don't trip it    
data$MDD[which(data$MSW3 == "Vulpes_rueppellii")] <- "Vulpes_rueppelli" 
data$MDD[which(data$MSW3 == "Leopardus_colocolo")] <- "Leopardus_colocola"
data$MDD[which(data$MSW3 == "Herpestes_brachyurus")] <- "Urva_brachyura"
data$MDD[which(data$MSW3 == "Herpestes_javanicus")] <- "Urva_javanica"
data$MDD[which(data$MSW3 == "Herpestes_semitorquatus")] <- "Urva_semitorquata"
data$MDD[which(data$MSW3 == "Herpestes_vitticollis")] <- "Urva_vitticolla"
data$MDD[which(data$MSW3 == "Aonyx_cinerea")] <- "Aonyx_cinereus"
data$MDD[which(data$MSW3 == "Triaenops_rufus")] <- "Triaenops_menamena"
data$MDD[which(data$MSW3 == "Sminthopsis_fuliginosus")] <- "Sminthopsis_griseoventer"
data$MDD[which(data$MSW3 == "Lepus_microtis")] <- "Lepus_victoriae"                # searched IUCN Red List data
data$MDD[which(data$MSW3 == "Heterocephalus_glaber")] <- "Heterocephalus_glaber"   # wrongly named glader in Burgin 2018
data$MDD[which(data$MSW3 == "Abrothrix_olivaceus")] <- "Abrothrix_olivacea"
data$MDD[which(data$MSW3 == "Phyllomys_blainvillii")] <- "Phyllomys_blainvillii"   # missing one 'l' in Burgin 2018  
data$MDD[which(data$MSW3 == "Hylopetes_lepidus")] <- "Hylopetes_sagitta"           # based on notes in Burgin 2018
data$MDD[which(data$MSW3 == "Callicebus_dubius")] <- "Plecturocebus_caligatus"
data$MDD[which(data$MSW3 == "Vampyressa_bidens")] <- "Vampyriscus_bidens"          # the loop fails due to >1 candidates in the family
data$MDD[which(data$MSW3 == "Cercopithecus_preussi")] <- "Allochrocebus_preussi"   # the loop fails due to >1 candidates in the family
}

# finding species-name matches in the family (for when the genus transfer notes fall short)
manual <- data$MSW3[which(is.na(data$MDD))]
for (i in 1:length(manual)){
    spp <- strsplit(manual[i], '_', fixed = TRUE)[[1]][2]           # get species name
    fam <- as.character(tax$SciName[which(tax$Family == toupper(data$Family[which(data$MSW3 == manual[i])]))]) # all sp. in the family
    # get all species names in the family
    confams <- as.character()
    for (j in 1:length(fam)){
        confams <- c(confams, strsplit(fam[j], '_', fixed = TRUE)[[1]][2])
    }
    if (length(which(confams == spp)) == 1){ # if exactly one match
        data$MDD[which(data$MSW3 == manual[i])] <- fam[which(confams == spp)] # this should be it
    } else {print(i)} # flag to screen for manual verification
}
# consistency test
length(which(is.na(data$MDD))) == 0 # nice


# * 2.3.3 Insert splits, remove duplicates (merges) -----------------------
# splits of species in the manual matching above
incumbent <- which(data$MDD == "Nyctophilus_major")
template <- data[incumbent,]
template$MDD <- "Nyctophilus_corbeni" # taxon split
insert <- template
template$MDD <- "Nyctophilus_sherrini" # taxon split
insert <- rbind(insert, template)
template$MDD <- "Nyctophilus_shirleyae" # taxon split
insert <- rbind(insert, template)

# merge into data
data <- rbind(data[1:incumbent,], insert, data[(incumbent+1):nrow(data),])

## and again for Marmosa
incumbent <- which(data$MDD == "Marmosa_isthmica")
template <- data[incumbent,]
template$MDD <- "Marmosa_germana" # genus transfer and split

# merge into data
data <- rbind(data[1:incumbent,], template, data[(incumbent+1):nrow(data),])

# I THINK THIS SHOULD BE WHERE THE BIG LOOP ENDS BECAUSE THE FOLLOWING SECTION MAY DIFFER BY TAXONOMY ##
########        ########        ########        ########        ########        ########        ########

# * 2.3.4 Updating AP following species merges ----------------------------
# removing all redundancies due to merging (>1 species from MSW3 become 1 MDD species)
# first deal with Bassaricyon_gabbii alone because its' the only species with more than 2 entries
AllRecs <- data[which(data$MDD == "Bassaricyon_gabbii"),]       # check that AP is consistent
data$MDD[which(data$MDD == "Bassaricyon_gabbii")[2:3]] <- NA    # flag and remove redundant entries
data <- data[-which(is.na(data$MDD)),]

dups <- data$MDD[which(duplicated(data$MDD))]                   # iterate over all spp with >1 entry
for (merged in dups){                                       
    AllRecs <- data[which(data$MDD == merged),]                 # take both entries for a species
    if (length(unique(AllRecs$AP))==1) {                        # if all AP entries are identical
        data$MDD[which(data$MDD == merged)[2]] <- NA            # flag redundant entries    
    } else {                                                    # if APs differ
        data$AP[which(data$MDD == merged)[1]] <- "Cathemeral"   # make one entry "Cathemeral"
        data$MDD[which(data$MDD == merged)[2]] <- NA            # flag the other
        print(data[which(data$MDD == merged),])
    }
}
data <- data[-which(is.na(data$MDD)),]      # remove redundant entries

write.csv(data, file="ActivityData_MDD_v1.csv", row.names = FALSE)


#  3. Updating to MDD v1.1 ------------------------------------------------
# (MDD data accessed 04/02/2020)

# * 3.1 Taxonomy ----------------------------------------------------------
tax <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MamTax2018_v1.1.csv", header = TRUE)
tax <- tax[,c(3,4,8:11,21,34,35)]

# * 3.2 Retain unchanged binomials ----------------------------------------
unchanged <- which(data$MSW3 %in% tax$SciName)
data$MDD[unchanged] <- data$MSW3[unchanged]

# consistency test
length(which(is.na(data$MDD))) == 245 # 90% of binomials remain unchanged

##  NOTE: The change of Pteropus_giganteus to Pteropus_medius is misrepresented in the data, so the number of 
# changed binomials should be 246.

# * 3.3 Converting MSW3 taxonomy to Burgin et al. 2018 (MDD v1)------------
## There are several types of changes: genus transfer, split/elevation to full species, spelling errors, historical errors

# * 3.3.1 Match species based on MDD name-change records ------------------
updated <- which(data$MSW3 %in% tax$IfTransfer_oldSciName) # this wrongly includes "Leptonycteris_yerbabuenae" which is in MSW3
for (i in 1:length(updated)){
    data$MDD[updated[i]] <- as.character(tax$SciName[which(tax$IfTransfer_oldSciName == data$MSW3[updated[i]])])
}
# consistency test
length(which(is.na(data$MDD))) == 129


# * 3.3.2 Manual matching for exceptions ---------------------------------
# done manually because these entries trip the following loop of automatic text search
{data$MDD[which(data$MSW3 == "Monodelphis_maraxina")] <- "Monodelphis_glirina"
    data$MDD[which(data$MSW3 == "Monodelphis_rubida")] <- "Monodelphis_americana" 
    data$MDD[which(data$MSW3 == "Capricornis_milneedwardsii")] <- "Capricornis_sumatraensis"
    data$MDD[which(data$MSW3 == "Prionailurus_iriomotensis")] <- "Prionailurus_bengalensis" # searched IUCN Red List data
    data$MDD[which(data$MSW3 == "Bassaricyon_beddardi")] <- "Bassaricyon_alleni"              
    data$MDD[which(data$MSW3 == "Pteropus_giganteus")] <- "Pteropus_medius" 
    data$MDD[which(data$MSW3 == "Hypsugo_bodenheimeri")] <- "Hypsugo_ariel"
    data$MDD[which(data$MSW3 == "Scotophilus_borbonicus")] <- "Scotophilus_trujilloi"
    data$MDD[which(data$MSW3 == "Monodelphis_sorex")] <- "Monodelphis_dimidiata"
    data$MDD[which(data$MSW3 == "Eulemur_albocollaris")] <- "Eulemur_cinereiceps" 
    data$MDD[which(data$MSW3 == "Cratogeomys_neglectus")] <- "Cratogeomys_fumosus"
    data$MDD[which(data$MSW3 == "Akodon_molinae")] <- "Akodon_dolores"
    data$MDD[which(data$MSW3 == "Orthogeomys_thaeleri")] <- "Heterogeomys_dariensis"
    data$MDD[which(data$MSW3 == "Microtus_bavaricus")] <- "Microtus_liechtensteini"
    data$MDD[which(data$MSW3 == "Orthogeomys_cuniculus")] <- "Orthogeomys_grandis"
    data$MDD[which(data$MSW3 == "Sorex_bairdi")] <- "Sorex_bairdii"
    data$MDD[which(data$MSW3 == "Proechimys_magdalenae")] <- "Proechimys_chrysaeolus"
    data$MDD[which(data$MSW3 == "Nyctophilus_timoriensis")] <- "Nyctophilus_major"  # taxon split
    data$MDD[which(data$MSW3 == "Micoureus_regina")] <- "Marmosa_isthmica"          # genus transfer and split
## all species in the next block are name changes that are not flagged as such MDD_v1.1, but some are noted in v1.31     
    data$MDD[which(data$MSW3 == "Sus_salvanius")] <- "Porcula_salvania"     # genus transfer and epithet change
    data$MDD[which(data$MSW3 == "Leopardus_jacobitus")] <- "Leopardus_jacobita"
    data$MDD[which(data$MSW3 == "Herpestes_edwardsi")] <- "Urva_edwardsii"  # new genus and epithet spelling
    data$MDD[which(data$MSW3 == "Physeter_catodon")] <-"Physeter_macrocephalus"
    data$MDD[which(data$MSW3 == "Neoromicia_nanus")] <-"Neoromicia_nana"    # noted as genus transfer from Pippistrellus
    data$MDD[which(data$MSW3 == "Neophascogale_lorentzi")] <-"Neophascogale_lorentzii"
    data$MDD[which(data$MSW3 == "Ningaui_yvonnae")] <-"Ningaui_yvonneae"
    data$MDD[which(data$MSW3 == "Galeopterus_variegates")] <- "Galeopterus_variegatus"
    data$MDD[which(data$MSW3 == "Micoureus_paraguayanus")] <- "Marmosa_paraguayana"
    data$MDD[which(data$MSW3 == "Callithrix_argentata")] <- "Mico_argentatus"       # elevated to full genus 
    data$MDD[which(data$MSW3 == "Callithrix_chrysoleuca")] <- "Mico_chrysoleucus"   # elevated to full genus 
    data$MDD[which(data$MSW3 == "Callithrix_humeralifera")] <- "Mico_humeralifer"   # elevated to full genus 
    data$MDD[which(data$MSW3 == "Callithrix_intermedia")] <- "Mico_intermedius"     # elevated to full genus 
    data$MDD[which(data$MSW3 == "Callithrix_nigriceps")] <- "Mico_nigriceps"        # elevated to full genus 
    data$MDD[which(data$MSW3 == "Galago_demidoff")] <- "Galagoides_demidovii"       # genus transfer 
    data$MDD[which(data$MSW3 == "Oryzomys_macconnelli")] <- "Euryoryzomys_macconnelli"  # located by tribe (2 hits on epithet in the family)
    data$MDD[which(data$MSW3 == "Oryzomys_melanotis")] <- "Handleyomys_melanotis"       # located by tribe (2 hits on epithet in the family)
    data$MDD[which(data$MSW3 == "Coendou_rothschildi")] <- "Coendou_quichua"            # notes in MDD v1.31
# technically there's no need to update the below species because these were merged into existing entries w same 
# ActPat but it's easier to update and delete than to manually bypass the loop every time one of them trips it
    data$MDD[which(data$MSW3 == "Alcelaphus_lichtensteinii")] <- "Alcelaphus_buselaphus"
    data$MDD[which(data$MSW3 == "Pseudois_schaeferi")] <- "Pseudois_nayaur"          
    data$MDD[which(data$MSW3 == "Galidictis_grandidieri")] <- "Galidictis_fasciata"
    data$MDD[which(data$MSW3 == "Conepatus_humboldtii")] <- "Conepatus_chinga"
    data$MDD[which(data$MSW3 == "Bassaricyon_lasius")] <- "Bassaricyon_gabbii"          
    data$MDD[which(data$MSW3 == "Bassaricyon_pauli")] <- "Bassaricyon_gabbii"
    data$MDD[which(data$MSW3 == "Chaetophractus_nationi")] <- "Chaetophractus_vellerosus"
    data$MDD[which(data$MSW3 == "Sminthopsis_aitkeni")] <- "Sminthopsis_fuliginosus"
    data$MDD[which(data$MSW3 == "Sminthopsis_griseoventer")] <- "Sminthopsis_fuliginosus"
    data$MDD[which(data$MSW3 == "Marmosops_cracens")] <- "Marmosops_fuscatus"
    data$MDD[which(data$MSW3 == "Marmosops_dorothea")] <- "Marmosops_noctivagus"
    data$MDD[which(data$MSW3 == "Monodelphis_theresa")] <- "Monodelphis_scalops"# searched IUCN Red List data
    data$MDD[which(data$MSW3 == "Trichosurus_arnhemensis")] <- "Trichosurus_vulpecula" 
    data$MDD[which(data$MSW3 == "Aotus_hershkovitzi")] <- "Aotus_lemurinus"
    data$MDD[which(data$MSW3 == "Alouatta_coibensis")] <- "Alouatta_palliata"
    data$MDD[which(data$MSW3 == "Cercopithecus_albogularis")] <- "Cercopithecus_mitis"
    data$MDD[which(data$MSW3 == "Callicebus_dubius")] <- "Plecturocebus_caligatus"
    data$MDD[which(data$MSW3 == "Lagidium_peruanum")] <- "Lagidium_viscacia"    # searched IUCN Red List data
    data$MDD[which(data$MSW3 == "Microtus_breweri")] <- "Microtus_pennsylvanicus"
    data$MDD[which(data$MSW3 == "Neotoma_martinensis")] <- "Neotoma_bryanti"
    data$MDD[which(data$MSW3 == "Oligoryzomys_delticola")] <- "Oligoryzomys_nigripes"
    data$MDD[which(data$MSW3 == "Dasyprocta_cristata")] <- "Dasyprocta_leporina"
    data$MDD[which(data$MSW3 == "Proechimys_poliopus")] <- "Proechimys_guairae"
    data$MDD[which(data$MSW3 == "Sphiggurus_villosus")] <- "Coendou_spinosus"   # genus transfer
    data$MDD[which(data$MSW3 == "Pseudomys_pilligaensis")] <- "Pseudomys_delicatulus"
    data$MDD[which(data$MSW3 == "Petinomys_sagitta")] <- "Hylopetes_sagitta"    # based on notes in Burgin 2018
# For some transfers the previous binomial is never mentioned in the new taxonomy table (e.g. Bison_bison). 
# SOLUTION: finding matches for sp. name in its family (that's what the loop below does). For families w multiple 
# matches I checked taxonomy and IUCN website and corrected manually. Manual corrections are placed before the 
# loop so they don't trip it    
    data$MDD[which(data$MSW3 == "Vulpes_rueppellii")] <- "Vulpes_rueppelli" 
    data$MDD[which(data$MSW3 == "Leopardus_colocolo")] <- "Leopardus_colocola"
    data$MDD[which(data$MSW3 == "Herpestes_brachyurus")] <- "Urva_brachyura"
    data$MDD[which(data$MSW3 == "Herpestes_javanicus")] <- "Urva_javanica"
    data$MDD[which(data$MSW3 == "Herpestes_semitorquatus")] <- "Urva_semitorquata"
    data$MDD[which(data$MSW3 == "Herpestes_vitticollis")] <- "Urva_vitticolla"
    data$MDD[which(data$MSW3 == "Aonyx_cinerea")] <- "Aonyx_cinereus"
    data$MDD[which(data$MSW3 == "Triaenops_rufus")] <- "Triaenops_menamena"
    data$MDD[which(data$MSW3 == "Sminthopsis_fuliginosus")] <- "Sminthopsis_griseoventer"
    data$MDD[which(data$MSW3 == "Lepus_microtis")] <- "Lepus_victoriae"                # searched IUCN Red List data
    data$MDD[which(data$MSW3 == "Heterocephalus_glaber")] <- "Heterocephalus_glaber"   # wrongly named glader in Burgin 2018
    data$MDD[which(data$MSW3 == "Abrothrix_olivaceus")] <- "Abrothrix_olivacea"
    data$MDD[which(data$MSW3 == "Phyllomys_blainvillii")] <- "Phyllomys_blainvillii"   # missing one 'l' in Burgin 2018  
    data$MDD[which(data$MSW3 == "Hylopetes_lepidus")] <- "Hylopetes_sagitta"           # based on notes in Burgin 2018
    data$MDD[which(data$MSW3 == "Callicebus_dubius")] <- "Plecturocebus_caligatus"
    data$MDD[which(data$MSW3 == "Vampyressa_bidens")] <- "Vampyriscus_bidens"          # the loop fails due to >1 candidates in the family
    data$MDD[which(data$MSW3 == "Cercopithecus_preussi")] <- "Allochrocebus_preussi"   # the loop fails due to >1 candidates in the family
}

# finding species-name matches in the family (for when the genus transfer notes fall short)
manual <- data$MSW3[which(is.na(data$MDD))]
for (i in 1:length(manual)){
    spp <- strsplit(manual[i], '_', fixed = TRUE)[[1]][2]           # get species name
    fam <- as.character(tax$SciName[which(tax$Family == toupper(data$Family[which(data$MSW3 == manual[i])]))]) # all sp. in the family
    # get all species names in the family
    confams <- as.character()
    for (j in 1:length(fam)){
        confams <- c(confams, strsplit(fam[j], '_', fixed = TRUE)[[1]][2])
    }
    if (length(which(confams == spp)) == 1){ # if exactly one match
        data$MDD[which(data$MSW3 == manual[i])] <- fam[which(confams == spp)] # this should be it
    } else {print(i)} # flag to screen for manual verification
}
# consistency test
length(which(is.na(data$MDD))) == 0 # nice


# * 3.3.3 Insert splits, remove duplicates (merges) -----------------------
# splits of species in the manual matching above
incumbent <- which(data$MDD == "Nyctophilus_major")
template <- data[incumbent,]
template$MDD <- "Nyctophilus_corbeni" # taxon split
insert <- template
template$MDD <- "Nyctophilus_sherrini" # taxon split
insert <- rbind(insert, template)
template$MDD <- "Nyctophilus_shirleyae" # taxon split
insert <- rbind(insert, template)

# merge into data
data <- rbind(data[1:incumbent,], insert, data[(incumbent+1):nrow(data),])

## and again for Marmosa
incumbent <- which(data$MDD == "Marmosa_isthmica")
template <- data[incumbent,]
template$MDD <- "Marmosa_germana" # genus transfer and split

# merge into data
data <- rbind(data[1:incumbent,], template, data[(incumbent+1):nrow(data),])


# * 3.3.4 Updating AP following species merges ----------------------------
# removing all redundancies due to merging (>1 species from MSW3 become 1 MDD species)
# first deal with Bassaricyon_gabbii alone because its' the only species with more than 2 entries
AllRecs <- data[which(data$MDD == "Bassaricyon_gabbii"),]       # check that AP is consistent
data$MDD[which(data$MDD == "Bassaricyon_gabbii")[2:3]] <- NA    # flag and remove redundant entries
data <- data[-which(is.na(data$MDD)),]

dups <- data$MDD[which(duplicated(data$MDD))]                   # iterate over all spp with >1 entry
for (merged in dups){                                       
    AllRecs <- data[which(data$MDD == merged),]                 # take both entries for a species
    if (length(unique(AllRecs$AP))==1) {                        # if all AP entries are identical
        data$MDD[which(data$MDD == merged)[2]] <- NA            # flag redundant entries    
    } else {                                                    # if APs differ
        data$AP[which(data$MDD == merged)[1]] <- "Cathemeral"   # make one entry "Cathemeral"
        data$MDD[which(data$MDD == merged)[2]] <- NA            # flag the other
        print(data[which(data$MDD == merged),])
    }
}
data <- data[-which(is.na(data$MDD)),]      # remove redundant entries

write.csv(data, file="ActivityData_MDD_v1.1.csv", row.names = FALSE)



#  4. Updating to MDD v1.31 ------------------------------------------------
# (MDD data accessed 08/01/2021, 23/02/2021)

# * 4.1 Taxonomy ----------------------------------------------------------
tax <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MDD_v1.31_6513species.csv") # taxonomy from Mammal Diversity Database [https://www.mammaldiversity.org/], last accessed 23/02/21
tax <- tax[,c(1,8,13,16,17,21,26,27,35:39)]


# * 4.2 Retain data w unchanged binomials ---------------------------------
data <- DATA    # start with original data table (section 1. in this file)
unchanged <- which(data$MSW3 %in% tax$sciName)
data$MDD[unchanged] <- data$MSW3[unchanged]
# consistency test
length(which(is.na(data$MDD))) == 301


# * 4.3 Converting MSW3 taxonomy to Burgin et al. 2018 --------------------
## There are several types of changes: genus transfer, split/elevation to full species, spelling errors, historical errors

# * 4.3.1 Match species based on MDD name-change records ------------------
for (i in 1:nrow(data)){
    if (is.na(data$MDD[i])){
        try(data$MDD[i] <- as.character(tax$sciName[which(tax$MSW3_sciName == data$MSW3[i])]), silent = TRUE)
    }
}
# consistency test
length(which(is.na(data$MDD))) == 45


# * 4.3.2 Manual matching for exceptions ----------------------------------
# done manually because automated text searches give false hits
# done manually because these entries trip the following loop of automatic text search
{data$MDD[which(data$MSW3 == "Monodelphis_maraxina")] <- "Monodelphis_glirina"
data$MDD[which(data$MSW3 == "Monodelphis_rubida")] <- "Monodelphis_americana" 
data$MDD[which(data$MSW3 == "Capricornis_milneedwardsii")] <- "Capricornis_sumatraensis"
data$MDD[which(data$MSW3 == "Prionailurus_iriomotensis")] <- "Prionailurus_bengalensis" # searched IUCN Red List data
data$MDD[which(data$MSW3 == "Bassaricyon_beddardi")] <- "Bassaricyon_alleni"              
data$MDD[which(data$MSW3 == "Pteropus_giganteus")] <- "Pteropus_medius" 
data$MDD[which(data$MSW3 == "Hypsugo_bodenheimeri")] <- "Hypsugo_ariel"
data$MDD[which(data$MSW3 == "Scotophilus_borbonicus")] <- "Scotophilus_trujilloi"
data$MDD[which(data$MSW3 == "Monodelphis_sorex")] <- "Monodelphis_dimidiata"
data$MDD[which(data$MSW3 == "Eulemur_albocollaris")] <- "Eulemur_cinereiceps" 
data$MDD[which(data$MSW3 == "Cratogeomys_neglectus")] <- "Cratogeomys_fumosus"
data$MDD[which(data$MSW3 == "Akodon_molinae")] <- "Akodon_dolores"
data$MDD[which(data$MSW3 == "Orthogeomys_thaeleri")] <- "Heterogeomys_dariensis"
data$MDD[which(data$MSW3 == "Microtus_bavaricus")] <- "Microtus_liechtensteini"
data$MDD[which(data$MSW3 == "Orthogeomys_cuniculus")] <- "Orthogeomys_grandis"
data$MDD[which(data$MSW3 == "Sorex_bairdi")] <- "Sorex_bairdii"
data$MDD[which(data$MSW3 == "Proechimys_magdalenae")] <- "Proechimys_chrysaeolus"
data$MDD[which(data$MSW3 == "Nyctophilus_timoriensis")] <- "Nyctophilus_major"  # taxon split
data$MDD[which(data$MSW3 == "Micoureus_regina")] <- "Marmosa_isthmica"          # genus transfer and split
## all species in the next block are name changes that are not flagged as such MDD_v1.1, but some are noted in v1.31     
data$MDD[which(data$MSW3 == "Sus_salvanius")] <- "Porcula_salvania"     # genus transfer and epithet change
data$MDD[which(data$MSW3 == "Leopardus_jacobitus")] <- "Leopardus_jacobita"
data$MDD[which(data$MSW3 == "Herpestes_edwardsi")] <- "Urva_edwardsii"  # new genus and epithet spelling
data$MDD[which(data$MSW3 == "Physeter_catodon")] <-"Physeter_macrocephalus"
data$MDD[which(data$MSW3 == "Neoromicia_nanus")] <-"Neoromicia_nana"    # noted as genus transfer from Pippistrellus
data$MDD[which(data$MSW3 == "Neophascogale_lorentzi")] <-"Neophascogale_lorentzii"
data$MDD[which(data$MSW3 == "Ningaui_yvonnae")] <-"Ningaui_yvonneae"
data$MDD[which(data$MSW3 == "Galeopterus_variegates")] <- "Galeopterus_variegatus"
data$MDD[which(data$MSW3 == "Micoureus_paraguayanus")] <- "Marmosa_paraguayana"
data$MDD[which(data$MSW3 == "Callithrix_argentata")] <- "Mico_argentatus"       # elevated to full genus 
data$MDD[which(data$MSW3 == "Callithrix_chrysoleuca")] <- "Mico_chrysoleucus"   # elevated to full genus 
data$MDD[which(data$MSW3 == "Callithrix_humeralifera")] <- "Mico_humeralifer"   # elevated to full genus 
data$MDD[which(data$MSW3 == "Callithrix_intermedia")] <- "Mico_intermedius"     # elevated to full genus 
data$MDD[which(data$MSW3 == "Callithrix_nigriceps")] <- "Mico_nigriceps"        # elevated to full genus 
data$MDD[which(data$MSW3 == "Galago_demidoff")] <- "Galagoides_demidovii"       # genus transfer 
data$MDD[which(data$MSW3 == "Oryzomys_macconnelli")] <- "Euryoryzomys_macconnelli"  # located by tribe (2 hits on epithet in the family)
data$MDD[which(data$MSW3 == "Oryzomys_melanotis")] <- "Handleyomys_melanotis"       # located by tribe (2 hits on epithet in the family)
data$MDD[which(data$MSW3 == "Coendou_rothschildi")] <- "Coendou_quichua"            # notes in MDD v1.31
# technically there's no need to update the below species because these were merged into existing entries w same 
# ActPat but it's easier to update and delete than to manually bypass the loop every time one of them trips it
data$MDD[which(data$MSW3 == "Alcelaphus_lichtensteinii")] <- "Alcelaphus_buselaphus"
data$MDD[which(data$MSW3 == "Pseudois_schaeferi")] <- "Pseudois_nayaur"          
data$MDD[which(data$MSW3 == "Galidictis_grandidieri")] <- "Galidictis_fasciata"
data$MDD[which(data$MSW3 == "Conepatus_humboldtii")] <- "Conepatus_chinga"
data$MDD[which(data$MSW3 == "Bassaricyon_lasius")] <- "Bassaricyon_gabbii"          
data$MDD[which(data$MSW3 == "Bassaricyon_pauli")] <- "Bassaricyon_gabbii"
data$MDD[which(data$MSW3 == "Chaetophractus_nationi")] <- "Chaetophractus_vellerosus"
data$MDD[which(data$MSW3 == "Sminthopsis_aitkeni")] <- "Sminthopsis_fuliginosus"
data$MDD[which(data$MSW3 == "Sminthopsis_griseoventer")] <- "Sminthopsis_fuliginosus"
data$MDD[which(data$MSW3 == "Marmosops_cracens")] <- "Marmosops_fuscatus"
data$MDD[which(data$MSW3 == "Marmosops_dorothea")] <- "Marmosops_noctivagus"
data$MDD[which(data$MSW3 == "Monodelphis_theresa")] <- "Monodelphis_scalops"# searched IUCN Red List data
data$MDD[which(data$MSW3 == "Trichosurus_arnhemensis")] <- "Trichosurus_vulpecula" 
data$MDD[which(data$MSW3 == "Aotus_hershkovitzi")] <- "Aotus_lemurinus"
data$MDD[which(data$MSW3 == "Alouatta_coibensis")] <- "Alouatta_palliata"
data$MDD[which(data$MSW3 == "Cercopithecus_albogularis")] <- "Cercopithecus_mitis"
data$MDD[which(data$MSW3 == "Callicebus_dubius")] <- "Plecturocebus_caligatus"
data$MDD[which(data$MSW3 == "Lagidium_peruanum")] <- "Lagidium_viscacia"    # searched IUCN Red List data
data$MDD[which(data$MSW3 == "Microtus_breweri")] <- "Microtus_pennsylvanicus"
data$MDD[which(data$MSW3 == "Neotoma_martinensis")] <- "Neotoma_bryanti"
data$MDD[which(data$MSW3 == "Oligoryzomys_delticola")] <- "Oligoryzomys_nigripes"
data$MDD[which(data$MSW3 == "Dasyprocta_cristata")] <- "Dasyprocta_leporina"
data$MDD[which(data$MSW3 == "Proechimys_poliopus")] <- "Proechimys_guairae"
data$MDD[which(data$MSW3 == "Sphiggurus_villosus")] <- "Coendou_spinosus"   # genus transfer
data$MDD[which(data$MSW3 == "Pseudomys_pilligaensis")] <- "Pseudomys_delicatulus"
data$MDD[which(data$MSW3 == "Petinomys_sagitta")] <- "Hylopetes_sagitta"    # based on notes in Burgin 2018
# For some transfers the previous binomial is never mentioned in the new taxonomy table (e.g. Bison_bison). 
# SOLUTION: finding matches for sp. name in its family (that's what the loop below does). For families w multiple 
# matches I checked taxonomy and IUCN website and corrected manually. Manual corrections are placed before the 
# loop so they don't trip it    
data$MDD[which(data$MSW3 == "Vulpes_rueppellii")] <- "Vulpes_rueppelli" 
data$MDD[which(data$MSW3 == "Leopardus_colocolo")] <- "Leopardus_colocola"
data$MDD[which(data$MSW3 == "Herpestes_brachyurus")] <- "Urva_brachyura"
data$MDD[which(data$MSW3 == "Herpestes_javanicus")] <- "Urva_javanica"
data$MDD[which(data$MSW3 == "Herpestes_semitorquatus")] <- "Urva_semitorquata"
data$MDD[which(data$MSW3 == "Herpestes_vitticollis")] <- "Urva_vitticolla"
data$MDD[which(data$MSW3 == "Aonyx_cinerea")] <- "Aonyx_cinereus"
data$MDD[which(data$MSW3 == "Triaenops_rufus")] <- "Triaenops_menamena"
data$MDD[which(data$MSW3 == "Sminthopsis_fuliginosus")] <- "Sminthopsis_griseoventer"
data$MDD[which(data$MSW3 == "Lepus_microtis")] <- "Lepus_victoriae"                # searched IUCN Red List data
data$MDD[which(data$MSW3 == "Heterocephalus_glaber")] <- "Heterocephalus_glaber"   # wrongly named glader in Burgin 2018
data$MDD[which(data$MSW3 == "Abrothrix_olivaceus")] <- "Abrothrix_olivacea"
data$MDD[which(data$MSW3 == "Phyllomys_blainvillii")] <- "Phyllomys_blainvillii"   # missing one 'l' in Burgin 2018  
data$MDD[which(data$MSW3 == "Hylopetes_lepidus")] <- "Hylopetes_sagitta"           # based on notes in Burgin 2018
data$MDD[which(data$MSW3 == "Callicebus_dubius")] <- "Plecturocebus_caligatus"
data$MDD[which(data$MSW3 == "Vampyressa_bidens")] <- "Vampyriscus_bidens"          # the loop fails due to >1 candidates in the family
data$MDD[which(data$MSW3 == "Cercopithecus_preussi")] <- "Allochrocebus_preussi"   # the loop fails due to >1 candidates in the family
}
# consistency test
length(which(is.na(data$MDD))) == 0 # good


# * 4.3.3 Insert splits, remove duplicates (merges) -----------------------
# splits of species in the last two lines of manual matching above
incumbent <- which(data$MDD == "Nyctophilus_major")
template <- data[incumbent,]
template$MDD <- "Nyctophilus_corbeni" # taxon split
insert <- template
template$MDD <- "Nyctophilus_sherrini" # taxon split
insert <- rbind(insert, template)
template$MDD <- "Nyctophilus_shirleyae" # taxon split
insert <- rbind(insert, template)

# merge into data
data <- rbind(data[1:incumbent,], insert, data[(incumbent+1):nrow(data),])

## and again for Marmosa
incumbent <- which(data$MDD == "Marmosa_isthmica")
template <- data[incumbent,]
template$MDD <- "Marmosa_germana" # genus transfer and split

# merge into data
data <- rbind(data[1:incumbent,], template, data[(incumbent+1):nrow(data),])

## removing all entries not updated to burgin 2018 - these are redundancies due to merging (>1 species from MSW3 become 1 MDD species)
data <- data[-which(is.na(data$MDD)),]  # empty


# * 4.3.4 Updating AP following species merges ----------------------------
# the remaining duplicate records are species that contain merged species with different APs so I made them cathemeral

dups <- data$MDD[which(duplicated(data$MDD))]                   # iterate over all spp with >1 entry
for (merged in dups){                                       
    AllRecs <- data[which(data$MDD == merged),]                 # take both entries for a species
    if (length(unique(AllRecs$AP))==1) {                        # if all AP entries are identical
        data$MDD[which(data$MDD == merged)[2]] <- NA            # flag redundant entries    
    } else {                                                    # if APs differ
        data$AP[which(data$MDD == merged)[1]] <- "Cathemeral"   # make one entry "Cathemeral"
        data$MDD[which(data$MDD == merged)[2]] <- NA            # flag the other
        print(data[which(data$MDD == merged),])
    }
}
data <- data[-which(is.na(data$MDD)),]      # remove redundant entries
# save 
write.csv(data, file="ActivityData_MDD_v1.31.csv", row.names = FALSE)



#  5. Matching names in data to phylogeny ---------------------------------
# (suprisingly this is NOT straightforward)
# * 5.1 Phylogeny (see also 1.4) ------------------------------------------
Btree <- read.tree("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre")
tree <- ladderize(Btree)

# plot tree for inspection if needed
pdf(file='tree.pdf', height=40, width=40)
plot(tree, type='fan', cex=0.15, label.offset = 0.05)
dev.off()

## random sample of 250 trees out of 10k from Upham 2019
indices <- sprintf("%04d", sample(0:9999, 250, replace=FALSE)) # sprintf makes all indices 4-digit long

#clades <- prop.part(tree)
#MCCtree <- maxCladeCred(tree, part = clades)

# * 5.2 Adjusting names in data to match phylogeny ------------------------

data <- read.csv(file="ActivityData_MDD_v1.csv")
# checks below use the "%in%" syntax becasue setdiff() returns unique entries only, messing up the totals 
length(which(!data$MDD %in% tree$tip.label)) == 153     # checking names match - they don't
length(which(!data$MSW3 %in% tree$tip.label)) == 136    # ...double checking - somehow MSW3 match better

# new columns for the names that are actually in the phylogeny, and where they're taken from
data$Phylo_name <- NA
data$Phylo_name_note <- NA

# Species names were updated in the data table's column 'Phylo_name': , 
# if it isn't but the MSW3 form is - MSW3 was retained, and 
for (i in 1:nrow(data)){ 
    if (as.character(data$MDD[i]) %in% tree$tip.label){ # if the MDD form is in the tree, it was retained
        data$Phylo_name[i] <- as.character(data$MDD[i])
        data$Phylo_name_note[i] <- "MDD_v1.1"
    } else if (as.character(data$MSW3[i]) %in% tree$tip.label){ # if MDD isn't but the MSW3 form is - MSW3 was retained
        data$Phylo_name[i] <- as.character(data$MSW3[i])
        data$Phylo_name_note[i] <- "MSW3"
    } else {                                            # if neither taxonomies match, manually verified form was inserted (see below)
        data$Phylo_name[i] <- "_"
        data$Phylo_name_note[i] <- "_"
    } 
}

## Correcting mismatches manually is the safest solution:

# find the species for which both MSW3 and MDD are not in the tree
new_out <- which(!data$MDD %in% tree$tip.label)
old_out <- which(!data$MSW3 %in% tree$tip.label)

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



### Updating activity data to the Burgin et al. (2018, J Zool) taxonomy. 
# This taxonomy is linked to the Mammal Diversity Database [https://www.mammaldiversity.org/] which is continually 
# updated, so newer versions of the taxonomy don't match the phylogeny of Upham et al. (2019, PLoS Biol)


#  1. Activity data -------------------------------------------------------
read.csv('Appendix_1-ActivityData.csv') -> DATA
DATA[,c(1:4)] -> DATA
DATA[,3] <- gsub(' ', '_', DATA[,3], ignore.case = FALSE)
DATA <- cbind(DATA[,c(1:3)], DATA[,c(3:4)])
colnames(DATA) <- c("Order", "Family", "MSW3", "MDD", "AP")
DATA$MDD <- NA
data <- DATA


#  2. Updating to MDD v1.0 ------------------------------------------------
# (MDD data accessed 07/01/2019)

# * 2.1 Taxonomy ----------------------------------------------------------
read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/asm-species-2019-01-07.csv", header = TRUE) -> TAX
TAX[-1,c(1,2,4,15,16,19,20)] -> tax 

# * 2.2 Retain unchanged binomials ----------------------------------------


# 3. Updating to MDD v1.31 ------------------------------------------------
# (MDD data accessed 08/01/2021, 23/02/2021)

# * 3.1 Taxonomy ----------------------------------------------------------
read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MDD_v1.31_6513species.csv") -> TAX   # taxonomy from Mammal Diversity Database [https://www.mammaldiversity.org/], last accessed 23/02/21
TAX[,c(1,8,13,16,17,21,26,27,35:39)] -> tax


# * 3.2 Retain data w unchanged binomials ---------------------------------
unchanged <- which(data$MSW3 %in% tax$sciName)
data$MDD[unchanged] <- data$MSW3[unchanged]
# consistency test
length(which(is.na(data$MDD))) == 301


# * 3.3 Converting MSW3 taxonomy to Burgin et al. 2018 --------------------
## There are several types of changes: genus transfer, split/elevation to full species, spelling errors, historical errors

# * 3.3.1 Match species based on MDD name-change records ------------------
for (i in 1:nrow(data)){
    if (is.na(data$MDD[i])){
        try(data$MDD[i] <- as.character(tax$sciName[which(tax$MSW3_sciName == data$MSW3[i])]), silent = TRUE)
    }
}
# consistency test
length(which(is.na(data$MDD))) == 45


#  * 3.3.2 Manual matching for exceptions ---------------------------------
# done manually because automated text searches give false hits
{data$MDD[which(data$MSW3 == "Capricornis_milneedwardsii")] <- "Capricornis_sumatraensis"
    data$MDD[which(data$MSW3 == "Prionailurus_iriomotensis")] <- "Prionailurus bengalensis" # searched IUCN Red List data
    data$MDD[which(data$MSW3 == "Bassaricyon_beddardi")] <- "Bassaricyon_alleni"              
    data$MDD[which(data$MSW3 == "Pteropus_giganteus")] <- "Pteropus_medius" 
    data$MDD[which(data$MSW3 == "Hypsugo_bodenheimeri")] <- "Hypsugo_ariel"
    data$MDD[which(data$MSW3 == "Scotophilus_borbonicus")] <- "Scotophilus_trujilloi"
    data$MDD[which(data$MSW3 == "Monodelphis_sorex")] <- "Monodelphis_dimidiata"
    data$MDD[which(data$MSW3 == "Eulemur_albocollaris")] <- "Eulemur_cinereiceps" 
    data$MDD[which(data$MSW3 == "Cratogeomys_neglectus")] <- "Cratogeomys_fumosus"
    data$MDD[which(data$MSW3 == "Monodelphis_maraxina")] <- "Monodelphis_glirina"
    data$MDD[which(data$MSW3 == "Monodelphis_rubida")] <- "Monodelphis_americana" 
    data$MDD[which(data$MSW3 == "Akodon_molinae")] <- "Akodon_dolores"
    data$MDD[which(data$MSW3 == "Orthogeomys_thaeleri")] <- "Heterogeomys_dariensis"
    data$MDD[which(data$MSW3 == "Microtus_bavaricus")] <- "Microtus_liechtensteini"
    data$MDD[which(data$MSW3 == "Orthogeomys_cuniculus")] <- "Orthogeomys_grandis"
    data$MDD[which(data$MSW3 == "Sorex_bairdi")] <- "Sorex_bairdii"
    data$MDD[which(data$MSW3 == "Proechimys_magdalenae")] <- "Proechimys_chrysaeolus"
    data$MDD[which(data$MSW3 == "Nyctophilus_timoriensis")] <- "Nyctophilus_major" # taxon split
    data$MDD[which(data$MSW3 == "Micoureus_regina")] <- "Marmosa_isthmica" # genus transfer and split
}   
# consistency test
length(which(is.na(data$MDD))) == 26


# * 3.3.4 Insert splits, remove duplicates (merges) -----------------------
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
data <- data[-which(is.na(data$MDD)),]


# * 3.3.5 Updating AP following species merges ----------------------------
# the remaining duplicate records are species that contain merged species with different APs so I made them cathemeral
data$MDD[which(duplicated(data$MDD))]
#[1] "Monodelphis_americana"  "Heterogeomys_dariensis"

# first remove excess records
rm1 <- which(data$MSW3 == "Monodelphis_rubida")
rm2 <- which(data$MSW3 == "Orthogeomys_thaeleri")
data <- data[-c(rm1, rm2),]
# then change the AP of the correspondigng record to cathemeral
rm1 <- which(data$MDD == "Monodelphis_americana")
rm2 <- which(data$MDD == "Heterogeomys_dariensis")
data$AP[c(rm1,rm2)] <- "Cathemeral"
# save 
write.csv(data, file="ActivityData_MDD_v1.31.csv", row.names = FALSE)

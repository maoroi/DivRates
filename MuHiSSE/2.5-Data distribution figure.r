### making the data distribution figure

library(ape)
library(phytools)
library(phangorn)
library(scales)

setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")

# 1. prepare the data ---------------------------------------------------------

## sections 2.1 and 4 from "2-Preparing for MuHiSSE.r"
# automatically correct tip names to match taxonomy - requires the data file loaded (based on section 1.4 below)
phy1 <- read.tree("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre")
tree <- ladderize(phy1)

data <- read.csv(file="ActivityData_MDD_v1_match.csv")
data$Phylo_name <- NA
data$Phylo_name_conforms_to <- NA

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


# 2. match activity data and tree ---------------------------------------------

unwantedspecies <- which(is.na(pmatch(tree$tip.label,as.character(data$Phylo_name))) == TRUE) 
datatree <- drop.tip(tree,unwantedspecies)
trimdat <- data[which(data$Phylo_name %in% datatree$tip.label),] 


# 3. prepare labels -----------------------------------------------------------


## activity pattern labels for tips
ACTLabel <- rep("magenta", length(datatree$tip.label))
names(ACTLabel) <- trimdat$Phylo_name
ACTLabel[which(trimdat$AP == "Nocturnal")] <- "dodgerblue3"
ACTLabel[which(trimdat$AP == "Cathemeral")] <- "#22dd11"
ACTLabel[which(trimdat$AP == "Diurnal")] <- "gold"

# branch colours by order
# get MRCAs
mar <- getMRCA(datatree, c('Petaurus_australis','Lestoros_inca'))    # marsupials
afr <- getMRCA(datatree, c('Dugong_dugon','Chrysochloris_asiatica')) # afrotheria
rod <- getMRCA(datatree, c('Castor_canadensis','Sciurus_niger'))     # Rodentia
cet <- getMRCA(datatree, c('Camelus_dromedarius','Axis_kuhlii'))     # Cetartiodactyla
pri <- getMRCA(datatree, c('Indri_indri','Cebus_albifrons'))         # Primates
chi <- getMRCA(datatree, c('Pteropus_alecto','Noctilio_leporinus'))  # Chiroptera
car <- getMRCA(datatree, c('Felis_manul','Nasua_nasua'))             # Carnivora
sor <- getMRCA(datatree, c('Solenodon_cubanus','Talpa_levantis'))    # Soricomorpha

# assign clade colours
orders <- c(mar,afr,rod,cet,pri,chi,car,sor)
colours <- c("deeppink1","#4d1d0ad0","darkorange2","goldenrod1","purple4","blue","gray60","green4")

EDGE <- as.matrix(datatree$edge)
EDGE <- cbind(EDGE, "black")
names(EDGE) <- names(datatree$edge.length)
# loop over edge from the ancestor node to all descendants to assign colour
for (k in 1:length(orders)) {
    AncNode <- orders[k]
    CladeSpp <- Descendants(datatree, AncNode, type="all")
    edges <- EDGE[which(EDGE[,2] %in% CladeSpp),]
    EDGE[which(EDGE[,2] == AncNode),3] <- as.character(colours[k]) ## colouring the order root branch
    for (i in 1:length(edges[,2])) {
        edges[i,3] <- as.character(colours[k])
        EDGE[which(EDGE[,2] == edges[i,2]),3] <- edges[i,3]
    }
} 

pdf(file="Figure 1 - Data distribution.pdf",width=40, height=40) 
plot(datatree, type = "fan", open.angle = 7, lwd = 0.4, label.offset = 0.4, cex = 0.3, show.tip.label = TRUE, 
     tip.color = "white", edge.color = EDGE[,3])
ring(20, datatree, style = "ring", offset = 1.5, col = ACTLabel[datatree$tip.label])
ring(-101.642, datatree, style = "ring", offset = -66.02, col = alpha("#67c5ca", 0.25))
axis(side = 1, at = c(17.662, 42.662, 67.662, 92.662, 117.662, 142.662), 
     labels = c(150, 125, 100, 75, 50, 25), tick = TRUE, pos = -0.3, cex.axis = 4, lwd = 3, lwd.ticks = 6, padj = 0.4) 
text(158, -5, "Mya", cex = 4.5)
dev.off()


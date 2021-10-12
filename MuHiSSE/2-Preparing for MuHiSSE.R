### This file contains data preparation steps for running MuHiSSE models ###

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
library(phangorn)
#> packageVersion('phangorn')
#[1] ‘2.5.5’
require(scales)
#> packageVersion('scales')
#[1] ‘1.0.0’

## reproducible evironemnt using {renv}: The packages used in this project are recorded into a lockfile, renv.lock.
## use renv::restore() to reinstall all of the packages as declared in the lockfile.
set.seed(88)

# * 1 Functions to iron out phylogeny -----------------------------------
# the Upham 2019 trees aren't precisely ultrametric, so I used this function from Liam Revell: 
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

# loading activity data here to enable the following function
read.csv("ActivityData_MDD_v1.csv") -> data   # file made in "1-Updating taxonomy.r"


# automatically correct tip names to match taxonomy - requires the data file loaded (based on section 1.4 below)
CorTax <- function(tree){
    tree$tip.label[which(tree$tip.label == "Equus_africanus")]  <- "Equus_asinus"   # asinus is domestic E. africanus
    tree$tip.label[which(tree$tip.label == "Pygeretmus_zhitkovi")]  <- "Pygeretmus_shitkovi"    # correct spelling
    tree$tip.label[which(tree$tip.label == "Aethomys_namaquensis")] <- "Micaelamys_namaquensis" # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Tamiops_macclellandii")]<- "Tamiops_mcclellandii"   # typo
    tree$tip.label[which(tree$tip.label == "Nesotragus_moschatus")] <- "Neotragus_moschatus"    # typo
    tree$tip.label[which(tree$tip.label == "Pseudalopex_culpaeus")] <- "Lycalopex_culpaeus"     # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Pseudalopex_griseus")]  <- "Lycalopex_griseus"      # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Pseudalopex_gymnocercus")]  <- "Lycalopex_gymnocercus"  # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Pseudalopex_gymnocercus")]  <- "Lycalopex_sechurae" # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Pseudalopex_vetulus")]  <- "Lycalopex_vetulus"      # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Lutra_maculicollis")]   <- "Hydrictis_maculicollis" # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Monodelphis_unistriatus")]  <- "Monodelphis_unistriata" # Latin grammar correction
    tree$tip.label[which(tree$tip.label == "Dactylonax_palpator")]  <- "Dactylopsila_palpator"  # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Elephantulus_revoilii")]<- "Elephantulus_revoili"   # spelling error
    tree$tip.label[which(tree$tip.label == "Zaglossus_bruijnii")]   <- "Zaglossus_bruijni"      # spelling error
    tree$tip.label[which(tree$tip.label == "Procolobus_badius")]    <- "Piliocolobus_badius"    # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Procolobus_kirkii")]    <- "Piliocolobus_kirkii"    # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Procolobus_preussi")]   <- "Piliocolobus_preussi"   # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Procolobus_pennantii")] <- "Piliocolobus_pennantii" # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Procolobus_rufomitratus")]  <- "Piliocolobus_rufomitratus"   # genus transfer v1
    tree <- drop.tip(tree, which(!tree$tip.label %in% data$Phylo_name))
}

# * 2.1 Phylogenetic data (see also 1.4) ----------------------------------
phy1 <- read.tree("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree0000.tre")
tree <- ladderize(phy1)

# plot tree for inspection if needed
#pdf(file='tree.pdf', height=40, width=40)
#plot(phy1, type='fan', cex=0.15, label.offset = 0.05)
#dev.off()

## random sample of trees out of 10k variants from Upham et al. 2019
indices <- sprintf("%04d", sample(0:9999, 5000, replace=FALSE)) # sprintf makes all samples indices 4-digit long
# make MCC tree
for (i in 1:length(indices)){
    phyfile <- paste0("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree",indices[i],".tre")
    treevar <- read.tree(phyfile)
    treevar <- CorTax(treevar)
    #write.tree(treevar, file=paste0('treevar',indices[i],'.nex'))
    if (i == 1) {trees <- treevar} else {trees <- ape:::c.phylo(trees, treevar)} # group trees to one object
}
# MCC tree made of 5000 variants
clades <- prop.part(trees)
MCCtree <- maxCladeCred(trees, part = clades)

is.ultrametric(MCCtree)
UltMCCtree <- force.ultrametric(MCCtree)
is.ultrametric(UltMCCtree)

# * 2.2 Sample tree variants for analysis ---------------------------------
## variants sampled out of the indices that were used for the MCC tree, to ensure consistency

vars <- sample(1:5000, 50, replace=FALSE)
for (j in 1:length(vars)){
    phyfile <- paste0("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree",indices[vars[j]],".tre")
    treevar <- read.tree(phyfile)
    treevar <- CorTax(treevar)
    write.tree(treevar, file=paste0("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Data files/tree variants/treevar",indices[vars[j]],".nex"))
}

# * 3 Activity data  ------------------------------------------------------

# making sure names in data and in tree match
# checks below use the "%in%" syntax becasue setdiff() returns unique entries only, messing up the totals 
length(which(!data$MDD %in% tree$tip.label)) == 153     # checking names match - they don't
length(which(!data$MSW3 %in% tree$tip.label)) == 136    # somehow MSW3 has fewer mismatches


# * 4 Matching names in the data to phylogeny ---------------------------
# new columns for the names that are actually in the phylogeny, and where they're taken from
data$Phylo_name <- NA
data$Phylo_name_conforms_to <- NA

# * 4.1 Matching binomials in phylogeny to taxonomy ---------------------
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


# * 4.2 Correcting mismatches manually ----------------------------------
# (it's the safest solution)

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
#write.csv(data, file="ActivityData_MDD_v1_match.csv", row.names = FALSE)


# * 5 Final data processing steps ---------------------------------------
# read full data table (to save repeating the previous steps)
data <- read.csv(file="ActivityData_MDD_v1_match.csv")
act <- data[which(data$Phylo_name %in% tree$tip.label),]           

# remove crepusculars
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

# prepare data for use on cluster
write.table(act3[,c(6,5)], file="APdata_2400spp.txt", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)

# force ultrametric - placing this step after trimming the tree to save processing time
is.ultrametric(tree)
tree <- force.ultrametric(tree, "nnls")
is.ultrametric(tree)
#write.tree(tree, file='UltMCC.nex')


# * 6 Transforming to binary data (optional) ----------------------------
### THIS COULD ALSO BE DONE BY SETTING PARAMETER IDENTITIES IN THE MODELS
# * 6.1 Any daytime activity (N <-> CD dichotomy) ----------------------


# * 6.2 Only daytime activity (NC <-> D dichotomy ----------------------


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





# # 7 Plotting activity data on trees -------------------------------------
phy1 <- read.tree("UltMCC_2424spp_4AP_5kvars.nex")
tree <- ladderize(phy1)
 
label <- character(length(tree$tip.label))
names(label) <- act[which(act$Phylo_name %in% tree$tip.label),6]
label[act$AP=="Nocturnal"] <- "dodgerblue3"
label[act$AP=="Cathemeral"] <- "green3"
label[act$AP=="Diurnal"] <- "gold"
label[act$AP=="Crepuscular"] <- "magenta2"
label[which(label[] == '')] <- "white" 



# marking main orders
{artio <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Artiodactyla')])))
    chiro <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Chiroptera')])))
    carni <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Carnivora')])))
    prim <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Primates')])))
    roden <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order == 'Rodentia')])))
    eulip <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order %in% c('Erinaceomorpha','Soricomorpha'))])))
    xenar <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order %in% c('Cingulata','Pholidota','Pilosa'))])))
    afro <- getMRCA(tree,c(which(tree$tip.label %in% act$Phylo_name[which(act$Order %in% c('Afrosoricida','Hyracoidea','Macroscelidea','Proboscidea','Sirenia','Tubulidentata'))])))
}
orders <- c(afro, xenar, eulip, chiro, roden, artio, prim, carni)
labels <- c("Afrotheria","Xenarthra","Eulipotyphla","Chiroptera","Rodentia","Artiodactyla","Primates","Carnivora")


pdf(file="DataMCC.pdf",width=40, height=40)
plot(tree, type = "fan", lwd = 0.4, label.offset = 0.4, cex = 0.3, show.tip.label = TRUE, 
     tip.color = "white")#, edge.color = BemEDGE[,3])
ring(20, tree, style = "ring", offset = 5, col = label[tree$tip.label])
# adding arcs for orders
for(i in 1:4){#length(orders)) {
    arc.cladelabels(tree, labels[i], orders[i], cex=6,
                    MoreArgs = list(mark.node=FALSE, lab.offset=25, wing.length=0.5))
}
#ring(, tree, style = "ring", offset = -65.5, col = alpha("gray50", 0.25))
dev.off()


##OPTIONAL:
# colour main order brnaches
colours <- c("brown","cyan","green4","blue","darkorange2","goldenrod1","purple4","gray64")
for (k in 1:length(orders)) {
    AncNode <- orders[k]
    CladeSpp <- Descendants(tree, AncNode, type="all")
    edges <- BemEDGE[which(BemEDGE[,2] %in% CladeSpp),]
    BemEDGE[which(BemEDGE[,2] == AncNode),3] <- as.character(colours[k]) ## colouring the order root branch
    for (i in 1:length(edges[,2])) {
        edges[i,3] <- as.character(colours[k])
        BemEDGE[which(BemEDGE[,2] == edges[i,2]),3] <- edges[i,3]
    }
} 


### Order-wise analyses (Artiodactyla, Carnivora, Primates, Rodentia, Terrestrial-only[no bats, whales, sirenia)

## consider also: all terrestrial mammals, Euarchontoglires, Laurasiatheria, Afrotheria

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
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")
act <- read.csv(file="ActivityData_MDD_v1_match.csv")

# 1. Prune tree and trim dataset ----------------------------------------------
trees <- c("0028", "1166", "1845", "2966", "3865", "4496", "5221", "5333", "5455", "6404", "7198", "9128")
orders <-  c("Art","Car","Pri","Rod")

for (k in 1:length(trees)){
    phyfile <- paste0("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree",tree_no[i],".tre")
    treevar <- read.tree(phyfile)
    treevar <- CorTax(treevar)
    treevar <- drop.tip(treevar, which(!treevar$tip.label %in% act3$Phylo_name))
    write.tree(treevar, file=paste0('tree',tree_no[i],'.nex'))
}

for (ord in orders) {
    ordat <- act[which(act$Order == ord),c(6,5)]
    # first loop over all trees to extract the clade 
    for (k in 1:length(trees)) {
        # adjust tree
        tfile <- paste0("treevar",trees[k],".nex")
        tree <- read.tree(tfile) %>%
                CorTax() %>%
                force.ultrametric()
        # trim data and tree to match perfectly
        ordat <- ordat[which(ordat$Phylo_name %in% tree$tip.label),]
        LCA <- getMRCA(tree, ordat$Phylo_name)
        trim <- extract.clade(tree, LCA)
        write.tree(trim, file = paste0(ord, trees[k], ".nex"))
    }
    # write the data file before looping further
    write.csv(ordat, file = paste0(ord,"MatchedData.csv"))
    
    # then carry on to generate analysis files
    dir.create(file.path("./Analyses/ByOrder/",ord), recursive = TRUE)
    
    for (n in 1:5){
        dir.create(file.path(paste0("./Analyses/ByOrder/",ord,"/vrMuHiSSE",n+1)))
        dir.create(file.path(paste0("./Analyses/ByOrder/",ord,"/vrCID",n+1)))
        for (k in trees){
            # MuHiSSE Shell file
            f <- file(paste0("./Analyses/ByOrder/",ord,"/vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_tree", k,".sh"), open = "wb")
            sink(file=f)
                cat("#!/bin/bash -l", sep="\n")
                cat("", sep="\n")
                cat("# state program running time", sep="\n")
                cat("/usr/bin/time --verbose", sep="\n")
                cat("", sep="\n")
                cat("# Request wallclock time", sep="\n")
                cat("#$ -l h_rt=71:50:0", sep="\n")
                cat("", sep="\n")
                cat("# Request RAM", sep="\n")
                cat("#$ -l mem=15G", sep="\n")
                cat("", sep="\n")
                cat("# Request TMPDIR space", sep="\n")
                cat("#$ -l tmpfs=10G", sep="\n")
                cat("", sep="\n")
                cat("# Set the name of the job", sep="\n")
                # variable line
                cat(paste0("#$ -N ", ord,"_vrMuHiSSE",n+1,"_tree",k), sep="\n")
                cat("", sep="\n")
                cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
                cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
                cat("# NOTE: this directory must exist.", sep="\n")
                # variable line
                cat(paste0("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR/",ord), sep="\n")
                cat("", sep="\n")
                cat("# Load the R module and run your R program", sep="\n")
                cat("module unload compilers mpi", sep="\n")
                cat("module load r/recommended", sep="\n")
                cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS", sep="\n")
                cat("", sep="\n")
                cat("# Copy all necessary files into TMPDIR", sep="\n")
                # variable line
                cat(paste0("cp ",ord, "_MatchedData.csv $TMPDIR"), sep="\n")
                # variable line
                cat(paste0("cp ", ord, k,".nex $TMPDIR"), sep="\n")
                # variable line
                cat(paste0("cp ", ord, "_vrMuHiSSE",n+1,"_tree", k,"Code.R $TMPDIR"), sep="\n")
                cat("", sep="\n")
                cat("# Your work should be done in $TMPDIR", sep="\n")
                cat("cd $TMPDIR", sep="\n")
                cat("", sep="\n")
                # variable line
                cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR/", ord, "/", ord,"_vrMuHiSSE",n+1,"_tree", k,"Code.R> ", ord, "_vrMuHiSSE",n+1,"_tree", k,"Code.out"), sep="\n")
                cat("", sep="\n")
                cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
                # variable line
                cat(paste0("tar -zcvf $HOME/Scratch/MuHiSSE/VR/", ord, "/files_from_job_$JOB_ID.tgz $TMPDIR"), sep="\n")
                cat("", sep="\n")
                # variable line
                cat(paste0("cp ", ord, "_vrMuHiSSE",n+1,"_tree", k,".RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR/", ord), sep="\n")
            sink()
            close(f)
            
            # MuHiSSE R file  (the messy tab-alignment here makes the resulting file tidier)
            sink(paste0("./Analyses/VR/VR_MuHiSSE/vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.R"))
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools")\n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3523062, 0.4836342, 0, 0.4083013) \n")
            cat(paste0("tree <- read.tree(\"tree", tree_no[k],".nex\") \n\n")) # load one tree variant
            cat('# format data
for (i in 1:length(tree$tip.label)){
tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
} \n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state, tree$tip.state)
states_trans <- states
for(i in 1:Ntip(tree)){
if(states[i,2] == 1){
    states_trans[i,2] = 0
    states_trans[i,3] = 0
}
if(states[i,2] == 2){
    states_trans[i,2] = 0
    states_trans[i,3] = 1
}
if(states[i,2] == 3){
    states_trans[i,2] = 1
    states_trans[i,3] = 1
}
} \n\n')
            cat("# MuHiSSE \n") 
            cat(paste0("n <- ", n+1, "\n")) # set the number of hidden states
            cat("set_to_0 <- 4*(1:n)-1
TO_rates <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24)
drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64) \n")
            
            cat('# div rates 
turnover <- TO_rates[1:(4*n)]
extinction.fraction <- rep(1,4*n) 
extinction.fraction[set_to_0] <- 0 \n')
            
            cat('# trans rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = TRUE)
trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)])   
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
           turnover=turnover, 
           eps=extinction.fraction, 
           trans.rate=trans.rate.mod,
           hidden.states=TRUE, 
           turnover.upper=1000, 
           sann = TRUE) \n')
            
            cat(paste0("saveRDS(tmp, file=\"vrMuHiSSE",n+1,"_tree", tree_no[k],".RDS\")"))
            sink()
            
            
            # CID Shell file
            f <- file(paste0("./Analyses/VR/VR_CID/vrCID",n+1,"/vrCID",n+1,"_tree", tree_no[k],".sh"), open = "wb")
            sink(file=f)
            cat("#!/bin/bash -l", sep="\n")
            cat("", sep="\n")
            cat("# state program running time", sep="\n")
            cat("/usr/bin/time --verbose", sep="\n")
            cat("", sep="\n")
            cat("# Request wallclock time", sep="\n")
            cat("#$ -l h_rt=71:50:0", sep="\n")
            cat("", sep="\n")
            cat("# Request RAM", sep="\n")
            cat("#$ -l mem=31G", sep="\n")
            cat("", sep="\n")
            cat("# Request TMPDIR space", sep="\n")
            cat("#$ -l tmpfs=10G", sep="\n")
            cat("", sep="\n")
            cat("# Set the name of the job", sep="\n")
            # variable line
            cat(paste0("#$ -N vrCID",n+1,"_tree",tree_no[k]), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/CID/VR", sep="\n")
            cat("", sep="\n")
            cat("# Load the R module and run your R program", sep="\n")
            cat("module unload compilers mpi", sep="\n")
            cat("module load r/recommended", sep="\n")
            cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS", sep="\n")
            cat("", sep="\n")
            cat("# Copy all necessary files into TMPDIR", sep="\n")
            cat("cp APdata_2400spp.txt $TMPDIR", sep="\n")
            # variable line
            cat(paste0("cp tree", tree_no[k],".nex $TMPDIR"), sep="\n")
            # variable line
            cat(paste0("cp vrCID",n+1,"_tree", tree_no[k],"Code.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/CID/VR/vrCID",n+1,"_tree", tree_no[k],"Code.R> vrCID",n+1,"_tree", tree_no[k],"Code.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/CID/VR/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp vrCID",n+1,"_tree", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/CID/VR"), sep="\n")
            sink()
            close(f)
            
            # CID R file (the messy tab-alignment here makes the resulting file tidier)
            sink(paste0("./Analyses/VR/VR_CID/vrCID",n+1,"/vrCID",n+1,"_tree", tree_no[k],"Code.R"))
            cat('# load required packages 
Packages <- c("hisse", "diversitree", "phytools")
lapply(Packages, library, character.only = TRUE) \n\n')
            cat('# load data and tree
act3 <- read.table(file="APdata_2400spp.txt")
freq <- c(0.3523062, 0.4836342, 0, 0.4083013) \n')
            cat(paste0("tree <- read.tree(\"tree", tree_no[k],".nex\") \n\n")) # load one tree variant
            cat('# format data
for (i in 1:length(tree$tip.label)){
tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
} \n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state, tree$tip.state)
states_trans <- states
for(i in 1:Ntip(tree)){
if(states[i,2] == 1){
    states_trans[i,2] = 0
    states_trans[i,3] = 0
}
if(states[i,2] == 2){
    states_trans[i,2] = 0
    states_trans[i,3] = 1
}
if(states[i,2] == 3){
    states_trans[i,2] = 1
    states_trans[i,3] = 1
}
} \n\n')
            cat("# CID \n") 
            cat(paste0("n <- ", n+1, "\n")) # set the number of hidden states
            cat("set_to_0 <- 4*(1:n)-1 \n\n")
            cat("# div rates 
for (z in 1:n) {reps <- rep(z, 4)
if (z == 1) {turnover <- reps} else {turnover <- c(turnover, reps)}}
extinction.fraction <- rep(1,4*n) 
turnover[set_to_0] <- extinction.fraction[set_to_0] <- 0 \n\n")
            cat('# transition rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, make.null = TRUE, cat.trans.vary = TRUE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))
if (n == 5){
trans.rate.mod[which(trans.rate.mod == 23)] <- 24
trans.rate.mod[which(trans.rate.mod == 21)] <- 22
trans.rate.mod[which(trans.rate.mod == 20)] <- 21
trans.rate.mod[which(trans.rate.mod == 18)] <- 19
trans.rate.mod[which(trans.rate.mod == 17)] <- 18
trans.rate.mod[which(trans.rate.mod == 16)] <- 17
trans.rate.mod[which(trans.rate.mod == 15)] <- 16

trans.rate.mod[17,] <- c(23,0,0,0, 20,0,0,0, 15,0,0,0, 8,0,0,0, NA,2,0,0)
trans.rate.mod[18,] <- c(0,23,0,0, 0,20,0,0, 0,15,0,0, 0,8,0,0, 1,NA,0,4)
trans.rate.mod[19,] <- c(0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0, 0,0,NA,0)
trans.rate.mod[20,] <- c(0,0,0,23, 0,0,0,20, 0,0,0,15, 0,0,0,8, 0,3,0,NA)
}
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
           turnover=turnover, 
           eps=extinction.fraction, 
           trans.rate=trans.rate.mod,
           hidden.states=TRUE, 
           turnover.upper=1000, 
           sann = TRUE) \n\n')
            cat(paste0("saveRDS(tmp, file=\"vrCID",n+1,"_tree", tree_no[k],".RDS\")"))
            sink()
        }
    }
}



# 2. calculate sampling fractions ---------------------------------------------


# 3. set up relevant models ---------------------------------------------------


# 4. generate input and shell files for cluster -------------------------------




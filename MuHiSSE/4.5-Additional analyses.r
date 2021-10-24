library(phytools)
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")

# 1. Preliminary data ---------------------------------------------------------
# (this is a randomisation section and, to save resources, doesn't have to be re-run every time) 
set.seed(88)

# additional variants to increase sampling in MuHiSSE3 and MuHiSSE4 models
used <- read.csv(file="Indices_of_24Semi-finalists.csv")
used <- sprintf("%04d",used$x)
to_use <- sprintf("%04d", sample(1:9999, 100, replace=FALSE)) # choose 48 trees using 4-digit format
# loop to ensure new sample doesn't introduce duplicates
while(length(which(to_use %in% used)) > 0) {
    to_add <- length(which(to_use %in% used))                           
    addsamp <- sprintf("%04d", sample(1:9999, to_add, replace=FALSE))   # additional sample 
    to_use <- c(to_use, addsamp)                        
    to_use <- unique(to_use[-which(to_use %in% used)])                  # remove duplicates
}
tree_no <- to_use

# loading activity data here to enable the following function
data <- read.csv("ActivityData_MDD_v1_match.csv")   # file made in "2-Preparing for MuHiSSE.r"

# function to correct binomials
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

# get trees and remove crepuscular tips
creps <- data$Phylo_name[which(data$AP == 'Crepuscular')]           
for (i in 1:length(tree_no)){
    phyfile <- paste0("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree",tree_no[i],".tre")
    treevar <- read.tree(phyfile)
    treevar <- CorTax(treevar)
    treevar <- drop.tip(treevar, creps)
    write.tree(treevar, file=paste0('tree',tree_no[i],'.nex'))
}

# get all additional tree names
files <- list.files(pattern = "tree....\\.nex$") # list all tree****.nex files
used <- paste0("tree",used,".nex")
files <- files[!files %in% used]
fnames <- 1:100
for (i in 1:length(files)){
    fnames[i] <- strsplit(files[i], "\\.")[[1]][1]
}


dir.create(file.path("./Analyses/Additional"))   # directory for additional analyses

for (n in 1:4){
    for (k in 1:length(tree_no)){
        # MuHiSSE Shell file
        f <- file(paste0("./Analyses/Additional/vrMuHiSSE",n+1,"_tree", tree_no[k],".sh"), open = "wb")
        sink(file=f)
        cat("#!/bin/bash -l", sep="\n")
        cat("", sep="\n")
        cat("# Request wallclock time", sep="\n")
        cat("#$ -l h_rt=71:58:0", sep="\n")
        cat("", sep="\n")
        cat("# Request RAM", sep="\n")
        cat("#$ -l mem=15G", sep="\n")
        cat("", sep="\n")
        cat("# Request TMPDIR space", sep="\n")
        cat("#$ -l tmpfs=8G", sep="\n")
        cat("", sep="\n")
        cat("# Set the name of the job", sep="\n")
        # variable line
        cat(paste0("#$ -N vrMuHiSSE",n+1,"_tree",tree_no[k]), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/Additional", sep="\n")
        cat("", sep="\n")
        cat("# Load the R module and run your R program", sep="\n")
        cat("module -f unload compilers mpi gcc-libs", sep="\n")
        cat("module load beta-modules", sep="\n")
        cat("module load r/r-4.0.2_bc-3.11", sep="\n")
        cat("export R_LIBS=/home/ucbtmao/r/r-4.0.2_bc-3.11:$R_LIBS", sep="\n")
        cat("", sep="\n")
        cat("# Copy all necessary files into TMPDIR", sep="\n")
        cat("cp APdata_2400spp.txt $TMPDIR", sep="\n")
        # variable line
        cat(paste0("cp tree", tree_no[k],".nex $TMPDIR"), sep="\n")
        # variable line
        cat(paste0("cp vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        cat("# state program running time", sep="\n")
        # variable line
        cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/Additional/vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.R> vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        cat("tar -zcvf $HOME/Scratch/MuHiSSE/Additional/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp vrMuHiSSE",n+1,"_tree", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE/Additional"), sep="\n")
        sink()
        close(f)
        
        # MuHiSSE R file  (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/Additional/vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.R"))
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools") \n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3523062, 0.4836342, 0, 0.4083013) \n")
            cat(paste0("tree <- read.tree(\"tree", tree_no[k],".nex\") \n\n")) # load one tree variant
            cat('# format data \n')
            cat('for (i in 1:length(tree$tip.label)){ 
            tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
            } \n\n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state, tree$tip.state) \n')
            cat('states_trans <- states \n')
            cat('for(i in 1:Ntip(tree)){
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
            cat("set_to_0 <- 4*(1:n)-1 \n")
            cat("TO_rates <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24) \n")
            cat("drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64) \n\n")
            cat('# div rates \n')
            cat('turnover <- TO_rates[1:(4*n)] \n')
            cat('extinction.fraction <- rep(1,4*n) \n')
            cat('extinction.fraction[set_to_0] <- 0 \n\n')
            cat('# trans rates \n')
            cat('trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = TRUE) \n')
            cat('trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)]) \n')
            cat('trans.rate.mod[,set_to_0] <- 0 \n')
            cat('diag(trans.rate.mod) <- NA \n\n')
            cat('tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
                           turnover=turnover, 
                           eps=extinction.fraction, 
                           trans.rate=trans.rate.mod,
                           hidden.states=TRUE, 
                           turnover.upper=1000, 
                           sann = TRUE) \n\n')
            cat(paste0("saveRDS(tmp, file=\"vrMuHiSSE",n+1,"_tree", tree_no[k],".RDS\")"))
        sink()
        
        
        # CID Shell file
        f <- file(paste0("./Analyses/Additional/vrCID",n+1,"_tree", tree_no[k],".sh"), open = "wb")
        sink(file=f)
            cat("#!/bin/bash -l", sep="\n")
            cat("", sep="\n")
            cat("# Request wallclock time", sep="\n")
            cat("#$ -l h_rt=71:58:0", sep="\n")
            cat("", sep="\n")
            cat("# Request RAM", sep="\n")
            cat("#$ -l mem=15G", sep="\n")
            cat("", sep="\n")
            cat("# Request TMPDIR space", sep="\n")
            cat("#$ -l tmpfs=8G", sep="\n")
            cat("", sep="\n")
            cat("# Set the name of the job", sep="\n")
            # variable line
            cat(paste0("#$ -N vrCID",n+1,"_tree",tree_no[k]), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/CID/Additional", sep="\n")
            cat("", sep="\n")
            cat("# Load the R module and run your R program", sep="\n")
            cat("module -f unload compilers mpi gcc-libs", sep="\n")
            cat("module load beta-modules", sep="\n")
            cat("module load r/r-4.0.2_bc-3.11", sep="\n")
            cat("export R_LIBS=/home/ucbtmao/r/r-4.0.2_bc-3.11:$R_LIBS", sep="\n")
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
            cat("# state program running time", sep="\n")
            # variable line
            cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/CID/Additional/vrCID",n+1,"_tree", tree_no[k],"Code.R> vrCID",n+1,"_tree", tree_no[k],"Code.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/CID/Additional/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp vrCID",n+1,"_tree", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/CID/Additional"), sep="\n")
        sink()
        close(f)
        
        # CID R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/Additional/vrCID",n+1,"_tree", tree_no[k],"Code.R"))
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools") \n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3523062, 0.4836342, 0, 0.4083013) \n")
            cat(paste0("tree <- read.tree(\"tree", tree_no[k],".nex\") \n\n")) # load one tree variant
            cat('# format data \n')
            cat('for (i in 1:length(tree$tip.label)){
            tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
            } \n\n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state, tree$tip.state) \n')
            cat('states_trans <- states \n')
            cat('for(i in 1:Ntip(tree)){
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
            cat("# div rates \n")
            cat("for (z in 1:n) {reps <- rep(z, 4) \n")
            cat("if (z == 1) {turnover <- reps} else {turnover <- c(turnover, reps)}} \n")
            cat("extinction.fraction <- rep(1,4*n) \n")
            cat("turnover[set_to_0] <- extinction.fraction[set_to_0] <- 0 \n\n")
            cat("# transition rates \n")
            cat("trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, make.null = TRUE, cat.trans.vary = TRUE) \n")
            cat("trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8)) \n")
            cat("trans.rate.mod[,set_to_0] <- 0 \n")
            cat("diag(trans.rate.mod) <- NA \n\n")
            cat("tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
                           turnover=turnover, 
                           eps=extinction.fraction, 
                           trans.rate=trans.rate.mod,
                           hidden.states=TRUE, 
                           turnover.upper=1000, 
                           sann = TRUE) \n\n")
            cat(paste0("saveRDS(tmp, file=\"vrCID",n+1,"_tree", tree_no[k],".RDS\")"))
        sink()
    }
}

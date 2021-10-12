### Automatic shell file and R script generation for running models on the cluster
library(phytools)

# 1. Preliminary data ---------------------------------------------------------
# (this is a randomisation section and, to save resources, doesn't have to be re-run every time) 
set.seed(88)

#indices <- read.csv(file="Indices_of_5k_vars_used_for_MCC.csv")
#to_use <- sample(1:5000, 24, replace=FALSE) # choose 24 of the 5k trees used in the MCC
#tree_no <- sprintf("%04d",indices[to_use,1])
#write.csv(tree_no, file="Indices_of_24Semi-finalists.csv")
#tree_no <- tree_no[c(2,4,6,8,10,12,14,16,18,20,22,24)] # to halve number of repeats (if needed)


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


# 2. Make files for MCC NULL and MuSSE models ---------------------------------

dir.create(file.path("./Analyses/MuHiSSE/MuSSE"))
dir.create(file.path("./Analyses/CID/null"))


## 2.1 Generate files for the MCC tree analyses ------------------------------

## MuHiSSE Shell file
# this format is used to avoid the DOS (^M) end-of-line which trips the cluster. End-of-line DOS characters can be 
# removed from a bash file by using dos2nix filename.sh in the command line (use cat -v filename.sh to see if those 
# characters are there at all), but I want to automate it rather than repeating a manual command for every job submitted
{f <- file(paste0("./Analyses/MuHiSSE/MuSSE/MuSSE_MCCtree.sh"), open = "wb")
sink(file=f)
    cat("#!/bin/bash -l", sep="\n")
    cat("", sep="\n")
    cat("# state program running time", sep="\n")
    cat("/usr/bin/time --verbose", sep="\n")
    cat("", sep="\n")
    cat("# Request wallclock time", sep="\n")
    cat("#$ -l h_rt=50:50:0", sep="\n")
    cat("", sep="\n")
    cat("# Request RAM", sep="\n")
    cat("#$ -l mem=15G", sep="\n")
    cat("", sep="\n")
    cat("# Request TMPDIR space", sep="\n")
    cat("#$ -l tmpfs=10G", sep="\n")
    cat("", sep="\n")
    cat("# Set the name of the job", sep="\n")
    # variable line
    cat("#$ -N MuSSE_MCCtree", sep="\n")
    cat("", sep="\n")
    cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
    cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
    cat("# NOTE: this directory must exist.", sep="\n")
    # variable line
    cat("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/", sep="\n")
    cat("", sep="\n")
    cat("# Load the R module and run your R program", sep="\n")
    cat("module unload compilers mpi", sep="\n")
    cat("module load r/recommended", sep="\n")
    cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS", sep="\n")
    cat("", sep="\n")
    cat("# Copy all necessary files into TMPDIR", sep="\n")
    cat("cp APdata_2400spp.txt $TMPDIR", sep="\n")
    # variable line
    cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
    # variable line
    cat("cp MuSSE_MCCtreeCode.R $TMPDIR", sep="\n")
    cat("", sep="\n")
    cat("# Your work should be done in $TMPDIR", sep="\n")
    cat("cd $TMPDIR", sep="\n")
    cat("", sep="\n")
    # variable line
    cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/MuSSE_MCCtreeCode.R> MuSSE_MCCtreeCode.out"), sep="\n")
    cat("", sep="\n")
    cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
    cat("tar -zcvf $HOME/Scratch/MuHiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
    cat("", sep="\n")
    # variable line
    cat(paste0("cp MuSSE_MCCtree.RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE"), sep="\n")
sink()
close(f)
    
    # MuHiSSE R file (the messy tab-alignment here makes the resulting file tidier)
sink("./Analyses/MuHiSSE/MuSSE/MuSSE_MCCtreeCode.R")
    cat("# load required packages \n")
    cat('Packages <- c("hisse", "diversitree", "phytools")\n')
    cat("lapply(Packages, library, character.only = TRUE) \n\n")
    cat("# load data and tree \n")
    cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
    cat("freq <- c(0.3523062, 0.4836342, 0, 0.4083013) \n")
    cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
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
    cat("# MuSSE simple \n") 
    cat('# div rates 
turnover <- c(1,2,0,3)          
extinction.fraction <- c(1,1,0,1) \n')
    cat('# trans rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0, include.diagonals = FALSE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8)) 
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
           turnover=turnover, 
           eps=extinction.fraction, 
           trans.rate=trans.rate.mod,
           hidden.states=FALSE, 
           sann = TRUE) \n')
    
    cat(paste0("saveRDS(tmp, file=\"MuSSE_MCCtree.RDS\")"))
sink()
    
    
    # CID Shell file
f <- file(paste0("./Analyses/CID/null/null_MCCtree.sh"), open = "wb")
sink(file=f)
    cat("#!/bin/bash -l", sep="\n")
    cat("", sep="\n")
    cat("# state program running time", sep="\n")
    cat("/usr/bin/time --verbose", sep="\n")
    cat("", sep="\n")
    cat("# Request wallclock time", sep="\n")
    cat("#$ -l h_rt=50:50:0", sep="\n")
    cat("", sep="\n")
    cat("# Request RAM", sep="\n")
    cat("#$ -l mem=15G", sep="\n")
    cat("", sep="\n")
    cat("# Request TMPDIR space", sep="\n")
    cat("#$ -l tmpfs=10G", sep="\n")
    cat("", sep="\n")
    cat("# Set the name of the job", sep="\n")
    # variable line
    cat(paste0("#$ -N null_MCCtree"), sep="\n")
    cat("", sep="\n")
    cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
    cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
    cat("# NOTE: this directory must exist.", sep="\n")
    # variable line
    cat("#$ -wd /lustre/scratch/scratch/ucbtmao/CID/", sep="\n")
    cat("", sep="\n")
    cat("# Load the R module and run your R program", sep="\n")
    cat("module unload compilers mpi", sep="\n")
    cat("module load r/recommended", sep="\n")
    cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS", sep="\n")
    cat("", sep="\n")
    cat("# Copy all necessary files into TMPDIR", sep="\n")
    cat("cp APdata_2400spp.txt $TMPDIR", sep="\n")
    # variable line
    cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
    # variable line
    cat(paste0("cp null_MCCtreeCode.R $TMPDIR"), sep="\n")
    cat("", sep="\n")
    cat("# Your work should be done in $TMPDIR", sep="\n")
    cat("cd $TMPDIR", sep="\n")
    cat("", sep="\n")
    # variable line
    cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/CID/null_MCCtreeCode.R> null_MCCtreeCode.out"), sep="\n")
    cat("", sep="\n")
    cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
    cat("tar -zcvf $HOME/Scratch/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
    cat("", sep="\n")
    # variable line
    cat(paste0("cp null_MCCtree.RDS /lustre/scratch/scratch/ucbtmao/CID"), sep="\n")
sink()
close(f)
    
    # CID R file (the messy tab-alignment here makes the resulting file tidier)
sink(paste0("./Analyses/CID/null/null_MCCtreeCode.R"))
    cat('# load required packages 
Packages <- c("hisse", "diversitree", "phytools")
lapply(Packages, library, character.only = TRUE) \n\n')
    cat('# load data and tree
act3 <- read.table(file="APdata_2400spp.txt")
freq <- c(0.3523062, 0.4836342, 0, 0.4083013) # calculated in file: "3-MuHiSSE code.r" \n')
    cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n")) # load one tree variant
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
    cat("# dull null \n") 
    cat("# div rates 
turnover <- c(1,1,0,1)
extinction.fraction <- c(1,1,0,1) \n\n")
    cat('# transition rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
           turnover=turnover, 
           eps=extinction.fraction, 
           trans.rate=trans.rate.mod,
           hidden.states=FALSE, 
           sann = TRUE) \n\n')
    cat(paste0("saveRDS(tmp, file=\"null_MCCtree.RDS\")"))
sink()
}

## 2.2 Generate files for null and MuSSE on the tree variants ------------------
for (k in 1:length(tree_no)){
    # MuHiSSE Shell file
    f <- file(paste0("./Analyses/MuHiSSE/MuSSE/MuSSE_tree", tree_no[k],".sh"), open = "wb")
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
        cat(paste0("#$ -N MuSSE_tree",tree_no[k]), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/", sep="\n")
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
        cat(paste0("cp MuSSE_tree", tree_no[k],"Code.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/MuSSE_tree", tree_no[k],"Code.R> MuSSE_tree", tree_no[k],"Code.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        cat("tar -zcvf $HOME/Scratch/MuHiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp MuSSE_tree", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE"), sep="\n")
    sink()
    close(f)
        
    # MuHiSSE R file  (the messy tab-alignment here makes the resulting file tidier)
    sink(paste0("./Analyses/MuHiSSE/MuSSE/MuSSE_tree", tree_no[k],"Code.R"))
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
        cat("# MuSSE simple \n") 
        cat('# div rates 
turnover <- c(1,2,0,3)          
extinction.fraction <- c(1,1,0,1) \n')
        cat('# trans rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0, include.diagonals = FALSE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))  
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.mod,
               hidden.states=FALSE, 
               sann = TRUE) \n')
        
        cat(paste0("saveRDS(tmp, file=\"MuSSE_tree", tree_no[k],".RDS\")"))
    sink()
        
        
    # CID Shell file
    f <- file(paste0("./Analyses/CID/null/null_tree", tree_no[k],".sh"), open = "wb")
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
        cat(paste0("#$ -N null_tree",tree_no[k]), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/CID/", sep="\n")
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
        cat(paste0("cp null_tree", tree_no[k],"Code.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/CID/null_tree", tree_no[k],"Code.R> null_tree", tree_no[k],"Code.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        cat("tar -zcvf $HOME/Scratch/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp null_tree", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/CID"), sep="\n")
    sink()
    close(f)
        
    # CID R file (the messy tab-alignment here makes the resulting file tidier)
    sink(paste0("./Analyses/CID/null/null_tree", tree_no[k],"Code.R"))
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
        cat("# dull null \n") 
        cat("# div rates 
turnover <- c(1,1,0,1)
extinction.fraction <- c(1,1,0,1) \n\n")
        cat('# transition rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.mod,
               hidden.states=FALSE, 
               sann = TRUE) \n\n')
        cat(paste0("saveRDS(tmp, file=\"null_tree", tree_no[k],".RDS\")"))
    sink()
}


# 3. Make files for all hidden states models ---------------------------------
dir.create(file.path("./Analyses/MuHiSSE"), recursive = TRUE)   # create 'Analyses' directory as well
dir.create(file.path("./Analyses/CID"))

for (n in 1:7){
    dir.create(file.path(paste0("./Analyses/MuHiSSE/MuHiSSE",n+1)))
    dir.create(file.path(paste0("./Analyses/CID/CID",n+1)))
    
    ## 3.1 MCC tree analyses ---------------------------------------------------
    
    ## MuHiSSE Shell file
    # this format is used to avoid the DOS (^M) end-of-line which trips the cluster. End-of-line DOS characters can be 
    # removed from a bash file by using dos2nix filename.sh in the command line (use cat -v filename.sh to see if those 
    # characters are there at all), but I want to automate it rather than repeating a manual command for every job submitted
    {f <- file(paste0("./Analyses/MuHiSSE/MuHiSSE",n+1,"/MuHiSSE",n+1,"_MCCtree.sh"), open = "wb")
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
            cat(paste0("#$ -N MuHiSSE",n+1,"_MCCtree_31MB"), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/", sep="\n")
            cat("", sep="\n")
            cat("# Load the R module and run your R program", sep="\n")
            cat("module unload compilers mpi", sep="\n")
            cat("module load r/recommended", sep="\n")
            cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS", sep="\n")
            cat("", sep="\n")
            cat("# Copy all necessary files into TMPDIR", sep="\n")
            cat("cp APdata_2400spp.txt $TMPDIR", sep="\n")
            # variable line
            cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
            # variable line
            cat(paste0("cp MuHiSSE",n+1,"_MCCtreeCode.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/MuHiSSE",n+1,"_MCCtreeCode.R> MuHiSSE",n+1,"_MCCtreeCode.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/MuHiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp MuHiSSE",n+1,"_MCCtree.RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE"), sep="\n")
        sink()
        close(f)
        
        # MuHiSSE R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/MuHiSSE/MuHiSSE",n+1,"/MuHiSSE",n+1,"_MCCtreeCode.R"))
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools")\n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3523062, 0.4836342, 0, 0.4083013) \n")
            cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
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
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = FALSE) # simplifying assumption: transitions between hidden states are equal 
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
        
            cat(paste0("saveRDS(tmp, file=\"MuHiSSE",n+1,"_MCCtree.RDS\")"))
        sink()
        
        
        # CID Shell file
        f <- file(paste0("./Analyses/CID/CID",n+1,"/CID",n+1,"_MCCtree.sh"), open = "wb")
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
            cat(paste0("#$ -N CID",n+1,"_MCCtree_31MB"), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/CID/", sep="\n")
            cat("", sep="\n")
            cat("# Load the R module and run your R program", sep="\n")
            cat("module unload compilers mpi", sep="\n")
            cat("module load r/recommended", sep="\n")
            cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS", sep="\n")
            cat("", sep="\n")
            cat("# Copy all necessary files into TMPDIR", sep="\n")
            cat("cp APdata_2400spp.txt $TMPDIR", sep="\n")
            # variable line
            cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
            # variable line
            cat(paste0("cp CID",n+1,"_MCCtreeCode.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/CID/CID",n+1,"_MCCtreeCode.R> CID",n+1,"_MCCtreeCode.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp CID",n+1,"_MCCtree.RDS /lustre/scratch/scratch/ucbtmao/CID"), sep="\n")
        sink()
        close(f)
        
        # CID R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/CID/CID",n+1,"/CID",n+1,"_MCCtreeCode.R"))
            cat('# load required packages 
Packages <- c("hisse", "diversitree", "phytools")
lapply(Packages, library, character.only = TRUE) \n\n')
            cat('# load data and tree
act3 <- read.table(file="APdata_2400spp.txt")
freq <- c(0.3523062, 0.4836342, 0, 0.4083013) # calculated in file: "3-MuHiSSE code.r" \n')
            cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n")) # load one tree variant
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
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, make.null = TRUE, cat.trans.vary = FALSE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))   
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.mod,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE) \n\n')
            cat(paste0("saveRDS(tmp, file=\"CID",n+1,"_MCCtree.RDS\")"))
        sink()
    }
    

    ## 3.2 Generate files for the tree variants -------------------------------
    
    for (k in 1:length(tree_no)){
        
        # MuHiSSE Shell file
        f <- file(paste0("./Analyses/MuHiSSE/MuHiSSE",n+1,"/MuHiSSE",n+1,"_tree", tree_no[k],".sh"), open = "wb")
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
            cat(paste0("#$ -N MuHiSSE",n+1,"_tree",tree_no[k]), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/", sep="\n")
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
            cat(paste0("cp MuHiSSE",n+1,"_tree", tree_no[k],"Code.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/MuHiSSE",n+1,"_tree", tree_no[k],"Code.R> MuHiSSE",n+1,"_tree", tree_no[k],"Code.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/MuHiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp MuHiSSE",n+1,"_tree", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE"), sep="\n")
        sink()
        close(f)
        
        # MuHiSSE R file  (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/MuHiSSE/MuHiSSE",n+1,"/MuHiSSE",n+1,"_tree", tree_no[k],"Code.R"))
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
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = FALSE) # simplifying assumption: transitions between hidden states are equal 
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
        
        cat(paste0("saveRDS(tmp, file=\"MuHiSSE",n+1,"_tree", tree_no[k],".RDS\")"))
        sink()
        
        
        # CID Shell file
        f <- file(paste0("./Analyses/CID/CID",n+1,"/CID",n+1,"_tree", tree_no[k],".sh"), open = "wb")
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
            cat(paste0("#$ -N CID",n+1,"_tree",tree_no[k]), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/CID/", sep="\n")
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
            cat(paste0("cp CID",n+1,"_tree", tree_no[k],"Code.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/CID/CID",n+1,"_tree", tree_no[k],"Code.R> CID",n+1,"_tree", tree_no[k],"Code.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp CID",n+1,"_tree", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/CID"), sep="\n")
        sink()
        close(f)
        
        # CID R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/CID/CID",n+1,"/CID",n+1,"_tree", tree_no[k],"Code.R"))
        cat('# load required packages 
Packages <- c("hisse", "diversitree", "phytools")
lapply(Packages, library, character.only = TRUE) \n\n')
        cat('# load data and tree
act3 <- read.table(file="APdata_2400spp.txt")
freq <- c(0.3523062, 0.4836342, 0, 0.4083013) # calculated in file: "3-MuHiSSE code.r" \n')
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
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, make.null = TRUE, cat.trans.vary = FALSE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))   
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.mod,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE) \n\n')
        cat(paste0("saveRDS(tmp, file=\"CID",n+1,"_tree", tree_no[k],".RDS\")"))
        sink()
    }
}  


# 4. Make VR files (between-states transition rates vary) ---------------------
dir.create(file.path("./Analyses/VR/VR_MuHiSSE"), recursive = TRUE)   # create 'Analyses' directory as well
dir.create(file.path("./Analyses/VR/VR_CID"))

for (n in 1:4){
    #dir.create(file.path(paste0("./Analyses/VR/VR_MuHiSSE/vrMuHiSSE",n+1)))
    #dir.create(file.path(paste0("./Analyses/VR/VR_CID/vrCID",n+1)))
    

    ## 4.1 Files for MCC analyses ---------------------------------------------

    ## MuHiSSE Shell file
    # this format is used to avoid the DOS (^M) end-of-line which trips the cluster. End-of-line DOS characters can be 
    # removed from a bash file by using dos2nix filename.sh in the command line (use cat -v filename.sh to see if those 
    # characters are there at all), but I want to automate it rather than repeating a manual command for every job submitted
    {f <- file(paste0("./Analyses/VR/VR_MuHiSSE/vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_MCCtree.sh"), open = "wb")
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
        cat(paste0("#$ -N vrMuHiSSE",n+1,"_MCCtree"), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR", sep="\n")
        cat("", sep="\n")
        cat("# Load the R module and run your R program", sep="\n")
        cat("module unload compilers mpi", sep="\n")
        cat("module load r/recommended", sep="\n")
        cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS", sep="\n")
        cat("", sep="\n")
        cat("# Copy all necessary files into TMPDIR", sep="\n")
        cat("cp APdata_2400spp.txt $TMPDIR", sep="\n")
        # variable line
        cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
        # variable line
        cat(paste0("cp vrMuHiSSE",n+1,"_MCCtreeCode.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR/vrMuHiSSE",n+1,"_MCCtreeCode.R> vrMuHiSSE",n+1,"_MCCtreeCode.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        cat("tar -zcvf $HOME/Scratch/MuHiSSE/VR/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp vrMuHiSSE",n+1,"_MCCtree.RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR"), sep="\n")
        sink()
        close(f)
        
        # MuHiSSE R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/VR/VR_MuHiSSE/vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_MCCtreeCode.R"))
        cat("# load required packages \n")
        cat('Packages <- c("hisse", "diversitree", "phytools")\n')
        cat("lapply(Packages, library, character.only = TRUE) \n\n")
        cat("# load data and tree \n")
        cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
        cat("freq <- c(0.3523062, 0.4836342, 0, 0.4083013) \n")
        cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
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
        
        cat(paste0("saveRDS(tmp, file=\"vrMuHiSSE",n+1,"_MCCtree.RDS\")"))
        sink()
        
        
        # CID Shell file
        f <- file(paste0("./Analyses/VR/VR_CID/vrCID",n+1,"/vrCID",n+1,"_MCCtree.sh"), open = "wb")
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
        cat(paste0("#$ -N vrCID",n+1,"_MCCtree"), sep="\n")
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
        cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
        # variable line
        cat(paste0("cp vrCID",n+1,"_MCCtreeCode.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/CID/VR/vrCID",n+1,"_MCCtreeCode.R> vrCID",n+1,"_MCCtreeCode.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        cat("tar -zcvf $HOME/Scratch/CID/VR/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp vrCID",n+1,"_MCCtree.RDS /lustre/scratch/scratch/ucbtmao/CID/VR"), sep="\n")
        sink()
        close(f)
        
        # CID R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/VR/VR_CID/vrCID",n+1,"/vrCID",n+1,"_MCCtreeCode.R"))
        cat('# load required packages 
Packages <- c("hisse", "diversitree", "phytools")
lapply(Packages, library, character.only = TRUE) \n\n')
        cat('# load data and tree
act3 <- read.table(file="APdata_2400spp.txt")
freq <- c(0.3523062, 0.4836342, 0, 0.4083013) \n')
        cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n")) # load one tree variant
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
#if (n == 5){ 
#    trans.rate.mod[which(trans.rate.mod == 23)] <- 24
#    trans.rate.mod[which(trans.rate.mod == 21)] <- 22
#    trans.rate.mod[which(trans.rate.mod == 20)] <- 21
#    trans.rate.mod[which(trans.rate.mod == 18)] <- 19
#    trans.rate.mod[which(trans.rate.mod == 17)] <- 18
#    trans.rate.mod[which(trans.rate.mod == 16)] <- 17
#    trans.rate.mod[which(trans.rate.mod == 15)] <- 16
#    
#    trans.rate.mod[17,] <- c(23,0,0,0, 20,0,0,0, 15,0,0,0, 8,0,0,0, NA,2,0,0)
#    trans.rate.mod[18,] <- c(0,23,0,0, 0,20,0,0, 0,15,0,0, 0,8,0,0, 1,NA,0,4)
#    trans.rate.mod[19,] <- c(0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0, 0,0,NA,0)
#    trans.rate.mod[20,] <- c(0,0,0,23, 0,0,0,20, 0,0,0,15, 0,0,0,8, 0,3,0,NA)
#}
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.mod,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE) \n\n')
        cat(paste0("saveRDS(tmp, file=\"vrCID",n+1,"_MCCtree.RDS\")"))
        sink()
    }

    ## 4.2 Generate files for 12 tree variants --------------------------------

    for (k in 1:length(tree_no)){
        # MuHiSSE Shell file
        f <- file(paste0("./Analyses/VR/VR_MuHiSSE/vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_tree", tree_no[k],".sh"), open = "wb")
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
            cat(paste0("#$ -N vrMuHiSSE",n+1,"_tree",tree_no[k]), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR", sep="\n")
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
            cat(paste0("cp vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR/vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.R> vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/MuHiSSE/VR/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp vrMuHiSSE",n+1,"_tree", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE/VR"), sep="\n")
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
#if (n == 5){
#    trans.rate.mod[which(trans.rate.mod == 23)] <- 24
#    trans.rate.mod[which(trans.rate.mod == 21)] <- 22
#    trans.rate.mod[which(trans.rate.mod == 20)] <- 21
#    trans.rate.mod[which(trans.rate.mod == 18)] <- 19
#    trans.rate.mod[which(trans.rate.mod == 17)] <- 18
#    trans.rate.mod[which(trans.rate.mod == 16)] <- 17
#    trans.rate.mod[which(trans.rate.mod == 15)] <- 16
#    
#    trans.rate.mod[17,] <- c(23,0,0,0, 20,0,0,0, 15,0,0,0, 8,0,0,0, NA,2,0,0)
#    trans.rate.mod[18,] <- c(0,23,0,0, 0,20,0,0, 0,15,0,0, 0,8,0,0, 1,NA,0,4)
#    trans.rate.mod[19,] <- c(0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0, 0,0,NA,0)
#    trans.rate.mod[20,] <- c(0,0,0,23, 0,0,0,20, 0,0,0,15, 0,0,0,8, 0,3,0,NA)
#}
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



# 5. Creating files for additional analyses -----------------------------------

}

##################################################################################################################
###  R template to replicate for running on the cluster  ####
#########################################################

# load required packages 
Packages <- c("hisse", "diversitree", "phytools")
lapply(Packages, library, character.only = TRUE)

# load data and tree
act3 <- read.table(file="APdata_2400spp.txt")
freq <- c(0.3523062, 0.4836342, 0, 0.4083013) # calculated in file: "3-MuHiSSE code.r" 
tree <- read.tree("UltMCC_2400spp_3AP_5kvars.nex")

# format data
for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
}
states <- data.frame(tree$tip.label, tree$tip.state, tree$tip.state)
states_trans <- states
for(i in 1:Ntip(tree)){
    if(states[i,2] == 1){
        states_trans[i,2] = 0
        states_trans[i,3] = 0
    }
    if(states[i,2] == 2){
        states_trans[i,2] = 0   # change to '1' to lump cath with diur (CD dataset)
        states_trans[i,3] = 1   # change to '0' to lump cath with noct (NC dataset)
    }
    if(states[i,2] == 3){
        states_trans[i,2] = 1
        states_trans[i,3] = 1
    }
}

###############
##  MuHiSSE  ##
###############
n <- 4                  
set_to_0 <- 4*(1:n)-1
TO_rates <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24)
drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)

# div rates 
turnover <- TO_rates[1:(4*n)]
extinction.fraction <- rep(1,4*n) 
extinction.fraction[set_to_0] <- 0

# trans rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, include.diagonals = FALSE, cat.trans.vary = FALSE)
trans.rate.mod <- ParDrop(trans.rate, drops[1:(4*n)])   
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.mod,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE)

saveRDS(tmp, file="MuHiSSE-8.RDS")

##### ##### 
##  CID  ##
##### #####

n <- 3                  
set_to_0 <- 4*(1:n)-1

# div rates 
for (z in 1:n) {reps <- rep(z, 4)
if (z == 1) {turnover <- reps} else {turnover <- c(turnover, reps)}}
extinction.fraction <- rep(1,4*n) 
turnover[set_to_0] <- extinction.fraction[set_to_0] <- 0

# transition rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, make.null = TRUE, cat.trans.vary = FALSE)
trans.rate.mod <- ParDrop(trans.rate, c(2,5,6,8))   
trans.rate.mod[,set_to_0] <- 0
diag(trans.rate.mod) <- NA
tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.mod,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE)

saveRDS(tmp, file="CID8.RDS")



###  E N D  ###
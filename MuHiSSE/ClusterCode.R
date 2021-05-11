



### Automating shell file generation
#  1. Preliminary data (this is a randomisation section and, to save resources, doesn't have to be re-run every time) 

#indices <- read.csv(file="Indices_of_5k_vars_used_for_MCC.csv")
#to_use <- sample(1:5000, 24, replace=FALSE) # choose 24 of the 5k trees used in the MCC
#tree_no <- sprintf("%04d",indices[to_use,1])
#tree_no <- tree_no[c(2,4,6,8,10,12,14,16,18,20,22,24)] # halving repeats because life is too short
#
## get trees and remove crepuscular tips
#for (i in 1:length(tree_no)){
#    phyfile <- paste0("C:/Users/Roi Maor/Desktop/New Mam Phylo/Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_tree",tree_no[i],".tre")
#    treevar <- read.tree(phyfile)
#    treevar <- CorTax(treevar)
#    treevar <- drop.tip(treevar, which(!treevar$tip.label %in% act3$Phylo_name))
#    write.tree(treevar, file=paste0('tree',tree_no[i],'.nex'))
#}

#  2. Make files
dir.create(file.path("./Analyses/MuHiSSE"), recursive = TRUE)   # create 'Analyses' directory as well
dir.create(file.path("./Analyses/CID"))

for (n in 1:7){
    dir.create(file.path(paste0("./Analyses/MuHiSSE/MuHiSSE",n+1)))
    dir.create(file.path(paste0("./Analyses/CID/CID",n+1)))
    
    ## 2.1 Generate files for the MCC tree analyses
    
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
    
    ## 2.2 Generate files for the 24 tree variants
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

##################################################################################################################
###  R code to run on the cluster  ####
#######################################

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
        states_trans[i,2] = 0
        states_trans[i,3] = 1
    }
    if(states[i,2] == 3){
        states_trans[i,2] = 1
        states_trans[i,3] = 1
    }
}

###############
##  MuHiSSE  ##
###############
n <- 8                  
set_to_0 <- 4*(1:n)-1
TO_rates <- c(1,2,0,3, 4,5,0,6, 7,8,0,9, 10,11,0,12, 13,14,0,15, 16,17,0,18, 19,20,0,21, 22,23,0,24)
drops <- c(2,5,6,8, 10,13,14,16, 18,21,22,24, 26,29,30,32, 34,37,38,40, 42,45,46,48, 50,53,54,56, 58,61,62,64)

# div rates 
turnover <- TO_rates[1:(4*n)]
extinction.fraction <- rep(1,4*n) 
extinction.fraction[set_to_0] <- 0

# trans rates
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
               sann = TRUE)

saveRDS(tmp, file="MuHiSSE-8.RDS")

##### ##### 
##  CID  ##
##### #####

n <- 8                  
set_to_0 <- 4*(1:n)-1

# div rates 
for (z in 1:n) {reps <- rep(z, 4)
if (z == 1) {turnover <- reps} else {turnover <- c(turnover, reps)}}
extinction.fraction <- rep(1,4*n) 
turnover[set_to_0] <- extinction.fraction[set_to_0] <- 0

# transition rates
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=n-1, make.null = TRUE, cat.trans.vary = TRUE)
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



#####################################################################################


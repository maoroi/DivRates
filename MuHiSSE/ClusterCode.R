
# load required packages 
Packages <- c("hisse", "diversitree", "phytools")
lapply(Packages, library, character.only = TRUE)

# load data and tree
act3 <- read.table(file="APdata_2400spp.txt")
f <- c(0.3523062, 0.4836342, 0, 0.4083013) # calculated in file: "3-MuHiSSE code.r" 
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

# MuHiSSE
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
tmp <- MuHiSSE(phy=tree, data=states_trans, f=f, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.mod,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE)

saveRDS(tmp, file="MuHiSSE-8.RDS")
#####################################################################################

### automating shell file generation
indices <- read.csv(file="Indices_of_5k_vars_used_for_MCC.csv")
to_use <- sample(1:5000, 24, replace=FALSE) # choose 24 of the 5k trees used in the MCC
tree_no <- sprintf("%04d",indices[to_use,1])

dir.create(file.path("./Analyses/MuHiSSE"), recursive = TRUE)
dir.create(file.path("./Analyses/CID"))

for (n in 1:7){
    dir.create(file.path(paste0("./Analyses/MuHiSSE/MuHiSSE-",n+1)))
    dir.create(file.path(paste0("./Analyses/CID/CID-",n+1)))
    for (k in 1:length(tree_no)){
        # MuHiSSE Shell file
        sink(paste0("./Analyses/MuHiSSE/MuHiSSE-",n+1,"/MuHiSSE-",n+1,"_treevar", tree_no[k],".sh"))
            cat("#!/bin/bash -l \n")
            cat(" \n")
            cat("# state program running time \n")
            cat("/usr/bin/time --verbose \n")
            cat(" \n")
            cat("# Request wallclock time \n")
            cat("#$ -l h_rt=71:30:0 \n")
            cat(" \n")
            cat("# Request RAM \n")
            cat("#$ -l mem=125G \n")
            cat(" \n")
            cat("# Request TMPDIR space \n")
            cat("#$ -l tmpfs=10G \n")
            cat(" \n")
            cat("# Set the name of the job \n")
            cat(paste0("#$ -N MuHiSSE",n+1,"_tree",tree_no[k],"_125MB \n"))
            cat(" \n")
            cat("# Set the working directory to somewhere in your scratch space. \n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME. \n")
            cat("# NOTE: this directory must exist. \n")
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/MuHiSSE/ \n")
            cat(" \n")
            cat("# Load the R module and run your R program \n")
            cat("module unload compilers mpi \n")
            cat("module load r/recommended \n")
            cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS \n")
            cat(" \n")
            cat("# Copy all necessary files into TMPDIR \n")
            cat("cp APdata_2400spp.txt $TMPDIR \n")
            cat(paste0("cp treevar", tree_no[k], ".nex $TMPDIR \n"))
            cat("cp MuHiSSECode.R $TMPDIR \n")
            cat(" \n")
            cat("# Your work should be done in $TMPDIR \n")
            cat("cd $TMPDIR \n")
            cat(" \n")
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/MuHiSSE/MuHiSSE",n+1,"Code.R> MuHiSSE",n+1,"Code.out \n"))
            cat(" \n")
            cat("# Preferably, tar-up (archive) all output files to transfer them back \n")
            cat("# to your space. This will include the R_output file above. \n")
            cat("tar -zcvf $HOME/Scratch/MuHiSSE/files_from_job_$JOB_ID.tgz $TMPDIR \n")
            cat(" \n")
            cat(paste0("cp MuHiSSE-",n+1,"_treevar", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/MuHiSSE \n"))
        sink()
        
        # CID Shell file
        sink(paste0("./Analyses/CID/CID-",n+1,"/CID-",n+1,"_treevar", tree_no[k],".sh"))
            cat("#!/bin/bash -l \n")
            cat(" \n")
            cat("# state program running time \n")
            cat("/usr/bin/time --verbose \n")
            cat(" \n")
            cat("# Request wallclock time \n")
            cat("#$ -l h_rt=71:30:0 \n")
            cat(" \n")
            cat("# Request RAM \n")
            cat("#$ -l mem=125G \n")
            cat(" \n")
            cat("# Request TMPDIR space \n")
            cat("#$ -l tmpfs=10G \n")
            cat(" \n")
            cat("# Set the name of the job \n")
            cat(paste0("#$ -N CID",n+1,"_tree",tree_no[k],"_125MB \n"))
            cat(" \n")
            cat("# Set the working directory to somewhere in your scratch space. \n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME. \n")
            cat("# NOTE: this directory must exist. \n")
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/CID/ \n")
            cat(" \n")
            cat("# Load the R module and run your R program \n")
            cat("module unload compilers mpi \n")
            cat("module load r/recommended \n")
            cat("export R_LIBS=/home/ucbtmao/R/R-3.6.0-OpenBLAS:$R_LIBS \n")
            cat(" \n")
            cat("# Copy all necessary files into TMPDIR \n")
            cat("cp APdata_2400spp.txt $TMPDIR \n")
            cat(paste0("cp treevar", tree_no[k], ".nex $TMPDIR \n"))
            cat(paste0("cp CID",n+1,"Code.R $TMPDIR \n"))
            cat(" \n")
            cat("# Your work should be done in $TMPDIR \n")
            cat("cd $TMPDIR \n")
            cat(" \n")
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/CID/CID",n+1,"Code.R> CID",n+1,"Code.out \n"))
            cat(" \n")
            cat("# tar-up (archive) all output files to transfer them back to your space. \n")
            cat("tar -zcvf $HOME/Scratch/CID/files_from_job_$JOB_ID.tgz $TMPDIR \n")
            cat(" \n")
            cat(paste0("cp CID-",n+1,"_treevar", tree_no[k],".RDS /lustre/scratch/scratch/ucbtmao/CID \n"))
        sink()
    }
}  
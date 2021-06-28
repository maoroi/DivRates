### Binary analyses

# 1. Prepare data for binary analyses -----------------------------------------
data <- read.csv('ActivityData_MDD_v1_match.csv')
# total count in the entire dataset (inc species not in the phylogeny)
TotN <- length(which(data$AP == 'Nocturnal'))
TotCD <- length(which(data$AP %in% c('Cathemeral','Diurnal')))
# proportions in data (sum up to 99%, crepuscular and others make the rest)
datN <- TotN / nrow(data)
datCD <- TotCD / nrow(data)

## 1.1 Assumptions ------------------------------------------------------------
tax <- read.csv("C:/Users/Roi Maor/Desktop/New Mam Phylo/Burgin taxonomy/MDD_v1_6495species_JMamm.csv", header = TRUE)

# unsampled bats presumed nocturnal (n = 1096 in v1; 1122 in v1.1) 
bats_unsamp <- length(which(tax$Order == 'CHIROPTERA')) - length(which(data$Order == 'Chiroptera')) - 6  # unsampled bats minus 6 extinct spp

# unsampled squirrels (27 flying squirrels presumed noct; 87 other sciurids presumed diur)
SQRLS <- tax[which(tax$Family == 'SCIURIDAE'),]                         # all sciurids (all extant)
SQRLS_DATA <- data[which(data$Family == 'Sciuridae'),]                  # sampled sciurids 
flying <- length(which(SQRLS$Tribe == 'PTEROMYINI'))                    # 57 flying squirrels in MDD
noct_sqrls <- length(which(SQRLS_DATA$AP == "Nocturnal"))               # flying squirels are the only nocturnal sciurids
FLY_SQRLS_unsamp <- flying - noct_sqrls
Diur_sqrls_unsamp <- nrow(SQRLS) - nrow(SQRLS_DATA) - FLY_SQRLS_unsamp  # unsampled scuirids minus unsampled noct sciurids

# Simiiformes presumed diurnal (n = 148, excludes Aotus) 
PRIMS <- tax[which(tax$Order == "PRIMATES"),]
PRIMS <- PRIMS[which(PRIMS$extinct. == 0),]                             # only extant species
simians <- c("Atelidae", "Cebidae", "Cercopithecidae", "Pitheciidae", "Hominidae", "Hylobatidae")
simian_count <- length(which(PRIMS$Family %in% toupper(simians)))       # all simians recognised in MDD
SIMIAN_DATA <- data[which(data$Family %in% c(simians, "Aotidae")),]     # Aotidae demoted into Cebidae (Aotinae) in MDD but data$Family relates to MSW3   
# Out of 11 species of Aotus in MDD, 7 have AP info (no assumptions about remaining 4 spp)
Aotus_unsamp <- length(which(tax$Subfamily == "AOTINAE")) - length(which(data$Family == 'Aotidae'))
simi_diur_unsamp <- simian_count - nrow(SIMIAN_DATA) - Aotus_unsamp          # unsmapled (diur) simians minus unsampled (noct) Aotinae

# total unsampled (Noct=1123, Diur=235)
Nunsamp <- bats_unsamp + FLY_SQRLS_unsamp
Dunsamp <- Diur_sqrls_unsamp + simi_diur_unsamp
rm(bats_unsamp, FLY_SQRLS_unsamp, Diur_sqrls_unsamp, simi_diur_unsamp, PRIMS, simians, 
   simian_count, SIMIAN_DATA, Aotus_unsamp, noct_sqrls, flying, SQRLS_DATA, SQRLS)
##      ** end assumptions **       ##

## 1.2 Estimate sampling fraction ---------------------------------------------
# species with no data or assumptions
all_extant <- length(unique(tax$SciName[which(tax$extinct. == 0)]))
unkn <- all_extant - (nrow(data) + Nunsamp + Dunsamp)                   # total number of extant species based on MDD v1
# expected AP distribution in unknown species if same proportion as in the data
expN <- datN * unkn
expCD <- datCD * unkn
# proportion sampled out of likely total
Nprop <- TotN / (TotN + Nunsamp + expN)
CDprop <- TotCD / (TotCD + Dunsamp + expCD)
#> c(Nprop, CDprop)
#[1] 0.3523062 0.4322946
freq <- c(Nprop, CDprop)


#### This code generates shell files and R scripts for running MuHiSSE and MuCID models
# on _binary_ AP data (Cath lumped with Noct or Diur)

dir.create(file.path("./Analyses/Binary/HiSSE"), recursive = TRUE)   # create 'Analyses' directory as well
dir.create(file.path("./Analyses/Binary/CID"))


# 2. MuSSE and CID null models ------------------------------------------------

for (n in 1:7) {
    #dir.create(file.path(paste0("./Analyses/Binary/HiSSE/bvrMuHiSSE",n+1)))
    #dir.create(file.path(paste0("./Analyses/Binary/CID/bvrCID",n+1)))
    
# 1.1 Generate null and MuSSE files for MCC tree ----------------------------- 
    {f <- file(paste0("./Analyses/Binary/HiSSE/BiSSE_MCCtree.sh"), open = "wb")
        sink(file=f)
            cat("#!/bin/bash -l", sep="\n")
            cat("", sep="\n")
            cat("# Request wallclock time", sep="\n")
            cat("#$ -l h_rt=55:50:0", sep="\n")
            cat("", sep="\n")
            cat("# Request RAM", sep="\n")
            cat("#$ -l mem=15G", sep="\n")
            cat("", sep="\n")
            cat("# Request TMPDIR space", sep="\n")
            cat("#$ -l tmpfs=10G", sep="\n")
            cat("", sep="\n")
            cat("# Set the name of the job", sep="\n")
            # variable line
            cat("#$ -N BiSSE_MCCtree", sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/HiSSE/", sep="\n")
            cat("", sep="\n")
            cat("# Load the R module and run your R program", sep="\n")
            cat("module -f unload compilers mpi gcc-libs", sep="\n")
            cat("module load beta-modules", sep="\n")
            cat("module load r/r-4.0.2_bc-3.11", sep="\n")
            cat("export R_LIBS=/home/ucbtmao/r/r-4.0.2_bc-3.11:$R_LIBS", sep="\n")
            cat("", sep="\n")
            cat("# Copy all necessary files into TMPDIR", sep="\n")
            cat("cp MDD_v1_6495species_JMamm.csv $TMPDIR", sep="\n")
            cat("cp ActivityData_MDD_v1_match.csv $TMPDIR", sep="\n")
            # variable line
            cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
            # variable line
            cat("cp BiSSE_MCCtreeCode.R $TMPDIR", sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            cat("# state program running time", sep="\n")
            # variable line
            cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/HiSSE/BiSSE_MCCtreeCode.R> BiSSE_MCCtreeCode.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/MuHiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp BiSSE_MCCtree.RDS /lustre/scratch/scratch/ucbtmao/Binary/HiSSE"), sep="\n")
        sink()
        close(f)
        
        # MuHiSSE R file (the messy tab-alignment here makes the resulting file tidier)
        sink("./Analyses/Binary/HiSSE/BiSSE_MCCtreeCode.R")
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools")\n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3523062, 0.4322946) \n")
            cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
            cat('# format data \n')
            cat('for (i in 1:length(tree$tip.label)){ \n')
            cat('tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])} \n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state) \n')
            cat('states_trans <- states \n')
            cat('states_trans[,2] <- states_trans[,2] - 1 \n')
            cat('states_trans[which(states_trans[,2] > 0),2] <- 1 # lump C + D together \n\n')
            cat("# BiSSE \n") 
            cat('# div rates \n')
            cat('turnover <- c(1,2) \n')      
            cat('extinction.fraction <- c(1,1) \n')
            cat('# trans rates \n')
            cat('trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits = 0) \n')

            cat('tmp <- hisse(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rates.bisse,
               hidden.states=FALSE, 
               sann = TRUE) \n')
            
            cat(paste0("saveRDS(tmp, file=\"BiSSE_MCCtree.RDS\")"))
        sink()
        
        
        # CID Shell file
        f <- file(paste0("./Analyses/Binary/CID/nullbin_MCCtree.sh"), open = "wb")
        sink(file=f)
            cat("#!/bin/bash -l", sep="\n")
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
            cat(paste0("#$ -N nullbin_MCCtree"), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/CID/", sep="\n")
            cat("", sep="\n")
            cat("module -f unload compilers mpi gcc-libs", sep="\n")
            cat("module load beta-modules", sep="\n")
            cat("module load r/r-4.0.2_bc-3.11", sep="\n")
            cat("export R_LIBS=/home/ucbtmao/r/r-4.0.2_bc-3.11:$R_LIBS", sep="\n")
            cat("", sep="\n")
            cat("# Copy all necessary files into TMPDIR", sep="\n")
            cat("cp APdata_2400spp.txt $TMPDIR", sep="\n")
            # variable line
            cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
            # variable line
            cat(paste0("cp nullbin_MCCtreeCode.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            cat("# state program running time", sep="\n")
            # variable line
            cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/CID/nullbin_MCCtreeCode.R> nullbin_MCCtreeCode.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp nullbin_MCCtree.RDS /lustre/scratch/scratch/ucbtmao/Binary/CID"), sep="\n")
        sink()
        close(f)
        
        # CID R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/Binary/CID/nullbin_MCCtreeCode.R"))
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools")\n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3523062, 0.4322946) \n")
            cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
            cat('# format data \n')
            cat('for (i in 1:length(tree$tip.label)){ \n')
            cat('tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])} \n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state) \n')
            cat('states_trans <- states \n')
            cat('states_trans[,2] <- states_trans[,2] - 1 \n')
            cat('states_trans[which(states_trans[,2] > 0),2] <- 1 # lump C + D together \n\n')
            cat('## dull null ## \n') 
            cat('# div rates \n')
            cat('turnover <- c(1,1) \n')
            cat('extinction.fraction <- c(1,1) \n\n')
            cat('# transition rates \n')
            cat('trans.rate <- TransMatMakerHiSSE(hidden.traits = 0) \n')
            cat('tmp <- MuHiSSE(phy=tree, data=states_trans, f=freq, 
                    turnover=turnover, eps=extinction.fraction, 
                    trans.rate=trans.rate, hidden.states=FALSE, 
                    sann = TRUE) \n\n')
            cat(paste0("saveRDS(tmp, file=\"nullbin_MCCtree.RDS\")"))
        sink()
    }


# 1.2 Generate files for null and MuSSE on the tree variants  ----------------

for (k in 1:length(tree_no)){
    # MuHiSSE Shell file
    f <- file(paste0("./Analyses/.../MuSSE_tree", tree_no[k],".sh"), open = "wb")
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
    sink(paste0("./Analyses/.../MuSSE_tree", tree_no[k],"Code.R"))
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
    f <- file(paste0("./Analyses/.../null_tree", tree_no[k],".sh"), open = "wb")
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
    sink(paste0("./Analyses/.../null_tree", tree_no[k],"Code.R"))
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



# 2. Make VR files (between-states transition rates vary) ---------------------
dir.create(file.path("./Analyses/VR/VR_MuHiSSE"), recursive = TRUE)   # create 'Analyses' directory as well
dir.create(file.path("./Analyses/VR/VR_CID"))

for (n in 1:4){
    #dir.create(file.path(paste0("./Analyses/VR/VR_MuHiSSE/vrMuHiSSE",n+1)))
    #dir.create(file.path(paste0("./Analyses/VR/VR_CID/vrCID",n+1)))
    
    
    # 2.1 Files for MCC analyses -------------------------------------------------
    
    ## MuHiSSE Shell file
    # this format is used to avoid the DOS (^M) end-of-line which trips the cluster. End-of-line DOS characters can be 
    # removed from a bash file by using dos2nix filename.sh in the command line (use cat -v filename.sh to see if those 
    # characters are there at all), but I want to automate it rather than repeating a manual command for every job submitted
    {f <- file(paste0("./Analyses/.../vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_MCCtree.sh"), open = "wb")
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
    sink(paste0("./Analyses/.../vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_MCCtreeCode.R"))
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
    f <- file(paste0("./Analyses/.../vrCID",n+1,"/vrCID",n+1,"_MCCtree.sh"), open = "wb")
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
    sink(paste0("./Analyses/.../vrCID",n+1,"/vrCID",n+1,"_MCCtreeCode.R"))
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
    
    # 2.2 Generate files for 12 tree variants ------------------------------------
    
    for (k in 1:length(tree_no)){
        # MuHiSSE Shell file
        f <- file(paste0("./Analyses/.../vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_tree", tree_no[k],".sh"), open = "wb")
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
        sink(paste0("./Analyses/.../vrMuHiSSE",n+1,"/vrMuHiSSE",n+1,"_tree", tree_no[k],"Code.R"))
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
        f <- file(paste0("./Analyses/.../vrCID",n+1,"/vrCID",n+1,"_tree", tree_no[k],".sh"), open = "wb")
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
        sink(paste0("./Analyses/.../vrCID",n+1,"/vrCID",n+1,"_tree", tree_no[k],"Code.R"))
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


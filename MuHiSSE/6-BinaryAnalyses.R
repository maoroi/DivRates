### Binary analyses
library(hisse)

setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")

# 1. Prepare data for binary analyses -----------------------------------------
data <- read.csv('ActivityData_MDD_v1_match.csv')

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
###    ** end assumptions **    ###


## 1.2 *** GROUPING N-CD *** ---------------------------------------------------
# total count in the entire dataset (inc species not in the phylogeny)
#TotN <- length(which(data$AP == 'Nocturnal'))
#TotCD <- length(which(data$AP %in% c('Cathemeral','Diurnal')))
#datN <- TotN / nrow(data)
#datCD <- TotCD / nrow(data)

## 1.3 *** ALTERNATIVE GROUPING: D-CN ------------------------------------------
# total count in the entire dataset (inc species not in the phylogeny)
TotCN <- length(which(data$AP %in% c('Cathemeral','Nocturnal')))
TotD <- length(which(data$AP == 'Diurnal'))
datCN <- TotCN / nrow(data)
datD <- TotD / nrow(data)

## 1.4 Estimate sampling fraction ---------------------------------------------
# species with no data or assumptions
all_extant <- length(unique(tax$SciName[which(tax$extinct. == 0)]))
unkn <- all_extant - (nrow(data) + Nunsamp + Dunsamp)                   # total number of extant species based on MDD v1
# expected AP distribution in unknown species if same proportion as in the data
expCN <- datCN * unkn
expD <- datD * unkn
# proportion sampled out of likely total
CNprop <- TotCN / (TotCN + Nunsamp + expCN)
Dprop <- TotD / (TotD + Dunsamp + expD)

#> c(Nprop, CDprop)
#[1] 0.3523062 0.4322946
freq <- c(Nprop, CDprop)

#> c(CNprop, Dprop)
#[1] 0.3714350 0.4083013
freq <- c(CNprop, Dprop)

# 2. BiSSE and CID null models ------------------------------------------------

#### This code generates shell files and R scripts for running HiSSE and MuCID models
# on _binary_ AP data (Cath lumped with Noct or Diur)

#dir.create(file.path("./Analyses/Binary/D-CD/HiSSE"), recursive = TRUE)   # create 'Analyses' directory as well
#dir.create(file.path("./Analyses/Binary/N-CD/CID"))

dir.create(file.path("./Analyses/Binary/D-CN/HiSSE"), recursive = TRUE)   # create 'Analyses' directory as well
dir.create(file.path("./Analyses/Binary/D-CN/CID"))

for (n in 1:7) {
    
    ## 1.1 Generate null and BiSSE files for MCC tree ------------------------- 
    {f <- file(paste0("./Analyses/Binary/D-CN/HiSSE/BiSSE_MCCtree_D-CN.sh"), open = "wb")
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
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE/", sep="\n")
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
        cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
        # variable line
        cat("cp BiSSE_MCCtreeCode_D-CN.R $TMPDIR", sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        cat("# state program running time", sep="\n")
        # variable line
        cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE/BiSSE_MCCtreeCode_D-CN.R> BiSSE_MCCtreeCode.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        # variable line
        cat("tar -zcvf $HOME/Scratch/Binary/D-CN/HiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp BiSSE_MCCtree_D-CN.RDS /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE"), sep="\n")
    sink()
    close(f)
    
    # HiSSE R file (the messy tab-alignment here makes the resulting file tidier)
    sink("./Analyses/Binary/D-CN/HiSSE/BiSSE_MCCtreeCode_D-CN.R")
        cat("# load required packages \n")
        cat('Packages <- c("hisse", "diversitree", "phytools")\n')
        cat("lapply(Packages, library, character.only = TRUE) \n\n")
        cat("# load data and tree \n")
        cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
        cat("freq <- c(0.3714350, 0.4083013) \n") # c(0.3714350, 0.4083013) is correct for D-CN
        cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
        cat('# format data \n')
        cat('for (i in 1:length(tree$tip.label)){ \n')
        cat('tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])} \n')
        cat('states <- data.frame(tree$tip.label, tree$tip.state) \n')
        cat('states_trans <- states \n')
        cat('states_trans[which(states_trans[,2] < 3),2] <- 0 \n')
        cat('states_trans[which(states_trans[,2] == 3),2] <- 1 \n\n')
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
        
        cat(paste0("saveRDS(tmp, file=\"BiSSE_MCCtree_D-CN.RDS\")"))
    sink()
}    
        
        # CID Shell file
        f <- file(paste0("./Analyses/Binary/D-CN/CID/nullbin_MCCtree_D-CN.sh"), open = "wb")
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
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID/", sep="\n")
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
            cat(paste0("cp nullbin_MCCtreeCode_D-CN.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            cat("# state program running time", sep="\n")
            # variable line
            cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID/nullbin_MCCtreeCode_D-CN.R> nullbin_MCCtreeCode.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            # variable line
            cat("tar -zcvf $HOME/Scratch/Binary/D-CN/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp nullbin_MCCtree_D-CN.RDS /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID"), sep="\n")
        sink()
        close(f)
        
        # CID R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/Binary/D-CN/CID/nullbin_MCCtreeCode_D-CN.R"))
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools")\n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3714350, 0.4083013) \n") # c(0.3714350, 0.4083013) is correct for D-CN
            cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
            cat('# format data \n')
            cat('for (i in 1:length(tree$tip.label)){ \n')
            cat('tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])} \n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state) \n')
            cat('states_trans <- states \n')
            cat('states_trans[which(states_trans[,2] < 3),2] <- 0 \n')
            cat('states_trans[which(states_trans[,2] == 3),2] <- 1 \n\n')
            cat('## dull null ## \n') 
            cat('# div rates \n')
            cat('turnover <- c(1,1) \n')
            cat('extinction.fraction <- c(1,1) \n\n')
            cat('# transition rates \n')
            cat('trans.rate <- TransMatMakerHiSSE(hidden.traits = 0) \n')
            cat('tmp <- hisse(phy=tree, data=states_trans, f=freq, 
                turnover=turnover, eps=extinction.fraction, 
                trans.rate=trans.rate, hidden.states=FALSE, 
                sann = TRUE) \n\n')
            cat(paste0("saveRDS(tmp, file=\"nullbin_MCCtree_D-CN.RDS\")"))
        sink()
}


    ## 1.2 Generate files for null and BiSSE on the tree variants  ------------
tree_no <- read.csv(file="Indices_of_24Semi-finalists.csv")
tree_no <- sprintf("%04d",tree_no[,2])

for (k in 1:length(tree_no)){
    # HiSSE Shell file
    f <- file(paste0("./Analyses/Binary/D-CN/HiSSE/BiSSE_tree", tree_no[k],"_D-CN.sh"), open = "wb")
    sink(file=f)
        cat("#!/bin/bash -l", sep="\n")
        cat("", sep="\n")
        cat("# Request wallclock time", sep="\n")
        cat("#$ -l h_rt=25:50:0", sep="\n")
        cat("", sep="\n")
        cat("# Request RAM", sep="\n")
        cat("#$ -l mem=15G", sep="\n")
        cat("", sep="\n")
        cat("# Request TMPDIR space", sep="\n")
        cat("#$ -l tmpfs=10G", sep="\n")
        cat("", sep="\n")
        cat("# Set the name of the job", sep="\n")
        # variable line
        cat(paste0("#$ -N BiSSE_tree",tree_no[k]), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE/", sep="\n")
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
        cat(paste0("cp BiSSE_tree", tree_no[k],"Code_D-CN.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        cat("# state program running time", sep="\n")
        # variable line
        cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE/BiSSE_tree", tree_no[k],"Code_D-CN.R> BiSSE_tree", tree_no[k],"Code.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        # variable line
        cat("tar -zcvf $HOME/Scratch/Binary/D-CN/HiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp BiSSE_tree", tree_no[k],"_D-CN.RDS /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE"), sep="\n")
    sink()
    close(f)
    
    # HiSSE R file  (the messy tab-alignment here makes the resulting file tidier)
    sink(paste0("./Analyses/Binary/D-CN/HiSSE/BiSSE_tree", tree_no[k],"Code_D-CN.R"))
        cat("# load required packages \n")
        cat('Packages <- c("hisse", "diversitree", "phytools")\n')
        cat("lapply(Packages, library, character.only = TRUE) \n\n")
        cat("# load data and tree \n")
        cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
        cat("freq <- c(0.3714350, 0.4083013) \n") # c(0.3714350, 0.4083013) is correct for D-CN
        cat(paste0("tree <- read.tree(\"tree", tree_no[k],".nex\") \n\n")) # load one tree variant
        cat('# format data
for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
} \n')
        cat('states <- data.frame(tree$tip.label, tree$tip.state)
states_trans <- states
states_trans[which(states_trans[,2] < 3),2] <- 0
states_trans[which(states_trans[,2] == 3),2] <- 1 \n\n')
    cat("# BiSSE simple \n") 
    cat('# div rates 
turnover <- c(1,2)          
extinction.fraction <- c(1,1) \n')
    cat('# trans rates
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=0)
tmp <- hisse(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.hisse,
               hidden.states=FALSE, 
               sann = TRUE) \n')
    
    cat(paste0("saveRDS(tmp, file=\"BiSSE_tree", tree_no[k],"_D-CN.RDS\")"))
    sink()
    
   
    # CID Shell file
    f <- file(paste0("./Analyses/Binary/D-CN/CID/nullbin_tree", tree_no[k],"_D-CN.sh"), open = "wb")
    sink(file=f)
        cat("#!/bin/bash -l", sep="\n")
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
        cat(paste0("#$ -N nullbin_tree",tree_no[k]), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID/", sep="\n")
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
        cat(paste0("cp nullbin_tree", tree_no[k],"Code_D-CN.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        cat("# state program running time", sep="\n")
        # variable line
        cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID/nullbin_tree", tree_no[k],"Code_D-CN.R> nullbin_tree", tree_no[k],"Code.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        # variable line
        cat("tar -zcvf $HOME/Scratch/Binary/D-CN/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp nullbin_tree", tree_no[k],"_D-CN.RDS /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID"), sep="\n")
    sink()
    close(f)
    
    # CID R file (the messy tab-alignment here makes the resulting file tidier)
    sink(paste0("./Analyses/Binary/D-CN/CID/nullbin_tree", tree_no[k],"Code_D-CN.R"))
        cat("# load required packages \n")
        cat('Packages <- c("hisse", "diversitree", "phytools")\n')
        cat("lapply(Packages, library, character.only = TRUE) \n\n")
        cat("# load data and tree \n")
        cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
        cat("freq <- c(0.3714350, 0.4083013) \n") # c(0.3714350, 0.4083013) is correct for D-CN
        cat(paste0("tree <- read.tree(\"tree", tree_no[k],".nex\") \n\n")) # load one tree variant
        cat('# format data
for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
} \n')
        cat('states <- data.frame(tree$tip.label, tree$tip.state)
states_trans <- states
states_trans[which(states_trans[,2] < 3),2] <- 0
states_trans[which(states_trans[,2] == 3),2] <- 1 \n\n')
        cat("# dull null \n") 
        cat("# div rates 
turnover <- c(1,1)
extinction.fraction <- c(1,1) \n\n")
        cat('# transition rates
trans.rate.cid <- TransMatMakerHiSSE(hidden.traits=0)
tmp <- hisse(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.cid,
               hidden.states=FALSE, 
               sann = TRUE) \n\n')
        cat(paste0("saveRDS(tmp, file=\"nullbin_tree", tree_no[k],"_D-CN.RDS\")"))
    sink()
} # optional: close looping over tree numbers (to yield BiSSE and null models only)

# 3. Make only VR files (between-states transition rates vary) ----------------

for (n in 1:5){
    
    #dir.create(file.path(paste0("./Analyses/Binary/HiSSE/vrHiSSE",n+1)))
    #dir.create(file.path(paste0("./Analyses/Binary/CID/bvrCID",n+1)))
    
    ## 3.1 Files for MCC analyses ---------------------------------------------
    
    ## HiSSE Shell file
    # this format is used to avoid the DOS (^M) end-of-line which trips the cluster. End-of-line DOS characters can be 
    # removed from a bash file by using dos2nix filename.sh in the command line (use cat -v filename.sh to see if those 
    # characters are there at all), but I want to automate it rather than repeating a manual command for every job submitted
    {f <- file(paste0("./Analyses/Binary/D-CN/HiSSE/vrHiSSE",n+1,"_MCCtree_D-CN.sh"), open = "wb")
    sink(file=f)
        cat("#!/bin/bash -l", sep="\n")
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
        cat(paste0("#$ -N vrHiSSE",n+1,"_MCCtree"), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE", sep="\n")
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
        cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
        # variable line
        cat(paste0("cp vrHiSSE",n+1,"_MCCtreeCode_D-CN.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE/vrHiSSE",n+1,"_MCCtreeCode_D-CN.R> vrHiSSE",n+1,"_MCCtreeCode.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        cat("tar -zcvf $HOME/Scratch/Binary/D-CN/HiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp vrHiSSE",n+1,"_MCCtree_D-CN.RDS /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE"), sep="\n")
    sink()
    close(f)
    
    # HiSSE R file (the messy tab-alignment here makes the resulting file tidier)
    sink(paste0("./Analyses/Binary/D-CN/HiSSE/vrHiSSE",n+1,"_MCCtreeCode_D-CN.R"))
        cat("# load required packages \n")
        cat('Packages <- c("hisse", "diversitree", "phytools")\n')
        cat("lapply(Packages, library, character.only = TRUE) \n\n")
        cat("# load data and tree \n")
        cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
        cat("freq <- c(0.3714350, 0.4083013) \n") # c(0.3714350, 0.4083013) is correct for D-CN
        cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
        cat('# format data
for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
} \n')
        cat('states <- data.frame(tree$tip.label, tree$tip.state)
states_trans <- states
states_trans[which(states_trans[,2] < 3),2] <- 0
states_trans[which(states_trans[,2] == 3),2] <- 1 \n\n')
        cat("# HiSSE \n") 
        cat(paste0("n <- ", n+1, "\n")) # set the number of hidden states
        cat('# div rates 
turnover <- 1:(2*n)
extinction.fraction <- rep(1,2*n) \n\n')
        cat('# trans rates
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits = n-1, cat.trans.vary = TRUE)
tmp <- hisse(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.hisse,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE) \n')
    
        cat(paste0("saveRDS(tmp, file=\"vrHiSSE",n+1,"_MCCtree_D-CN.RDS\")"))
    sink()
    }
    
    # CID Shell file
    f <- file(paste0("./Analyses/Binary/D-CN/CID/bvrCID",n+1,"_MCCtree_D-CN.sh"), open = "wb")
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
        cat(paste0("#$ -N bvrCID",n+1,"_MCCtree"), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID", sep="\n")
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
        cat("cp UltMCC_2400spp_3AP_5kvars.nex $TMPDIR", sep="\n")
        # variable line
        cat(paste0("cp bvrCID",n+1,"_MCCtreeCode_D-CN.R $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID/bvrCID",n+1,"_MCCtreeCode_D-CN.R> bvrCID",n+1,"_MCCtreeCode.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
        # variable line
        cat("tar -zcvf $HOME/Scratch/Binary/D-CN/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp bvrCID",n+1,"_MCCtree_D-CN.RDS /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID"), sep="\n")
    sink()
    close(f)
    
    # CID R file (the messy tab-alignment here makes the resulting file tidier)
    sink(paste0("./Analyses/Binary/D-CN/CID/bvrCID",n+1,"_MCCtreeCode_D-CN.R"))
        cat("# load required packages \n")
        cat('Packages <- c("hisse", "diversitree", "phytools")\n')
        cat("lapply(Packages, library, character.only = TRUE) \n\n")
        cat("# load data and tree \n")
        cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
        cat("freq <- c(0.3714350, 0.4083013) \n") # c(0.3714350, 0.4083013) is correct for D-CN
        cat(paste0("tree <- read.tree(\"UltMCC_2400spp_3AP_5kvars.nex\") \n\n"))
        cat('# format data
for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
} \n')
        cat('states <- data.frame(tree$tip.label, tree$tip.state)
states_trans <- states
states_trans[which(states_trans[,2] < 3),2] <- 0
states_trans[which(states_trans[,2] == 3),2] <- 1 \n\n')
        cat("# bvrCID \n") 
        cat(paste0("n <- ", n+1, "\n")) # set the number of hidden states
        cat("# div rates 
for (z in 1:n) {reps <- rep(z, 2)
if (z == 1) {turnover <- reps} else {turnover <- c(turnover, reps)}}
extinction.fraction <- rep(1, 2*n) \n\n")
        cat('# transition rates
trans.rate.cid <- TransMatMakerHiSSE(hidden.traits = n-1, make.null = TRUE, cat.trans.vary = TRUE)
tmp <- hisse(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.cid,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE) \n\n')
        cat(paste0("saveRDS(tmp, file=\"bvrCID",n+1,"_MCCtree_D-CN.RDS\")"))
    sink()

    
    # 3.2 Generate files for tree variants ------------------------------------
    
    for (k in 1:length(tree_no)){
        # HiSSE Shell file
        f <- file(paste0("./Analyses/Binary/D-CN/HiSSE/vrHiSSE",n+1,"_tree", tree_no[k],"_D-CN.sh"), open = "wb")
        sink(file=f)
            cat("#!/bin/bash -l", sep="\n")
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
            cat(paste0("#$ -N vrHiSSE",n+1,"_tree",tree_no[k]), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE", sep="\n")
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
            cat(paste0("cp vrHiSSE",n+1,"_tree", tree_no[k],"Code_D-CN.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Your work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE/vrHiSSE",n+1,"_tree", tree_no[k],"Code_D-CN.R> vrHiSSE",n+1,"_tree", tree_no[k],"Code.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/Binary/D-CN/HiSSE/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp vrHiSSE",n+1,"_tree", tree_no[k],"_D-CN.RDS /lustre/scratch/scratch/ucbtmao/Binary/D-CN/HiSSE"), sep="\n")
        sink()
        close(f)
        
        # HiSSE R file  (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/Binary/D-CN/HiSSE/vrHiSSE",n+1,"_tree", tree_no[k],"Code_D-CN.R"))
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools")\n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3714350, 0.4083013) \n") # c(0.3714350, 0.4083013) is correct for D-CN
            cat(paste0("tree <- read.tree(\"tree", tree_no[k],".nex\") \n\n")) # load one tree variant
            cat('# format data
for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
} \n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state)
states_trans <- states
states_trans[which(states_trans[,2] < 3),2] <- 0
states_trans[which(states_trans[,2] == 3),2] <- 1 \n\n')
            cat("# vrHiSSE \n") 
            cat(paste0("n <- ", n+1, "\n")) # set the number of hidden states
            cat('# div rates 
turnover <- 1:(2*n)
extinction.fraction <- rep(1,2*n) \n')
        
            cat('# trans rates
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits = n-1, cat.trans.vary = TRUE)
tmp <- hisse(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.hisse,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE) \n')
        
            cat(paste0("saveRDS(tmp, file=\"vrHiSSE",n+1,"_tree", tree_no[k],"_D-CN.RDS\")"))
        sink()
        
        
        # CID Shell file
        f <- file(paste0("./Analyses/Binary/D-CN/CID/bvrCID",n+1,"_tree", tree_no[k],"_D-CN.sh"), open = "wb")
        sink(file=f)
            cat("#!/bin/bash -l", sep="\n")
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
            cat(paste0("#$ -N bvrCID",n+1,"_tree",tree_no[k]), sep="\n")
            cat("", sep="\n")
            cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
            cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
            cat("# NOTE: this directory must exist.", sep="\n")
            # variable line
            cat("#$ -wd /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID", sep="\n")
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
            cat(paste0("cp bvrCID",n+1,"_tree", tree_no[k],"Code_D-CN.R $TMPDIR"), sep="\n")
            cat("", sep="\n")
            cat("# Work should be done in $TMPDIR", sep="\n")
            cat("cd $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("R --no-save < /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID/bvrCID",n+1,"_tree", tree_no[k],"Code_D-CN.R> bvrCID",n+1,"_tree", tree_no[k],"Code.out"), sep="\n")
            cat("", sep="\n")
            cat("# tar-up (archive) all output files to transfer them back to your space.", sep="\n")
            cat("tar -zcvf $HOME/Scratch/Binary/D-CN/CID/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
            cat("", sep="\n")
            # variable line
            cat(paste0("cp bvrCID",n+1,"_tree", tree_no[k],"_D-CN.RDS /lustre/scratch/scratch/ucbtmao/Binary/D-CN/CID"), sep="\n")
        sink()
        close(f)
        
        # CID R file (the messy tab-alignment here makes the resulting file tidier)
        sink(paste0("./Analyses/Binary/D-CN/CID/bvrCID",n+1,"_tree", tree_no[k],"Code_D-CN.R"))
            cat("# load required packages \n")
            cat('Packages <- c("hisse", "diversitree", "phytools")\n')
            cat("lapply(Packages, library, character.only = TRUE) \n\n")
            cat("# load data and tree \n")
            cat("act3 <- read.table(file=\"APdata_2400spp.txt\") \n")
            cat("freq <- c(0.3714350, 0.4083013) \n") # c(0.3714350, 0.4083013) is correct for D-CN
            cat(paste0("tree <- read.tree(\"tree", tree_no[k],".nex\") \n\n")) # load one tree variant
            cat('# format data
for (i in 1:length(tree$tip.label)){
    tree$tip.state[i] <- as.numeric(act3$AP[which(act3$Phylo_name == tree$tip.label[i])])
} \n')
            cat('states <- data.frame(tree$tip.label, tree$tip.state)
states_trans <- states
states_trans[which(states_trans[,2] < 3),2] <- 0
states_trans[which(states_trans[,2] == 3),2] <- 1 \n\n')
            cat("# bvrCID \n") 
            cat(paste0("n <- ", n+1, "\n")) # set the number of hidden states
            cat("# div rates 
for (z in 1:n) {reps <- rep(z, 2)
if (z == 1) {turnover <- reps} else {turnover <- c(turnover, reps)}}
extinction.fraction <- rep(1, 2*n) \n\n")
            cat('# transition rates
trans.rate.cid <- TransMatMakerHiSSE(hidden.traits = n-1, make.null = TRUE, cat.trans.vary = TRUE)
tmp <- hisse(phy=tree, data=states_trans, f=freq, 
               turnover=turnover, 
               eps=extinction.fraction, 
               trans.rate=trans.rate.cid,
               hidden.states=TRUE, 
               turnover.upper=1000, 
               sann = TRUE) \n\n')
            cat(paste0("saveRDS(tmp, file=\"bvrCID",n+1,"_tree", tree_no[k],"_D-CN.RDS\")"))
        sink()
    }
} 


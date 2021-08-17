setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")

# files for estimating parameter support region
dir.create(file.path("./Analyses/SupRegInput"))

# list of best supported models
flist <- c("vrMuHiSSE3_tree0028.RDS", "vrMuHiSSE3_tree0612.RDS", "vrMuHiSSE4_tree0677.RDS",
           "vrMuHiSSE3_tree1030.RDS", "vrMuHiSSE3_tree1166.RDS", "vrMuHiSSE3_tree1774.RDS", "vrMuHiSSE3_tree1845.RDS",
           "vrMuHiSSE3_tree2966.RDS", "MuHiSSE3_tree3865.RDS", "vrMuHiSSE3_tree3865.RDS", "vrMuHiSSE3_tree4496.RDS",
           "vrMuHiSSE3_tree5024.RDS", "vrMuHiSSE3_tree5221.RDS", "MuHiSSE4_tree5333.RDS", "vrMuHiSSE4_tree5455.RDS",
           "vrMuHiSSE3_tree5684.RDS", "vrMuHiSSE3_tree5844.RDS", "vrMuHiSSE4_tree6259.RDS", "vrMuHiSSE3_tree6375.RDS",
           "vrMuHiSSE4_tree6375.RDS", "vrMuHiSSE3_tree6404.RDS", "vrMuHiSSE3_tree6668.RDS", "vrMuHiSSE3_tree7198.RDS",
           "vrMuHiSSE3_tree7306.RDS", "vrMuHiSSE3_tree8981.RDS", "vrMuHiSSE3_tree9128.RDS")
for(i in 1:length(flist)){
    # extract model name
    mod <- str_remove(flist[i], ".RDS")
    
    ## Shell file
    # this format is used to avoid the DOS (^M) end-of-line which trips the cluster. 
    # End-of-line DOS characters can also be removed from a bash file manually by using 
    # dos2nix filename.sh in the command line (use cat -v filename.sh to see if those 
    # characters are there at all).
    f <- file(paste0("./Analyses/SupRegInput/Pars",mod,".sh"), open = "wb")
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
        cat("#$ -l tmpfs=4G", sep="\n")
        cat("", sep="\n")
        cat("# Set the name of the job", sep="\n")
        # variable line
        cat(paste0("#$ -N Pars",mod), sep="\n")
        cat("", sep="\n")
        cat("# Set the working directory to somewhere in your scratch space.", sep="\n")
        cat("# This is a necessary step as compute nodes cannot write to $HOME.", sep="\n")
        cat("# NOTE: this directory must exist.", sep="\n")
        # variable line
        cat("#$ -wd /lustre/scratch/scratch/ucbtmao/ParEst/", sep="\n")
        cat("", sep="\n")
        cat("# Load the R module and run your R program", sep="\n")
        cat("module -f unload compilers mpi gcc-libs", sep="\n")
        cat("module load beta-modules", sep="\n")
        cat("module load r/r-4.0.2_bc-3.11", sep="\n")
        cat("export R_LIBS=/home/ucbtmao/r/r-4.0.2_bc-3.11:$R_LIBS", sep="\n")
        cat("", sep="\n")
        cat("# Copy all necessary files into TMPDIR", sep="\n")
        # variable line
        cat(paste0("cp ", mod, ".RDS $TMPDIR"), sep="\n")
        cat("", sep="\n")
        cat("# Your work should be done in $TMPDIR", sep="\n")
        cat("cd $TMPDIR", sep="\n")
        cat("", sep="\n")
        cat("# state program running time", sep="\n")
        # variable line
        cat(paste0("/usr/bin/time --verbose R --no-save < /lustre/scratch/scratch/ucbtmao/ParEst/Pars",mod,"Code.R> Pars",mod,"Code.out"), sep="\n")
        cat("", sep="\n")
        cat("# tar-up (archive) all output files to transfer back to own space", sep="\n")
        cat("tar -zcvf $HOME/Scratch/ParEst/files_from_job_$JOB_ID.tgz $TMPDIR", sep="\n")
        cat("", sep="\n")
        # variable line
        cat(paste0("cp Pars",mod,".RDS /lustre/scratch/scratch/ucbtmao/ParEst"), sep="\n")
    sink()
    close(f)    
            
            
    ## Generate R file
    f <- file(paste0("./Analyses/SupRegInput/Pars",mod,"Code.R"), open = "wb")
    sink(file=f)
        cat("# load required packages \n")
        cat("library(hisse)\n")
        cat("# load data \n")
        cat(paste0("model <- \"",mod,".RDS\" \n\n"))
        cat("fit_obj <- readRDS(model) \n")
        cat("# support region for param estimates \n")
        cat("supreg <- SupportRegionMuHiSSE(fit_obj, n.points=1000, scale.int=0.1, desired.delta=2, min.number.points=10, verbose=FALSE) \n")
        cat(paste0("saveRDS(supreg, file=\"Pars",mod,".RDS\")"))
    sink()
    close(f)  
}


####################################################################################
#####  T E M P L A T E  ######
##############################
library(hisse)

model <- "vrMuHiSSE3_tree6404.RDS"
fit_obj <- readRDS(model)

# support region for param estimates
supreg <- SupportRegionMuHiSSE(fit_obj, n.points=1000, scale.int=0.1, desired.delta=2, min.number.points=10, verbose=FALSE)
saveRDS(supreg, file="ParsvrMuHiSSE3_tree6404.RDS")


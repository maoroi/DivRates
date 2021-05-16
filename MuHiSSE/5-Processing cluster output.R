### Processing output from MuHiSSE models (run on cluster)
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster")

extRDS <- function(file){
    tmp <- readRDS(file)
    res <- tmp[c(1:3)]
}

allmods <- list.files(pattern = ".RDS")
modat <- data.frame(t(sapply(allmods, extRDS)))
score <- do.call(cbind, lapply(modat[,1:3], as.numeric))
modat$type <- NA
modat$states <- NA
modat$tree <- NA

#trees <- c("MCC", "0028", "1166", "1845", "2966", "3865", "4496", "5221", "5333", "5455", "6404", "7198", "9128")

# extracting model parameters from file name
for (i in 1:nrow(modat)){
    modat$type[i] <- gsub("[0-9]", "", strsplit(rownames(modat[i,]),"_")[[1]][1])       # model type
    modat$states[i] <- gsub("[^0-9]", "", strsplit(rownames(modat[i,]),"_")[[1]][1])    # no. of hidden states
    if (i %% 13 == 1) {
        modat$tree[i] <- "MCC"                                                          # tree variant
    } else {
        modat$tree[i] <- gsub("[^0-9]", "", strsplit(rownames(modat[i,]),"_")[[1]][2])
    }
} 

rownames(modat) <- NULL
write.csv(as.data.frame(modat), file="modat.csv")





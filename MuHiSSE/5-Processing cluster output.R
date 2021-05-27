### Processing output from MuHiSSE models (run on cluster)
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster")

extRDS <- function(file){
    tmp <- readRDS(file)
    res <- tmp[c(1:3)]
}

# extract the first 3 elements of each model into a data frame
allmods <- list.files(pattern = ".RDS")
modat <- data.frame(t(sapply(allmods, extRDS)))
score <- as.data.frame(do.call(cbind, lapply(modat[,1:3], as.numeric)))
rownames(score) <- rownames(modat)
score$type <- NA
score$states <- NA
score$tree <- NA

#trees <- c("MCC", "0028", "1166", "1845", "2966", "3865", "4496", "5221", "5333", "5455", "6404", "7198", "9128")

# extracting model parameters from file name
for (i in 1:nrow(score)){
    score$type[i] <- gsub("[0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1])       # model type
    score$states[i] <- gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1])    # no. of hidden states
    if (strsplit(rownames(score[i,]),"_")[[1]][2] == "MCCtree.RDS") {
        score$tree[i] <- "MCC"                                                          # tree variant
    } else {
        score$tree[i] <- gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][2])
    }
} 

score[which(score$type %in% c('MuSSE','null')),5] <- 1  # fill in info for single state models 
#rownames(score) <- NULL
#write.csv(score, file="score.csv")

### plotting likelihood as function of tree variant and model type
p <- ggplot(score) +
    geom_point(aes(tree, loglik, size = states, colour = type))
p



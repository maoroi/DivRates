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

## calculating mean values for tree variants
Mscore <- score[1:48,]
Mscore[1:24,] <- score[which(score$tree == "MCC"),]
rownames(Mscore) <- NULL
rownames(Mscore)[1:24] <- rownames(score)[which(score$tree == "MCC")]

for (j in 2:8){
    for (k in score$type){} 
    
### plotting likelihood as function of tree variant and model type
theme_set(theme_light(base_family = "Poppins"))
p <- ggplot(score[which(score$states != 1),]) +
    geom_point(aes(x=tree, loglik, colour = states)) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    facet_wrap(~type) + 
    scale_color_viridis_d(option = "viridis", 
                          alpha = .9,
                          begin = 0,
                          end = .9,
                          direction = -1)
p
    
 
d <- ggplot(score, aes(x = states, y = loglik, color = type)) #+
    scale_color_viridis_d(option = "viridis", 
                          alpha = .9,
                          begin = 0,
                          end = .9,
                          direction = -1)

d + geom_point()

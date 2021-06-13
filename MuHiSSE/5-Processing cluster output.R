### Processing output from MuHiSSE models (run on cluster)

library("tidyverse") 

#* Function definitions ------------------------------------------------------

extRDS <- function(file){
    tmp <- readRDS(file)
    res <- tmp[c(1:3)]
}

extRDS_par <- function(file){
    tmp <- readRDS(file)
    res <- tmp$solution
}
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster")

# extract the first 3 elements of each model into a data frame
allmods <- list.files(pattern = ".RDS")


# 1. Model fit ----------------------------------------------------------------

modat <- data.frame(t(sapply(allmods, extRDS)))
score <- as.data.frame(do.call(cbind, lapply(modat[,1:3], as.numeric)))
rownames(score) <- rownames(modat)
score$type <- NA
score$states <- NA
score$npar <- NA
score$tree <- NA

# extracting model parameters from file name
for (i in 1:nrow(score)){
    score$type[i] <- gsub("[0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1])       # model type
    score$states[i] <- as.numeric(gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1]))    # no. of hidden states
    if (score$type[i] == "CID") { 
        score$npar[i] <- score$states[i]+6
    } else if (score$type[i] == "MuHiSSE") { 
        score$npar[i] <- 7*score$states[i]+2
    } else if (score$type[i] == "vrCID") { 
        score$npar[i] <- score$states[i]*(score$states[i]+4)+1
    } else if (score$type[i] == "vrMuHiSSE") { 
        score$npar[i] <- score$states[i]*(score$states[i]+6)+1
    } 
    
    if (strsplit(rownames(score[i,]),"_")[[1]][2] == "MCCtree.RDS") {
        score$tree[i] <- "MCC"  # tree variant
    } else {
        score$tree[i] <- gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][2])}
} 
# adjust info for single state models 
score[which(score$type %in% c('MuSSE','null')),5] <- 1  
score$npar[which(score$type == 'MuSSE')] <- 8
score$type[which(score$type == 'MuSSE')] <- 'MuHiSSE'
score$npar[which(score$type == 'null')] <- 6
score$type[which(score$type == 'null')] <- 'CID'
#write.csv(score, file="score.csv", row.names = FALSE)

ord <- score[order(score$AICc),]
ord <- ord[,c("tree","type","states","loglik","npar","AIC","AICc")]
ordbyt <- data.frame()
# calculate relative model support for each tree
for (z in unique(ord$tree)) {
    dat <- ord[which(ord$tree == z),]
    dat <- dat[order(dat$AICc),]
    dat$diff_AIC <- dat$AICc - min(dat$AICc)
    dat$rel_lik <- exp(-0.5 * dat$diff_AIC)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik)
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_AIC, decreasing = FALSE),])
}
ordbyt <- ordbyt[order(ordbyt$tree, ordbyt$AICc),]
#write.csv(ordbyt, file="Model fit by tree.csv")


# * Diagnostic plots ----------------------------------------------------------

## plotting likelihood as function of tree variant and model type
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")

# this line orders the trees in the plot by increasing loglik (by eye)
#score$tree = factor(score$tree, levels=c("3865","5333","2966","5455","MCC","0028","5221","4496","1845","9128","7198","1166","6404"))

df1 <- score[score$states == 1,]
pdf(file = "likelihood by model type.pdf",width = 10, height = 10)  #png(file = "likelihood by model type.png",width = 7, height = 8.5, units = "in", res=800)
ggplot(score) +
    geom_point(aes(x=tree, loglik, colour = as.factor(states)), alpha= .8, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    geom_point(data = df1, aes(tree, loglik), colour = 'grey20', shape = 1, size = 2) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    facet_wrap(~type) + #, nrow = 1) + 
    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)
dev.off()

pdf(file = "likelihood by no. of states.pdf", width = 7, height = 8.5)
ggplot(score) +
    geom_point(aes(x=states, loglik, colour = reorder(tree, loglik, FUN = mean)), alpha= .5, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    geom_point(data = score[score$tree == "MCC",], aes(states, loglik), colour = 'grey20', shape = 3, size = 2) + 
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "Tree variant") +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .9, begin = 0, end = .9, direction = -1)
dev.off()

# comparing likelihood of all models, divided to model type and no. of states, grouped by tree
pdf(file="likelihood by states_MCC marked.pdf", width = 10, height = 8)
df1 <- score[score$states != 1,]
df2 <- subset(df1, df1$tree =="MCC")    # get only MCC models separately

ggplot(df1, aes(x = reorder(tree, loglik, FUN = max), y = loglik, colour = as.factor(states))) +
    geom_point(alpha= .8, size = 2) + 
    #geom_point(data = df2, colour = 'white', shape = 3, size = 2) + 
    geom_point(data = df2, aes(tree, loglik), colour = 'grey20', shape = 1, size = 3) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)
dev.off()

# comparing model type performance
pdf(file="model type performance.pdf", width = 10, height = 8)
ggplot(ordbyt) +
    geom_point(aes(x=states, loglik, colour = type), alpha= .7, size = 3) +  
    theme_light() +
    facet_wrap(~tree, nrow = 3) + 
    scale_color_viridis_d(option = "viridis", #inferno, turbo, mako
                          begin = 0, end = .9, direction = -1)
dev.off()

# comparing likelihood of all models, divided to trees and no. of states, grouped by model type
#pdf(file="likelihood by trees.pdf", width = 10, height = 8)
#ggplot(ordbyt) +
#    geom_point(aes(x=type, loglik, colour = factor(states)), alpha= .8, size = 2) +
#    theme_light() +
#    theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1)) +
#    facet_wrap(~tree, nrow = 3) + 
#    scale_color_viridis_d(option = "viridis", 
#                          begin = 0,
#                          end = .9,
#                          direction = -1)
#dev.off()

# 2. Parameter estimates ------------------------------------------------------
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster")

pardat <- data.frame(t(sapply(allmods, extRDS_par)))
params <- pardat[which(is.na(str_extract(names(pardat), ".10.")))]      # remove state "10" (disabled in all models)
params <- params[which(is.na(str_extract(names(params), ".00..11.|.11..00.")))]  # remove 'diagonal' change (disabled in all)


# 2.1 calculate lambda and mu from composite rates ----------------------------

#   reminder: turnover = spec+ext; netDiv = spec-ext; extFrac = ext/spec 
#   symbols: turnover-tau; extFrac-eps; spec_rate-lambda; ext_rate-mu; 
#   Developing two equations w two unknowns on paper, got this:
#
#  // lambda = tau/(1+eps) \\ 
#  \\ mu = tau-lambda      //


div <- params           # create a mirror table to write lambda and mu values into
div[,which(is.na(str_extract(names(params), "q.")))] <- NA  # remove composite rates
# iterate over all models and calculate lambda and mu for each
for (i in colnames(params)) {
    if (is.na(str_extract(i, "turn.")) == FALSE) {
        index <- strsplit(i, "over")[[1]][2]
        corr_rate <- which(colnames(params) == paste0("eps", index))
        for (j in rownames(params)) {
            div[j,i] <- params[j,i] / (1 + params[j,corr_rate])
        }
        colnames(div)[which(colnames(div) == i)] <- paste0("lambda",index)
    } else if (is.na(str_extract(i, "eps.")) == FALSE) {
        index <- strsplit(i, "eps")[[1]][2]
        corr_rate <- which(colnames(div) == paste0("lambda", index))
        for (j in rownames(params)) {
            div[j,i] <- params[j,i] * div[j,corr_rate]
        }
        colnames(div)[which(colnames(div) == i)] <- paste0("mu",index)
    }
}
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")
#write.csv(div, file="Rates of Evolution.csv", row.names = TRUE)

# combine fit and params estimates to one table
results <- ordbyt
filler <- matrix(NA, nrow = nrow(ordbyt), ncol = ncol(div))
colnames(filler) <- colnames(div)
results <- cbind(ordbyt, filler)
for (i in 1:nrow(ordbyt)) {
    results[i,] <- cbind(ordbyt[i,], div[which(rownames(div) == rownames(ordbyt)[i]),])
}
#write.csv(results, file = "Model estimates exploratory.csv")


# 2.2 Parameter estimates for leading models ----------------------------------

# take all supported models
leading <- results[which(rownames(results) %in% rownames(ordbyt)[which(ordbyt$diff_AIC < 6)]),]
means <- colMeans(leading[-which(leading$tree == "MCC"),10:ncol(leading)])
means <- c(rep(NA, 9), means)
leading <- rbind(leading, means)
#write.csv(leading, file="Best supported models.csv")

est <-matrix(data=NA, ncol = 3, nrow = 5)
colnames(est) <- c("Noct","Cath","Diur")
rownames(est) <- c("Noct","Cath","Diur","Spec","Ext")
est[c(3,11)] <- 0
#
#
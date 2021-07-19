### Processing output from MuHiSSE models (run on cluster)

library("tidyverse")
library("hisse")

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
score$BIC <- NA
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
    
    # tree variant index
    if (strsplit(rownames(score[i,]),"_")[[1]][2] == "MCCtree.RDS") {
        score$tree[i] <- "MCC"  # tree variant
    } else {
        score$tree[i] <- gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][2])}
} 
# adjust info for single state models 
score$states[which(score$type %in% c('MuSSE','null'))] <- 1  
score$npar[which(score$type == 'MuSSE')] <- 8
score$type[which(score$type == 'MuSSE')] <- 'MuHiSSE'
score$npar[which(score$type == 'null')] <- 6
score$type[which(score$type == 'null')] <- 'CID'

# calculate BIC [the number of data points (n=2400) should change in clade-wise models
for (i in 1:nrow(score)){
    score$BIC[i] <- -2*score$loglik[i] + log(2400) * score$npar[i]}


tert <- score[which(score$type %in% c("CID","MuHiSSE","vrCID","vrMuHiSSE")),]

## model selection by AICc and BIC
ordA <- tert[order(tert$AICc),]
ordA <- ordA[,c("tree","type","states","loglik","npar","AIC","AICc","BIC")]
ordB <- tert[order(tert$BIC),]
ordB <- ordB[,c("tree","type","states","loglik","npar","AIC","AICc","BIC")]

ordbyt <- data.frame()
# calculate delta AICc/BIC for each tree
for (z in unique(ordA$tree)) {
    dat <- ordA[which(ordA$tree == z),]
    dat <- dat[order(dat$AICc),]
    dat$diff_AIC <- dat$AICc - min(dat$AICc)
    dat$rel_lik <- exp(-0.5 * dat$diff_AIC)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik) # this is identical to rel_lik on most trees
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_AIC, decreasing = FALSE),])
}
ordbytA <- ordbyt[order(ordbyt$tree, ordbyt$AICc),]

# and again for BIC:
ordbyt <- data.frame()
# calculate relative model support for each tree
for (z in unique(ordB$tree)) {
    dat <- ordB[which(ordB$tree == z),]
    dat <- dat[order(dat$BIC),]
    dat$diff_BIC <- dat$BIC - min(dat$BIC)
    dat$rel_lik <- exp(-0.5 * dat$diff_BIC)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik) # this is identical to rel_lik on most trees
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_BIC, decreasing = FALSE),])
}
ordbytB<- ordbyt[order(ordbyt$tree, ordbyt$BIC),]


# * Summarise support (plots) ----------------------------------------------------------

setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")
#write.csv(score, file="score.csv", row.names = FALSE)
#write.csv(ordbyt, file="Model fit by tree.csv")

## plotting likelihood as function of tree variant and model type

# comparing all models by BIC, divided to model type and no. of states, grouped by tree
pdf(file="Model support by states_MCC marked.pdf", width = 10, height = 8)
df1 <- tert[tert$states != 1,]
df2 <- subset(df1, df1$tree =="MCC")    # get only MCC models separately
#par(mfrow=c(2,1))
ggplot(df1, aes(x = reorder(tree, AICc, FUN = min), y = AICc, colour = as.factor(states))) +
    geom_point(alpha= .5, size = 2.5) + 
    #geom_point(data = df2, colour = 'white', shape = 3, size = 2) + 
    geom_point(data = df2, aes(tree, AICc), colour = 'red', shape = 1, size = 3) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)

ggplot(df1, aes(x = reorder(tree, BIC, FUN = min), y = BIC, colour = as.factor(states))) +
    geom_point(alpha= .5, size = 2.5) + 
    #geom_point(data = df2, colour = 'white', shape = 3, size = 2) + 
    geom_point(data = df2, aes(tree, BIC), colour = 'red', shape = 1, size = 3) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    scale_color_viridis_d(option = "inferno", begin = 0, end = .9)#
dev.off()

df3 <- tert[tert$states == 1,]
pdf(file = "Model support by model type.pdf",width = 10, height = 10)  #png(file = "likelihood by model type.png",width = 7, height = 8.5, units = "in", res=800)
ggplot(tert) +
    geom_point(aes(x=tree, AICc, colour = as.factor(states)), alpha= .8, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    geom_point(data = df3, aes(tree, AICc), colour = 'grey20', shape = 1, size = 2) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    facet_wrap(~type) + #, nrow = 1) + 
    scale_color_viridis_d(option = "viridis", begin = .1, end = .9, direction = -1)

ggplot(tert) +
    geom_point(aes(x=tree, BIC, colour = as.factor(states)), alpha= .8, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    geom_point(data = df3, aes(tree, BIC), colour = 'grey20', shape = 1, size = 2) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    facet_wrap(~type) + #, nrow = 1) + 
    scale_color_viridis_d(option = "inferno", begin = .1, end = .9, direction = -1)
dev.off()

pdf(file = "Model support by no. of states.pdf", width = 7, height = 8.5)
## code for this figure adapted from: https://ourcodingclub.github.io/tutorials/dataviz-beautification-synthesis/#distributions 
# This code loads the function in the working environment
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

vars <- tert[-which(tert$tree == "MCC"),]
vrci <- vars[which(vars$type == "vrCID"),]
vrsse <- vars[(which(vars$type == "vrMuHiSSE")),]

library(ggridges)
ggplot(vrsse, aes(x = as.factor(states), y = BIC, fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.7) +
    geom_point(aes(y = BIC), colour = 'grey30', position = position_jitter(width = .1), shape = 21, size = 2, alpha = 0.8) +
    #theme_light() +
    scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    labs(y = "BIC", x = "States") +
    # Removing legends
    guides(fill = FALSE, color = FALSE) 


ggplot(vrci, aes(x = as.factor(states), y = BIC, fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.7) +
    geom_point(aes(y = BIC), colour = 'grey30', position = position_jitter(width = .1), shape = 21, size = 2, alpha = 0.8) +
    theme_light() + 
    scale_fill_viridis_d(begin = .3, end = .9, direction = -1, alpha=0.8) +
    labs(y = "BIC", x = "States") +
    # Removing legends
    guides(fill = FALSE, color = FALSE) +
    # Setting the limits of the y axis
    scale_y_continuous(limits = c(17400, 18600))# +
    # Picking nicer colours
    #scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    #scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C"))
    
ggplot(tert) +
    geom_point(aes(x=states, AICc, colour = reorder(tree, AICc, FUN = mean)), alpha= .5, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    geom_point(data = tert[tert$tree == "MCC",], aes(states, AICc), colour = 'grey20', shape = 3, size = 2) + 
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "Tree variant") +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .9, begin = 0, end = .9, direction = -1)

ggplot(tert) +
    geom_point(aes(x=states, BIC, colour = reorder(tree, BIC, FUN = mean)), alpha= .5, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    geom_point(data = tert[tert$tree == "MCC",], aes(states, BIC), colour = 'grey20', shape = 3, size = 2) + 
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "Tree variant") +
    facet_wrap(~type,) #+ 
    scale_color_viridis_d(option = "inferno", alpha = .9, begin = 0, end = .9, direction = -1)
dev.off()

# comparing likelihood of all models, divided to model type and no. of states, grouped by tree
pdf(file="likelihood by states_MCC marked.pdf", width = 10, height = 8)
df1 <- tert[tert$states != 1,]
df2 <- subset(df1, df1$tree =="MCC")    # get only MCC models separately

ggplot(df1, aes(x = reorder(tree, loglik, FUN = max), y = loglik, colour = as.factor(states))) +
    geom_point(alpha= .8, size = 2.5) + 
    #geom_point(data = df2, colour = 'white', shape = 3, size = 2) + 
    geom_point(data = df2, aes(tree, loglik), colour = 'red', shape = 1, size = 3) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)
dev.off()

# comparing model type performance
pdf(file="model type performance.pdf", width = 10, height = 8)
ggplot(ordbytA) +
    geom_point(aes(x=states, AICc, colour = type), alpha= .5, size = 3) +  
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "viridis", #inferno, turbo, mako
                          begin = 0, end = .9, direction = -1)
ggplot(ordbytB) +
    geom_point(aes(x=states, BIC, colour = type), alpha= .5, size = 3) +  
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "inferno", #inferno, turbo, mako
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


# 2.1 calculate s and mu from composite rates ----------------------------

#   reminder: turnover = spec+ext; netDiv = spec-ext; extFrac = ext/spec 
#   symbols: turnover-tau; extFrac-eps; spec_rate-s; ext_rate-mu; netDiv-lambda;
#   Developing two equations w two unknowns on paper, got this:
#
#  // s = tau/(1+eps) \\ 
#  \\ mu = tau-s      //

div <- params           # create a mirror table to write s and mu values into
div[,which(is.na(str_extract(names(params), "q.")))] <- NA  # remove composite rates
# iterate over all models and calculate s and mu for each
for (i in colnames(params)) {
    if (is.na(str_extract(i, "turn.")) == FALSE) {
        index <- strsplit(i, "over")[[1]][2]
        corr_rate <- which(colnames(params) == paste0("eps", index))
        for (j in rownames(params)) {
            div[j,i] <- params[j,i] / (1 + params[j,corr_rate])
        }
        colnames(div)[which(colnames(div) == i)] <- paste0("s",index)
    } else if (is.na(str_extract(i, "eps.")) == FALSE) {
        index <- strsplit(i, "eps")[[1]][2]
        corr_rate <- which(colnames(div) == paste0("s", index))
        for (j in rownames(params)) {
            div[j,i] <- params[j,i] * div[j,corr_rate]
        }
        colnames(div)[which(colnames(div) == i)] <- paste0("mu",index)
    }
}
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")
#write.csv(div, file="Rates of Evolution.csv", row.names = TRUE)

# combine fit and params estimates to one table
results <- ordbytB
filler <- matrix(NA, nrow = nrow(ordbytB), ncol = ncol(div))
colnames(filler) <- colnames(div)
results <- cbind(ordbytB, filler)
for (i in 1:nrow(ordbytB)) {
    results[i,] <- cbind(ordbytB[i,], div[which(rownames(div) == rownames(ordbytB)[i]),])
}
#write.csv(results, file = "Model estimates exploratory.csv")


# 2.2 Parameter estimates for leading models ----------------------------------

# take all supported models
leading <- results[which(rownames(results) %in% rownames(ordbytB)[which(ordbytB$diff_BIC < 6)]),]
#leading <- leading[order(leading$states),]
# calculate trait-dependent netDiv values per state
for (states in 1:5) {
    for (patt in c("00","01","11")) {
        s_rate <- paste0("s",patt,LETTERS[states])
        mu_rate <- paste0("mu",patt,LETTERS[states])
        leading$new <- leading[,which(colnames(leading) == s_rate)] - leading[,which(colnames(leading) == mu_rate)]
        if (states == 5) {
            leading$new[which(leading$new == 0)] <- NA
        }
    colnames(leading)[ncol(leading)] <- paste0("lambda",patt,LETTERS[states])
    }
}

means <- colMeans(leading[-which(leading$tree == "MCC"),10:ncol(leading)])
means <- c(rep(NA, 9), means)
leading <- rbind(leading, means)
rownames(leading)[which(is.na(leading$tree))] <- "mean_24trees"
# remove rates not considered in the model (states F, G, H)
leading <- leading[which(is.na(str_extract(colnames(leading), "..F|..G|..H")))]  
#write.csv(leading, file="Best supported models.csv")


#
#

# 3. Follow up analyses -------------------------------------------------------

# 3.1 Reconstruction based on fitted models -----------------------------------
TODO: extract parameters from estimated models
reco <- MarginReconMuHiSSE(tree, act3, f = freq, pars, hidden.states=1, 
                           condition.on.survival=TRUE, root.type="madfitz", 
                           root.p=NULL, AIC=NULL, get.tips.only=FALSE, 
                           verbose=TRUE, n.cores=NULL, dt.threads=1)



# 5. Repeat for binary models -------------------------------------------------



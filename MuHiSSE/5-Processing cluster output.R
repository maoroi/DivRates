### Processing output from MuHiSSE models (run on cluster)

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

score[which(score$type %in% c('MuSSE','null')),5] <- 1  # fill in info for single state models 
score$npar[which(score$type == 'MuSSE')] <- 8
score$type[which(score$type == 'MuSSE')] <- 'MuHiSSE'
score$npar[which(score$type == 'null')] <- 6
score$type[which(score$type == 'null')] <- 'CID'


#rownames(score) <- NULL
#write.csv(score, file="score.csv")
ord <- score[order(score$AICc),]
ord <- ord[,c("tree","type","states","loglik","npar","AIC","AICc")]
ordbyt <- data.frame()
for (z in unique(ord$tree)) {
    dat <- ord[ord$tree == z,]
    dat <- dat[order(dat$AICc),]
    dat$diff_AIC <- NA
    dat$rel_lik <- NA
    for (i in 1:nrow(dat)) {
        dat$diff_AIC[i] <- dat$AICc[i] - min(dat$AICc)
        dat$rel_lik[i] <- exp(-0.5 * dat$diff_AIC[i])
    }
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_AIC, decreasing = FALSE),])
}
ordbyt <- ord[order(ord$tree, ord$AICc),]
#write.csv(ord, file="score_ordered.csv")


# * Diagnostic plots ----------------------------------------------------------

## plotting likelihood as function of tree variant and model type
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")

# this line orders the trees in the plot by increasing loglik (by eye)
#score$tree = factor(score$tree, levels=c("3865","5333","2966","5455","MCC","0028","5221","4496","1845","9128","7198","1166","6404"))

df1 <- score[score$states == 1,]
pdf(file = "likelihood by model type.pdf",width = 10, height = 10)  #png(file = "likelihood by model type.png",width = 7, height = 8.5, units = "in", res=800)
ggplot(score) +
    geom_point(aes(x=tree, loglik, colour = states), alpha= .8, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    geom_point(data = df1, aes(tree, loglik), colour = 'grey20', shape = 1, size = 2) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    facet_wrap(~type) + #, nrow = 1) + 
    scale_color_viridis_d(option = "viridis", 
                          begin = 0,
                          end = .9,
                          direction = -1)
dev.off()
 

pdf(file = "likelihood by no. of states.pdf", width = 7, height = 8.5)
ggplot(score[which(score$states != 1),]) +
    geom_point(aes(x=states, loglik, colour = tree), alpha= .5, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "viridis", 
                          alpha = .9,
                          begin = 0,
                          end = .9,
                          direction = -1)
dev.off()


# comparing likelihood of all models, divided to model type and no. of states, grouped by tree
pdf(file="likelihood by states_MCC marked.pdf", width = 10, height = 8)
df2 <- subset(score, score$tree =="MCC")    # get only MCC models separately

ggplot(score, aes(x = tree, y = loglik, colour = states)) +
    geom_point(alpha= .8, size = 2) + 
    geom_point(data = df2, colour = 'white', shape = 3, size = 2) + 
    geom_point(data = df1, aes(tree, loglik), colour = 'grey20', shape = 1, size = 2) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    scale_color_viridis_d(option = "inferno", 
                          begin = 0.1,
                          end = .9,
                          direction = -1)
dev.off()


# comparing likelihood of all models, divided to trees and no. of states, grouped by model type
pdf(file="likelihood by trees.pdf", width = 10, height = 8)
ggplot(ordbyt) +
    geom_point(aes(x=type, loglik, colour = factor(states)), alpha= .8, size = 2) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1)) +
    facet_wrap(~tree, nrow = 3) + 
    scale_color_viridis_d(option = "viridis", 
                          begin = 0,
                          end = .9,
                          direction = -1)
dev.off()

# comparing model type performance
pdf(file="model type performance.pdf", width = 10, height = 8)
ggplot(ordbyt) +
    geom_point(aes(x=states, loglik, colour = type), alpha= .7, size = 3) +  
    theme_light() +
    facet_wrap(~tree) + 
    scale_color_viridis_d(option = "viridis", 
                          begin = 0,
                          end = .9,
                          direction = -1)
dev.off()


# 2. Parameter estimates ------------------------------------------------------

pardat <- data.frame(t(sapply(allmods, extRDS_par)))
params <- pardat[which(is.na(str_extract(names(pardat), ".10.")))]      # remove state "10" (disabled in all models)
params <- params[which(is.na(str_extract(names(params), ".00..11.|.11..00.")))]  # remove 'diagonal' change (disabled in all)

# reminder: turnover = spec+ext; netDiv = spec-ext; extFrac = ext/spec 

lambda = TO-myu = TO-(extFrac * lambda) = TO- 
mu = extFrac*lambda = 

    
    
for (stat in as.numeric(score$states)) {
    
}
#
#
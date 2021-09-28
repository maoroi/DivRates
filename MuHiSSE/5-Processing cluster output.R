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
score$BIC <- NA
score$type <- NA
score$states <- NA
score$npar <- NA
score$tree <- NA


## 1.1 Extract model parameters ------------------------------------------------
# extracting model parameters from file name
for (i in 1:nrow(score)){
    score$type[i] <- gsub("[0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1])       # model type
    score$states[i] <- as.numeric(gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1]))    # no. of hidden states
    if (score$type[i] == "CID") { 
        score$npar[i] <- score$states[i]+6  # no. of parameters calculated in "Evolutionary models.xlsx"
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


## 1.2 Model selection ---------------------------------------------------------

# model selection by AICc and BIC
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
    dat$weight <- dat$rel_lik / sum(dat$rel_lik)
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_BIC, decreasing = FALSE),])
}
ordbytB<- ordbyt[order(ordbyt$tree, ordbyt$BIC),]


## 1.3 Support summary plots --------------------------------------------------

setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Output")
#write.csv(score, file="score.csv", row.names = FALSE)
#write.csv(ordbyt, file="Model fit by tree.csv")

## plotting likelihood as function of tree variant and model type

# comparing all models by BIC, divided to model type and no. of states, grouped by tree
pdf(file="Model support ordered_MCC marked.pdf", width = 10, height = 8)
df1 <- tert[tert$states != 1,]
df2 <- subset(df1, df1$tree =="MCC")    # get only MCC models separately
#par(mfrow=c(2,1))
ggplot(df1, aes(x = reorder(tree, AICc, FUN = min), y = AICc, colour = as.factor(states))) +
    geom_point(alpha= .75, size = 3) + 
    geom_point(data = df2, aes(tree, AICc), colour = 'red', shape = 1, size = 3) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("Tree variant (reordered for visualisation)") +
    scale_color_viridis_d(option = "viridis", begin = .15, end = 1, alpha = 1, direction = -1) +
    scale_fill_viridis_d(option = "inferno", begin = .1, end = .9, alpha = 0.5, direction = -1)

ggplot(df1, aes(x = reorder(tree, BIC, FUN = min), y = BIC, colour = as.factor(states))) +
    geom_point(alpha= .75, size = 3) + 
    geom_point(data = df2, aes(tree, BIC), colour = 'red', shape = 1, size = 3) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("Tree variant (reordered for visualisation)") +
    scale_color_viridis_d(option = "inferno", begin = .1, end = .9, alpha = 1) +
    scale_fill_viridis_d(option = "inferno", begin = .1, end = .9, alpha = 0.5)
dev.off()

df3 <- tert[tert$states == 1,]
pdf(file = "Model support unordered.pdf",width = 10, height = 10)  #png(file = "likelihood by model type.png",width = 7, height = 8.5, units = "in", res=800)
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

#pdf(file = "Model support by no. of states.pdf", width = 7, height = 8.5)
#library(ggridges)

## code for this figure adapted from: https://ourcodingclub.github.io/tutorials/dataviz-beautification-synthesis/#distributions 
# This code loads the function in the working environment
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
# theme from: https://neuroconscience.wordpress.com/2018/03/15/introducing-raincloud-plots/
theme_niwot <- function(){
    theme_bw() +
        theme(text = element_text(family = "Helvetica Light"),
              axis.text = element_text(size = 16),
              axis.title = element_text(size = 18),
              axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"),
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_blank(),
              plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
              plot.title = element_text(size = 18, vjust = 1, hjust = 0),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              legend.position = c(0.95, 0.15),
              legend.key = element_blank(),
              legend.background = element_rect(color = "black",
                                               fill = "transparent",
                                               size = 2, linetype = "blank"))
}
    
raincloud_theme = theme(
    text = element_text(size = 10),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title=element_text(size=16),
    legend.text=element_text(size=16),
    legend.position = "right",
    plot.title = element_text(lineheight=.8, face="bold", size = 16),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

vars <- tert[which(tert$tree != "MCC"),]
vr <- vars[which(vars$type %in% c("vrCID", "vrMuHiSSE")),]
vr <- vr[-which(vr$states == 6),]

#vrci <- vars[which(vars$type == "vrCID"),]
#vrsse <- vars[(which(vars$type == "vrMuHiSSE")),]
#
#ggplot(vrsse, aes(x = as.factor(states), y = BIC, fill = as.factor(states))) +
#    #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#    geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.7) +
#    geom_point(aes(y = BIC), colour = 'grey30', position = position_jitter(width = .1), shape = 21, size = 2, alpha = 0.8) +
#    #theme_light() +
#    scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
#    scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
#    labs(y = "BIC", x = "States") +
#    # Removing legends
#    guides(fill = FALSE, color = FALSE) 
#
#ggplot(vrci, aes(x = as.factor(states), y = BIC, fill = as.factor(states))) +
#    #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
#    geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.7) +
#    geom_point(aes(y = BIC), colour = 'grey30', position = position_jitter(width = .1), shape = 21, size = 2, alpha = 0.8) +
#    theme_light() + 
#    scale_fill_viridis_d(begin = .3, end = .9, direction = -1, alpha=0.8) +
#    labs(y = "BIC", x = "States") +
#    # Removing legends
#    guides(fill = FALSE, color = FALSE) +
#    # Setting the limits of the y axis
#    scale_y_continuous(limits = c(17400, 18600)) +
#    # Picking nicer colours
#    scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","2CD11E")) +
#    scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","2CD11E")) 
   
# plotting support by no. states, aggregaetd over all variants
pdf(file = "Model support VR aggregated.pdf", width = 8.5, height = 6)
ggplot(data = vr, aes(y = AICc, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = AICc, x = as.factor(states), color = as.factor(states)), 
               position = position_jitter(width = .1), size = 3, alpha = 0.6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    guides(fill = FALSE, color = FALSE) + # don't show a legend
    #labs(color = "States", fill = "States") + # in case a legend is needed
    xlab("No. of states") +
    facet_wrap(~type,) +
    expand_limits(x = 5.25) +
    #scale_color_brewer(palette = "Spectral") +
    #scale_fill_brewer(palette = "Spectral") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .95, alpha = 0.8, direction = -1) +
    scale_color_viridis_d(option = "viridis", begin = 0.1, end = .95, alpha = 1, direction = -1) +
    #scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6E")) +
    #scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6E")) +
    theme_minimal(base_size = 12) 

ggplot(data = vr, aes(y = BIC, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = BIC, x = as.factor(states), color = as.factor(states)), 
               position = position_jitter(width = .1), size = 3, alpha = 0.6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    guides(fill = FALSE, color = FALSE) + # don't show a legend
    #labs(color = "States", fill = "States") + # in case a legend is needed
    xlab("No. of states") +
    facet_wrap(~type,) +
    expand_limits(x = 5.25) +
    #scale_color_brewer(palette = "Spectral") +
    #scale_fill_brewer(palette = "Spectral") +
    scale_fill_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 0.8, direction = -1) +
    scale_color_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 1, direction = -1) +
    #scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6E")) + #, "2CD11E")) +
    #scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6E")) + #, "2CD11E")) +
    theme_minimal(base_size = 12) 
dev.off()

# plotting support by no. states with trendlines to show that the pattern is uniform accross trees
pdf(file = "Model support by states_tree-trendlines.pdf", width = 7, height = 8.5)
ggplot(vars,aes(x=states, AICc, colour = reorder(tree, AICc, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    # use the below if need to plot MCC along with tree variants
    #geom_point(data = tert[tert$tree == "MCC",], aes(states, AICc), colour = 'grey10', shape = 4, size = 2) + 
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "Tree variant") +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .5, direction = -1)

ggplot(vars,aes(x=states, BIC, colour = reorder(tree, BIC, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    # use the below if need to plot MCC along with tree variants
    #geom_point(data = tert[tert$tree == "MCC",], aes(states, BIC), colour = 'grey10', shape = 4, size = 2) + 
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "Tree variant") +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .5, direction = -1)
dev.off()

# comparing likelihood of all models, divided to model type and no. of states, grouped by tree
#pdf(file="likelihood by states_MCC marked.pdf", width = 10, height = 8)
#df1 <- tert[tert$states != 1,]
#df2 <- subset(df1, df1$tree =="MCC")    # get only MCC models separately
#
#ggplot(df1, aes(x = reorder(tree, loglik, FUN = max), y = loglik, colour = as.factor(states))) +
#    geom_point(alpha= .8, size = 2.5) + 
#    #geom_point(data = df2, colour = 'white', shape = 3, size = 2) + 
#    geom_point(data = df2, aes(tree, loglik), colour = 'red', shape = 1, size = 3) + 
#    facet_wrap(~type) + 
#    theme_light() +
#    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
#    labs(color = "States") +
#    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)
#dev.off()

# comparing model type performance
pdf(file="model type performance by tree.pdf", width = 10, height = 8)
ggplot(ordbytA, aes(x=states, AICc, colour = type, group=type)) +
    geom_point(alpha= .5, size = 3) + 
    geom_line(size = 1.5, alpha = 0.5) +
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "viridis", #inferno, turbo, mako
                          begin = 0, end = .95, direction = -1)
ggplot(ordbytB, aes(x=states, BIC, colour = type, group=type)) +
    geom_point(alpha= .5, size = 3) +  
    geom_line(size = 1.5, alpha = 0.5) +
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


## 2.1 calculate s and mu from composite rates --------------------------------

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
#setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Output")
#write.csv(div, file="Raw Rates of Evolution.csv", row.names = TRUE)


## 2.2 Obtain parameter estimates from averaging supported models -------------

# combine fit and params estimates to one table
results <- ordbytB
filler <- matrix(NA, nrow = nrow(ordbytB), ncol = ncol(div))
colnames(filler) <- colnames(div)
results <- cbind(ordbytB, filler)
for (i in 1:nrow(ordbytB)) {
    results[i,] <- cbind(ordbytB[i,], div[which(rownames(div) == rownames(ordbytB)[i]),])
}
#write.csv(results, file = "Model estimates exploratory.csv")


### 2.2.1 Calculate diversification rates (lambda) ----------------------------

# take all BIC supported models (start from the previous step to repeat this process for AIC)
leading <- results[which(rownames(results) %in% rownames(ordbytB)[which(ordbytB$diff_BIC < 6)]),]
#leading <- leading[order(leading$states),]
# calculate trait-dependent netDiv values per state
for (state in 1:max(leading$states)) {
    for (patt in c("00","01","11")) {
        s_rate <- paste0("s", patt, LETTERS[state])
        mu_rate <- paste0("mu", patt, LETTERS[state])
        # claculate diversification rate
        leading$new <- leading[,which(colnames(leading) == s_rate)] - leading[,which(colnames(leading) == mu_rate)]
        # ensure that rates not in the model are marked NA instead of 0 (to exclude them from downstream calculations)
        if (state > min(leading$states)) { 
            leading$new[which(leading$states < state)] <- NA
        }
    # rename new column    
    colnames(leading)[ncol(leading)] <- paste0("lambda", patt, LETTERS[state])
    }
}


### 2.2.2 model weighting -----------------------------------------------------
## I did not average separately competing models of the same tree variant because they 
## indicate different number of parameters so not comparable. 

# weight estimated values by model support for all tree variants (not for MCC)
mod_avg <- leading
for (i in 1:length(which(mod_avg$tree != "MCC"))) {
    mod_avg[i,12:ncol(mod_avg)] <- mod_avg$weight[i] * mod_avg[i,12:ncol(mod_avg)]
} 

# ** the below averaging is likely wrong so commented out **
#means <- colSums(mod_avg[which(mod_avg$tree != "MCC"),12:ncol(mod_avg)]) / 24 # mean rates
#means <- c(rep(NA, 11), means)
#mod_avg <- rbind(mod_avg, means)
#rownames(mod_avg)[which(is.na(mod_avg$tree))] <- "averaged_24trees"
## remove rates not considered in the model (states F, G, H)
#mod_avg <- mod_avg[which(is.na(str_extract(colnames(mod_avg), "..E|..F|..G|..H")))]  
# ** end dodgy averaging procedure **

#write.csv(mod_avg, file="Best supported models.csv")


## 2.3 Rate estimates ---------------------------------------------------------

## split into evolutionary rates and character transition rates
evol <- mod_avg[,-c(1:11)]
evol <- evol[which(is.na(str_extract(colnames(evol), "q")))] # all rates that are NOT transition rates

tran <- mod_avg[,-c(1:11)]
tran <- tran[which(!is.na(str_extract(colnames(tran), "q")))] # transition rates only


### 2.3.1 Evolutionary rates --------------------------------------------------
evol$model <- str_replace(rownames(evol), ".RDS", "") # remove file suffix from model names

# melt into long form removing the (wrongly) averaged values
evo <- as_tibble(evol, subset = rownames(evol) != "averaged_24trees") %>%
    pivot_longer(!model, names_to = "rate", values_to = "value") 

{type <- modtype <- rtype <- AP <- state <- character()

type[which(!is.na(str_extract(evo$rate, "s")))] <- "spec"
type[which(!is.na(str_extract(evo$rate, "mu")))] <- "ext"
type[which(!is.na(str_extract(evo$rate, "lambda")))] <- "net.div"

for (i in 1:nrow(evo)){
    modtype[i] <- str_extract(str_split(evo$model[i], "_")[[1]][1], "[0-9]")
    rtype[i] <- str_sub(evo$rate[i], start = 1L, end = -2L)}

AP[which(str_sub(evo$rate, start = -3L, end = -2L) == "00")] <- "Noct"
AP[which(str_sub(evo$rate, start = -3L, end = -2L) == "01")] <- "Cath"
AP[which(str_sub(evo$rate, start = -3L, end = -2L) == "11")] <- "Diur"

state <- str_sub(evo$rate, start = -1L, end = -1L)
}

evo <- cbind(evo,type, modtype, rtype, AP, state)
####

## removing transitions between states that do not exist in the model
tri_state <- which(str_sub(evo$model, start = -10L, end = -10L) == 3) # 3-state models
quad_state <- which(str_sub(evo$model, start = -10L, end = -10L) == 4)# 4-state models
hi_states <- which(evo$state %in% c("D","E","F","G","H"))   # transitions outside states A-C
vhi_states <- which(evo$state %in% c("E","F","G","H"))  # transitions outside states A-D
remove <- c(intersect(tri_state, hi_states), intersect(quad_state, vhi_states))
evo <- evo[-remove,]
rm(tri_state, quad_state, hi_states, vhi_states, remove)

# plot rate estimates
# This code loads the function in the working environment
#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Output")

tb1 <- evo[which(evo$model != "vrMuHiSSE4_MCCtree"),] # exclude MCCtree
tb2 <- tb1[which(tb1$modtype == 3),] # only 3-state models
tb3 <- tb1[which(tb1$modtype == 4),] # only 4-state models

#png(file="Evolutionary rates estimate.png", width = 1600, height = 1200)
pdf(file="Evolutionary rates estimate.pdf", width = 14, height = 8)
ggplot(tb1, aes(x = rtype, y = value)) +
    geom_flat_violin(aes(fill = AP), position = position_nudge(x = .15, y = 0), alpha = .8) +
    geom_point(aes(fill = AP), position = position_jitter(width = .1), shape = 21, size = 3, alpha = 0.6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.7) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    facet_wrap(~modtype) + 
    scale_fill_manual(values = cols <- c("green3","gold2","blue")) +
    # the below line is used for coloring boxplots by AP
    #scale_colour_manual(values = c("green3","gold2","blue","blue")) + # no idea why 4 colours are needed, but 3 don't work
    guides(fill = FALSE, color = FALSE)

ggplot(tb1, aes(x = rtype, y = value)) +
    geom_flat_violin(aes(fill = AP), position = position_nudge(x = .15, y = 0), alpha = .8) +
    geom_point(aes(fill = AP), position = position_jitter(width = .1), shape = 21, size = 3, alpha = 0.6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.7) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    facet_wrap(~modtype) + 
    scale_fill_manual(values = cols <- c("green3","gold2","blue")) +
    guides(fill = FALSE, color = FALSE)
dev.off()
    
tb1 %>% group_by(state) %>% tally()


### 2.3.2 Transition rates ----------------------------------------------------
tran$model <- str_replace(rownames(tran), ".RDS", "") # remove file suffix from model names

# melt into long form removing the (wrongly) averaged values
Trates <- as_tibble(tran, subset = rownames(tran) != "averaged_24trees") %>%
    pivot_longer(!model, names_to = "rate", values_to = "value") 

# removing transitions between states that do not exist in the model
tri_state <- which(str_sub(Trates$model, start = -10L, end = -10L) == 3) # 3-state models
quad_state <- which(str_sub(Trates$model, start = -10L, end = -10L) == 4)# 4-state models
hi_states <- which(!is.na(str_extract(Trates$rate, "[D-H]")))   # transitions outside states A-C
vhi_states <- which(!is.na(str_extract(Trates$rate, "[E-H]")))  # transitions outside states A-D
remove <- c(intersect(tri_state, hi_states), intersect(quad_state, vhi_states))
Trates <- Trates[-remove,]

# adding useful information
modtype <- rtype <- rclass <- character()
for (i in 1:nrow(Trates)) {
    prefix <- str_split(Trates$model[i], "Mu")[[1]][1]
    num <- str_extract(str_split(Trates$model[i], "_")[[1]][1], "[0-9]") # get # of states
    modtype[i] <- paste0(prefix,num)
    rtype[i] <- gsub("[A-H]_", "->", str_sub(Trates$rate[i], start = 2L, end = -2L)) # transition type
    if (rtype[i] %in% c("00->00","01->01","11->11")) {
        rclass[i] <- "hidden change"
    } else {
        rclass[i] <- "character change"
    }
}
Trates <- cbind(Trates, modtype, rtype, rclass)
Trates$rtype <- sapply(Trates$rtype, gsub, pattern = "00", replacement = "N")
Trates$rtype <- sapply(Trates$rtype, gsub, pattern = "01", replacement = "C") 
Trates$rtype <- sapply(Trates$rtype, gsub, pattern = "11", replacement = "D")
rm(modtype, rtype, rclass, tri_state, quad_state, hi_states, vhi_states, remove)

tbr1 <- Trates[which(Trates$model != "vrMuHiSSE4_MCCtree"),] # exclude MCCtree
tbr2 <- tbr1[which(tbr1$rclass == "character change"),]        # exclude hidden rates
tbr3 <- tbr2[which(tbr2$modtype == 3),]
tbr4 <- tbr2[which(tbr2$modtype == 4),]
tbr5 <- tbr2
tbr5[which(tbr5$value < 10^-6),3] <- 10^-6

# plot transition rate estimates
pdf(file="Transition rates estimate.pdf", width = 10, height = 7)
ggplot(transform(tbr2, rtype=factor(rtype, levels=c("N->C","C->N","C->D","D->C"))),
       aes(x = rtype, y = log(value), fill = as.factor(rtype))) +
    #ggplot(tbr2, aes(x = rtype, y = log(value), fill = as.factor(rtype))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(colour = "grey20", position = position_jitter(width = .1), shape = 21, size = 1.5, alpha = 0.8) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    geom_hline(aes(yintercept = log(10^-6), colour = 'red'), linetype = 2) + # this marks the effective rate of 0, because 10^-6 transitions per lineage per million years is 1 transition per 10000 lineages per 100 million years, which is unlikely to be detected when I examine <2500 lineages over 160 million years
    theme_light() +
    #facet_wrap(~modtype) +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1)) +
    scale_fill_manual(values = cols <- c("#5A4A6F","#E47250","gold1","dodgerblue3")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    labs(y = "log(transition rate)", x = "Transition type (from -> to)") +
    guides(fill = FALSE, color = FALSE) #+
    # Setting the limits of the y axis
    #scale_y_continuous(limits = c(0, .03)) 
dev.off()

####################
** TODO:
    looking what happens when I round down all rates smaller than 10^-6 to zero
###################


# looking what happens when I omit worst performing models (based on BIC)
mods <- unique(Trates$model)
for (mod in unique(Trates$model)){
    if(str_sub(mod, start = -4L, end = -1L) %in% c("8981","5221","0028","5455","2966","5333","3865")){
        mods <- mods[-which(mods == mod)]
    }
}
NEWrates <- Trates[which(Trates$model %in% mods),]
nr1 <- NEWrates[which(NEWrates$model != "vrMuHiSSE4_MCCtree"),] # exclude MCCtree
nr2 <- nr1[which(nr1$rclass == "character change"),]  

TODO:
    ##  *** try MarginReconMuHiSSE() and use the output for 
    ##      GetModelAveRates()
    ##      for plot() (plot.muhisse.states)
    ##      and for SupportRegionMuHiSSE()
    
    
## summary of estimated rates
library(data.table)
rates <- as.data.table(mod_avg[,-c(1:11)])
mods <- NA
for (i in 1:nrow(mod_avg)){mods <- c(mods, strsplit(rownames(mod_avg)[i],"[.]")[[1]][1])}
#names(mods) <- "model"
div_rate <- cbind(mods[-which(is.na(mods))], rates)
colnames(div_rate)[1] <- "model"
# convert to long form 
all_rates <- melt(div_rate, id.vars = c("model","npar"), measure.vars = colnames(div_rate)[3:ncol(div_rate)], variable.name = "rate")
# aggregate rates over states by removing state name (capital letter) from rate name
#colnames(div_rate) <- str_sub(colnames(div_rate),1, nchar(colnames(div_rate))-1)




# 3. Follow up analyses -------------------------------------------------------
library("hisse")
# 3.1 Reconstruction based on fitted models -----------------------------------

mod <- readRDS(file="MuHiSSE2_tree6404.RDS")
parest <- SupportRegionMuHiSSE(mod, n.points=1000, scale.int=0.1, desired.delta=2, min.number.points=10, verbose=TRUE)

reco <- MarginReconMuHiSSE(tree, act3, f = freq, pars, hidden.states=1, 
                           condition.on.survival=TRUE, root.type="madfitz", 
                           root.p=NULL, AIC=NULL, get.tips.only=FALSE, 
                           verbose=TRUE, n.cores=NULL, dt.threads=1)



# 5. Processing binary models -------------------------------------------------

## 5.1 Any daytime activity ---------------------------------------------------
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster/Binary/N-CD")

## repeating stages 1-2 for binary models:

allmods <- list.files(pattern = ".RDS")

### 5.1.1 Model fit -----------------------------------------------------------

modat <- data.frame(t(sapply(allmods, extRDS)))
score <- as.data.frame(do.call(cbind, lapply(modat[,1:3], as.numeric)))
rownames(score) <- rownames(modat)
score$BIC <- NA
score$type <- NA
score$states <- NA
score$npar <- NA
score$tree <- NA

# extract model parameters from file name
for (i in 1:nrow(score)){
    score$type[i] <- gsub("[0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1])       # model type
    score$states[i] <- as.numeric(gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1]))    # no. of hidden states
    if (score$type[i] == "nullbin") { 
        score$type[i] <- "bvrCID"
        score$states[i] <- 1
        score$npar[i] <- 4              # no. of parameters calculated in "Evolutionary models.xlsx"
    } else if (score$type[i] == "BiSSE") { 
        score$type[i] <- "vrHiSSE"
        score$states[i] <- 1
        score$npar[i] <- 5
    } else if (score$type[i] == "bvrCID") { 
        score$npar[i] <- (score$states[i]+1)^2
    } else if (score$type[i] == "vrHiSSE") { 
        score$npar[i] <- score$states[i]*(score$states[i]+3)+1
    } 
    # tree variant index
    if (strsplit(rownames(score[i,]),"_")[[1]][2] == "MCCtree.RDS") {
        score$tree[i] <- "MCC"
    } else {
        score$tree[i] <- gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][2])}
} 

# calculate BIC [the number of data points (n=2400) should change in clade-wise models
for (i in 1:nrow(score)){
    score$BIC[i] <- -2*score$loglik[i] + log(2400) * score$npar[i]}


## model selection by AICc and BIC
binordA <- score[order(score$AICc),]
binordA <- binordA[,c("tree","type","states","loglik","npar","AIC","AICc","BIC")]
binordB <- score[order(score$BIC),]
binordB <- binordB[,c("tree","type","states","loglik","npar","AIC","AICc","BIC")]

ordbyt <- data.frame()
# calculate delta AICc/BIC for each tree
for (z in unique(binordA$tree)) {
    dat <- binordA[which(binordA$tree == z),]
    dat <- dat[order(dat$AICc),]
    dat$diff_AIC <- dat$AICc - min(dat$AICc)
    dat$rel_lik <- exp(-0.5 * dat$diff_AIC)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik) # this is identical to rel_lik on most trees
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_AIC, decreasing = FALSE),])
}
binordbytA <- ordbyt[order(ordbyt$tree, ordbyt$AICc),]

# and again for BIC:
ordbyt <- data.frame()
# calculate relative model support for each tree
for (z in unique(binordB$tree)) {
    dat <- binordB[which(binordB$tree == z),]
    dat <- dat[order(dat$BIC),]
    dat$diff_BIC <- dat$BIC - min(dat$BIC)
    dat$rel_lik <- exp(-0.5 * dat$diff_BIC)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik)
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_BIC, decreasing = FALSE),])
}
binordbytB<- ordbyt[order(ordbyt$tree, ordbyt$BIC),]


### 5.1.2 Support summary plots -----------------------------------------------
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Output/N-CD")
#write.csv(score, file="score.csv", row.names = FALSE)
#write.csv(ordbyt, file="Model fit by tree.csv")

## plotting likelihood as function of tree variant and model type

# comparing all models by BIC, divided to model type and no. of states, grouped by tree
pdf(file="Model support ordered_MCC marked.pdf", width = 10, height = 8)
#df1 <- score[score$states != 1,]
{df2 <- subset(score, score$tree =="MCC")    # get only MCC models separately
#par(mfrow=c(2,1))
ggplot(score, aes(x = reorder(tree, AICc, FUN = max), y = AICc, colour = as.factor(states))) +
    geom_point(alpha= .5, size = 2.5) + 
    geom_point(data = df2, aes(tree, AICc), colour = 'red', shape = 1, size = 3) +
    scale_y_reverse() + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("tree variant") +
    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)

ggplot(score, aes(x = reorder(tree, BIC, FUN = max), y = BIC, colour = as.factor(states))) +
    geom_point(alpha= .5, size = 2.5) + 
    geom_point(data = df2, aes(tree, BIC), colour = 'red', shape = 1, size = 3) + 
    scale_y_reverse() +
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("tree variant") +
    scale_color_viridis_d(option = "inferno", begin = 0, end = .9)}
dev.off()

#pdf(file = "Model support by no. of states.pdf", width = 7, height = 8.5)
#library(ggridges)

# plotting support by no. states, aggregaetd over all variants
vars <- score[-which(score$tree == "MCC"),]
pdf(file = "Model support aggregated.pdf", width = 8.5, height = 6)
{ggplot(data = vars, aes(y = AICc, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.5) +
    geom_point(aes(y = AICc, x = as.factor(states), color = as.factor(states)), 
               position = position_jitter(width = .1), size = 2, alpha = 0.6) +
    xlab("No. of states") +
    facet_wrap(~type,) +
    expand_limits(x = 5.25) +
    guides(fill = FALSE, color = FALSE) + # don't show a legend
    #scale_color_brewer(palette = "Spectral") +
    #scale_fill_brewer(palette = "Spectral") +
    scale_fill_viridis_d(option = "viridis", begin = .3, end = .9, alpha=0.8, direction = -1) +
    scale_color_viridis_d(option = "viridis", begin = .3, end = .9, alpha=0.8, direction = -1) +
    #scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6E")) +
    #scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6E")) +
    theme_minimal() 

ggplot(data = vars, aes(y = BIC, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.5) +
    geom_point(aes(y = BIC, x = as.factor(states), color = as.factor(states)), 
               position = position_jitter(width = .1), size = 2, alpha = 0.6) +
    xlab("No. of states") +
    facet_wrap(~type,) +
    expand_limits(x = 5.25) +
    guides(fill = FALSE, color = FALSE) + # don't show a legend
    scale_fill_viridis_d(option = "inferno", begin = .3, end = .9, alpha=0.8, direction = -1) +
    scale_color_viridis_d(option = "inferno", begin = .3, end = .9, alpha=0.8, direction = -1) +
    theme_minimal()}
dev.off()

# plotting support by no. states with trendlines to show that the pattern is uniform accross trees
pdf(file = "Model support by state_trendlines.pdf", width = 7, height = 8.5)
{ggplot(vars,aes(x=states, AICc, colour = reorder(tree, AICc, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    # use the below if need to plot MCC along with tree variants
    #geom_point(data = tert[tert$tree == "MCC",], aes(states, AICc), colour = 'grey10', shape = 4, size = 2) + 
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    theme_minimal() +
    guides(colour=FALSE) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .5, direction = -1)

ggplot(vars,aes(x=states, BIC, colour = reorder(tree, BIC, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    theme_minimal() +
    guides(colour=FALSE) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .5, direction = -1)}
dev.off()

# comparing likelihood of all models, divided to model type and no. of states, grouped by tree
#pdf(file="likelihood by states_MCC marked.pdf", width = 10, height = 8)
#df1 <- tert[tert$states != 1,]
#df2 <- subset(df1, df1$tree =="MCC")    # get only MCC models separately
#
#ggplot(df1, aes(x = reorder(tree, loglik, FUN = max), y = loglik, colour = as.factor(states))) +
#    geom_point(alpha= .8, size = 2.5) + 
#    #geom_point(data = df2, colour = 'white', shape = 3, size = 2) + 
#    geom_point(data = df2, aes(tree, loglik), colour = 'red', shape = 1, size = 3) + 
#    facet_wrap(~type) + 
#    theme_light() +
#    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
#    labs(color = "States") +
#    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)
#dev.off()

# comparing model type performance
pdf(file="model type performance by tree.pdf", width = 10, height = 8)
{ggplot(binordbytA, aes(x=states, AICc, colour = type, group=type)) +
    geom_point(alpha= .8, size = 3) + 
    geom_line(size = 1.5, alpha = 0.8) +
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "viridis", #inferno, turbo, mako
                          begin = .2, end = .8, direction = -1)
ggplot(binordbytB, aes(x=states, BIC, colour = type, group=type)) +
    geom_point(alpha= .8, size = 3) +  
    geom_line(size = 1.5, alpha = 0.8) +
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "inferno", #inferno, turbo, mako
                          begin = .2, end = .8, direction = -1)}
dev.off()


## 5.2 Strict diurnality only --------------------------------------------------
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster/Binary/D-CN")

## repeating stages 1-2 for binary models:

allmods <- list.files(pattern = ".RDS")

### 5.2.1 Model fit -----------------------------------------------------------

modat <- data.frame(t(sapply(allmods, extRDS)))
score <- as.data.frame(do.call(cbind, lapply(modat[,1:3], as.numeric)))
rownames(score) <- rownames(modat)
score$BIC <- NA
score$type <- NA
score$states <- NA
score$npar <- NA
score$tree <- NA

# extract model parameters from file name
for (i in 1:nrow(score)){
    score$type[i] <- gsub("[0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1])       # model type
    score$states[i] <- as.numeric(gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][1]))    # no. of hidden states
    if (score$type[i] == "nullbin") { 
        score$type[i] <- "bvrCID"
        score$states[i] <- 1
        score$npar[i] <- 4              # no. of parameters calculated in "Evolutionary models.xlsx"
    } else if (score$type[i] == "BiSSE") { 
        score$type[i] <- "vrHiSSE"
        score$states[i] <- 1
        score$npar[i] <- 5
    } else if (score$type[i] == "bvrCID") { 
        score$npar[i] <- (score$states[i]+1)^2
    } else if (score$type[i] == "vrHiSSE") { 
        score$npar[i] <- score$states[i]*(score$states[i]+3)+1
    } 
    # tree variant index
    if (strsplit(rownames(score[i,]),"_")[[1]][2] == "MCCtree") {
        score$tree[i] <- "MCC"
    } else {
        score$tree[i] <- gsub("[^0-9]", "", strsplit(rownames(score[i,]),"_")[[1]][2])}
} 

# calculate BIC [the number of data points (n=2400) should change in clade-wise models
for (i in 1:nrow(score)){
    score$BIC[i] <- -2*score$loglik[i] + log(2400) * score$npar[i]}


## model selection by AICc and BIC
binordA <- score[order(score$AICc),]
binordA <- binordA[,c("tree","type","states","loglik","npar","AIC","AICc","BIC")]
binordB <- score[order(score$BIC),]
binordB <- binordB[,c("tree","type","states","loglik","npar","AIC","AICc","BIC")]

ordbyt <- data.frame()
# calculate delta AICc/BIC for each tree
for (z in unique(binordA$tree)) {
    dat <- binordA[which(binordA$tree == z),]
    dat <- dat[order(dat$AICc),]
    dat$diff_AIC <- dat$AICc - min(dat$AICc)
    dat$rel_lik <- exp(-0.5 * dat$diff_AIC)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik) # this is identical to rel_lik on most trees
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_AIC, decreasing = FALSE),])
}
binordbytA <- ordbyt[order(ordbyt$tree, ordbyt$AICc),]

# and again for BIC:
ordbyt <- data.frame()
# calculate relative model support for each tree
for (z in unique(binordB$tree)) {
    dat <- binordB[which(binordB$tree == z),]
    dat <- dat[order(dat$BIC),]
    dat$diff_BIC <- dat$BIC - min(dat$BIC)
    dat$rel_lik <- exp(-0.5 * dat$diff_BIC)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik)
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_BIC, decreasing = FALSE),])
}
binordbytB<- ordbyt[order(ordbyt$tree, ordbyt$BIC),]


### 5.2.2 Support summary plots -----------------------------------------------
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Output/D-CN")
#write.csv(score, file="score_D-CN.csv", row.names = FALSE)
#write.csv(binordbytA, file="Model AICc fit by tree.csv")
#write.csv(binordbytB, file="Model BIC fit by tree.csv")

## plotting likelihood as function of tree variant and model type

# comparing all models by BIC, divided to model type and no. of states, grouped by tree
pdf(file="Model support ordered_MCC marked.pdf", width = 10, height = 8)
#df1 <- score[score$states != 1,]

df2 <- subset(score, score$tree =="MCC")    # get only MCC models separately
#par(mfrow=c(2,1))
ggplot(score, aes(x = reorder(tree, AICc, FUN = max), y = AICc, colour = as.factor(states))) +
    geom_point(alpha= .5, size = 2.5) + 
    geom_point(data = df2, aes(tree, AICc), colour = 'red', shape = 1, size = 3) +
    scale_y_reverse() + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("tree variant") +
    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)

ggplot(score, aes(x = reorder(tree, BIC, FUN = max), y = BIC, colour = as.factor(states))) +
    geom_point(alpha= .5, size = 2.5) + 
    geom_point(data = df2, aes(tree, BIC), colour = 'red', shape = 1, size = 3) + 
    scale_y_reverse() +
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("tree variant") +
    scale_color_viridis_d(option = "inferno", begin = 0, end = .9)
dev.off()

#pdf(file = "Model support by no. of states.pdf", width = 7, height = 8.5)
#library(ggridges)

## code for this figure adapted from: https://ourcodingclub.github.io/tutorials/dataviz-beautification-synthesis/#distributions 
# This code loads the function in the working environment
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
# theme from: https://neuroconscience.wordpress.com/2018/03/15/introducing-raincloud-plots/

raincloud_theme = theme(
    text = element_text(size = 10),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title=element_text(size=16),
    legend.text=element_text(size=16),
    legend.position = "right",
    plot.title = element_text(lineheight=.8, face="bold", size = 16),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

# plotting support by no. states, aggregaetd over all variants
vars <- score[-which(score$tree == "MCC"),]
pdf(file = "Model support aggregated.pdf", width = 8.5, height = 6)
ggplot(data = vars, aes(y = AICc, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.5) +
    geom_point(aes(y = AICc, x = as.factor(states), color = as.factor(states)), 
               position = position_jitter(width = .1), size = 2, alpha = 0.6) +
    xlab("No. of states") +
    facet_wrap(~type,) +
    expand_limits(x = 5.25) +
    guides(fill = FALSE, color = FALSE) + # don't show a legend
    #scale_color_brewer(palette = "Spectral") +
    #scale_fill_brewer(palette = "Spectral") +
    scale_fill_viridis_d(option = "viridis", begin = .3, end = .9, alpha=0.8, direction = -1) +
    scale_color_viridis_d(option = "viridis", begin = .3, end = .9, alpha=0.8, direction = -1) +
    #scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6E")) +
    #scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6E")) +
    theme_minimal() 

ggplot(data = vars, aes(y = BIC, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.5) +
    geom_point(aes(y = BIC, x = as.factor(states), color = as.factor(states)), 
               position = position_jitter(width = .1), size = 2, alpha = 0.6) +
    xlab("No. of states") +
    facet_wrap(~type,) +
    expand_limits(x = 5.25) +
    guides(fill = FALSE, color = FALSE) + # don't show a legend
    scale_fill_viridis_d(option = "inferno", begin = .3, end = .9, alpha=0.8, direction = -1) +
    scale_color_viridis_d(option = "inferno", begin = .3, end = .9, alpha=0.8, direction = -1) +
    theme_minimal()
dev.off()

# plotting support by no. states with trendlines to show that the pattern is uniform accross trees
pdf(file = "Model support by state_trendlines.pdf", width = 7, height = 8.5)
ggplot(vars,aes(x=states, AICc, colour = reorder(tree, AICc, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    # use the below if need to plot MCC along with tree variants
    #geom_point(data = tert[tert$tree == "MCC",], aes(states, AICc), colour = 'grey10', shape = 4, size = 2) + 
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    theme_minimal() +
    guides(colour=FALSE) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .5, direction = -1)

ggplot(vars,aes(x=states, BIC, colour = reorder(tree, BIC, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    theme_minimal() +
    guides(colour=FALSE) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .5, direction = -1)
dev.off()

# comparing likelihood of all models, divided to model type and no. of states, grouped by tree
#pdf(file="likelihood by states_MCC marked.pdf", width = 10, height = 8)
#df1 <- tert[tert$states != 1,]
#df2 <- subset(df1, df1$tree =="MCC")    # get only MCC models separately
#
#ggplot(df1, aes(x = reorder(tree, loglik, FUN = max), y = loglik, colour = as.factor(states))) +
#    geom_point(alpha= .8, size = 2.5) + 
#    #geom_point(data = df2, colour = 'white', shape = 3, size = 2) + 
#    geom_point(data = df2, aes(tree, loglik), colour = 'red', shape = 1, size = 3) + 
#    facet_wrap(~type) + 
#    theme_light() +
#    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
#    labs(color = "States") +
#    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)
#dev.off()

# comparing model type performance
pdf(file="model type performance by tree.pdf", width = 10, height = 8)
ggplot(binordbytA, aes(x=states, AICc, colour = type, group=type)) +
    geom_point(alpha= .8, size = 3) + 
    geom_line(size = 1.5, alpha = 0.8) +
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "viridis", #inferno, turbo, mako
                          begin = .2, end = .8, direction = -1)
ggplot(binordbytB, aes(x=states, BIC, colour = type, group=type)) +
    geom_point(alpha= .8, size = 3) +  
    geom_line(size = 1.5, alpha = 0.8) +
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "inferno", #inferno, turbo, mako
                          begin = .2, end = .8, direction = -1)
dev.off()

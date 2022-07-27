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

# place all .RDS files in one folder to eliminate this duplicity
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster")
#setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/Additional")

# extract the first 3 elements of each model into a data frame
allmods <- list.files(pattern = ".RDS")


# 1. Model fit ----------------------------------------------------------------

modat <- data.frame(t(sapply(allmods, extRDS)))
score <- as.data.frame(do.call(cbind, lapply(modat[,1:3], as.numeric)))
rownames(score) <- rownames(modat)
score$tree <- score$npar <- score$states <- score$type <- score$BIC <- NA


## 1.1 Extract model parameters -----------------------------------------------
# extracting model configuration from file name
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

# adjust info for single state models (MuSSE and CID-1)
score$states[which(score$type %in% c('MuSSE','null'))] <- 1  
score$npar[which(score$type == 'MuSSE')] <- 8
score$type[which(score$type == 'MuSSE')] <- 'MuHiSSE'
score$npar[which(score$type == 'null')] <- 6
score$type[which(score$type == 'null')] <- 'CID'

# calculate BIC [the number of data points (n=2400) should change in clade-wise models
for (i in 1:nrow(score)){
    score$BIC[i] <- -2*score$loglik[i] + log(2400) * score$npar[i]}

tert <- score[which(score$type %in% c("CID","MuHiSSE","vrCID","vrMuHiSSE")),]

# finding which trees don't have MuHiSSE4 results (timed out on cluster)
#Mu4 <- score[which(score$type == "vrMuHiSSE" & score$states == 4),]
#Mu3 <- score[which(score$type == "vrMuHiSSE" & score$states == 3),]
#Mu3$tree[which(!Mu3$tree %in% Mu4$tree)]
#[1] "4420" "6747" "8394"

## 1.2 Model selection ---------------------------------------------------------

ord <- tert

ordbyt <- data.frame()
# calculate relative model support (delta AICc) for each tree
for (z in unique(ord$tree)) {
    dat <- ord[which(ord$tree == z),]
    #dat <- dat[order(dat$AICc),]
    dat$diff_AICc <- dat$AICc - min(dat$AICc)
    dat$rel_lik <- exp(-0.5 * dat$diff_AICc)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik) # this is identical to rel_lik on most trees
    dat$diff_AICc <- dat$AICc - min(dat$AICc)
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_AICc, decreasing = FALSE),])
}
ordbytA <- ordbyt[order(ordbyt$tree, ordbyt$AICc),]

# and again for BIC
ordbyt <- data.frame()
for (z in unique(ord$tree)) {
    dat <- ord[which(ord$tree == z),]
    dat <- dat[order(dat$BIC),]
    dat$diff_BIC <- dat$BIC - min(dat$BIC)
    dat$rel_lik <- exp(-0.5 * dat$diff_BIC)
    dat$weight <- dat$rel_lik / sum(dat$rel_lik)
    ordbyt <- rbind(ordbyt, dat[order(dat$diff_BIC, decreasing = FALSE),])
}
ordbytB <- ordbyt[order(ordbyt$tree, ordbyt$BIC),]


## 1.3 Support summary plots --------------------------------------------------

setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Output")
#write.csv(score, file="score.csv", row.names = FALSE)
#write.csv(ordbyt, file="Model fit by tree.csv")

## plotting likelihood as function of tree variant and model type

# comparing all models by BIC, divided to model type and no. of states, grouped by tree
pdf(file="1b- Model support all_ordered_MCC marked_01JUN.pdf", width = 15, height = 10)
df1 <- tert[tert$states != 1,]
df2 <- subset(df1, df1$tree =="MCC")    # get only MCC models separately
#par(mfrow=c(2,1))
ggplot(df1, aes(x = reorder(tree, AICc, FUN = min), y = AICc, colour = as.factor(states))) +
    geom_point(alpha= .75, size = 3) + 
    geom_point(data = df2, aes(tree, AICc), colour = 'grey45', shape = 1, size = 3) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 5)) +
    labs(color = "States") +
    xlab("Tree variant (reordered for visualisation)") +
    scale_color_viridis_d(option = "viridis", begin = .15, end = 1, alpha = 1, direction = -1) +
    scale_fill_viridis_d(option = "viridis", begin = .1, end = .9, alpha = 0.5, direction = -1)

ggplot(df1, aes(x = reorder(tree, BIC, FUN = min), y = BIC, colour = as.factor(states))) +
    geom_point(alpha= .75, size = 3) + 
    geom_point(data = df2, aes(tree, BIC), colour = "grey75", shape = 1, size = 3) + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 5)) +
    labs(color = "States") +
    xlab("Tree variant (reordered for visualisation)") +
    scale_color_viridis_d(option = "inferno", begin = .1, end = .9, alpha = 1) +
    scale_fill_viridis_d(option = "inferno", begin = .1, end = .9, alpha = 0.5) #+ 
    # add arrow to mark MCC tree
    #geom_segment(aes(x = "MCC", y = 17300, xend = "MCC", yend = 17700), colour = 1,
    #             arrow = arrow(length = unit(0.5, "cm"), type = "open"))
dev.off()

# same plot as above but trees not ordered
#df3 <- tert[tert$states == 1,]
#pdf(file = "Model support unordered.pdf",width = 10, height = 10)  #png(file = "likelihood by model type.png",width = 7, height = 8.5, units = "in", res=800)
#ggplot(tert) +
#    geom_point(aes(x=tree, AICc, colour = as.factor(states)), alpha= .8, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
#    geom_point(data = df3, aes(tree, AICc), colour = 'grey20', shape = 1, size = 2) + 
#    theme_light() +
#    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
#    labs(color = "States") +
#    facet_wrap(~type) + #, nrow = 1) + 
#    scale_color_viridis_d(option = "viridis", begin = .1, end = .9, direction = -1)
#
#ggplot(tert) +
#    geom_point(aes(x=tree, BIC, colour = as.factor(states)), alpha= .8, size = 2) + # use: "reorder(tree, -loglik, FUN = mean)" instead of "tree" to order trees by loglik
#    geom_point(data = df3, aes(tree, BIC), colour = 'grey20', shape = 1, size = 2) + 
#    theme_light() +
#    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
#    labs(color = "States") +
#    facet_wrap(~type) + #, nrow = 1) + 
#    scale_color_viridis_d(option = "inferno", begin = .1, end = .9, direction = -1)
#dev.off()

#pdf(file = "Model support by no. of states.pdf", width = 7, height = 8.5)
library(ggridges)

## code for this figure adapted from: https://ourcodingclub.github.io/tutorials/dataviz-beautification-synthesis/#distributions 
# This code loads the function in the working environment
#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
# theme from: https://neuroconscience.wordpress.com/2018/03/15/introducing-raincloud-plots/
#theme_niwot <- function(){
#    theme_bw() +
#        theme(text = element_text(family = "Helvetica Light"),
#              axis.text = element_text(size = 16),
#              axis.title = element_text(size = 18),
#              axis.line.x = element_line(color="black"),
#              axis.line.y = element_line(color="black"),
#              panel.border = element_blank(),
#              panel.grid.major.x = element_blank(),
#              panel.grid.minor.x = element_blank(),
#              panel.grid.minor.y = element_blank(),
#              panel.grid.major.y = element_blank(),
#              plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
#              plot.title = element_text(size = 18, vjust = 1, hjust = 0),
#              legend.text = element_text(size = 12),
#              legend.title = element_blank(),
#              legend.position = c(0.95, 0.15),
#              legend.key = element_blank(),
#              legend.background = element_rect(color = "black",
#                                               fill = "transparent",
#                                               size = 2, linetype = "blank"))
#}
#    
#raincloud_theme = theme(
#    text = element_text(size = 10),
#    axis.title.x = element_text(size = 16),
#    axis.title.y = element_text(size = 16),
#    axis.text = element_text(size = 14),
#    axis.text.x = element_text(angle = 45, vjust = 0.5),
#    legend.title=element_text(size=16),
#    legend.text=element_text(size=16),
#    legend.position = "right",
#    plot.title = element_text(lineheight=.8, face="bold", size = 16),
#    panel.border = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.grid.major = element_blank(),
#    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
#    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

## plotting support by no. states, aggregaetd over all variants
varsA <- ordbytA[which(ordbytA$tree != "MCC"),]
varsB <- ordbytB[which(ordbytB$tree != "MCC"),]
vrA <- varsA[which(varsA$type %in% c("vrCID", "vrMuHiSSE")),]
vrB <- varsB[which(varsB$type %in% c("vrCID", "vrMuHiSSE")),]

# the zero values mess with the distribution curve so I plot them afterwards
vrB2 <- filter(varsB, diff_BIC != 0)
vrB4 <- filter(varsB, diff_BIC == 0)

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
   
pdf(file = "1- VR relative model support by AICc_light_01JUN.pdf", width = 12, height = 8)
ggplot(data = vrA, aes(y = diff_AICc, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +   #  adjust = 0.5 controls curve smoothing
    geom_point(aes(y = diff_AICc, x = as.factor(states), color = as.factor(states)), 
               position = position_jitter(width = .1), size = 2, alpha = 0.6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    guides(fill = "none", color = "none") + # don't show a legend
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
    theme_light(base_size = 24) 
dev.off()

## alternative plotting function because geom_flat_violin() was having issues in the panel of best supported models
pdf(file = "2- Relative model support by BIC_light_30JUN.pdf", width = 24, height = 8)
ggplot(varsB, aes(x = diff_BIC, y = as.factor(states), fill = as.factor(states))) +
    geom_density_ridges(position = position_nudge(y = 0.14), scale = 0.7, 
                        rel_min_height = 0.01, bandwidth = 7) +
    geom_point(aes(), position = position_jitter(height = .05), size = .85, alpha = 0.4, shape = 16) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.7) +
    guides(fill = "none", color = "none") + # don't show a legend
    scale_fill_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 0.8, direction = -1) +
    scale_color_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 1, direction = -1) +
    facet_wrap(~type, nrow = 1) + 
    ylab("No. of states") +
    coord_flip() +
    theme_light(base_size = 28)
dev.off()

# No. of estimated parameters instead of No. of states
pdf(file = "2.1- Relative model support by Est_Params_30JUN.pdf", width = 20, height = 18)
ggplot(varsB, aes(x = diff_BIC, y = npar, fill = as.factor(states))) +
    geom_density_ridges(position = position_nudge(y = 0.05), scale = 1.2, 
                        rel_min_height = 0.01, bandwidth = 8) +
    geom_point(position = position_jitter(height = .03), size = 0.7, alpha = 0.4, shape = 16) +
    guides(fill = "none", color = "none") + # don't show a legend
    scale_fill_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 0.8, direction = -1) +
    scale_color_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 1, direction = -1) +
    facet_wrap(~type, nrow = 4) + 
    coord_flip() +
    ylab("No. of estimated parameters") +
    theme_light(base_size = 26)
dev.off()

pdf(file = "2a- Relative model support by BIC_light_07JUN.pdf", width = 12, height = 8)
ggplot(data = varsB, aes(x = diff_BIC, y = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +  ##  adjust = 0.5 controls curve smoothing
    #geom_density_ridges(aes(x = diff_BIC, y = as.factor(states), colour = as.factor(states)),
    #                    jittered_points = TRUE, 
    #                    position = position_raincloud(adjust_vlines = TRUE, height = -0.1),
    #                    rel_min_height = 0.01, scale = 0.7, bandwidth = 7, size = 0.1, alpha = 0.7) +
    geom_point(aes(x = diff_BIC, y = as.factor(states), color = "grey20"), 
               position = position_jitter(height = .05), size = 1.2, shape = 21, alpha = 0.6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    guides(fill = "none", color = "none") + # don't show a legend
    ylab("No. of states") +
    xlab("Difference in BIC") +
    facet_wrap(~type) +
    coord_flip() +
    expand_limits(x = 5.25) +
    scale_fill_viridis_d(option = "inferno", begin = 0.1, end = .95, alpha = 0.8, direction = -1) +
    scale_color_viridis_d(option = "inferno", begin = 0.1, end = .95, alpha = 1, direction = -1) +
    theme_light(base_size = 24) 
dev.off()

# models with diff_BIC = 0 not included in density curves 
pdf(file = "2b- Relative model support by BIC_light_07JUN.pdf", width = 12, height = 9)
p <- ggplot(data = vrB2, aes(y = diff_BIC, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) + #  adjust = 0.5 controls curve smoothing
    geom_point(aes(y = diff_BIC, x = as.factor(states), fill = as.factor(states)), 
               position = position_jitter(width = .03), size = 1, alpha = 0.5, shape = 21) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    guides(fill = "none", color = "none") + # don't show a legend
    xlab("No. of states") +
    #facet_wrap(~type, scales = "fixed", nrow = 2) +
    facet_grid(. ~ type, scales = "free", space='free') +
    expand_limits(x = 5.25) +
    scale_fill_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 0.8, direction = -1) +
    scale_color_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 1, direction = -1) +
    theme_light(base_size = 24)
p +  geom_point(data = vrB4, aes(y = diff_BIC, x = as.factor(states), color = as.factor(states)),
                position = position_jitter(width = .1), size = 1.5, alpha = 0.6, col = "#00BDFF") 
dev.off()



# plotting support by no. states with trendlines to show that the pattern is uniform accross trees
pdf(file = "1c- Model support by states tree-wise trendlines_01_JUN.pdf", width = 7, height = 8.5)
ggplot(varsA,aes(x=states, AICc, colour = reorder(tree, AICc, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    # use the below if need to plot MCC along with tree variants
    #geom_point(data = tert[tert$tree == "MCC",], aes(states, AICc), colour = 'grey10', shape = 4, size = 2) + 
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    guides(fill = "none", color = "none") + # don't show a legend
    #labs(color = "Tree variant") +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .5, direction = -1)

ggplot(varsB,aes(x=states, BIC, colour = reorder(tree, BIC, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    # use the below if need to plot MCC along with tree variants
    #geom_point(data = tert[tert$tree == "MCC",], aes(states, BIC), colour = 'grey10', shape = 4, size = 2) + 
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    guides(fill = "none", color = "none") + # don't show a legend
    #labs(color = "Tree variant") +
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
pdf(file="2- model type performance by tree_AICc_01JUN.pdf", width = 15, height = 10)
ggplot(varsA, aes(x=states, AICc, colour = type, group=type)) +
    geom_point(alpha= .5, size = 3) + 
    geom_line(size = 1.2, alpha = 0.85) +
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "viridis", #inferno, turbo, mako
                          begin = 0, end = .95, direction = -1)
dev.off()

pdf(file="2- model type performance by tree_BIC_01JUN.pdf", width = 15, height = 10)
#pdf(file="FigS1 - model type performance by tree_BIC_04JUN.pdf", width = 15, height = 10)
ggplot(varsB, aes(x=states, BIC, colour = type, group=type)) +
    geom_point(alpha= .5, size = 3) +  
    geom_line(size = 1.2, alpha = 0.85) +
    theme_light() +
    facet_wrap(~tree, nrow = 5) + 
    scale_color_viridis_d(option = "inferno", #inferno, turbo, mako
                          begin = 0, end = .9, direction = -1)
dev.off()


# same as above but trees ordered by BIC (this code works, but tree labels are gone)
srd <- ordbytB[order(ordbytB$BIC),] 
sru <- unique(srd$tree)
names(sru) <- sru
#   for(i in 1:nrow(ordbytB)){ordbytB$ordBIC[i] <- which(sru == ordbytB$tree[i])}
#   
#   ggplot(ordbytB, aes(x=states, BIC, colour = type, group=type)) +
#       geom_point(alpha= .5, size = 3) +  
#       geom_line(size = 1.5, alpha = 0.5) +
#       theme_light() +
#       facet_grid(ordBIC ~ ., labeller = labeller(sru)) +
#       facet_wrap(~ordBIC, nrow = 6) + 
#       scale_color_viridis_d(option = "inferno", #inferno, turbo, mako
#                             begin = 0, end = .9, direction = -1)


## quick tally - which trees support MuHiSSE3 better than MuHiSSE2 and MuHiSSE4
## all supported models (BIC <6)
trees <- ordbytB[which(ordbytB$diff_BIC < 6),]
trees <- trees[-which(trees$tree == "MCC"),]
trees$best <- paste0(trees$type, trees$states)

## only one best model per tree
#   trees <- data.frame(sru)
#   trees$best <- NA
#   trees <- trees[-which(trees$sru == "MCC"),]
#   # make a table of all the best models
#   for(t in unique(ordbytB$tree)){
#       mods <- ordbytB[which(ordbytB$tree == t),]
#       best <- mods[which(mods$BIC == min(mods$BIC)),]
#       trees$best[which(trees$sru == t)] <- paste0(best$type,best$states)
#   }

best <- table(trees$best)

pdf(file="Figure 2- best configuration by BIC_04JUL.pdf", width = 7, height = 4)
# Plot, and store x-coordinates of bars in xx
xx <- barplot(best, xaxt = 'n', xlab = '', width = 0.85, ylim = c(0, max(best)+15),
              ylab = "No. of trees") #main = "Support for model configuration",
# Add text at top of bars
text(x = xx, y = best, label = best, col = "darkred", pos = 3, font = 2)
# Add x-axis labels 
axis(1, at = xx, labels = names(best), tick = FALSE, line = -0.5, cex = 0.5)
dev.off()

# 2. Parameter estimates ------------------------------------------------------
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster")

pardat <- data.frame(t(sapply(allmods, extRDS_par)))
params <- pardat[which(is.na(str_extract(names(pardat), ".10.")))]      # remove state "10" (disabled in all models)
params <- params[which(is.na(str_extract(names(params), ".00..11.|.11..00.")))]  # remove 'diagonal' change (disabled in all)


## 2.1 calculate s and mu from composite rates --------------------------------

#   reminder: turnover = spec+ext; netDiv = spec-ext; extFrac = ext/spec 
#   symbols: tau-turnover; eps-extFrac; s-speciation; mu-extinction; lambda-netDiv;
#   Developing two equations w two unknowns on paper, got this:
#
#  // s = tau/(1+eps) \\        # speciation rate
#  \\ mu = tau-s      //        # extinction rate

div <- params           # create a mirror table to write s and mu values into
div[,which(is.na(str_extract(names(params), "q.")))] <- NA  # remove transition rates
# iterate over all models and calculate s and mu for each
for (i in colnames(params)) {                       # iterate over each rate column
    if (is.na(str_extract(i, "turn.")) == FALSE) {  # if it's a turnover rate
        index <- strsplit(i, "over")[[1]][2]        # get the rate name
        corr_rate <- which(colnames(params) == paste0("eps", index)) # get corresponding epsilon rate
        for (j in rownames(params)) {               # calculate speciation rate    
            div[j,i] <- params[j,i] / (1 + params[j,corr_rate]) 
        }
        colnames(div)[which(colnames(div) == i)] <- paste0("s",index) # update name to 'speciation'
    } else if (is.na(str_extract(i, "eps.")) == FALSE) { # if it's epsilon column
        index <- strsplit(i, "eps")[[1]][2]         # get the rate name
        corr_rate <- which(colnames(div) == paste0("s", index)) # get corresponding speciation rate
        for (j in rownames(params)) {               # calculate extinction rate
            div[j,i] <- params[j,i] * div[j,corr_rate]  
        }
        colnames(div)[which(colnames(div) == i)] <- paste0("mu",index) # update name to 'extinction'
    }
}
#setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Output")
#write.csv(div, file="Raw Rates of Evolution_125trees.csv", row.names = TRUE)


## 2.2 Parameter estimates from averaging supported models --------------------

# combine fit and params estimates to one table
results <- ordbytB
filler <- matrix(NA, nrow = nrow(ordbytB), ncol = ncol(div))
colnames(filler) <- colnames(div)
#div_s <- div[which(rownames(div) %in% rownames(ordbytB)),]
results <- cbind(ordbytB, filler)
for (i in 1:nrow(div)) {
    r <- which(rownames(ordbytB) == rownames(div)[i])
    results[r,] <- cbind(ordbytB[r,], div[i,])
}
#write.csv(results, file = "Model estimates exploratory.csv")
#write.csv(results, file = "Parameter estimates_125trees.csv")


### 2.2.1 Calculate diversification rates (lambda) ----------------------------

# take all BIC supported models (start from the previous step to repeat this process with AIC)
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
## I did not average competing models of the same tree variant because they are 
## mutually exclusive so not comparable. 

# get all trees with equivocal results
#equi <- leading[which(leading$tree %in% leading$tree[duplicated(leading$tree)]),]

# weighing estimated values by model support for all tree variants (not for MCC)
#mod_avg <- leading[which(leading$weight > 0.6),] # remove less-well supported trees
#mod_avg <- leading[which(leading$tree != "MCC"),]
#    
#for (i in 1:length(mod_avg$tree)) {
#    mod_avg[i,12:ncol(mod_avg)] <- mod_avg$weight[i] * mod_avg[i,12:ncol(mod_avg)]
#} 

# ** the below averaging is wrong and was NOT used but kept for generating ideas **
#means <- colSums(mod_avg[which(mod_avg$tree != "MCC"),12:ncol(mod_avg)]) / 24 # mean rates
#means <- c(rep(NA, 11), means)
#mod_avg <- rbind(mod_avg, means)
#rownames(mod_avg)[which(is.na(mod_avg$tree))] <- "averaged_24trees"
## remove rates not considered in the model (states F, G, H)
#mod_avg <- mod_avg[which(is.na(str_extract(colnames(mod_avg), "..E|..F|..G|..H")))]  
# ** end dodgy averaging procedure **

#write.csv(mod_avg, file="Best supported models.csv")


## 2.3 Rate estimates ---------------------------------------------------------

non_rate <- c("loglik","AIC","AICc","BIC","type","states","npar","tree","diff_BIC","rel_lik","weight","ordBIC")
## split into evolutionary rates and character transition rates
tran <- leading[,setdiff(colnames(leading), non_rate)]          # raw rate estimates
#tran <- mod_avg[,-c(1:11)]                                     # weighted estimates
tran <- tran[which(!is.na(str_extract(colnames(tran), "q")))]   # transition rates only

evol <- leading[,setdiff(colnames(leading), non_rate)]          # raw estimates
# evol <- mod_avg[,-c(1:11)]                                    # weigthed estimates
evol <- evol[which(is.na(str_extract(colnames(evol), "q")))]    # rates that are NOT transition rates


### 2.3.1 Evolutionary rates --------------------------------------------------
evol$model <- str_replace(rownames(evol), ".RDS", "") # remove file suffix from model names

# find and remove model/s based on MCC tree
for(i in 1:nrow(evol)){if(strsplit(evol$model[i], "_")[[1]][2] == "MCCtree"){mcc <- i}}
evol <- evol[-mcc,]


#### 2.3.1.1 Individual evolutionary rates ------------------------------------

evo <- pivot_longer(evol, cols=1:ncol(evol)-1, names_to = "rate", values_to = "value") %>%
    drop_na()

# additional info ('metadata') on each rate
type <- modtype <- process <- AP <- state <- character()

type[which(!is.na(str_extract(evo$rate, "s")))] <- "spec"
type[which(!is.na(str_extract(evo$rate, "mu")))] <- "ext"
type[which(!is.na(str_extract(evo$rate, "lambda")))] <- "net.div"

for (i in 1:nrow(evo)){
    modtype[i] <- str_extract(str_split(evo$model[i], "_")[[1]][1], "[0-9]")
    process[i] <- str_sub(evo$rate[i], start = 1L, end = -2L)}

AP[which(str_sub(evo$rate, start = -3L, end = -2L) == "00")] <- "Noct"
AP[which(str_sub(evo$rate, start = -3L, end = -2L) == "01")] <- "Cath"
AP[which(str_sub(evo$rate, start = -3L, end = -2L) == "11")] <- "Diur"

state <- str_sub(evo$rate, start = -1L, end = -1L)


evo <- cbind(evo, type, modtype, process, AP, state)
rm(type, modtype, process, AP, state)
####

## removing transitions between states that do not exist in the model
bi_state <- which(str_sub(evo$model, start = -10L, end = -10L) == 2)    # 2-state models
tri_state <- which(str_sub(evo$model, start = -10L, end = -10L) == 3)   # 3-state models
quad_state <- which(str_sub(evo$model, start = -10L, end = -10L) == 4)  # 4-state models
mid_states <- which(evo$state %in% LETTERS[3:8])        # transitions outside states A-B
hi_states <- which(evo$state %in% LETTERS[4:8])         # transitions outside states A-C
vhi_states <- which(evo$state %in% LETTERS[5:8])        # transitions outside states A-D
remove <- c(intersect(bi_state, mid_states), intersect(tri_state, hi_states), intersect(quad_state, vhi_states))
evo <- evo[-remove,]
rm(bi_state, tri_state, quad_state, mid_states, hi_states, vhi_states, remove)

# plot rate estimates
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Output")

# unaggregated rates
tb1 <- evo
tb2 <- evo[which(evo$modtype == 3),] # only 3-state models
tb3 <- evo[which(evo$modtype == 4),] # only 4-state models

#### removing outliers 
tb1[which(tb1$value > 1.5),] # identify trees w exceptionally high rate estimates
# find indices of all rates based on the problematic tree
d <- numeric()
for(i in 1:nrow(evo)){if(strsplit(evo$model[i], "_")[[1]][2] == "tree1803"){d <- c(d,i)}}
# plot all rates and highlight estimates from outlier-bearing models
#ggplot(tb1, aes(x = rtype, y = value)) +
#    geom_flat_violin(aes(fill = AP), position = position_nudge(x = .15, y = 0), alpha = .8) +
#    geom_point(aes(fill = AP), position = position_jitter(width = .1), shape = 21, size = 2, alpha = 0.2) + #
#    geom_boxplot(aes(fill = AP), width = .15, outlier.shape = NA, alpha = 0.6) +
#    theme_light() +
#    theme(axis.text.x = element_text(vjust = 1, hjust = 0.5)) + #angle = 80, 
#    #ylim(c(0,1)) +
#    #facet_wrap(~type) + 
#    scale_fill_manual(values = cols <- c("#33DD00","#FFCC00","#2233FF"))+#,"grey")) +
#    # the below line should be used for coloring boxplots by AP with the argument 'colour='
#    #scale_colour_manual(values = c("green3","gold2","dodgerblue3"))+#,"blue")) +
#    guides(fill = "none", color = "none") +
#    geom_point(data = tb1[d,], aes(x = rtype, y = value), colour = "red", size = 2)
tb1 <- tb1[-d,]

pdf(file="4b- Evol rates est_no_outlier_01JUN.pdf", width = 14, height = 8)
# reorder to put speciation before extinction
tb1$process <- factor(tb1$process, levels = c("s00","s01","s11","mu00","mu01","mu11","lambda00","lambda01","lambda11"))
ggplot(tb1, aes(x = process, y = value)) +
    #geom_density_ridges(position = position_nudge(y = 0.14), scale = 0.7, 
    #                    rel_min_height = 0.01, bandwidth = 7) +
    geom_flat_violin(aes(fill = AP), position = position_nudge(x = .15, y = 0))+#, alpha = .8) +
    geom_point(aes(fill = AP), position = position_jitter(width = .05), shape = 21, size = 1.3, alpha = 0.3) + #
    geom_boxplot(aes(fill = AP), width = .15, outlier.shape = NA, alpha = 0.6) +
    theme_light() +
    theme(axis.text.x = element_text(vjust = 1, hjust = 0.5)) + #angle = 80, 
    scale_fill_manual(values = cols <- c("#33DD00","#FFCC00","#2233FF"))+#,"grey")) +
    # the below line should be used for coloring boxplots by AP with the argument 'colour='
    #scale_colour_manual(values = c("green3","gold2","dodgerblue3"))+#,"blue")) +
    guides(fill = "none", color = "none") +
    ylab("rate") +
    geom_rect(aes(xmin = 7 - 0.3, xmax = 9 + 0.55, ymin = 0 - 0.05, ymax = max(value, na.rm=TRUE) + 0.1),
              color = "#FF7799", fill = "transparent", size = 1.3) 
dev.off()
    
tb1 %>% group_by(state) %>% tally()


#### 2.3.1.2 Aggregated evolutionary rates (by AP) ----------------------------

## average rate estimates over states for each activity pattern
# make cols for aggergated rates
evol$div_N <- evol$ext_N <- evol$spec_N <- as.numeric(rep("", nrow(evol))) 
evol$div_C <- evol$ext_C <- evol$spec_C <- as.numeric(rep("", nrow(evol))) 
evol$div_D <- evol$ext_D <- evol$spec_D <- as.numeric(rep("", nrow(evol))) 

for(i in 1:nrow(evol)){
    # find number of states in the model
    nstate <- str_split(rownames(evol)[i], "_")[[1]][1] %>%
        str_sub(-1L,-1L) %>%
        as.numeric()
    # exclude columns of states not in the model (use clone to notmodify original data frame)
    clone <- evol[i,]
    dummies <- which(str_sub(colnames(clone), -1L, -1L) %in% LETTERS[(nstate+1):8])
    clone[,dummies] <- NA 
    
    ## calculate mean of relevant rate estimates
    for(stt in c("00","01","11")){
        if(stt == "00") AP <- "N"
        if(stt == "01") AP <- "C"
        if(stt == "11") AP <- "D"
        
        # isolate relevant rates
        spec_cols <- paste0("s",stt,LETTERS[1:nstate])
        ext_cols <- paste0("mu",stt,LETTERS[1:nstate])
        div_cols <- paste0("lambda",stt,LETTERS[1:nstate])
        
        # locate correct column to write mean rate into
        s_agg <- which(colnames(evol) == paste0("spec_",AP))
        e_agg <- which(colnames(evol) == paste0("ext_",AP))
        d_agg <- which(colnames(evol) == paste0("div_",AP))
        
        # calculate means and write to main table
        evol[i,s_agg] <- mean(as.numeric(clone[,which(colnames(clone) %in% spec_cols)]))
        evol[i,e_agg] <- mean(as.numeric(clone[,which(colnames(clone) %in% ext_cols)]))
        evol[i,d_agg] <- mean(as.numeric(clone[,which(colnames(clone) %in% div_cols)]))
    }
}
#store in separate object
delim <- which(colnames(evol) == "model")
mean_rates <- evol[,delim:ncol(evol)]
evol <- evol[,1:delim]

# melt into long form
agg <- pivot_longer(mean_rates, cols=2:ncol(mean_rates), names_to = "mean_rate", values_to = "value") %>%
    drop_na()

## remove outliers and plot
d <- numeric()
for(i in 1:nrow(agg)){if(strsplit(agg$model[i], "_")[[1]][2] == "tree1803"){d <- c(d,i)}}
agg_NOL <- agg[-d,]
agg_NOL$type <- agg_NOL$AP <- NA
agg_NOL$AP <- str_sub(agg_NOL$mean_rate, -1L, -1L)
agg_NOL$AP <- gsub("N","Noct", agg_NOL$AP)
agg_NOL$AP <- gsub("C","Cath", agg_NOL$AP)
agg_NOL$AP <- gsub("D","Diur", agg_NOL$AP)
agg_NOL$type <- str_sub(agg_NOL$mean_rate, 1L, -3L)
agg_NOL$type <- gsub("div","Diversification", agg_NOL$type)
agg_NOL$type <- gsub("ext","Extinction", agg_NOL$type)
agg_NOL$type <- gsub("spec","Speciation", agg_NOL$type)

# calculate median values for plot labels
agg_rate <- c("spec_N","spec_C","spec_D","ext_N","ext_C","ext_D","div_N","div_C","div_D")
agg_medians <- as.data.frame(agg_rate)
for(i in unique(agg_NOL$mean_rate)){
    agg_medians$median[which(agg_medians$agg_rate == i)] <- median(agg_NOL$value[which(agg_NOL$mean_rate == i)])
}


pdf(file="Figure 3 - Aggregated evol rates est_no_outlier_19JUL.pdf", width = 14, height = 8)
agg_NOL %>%
    mutate(mean_rate = factor(mean_rate, levels=c("spec_N", "spec_C", "spec_D", "ext_N", "ext_C", "ext_D", "div_N", "div_C", "div_D"))) %>%
    ggplot(aes(x = mean_rate, y = value)) +
        geom_flat_violin(aes(fill = AP), position = position_nudge(x = .15, y = 0), alpha = .8) +
        geom_point(aes(fill = AP), position = position_jitter(width = .03), shape = 21, size = 1.5, alpha = 0.4) + 
        geom_boxplot(aes(fill = AP), width = .15, outlier.shape = NA, alpha = 0.6) +
        theme_light(base_size = 24) +
        theme(axis.text.x = element_text(vjust = 1, hjust = 0.5)) + #angle = 80, 
        scale_fill_manual(values = cols <- c("#33DD00","#FFCC00","#2233FF")) +
        # the below line should be used for coloring boxplots by AP with the argument 'colour='
        scale_colour_manual(values = c("green3","gold2","dodgerblue3"))+#,"blue")) +
        guides(fill = "none", color = "none") +
        #facet_wrap(~type) +
        geom_label(label = round(agg_medians$median,3), data = agg_medians,
                   x = -0.3+c(1:9), y = agg_medians$median,
                   label.size = 0.07, label.padding = unit(0.2, "lines"), # Rectangle size around label
                   color = "black", fill=c("#6677FF","#33DD00","#FFCC00","#6677FF","#33DD00",
                                           "#FFCC00","#6677FF","#33DD00","#FFCC00"), alpha=0.6) +
        # surrounding rectangle
        #geom_rect(aes(xmin = 7 - 0.3, xmax = 9 + 0.55, ymin = 0 - 0.05, ymax = max(value, na.rm=TRUE) + 0.1),
        #          color = "#FF7799", fill = "transparent", size = 0.8)
        geom_rect(aes(xmin = 7 - 0.4, xmax = 9 + 0.55, ymin = -0.01, ymax = 0.01),
                  fill = "#FF7799", size = 0.8)
dev.off()


#### 2.3.1.3 Within-state relative rates --------------------------------------

## in each state, show the difference between nocturnal rates and diurnal or cathemeral
rel_rates <- evol[,-c(25:48)] ## exclude columns of states not in any model 

## calculate rate differences
for(z in 0:11){ # iterate over the 12 sets of rates(3 evol. rate types, up to 4 hidden states), 3 columns in each set (N,C, and D)
    # in each set of 3 columns the first one is the nocturnal baseline
    baseline <- rel_rates[,3*z+1]
    # caluclate differences for all columns
    rel_rates[,3*z+1] <- rel_rates[,3*z+1] - baseline
    rel_rates[,3*z+2] <- rel_rates[,3*z+2] - baseline
    rel_rates[,3*z+3] <- rel_rates[,3*z+3] - baseline
}
    
## average rate estimates over states for each activity pattern
# make cols for aggergated rates
rel_rates$div_N <- rel_rates$ext_N <- rel_rates$spec_N <- as.numeric(rep("", nrow(rel_rates))) 
rel_rates$div_C <- rel_rates$ext_C <- rel_rates$spec_C <- as.numeric(rep("", nrow(rel_rates))) 
rel_rates$div_D <- rel_rates$ext_D <- rel_rates$spec_D <- as.numeric(rep("", nrow(rel_rates))) 

# aggregate rate differences for each hidden state
for(i in 1:nrow(rel_rates)){
    # find number of states in the model
    nstate <- as.numeric(str_sub(rownames(rel_rates)[4], -14L,-14L))
        
    # alternative that was used in figures:
        #str_split(rownames(rel_rates)[i], "_")[[1]][1] %>%
        #str_sub(-1L,-1L) %>%
        #as.numeric()
    
    # exclude columns of states not in the model (use clone to notmodify original data frame)
    clone <- rel_rates[i,]
    dummies <- which(str_sub(colnames(clone), -1L, -1L) %in% LETTERS[(nstate+1):8])
    clone[,dummies] <- NA 
    
    ## calculate mean of relevant rate estimates
    for(chr in c("00","01","11")){
        if(chr == "00") AP <- "N"
        if(chr == "01") AP <- "C"
        if(chr == "11") AP <- "D"
        
        # isolate relevant rates
        spec_cols <- paste0("s",chr,LETTERS[1:nstate])
        ext_cols <- paste0("mu",chr,LETTERS[1:nstate])
        div_cols <- paste0("lambda",chr,LETTERS[1:nstate])
        
        # locate correct column to write mean rate into
        s_agg <- which(colnames(rel_rates) == paste0("spec_",AP))
        e_agg <- which(colnames(rel_rates) == paste0("ext_",AP))
        d_agg <- which(colnames(rel_rates) == paste0("div_",AP))
        
        # calculate means and write to main table
        rel_rates[i,s_agg] <- mean(as.numeric(clone[,which(colnames(clone) %in% spec_cols)]))
        rel_rates[i,e_agg] <- mean(as.numeric(clone[,which(colnames(clone) %in% ext_cols)]))
        rel_rates[i,d_agg] <- mean(as.numeric(clone[,which(colnames(clone) %in% div_cols)]))
    }
}

#store in separate object
delim <- which(colnames(rel_rates) == "model")
rate_diffs <- rel_rates[,1:delim]
relat <- rel_rates[,delim:ncol(rel_rates)]

# melt into long form
relat <- pivot_longer(relat, cols=2:ncol(relat), names_to = "relative_rate", values_to = "value") %>%
    drop_na()

relat$type <- relat$AP <- NA
relat$type <- str_sub(relat$relative_rate, 1L, -3L)
relat$AP <- str_sub(relat$relative_rate, -1L, -1L)
relat$AP <- gsub("N","Noct", relat$AP)
relat$AP <- gsub("C","Cath", relat$AP)
relat$AP <- gsub("D","Diur", relat$AP)


## remove outliers and plot
d <- numeric()
for(i in 1:nrow(relat)){if(strsplit(relat$model[i], "_")[[1]][2] == "tree1803"){d <- c(d,i)}}
relat_NOL <- relat[-d,]

# calculate median estimates for labels
rel_rate <- c("spec_N","spec_C","spec_D","ext_N","ext_C","ext_D","div_N","div_C","div_D")
rel_medians <- as.data.frame(rel_rate)
for(i in unique(relat_NOL$relative_rate)){
    rel_medians$median[which(rel_medians$rel_rate == i)] <- median(relat_NOL$value[which(relat_NOL$relative_rate == i)])
}


pdf(file="Figure 4 - Evol rates relative to NOCT_19JUL.pdf", width = 14, height = 8)
relat_NOL %>%
    mutate(relative_rate = factor(relative_rate, levels=c("spec_N", "spec_C", "spec_D", "ext_N", 
                                                          "ext_C", "ext_D", "div_N", "div_C", "div_D"))) %>%
    ggplot(aes(x = relative_rate, y = value)) +
    geom_flat_violin(aes(fill = AP), position = position_nudge(x = .15, y = 0), alpha = .8) +
    geom_point(aes(fill = AP), position = position_jitter(width = .03), shape = 21, 
               size = 1.5, alpha = 0.4) + 
    geom_boxplot(aes(fill = AP), width = .15, outlier.shape = NA, alpha = 0.6) +
    #facet_wrap(~type) +
    geom_label(label = round(rel_medians$median,3), data = rel_medians,
               x = c(1:9)-c(0.2,0.3,0.3,0.2,0.3,0.3,0.2,0.3,0.3), y = rel_medians$median,
               label.size = 0.08, label.padding = unit(0.2, "lines"), # Rectangle size around label
               color = "black", fill=c("#2233FF","#33DD00","#FFCC00","#2233FF","#33DD00",
                                       "#FFCC00","#2233FF","#33DD00","#FFCC00"), alpha=0.4) +
    theme_light(base_size = 24) +
    theme(axis.text.x = element_text(vjust = 1, hjust = 0.5)) + #, angle = 80)) + # 
    scale_fill_manual(values = cols <- c("#33DD00","#FFCC00","#2233FF")) +
    guides(fill = "none", color = "none") +
    geom_rect(aes(xmin = 7 - 0.4, xmax = 9 + 0.55, ymin = -0.2 - 0.007, ymax = -0.2 + 0.007),
              fill = "#FF7799", size = 0.8)
    #geom_rect(aes(xmin = 7 - 0.3, xmax = 9 + 0.65, ymin = -0.1 - 0.05, ymax = max(value, na.rm=TRUE) + 0.05),
    #          color = "#FF7799", fill = "transparent", size = 1.3)
dev.off()


### 2.3.2 Transition rates ----------------------------------------------------

tran$model <- str_replace(rownames(tran), ".RDS", "") # remove file suffix from model names

# find and then remove model/s based on MCC tree
for(i in 1:nrow(tran)){if(strsplit(tran$model[i], "_")[[1]][2] == "MCCtree"){drop_mcc <- i}}
tran <- tran[-drop_mcc,]

# melt rates into long form using the 'model' column as denominator
Trates <- pivot_longer(tran, cols=1:ncol(tran)-1, names_to = "rate", values_to = "value") %>%
    drop_na()

# removing transitions between states that do not exist in the model
bi_state <- which(str_sub(Trates$model, start = -10L, end = -10L) == 2) # 2-state models
tri_state <- which(str_sub(Trates$model, start = -10L, end = -10L) == 3) # 3-state models
hi_states <- which(!is.na(str_extract(Trates$rate, "C")))   # transitions to/from state C
vhi_states <-  which(!is.na(str_extract(Trates$rate, "D")))   # transitions to/from state D
nonexist <- which(!is.na(str_extract(Trates$rate, "[E-H]")))  # transitions outside states A-D
remove <- c(intersect(bi_state, hi_states), intersect(tri_state, vhi_states), nonexist)
Trates <- Trates[-remove,]

# adding useful information
modtype <- rtype <- rclass <- character()
for (i in 1:nrow(Trates)) {
    prefix <- str_split(Trates$model[i], "Mu")[[1]][1]
    num <- str_extract(str_split(Trates$model[i], "_")[[1]][1], "[0-9]") # get no. of states
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
rm(modtype, rtype, rclass, bi_state, tri_state, hi_states, vhi_states, nonexist, remove)

tbr1 <- Trates
tbr2 <- tbr1[which(tbr1$rclass == "character change"),]        # exclude hidden rates
tbr3 <- tbr2[which(tbr2$value < 50),]
#tbr4 <- tbr2[which(tbr2$modtype == 4),]

## removing outlier tree - 
#d <- numeric()
#for(i in 1:nrow(Trates)){if(strsplit(Trates$model[i], "_")[[1]][2] == "tree1803"){d <- c(d,i)}}
#tbr5 <- tbr2[-d,]


## get average transition rates over all hidden states for each model
# keep only character trasition rates
CTR <- Trates[which(Trates$rclass == "character change"),]

# make a table for output
avtran <- as.data.frame(matrix(data=NA, nrow=524, ncol=4))
colnames(avtran) <- colnames(CTR)[1:4]

# iterate over models to calculate mean transition rates
for(i in 1:131){
    evmod <- unique(CTR$model)[i]
    # find relevant rows for each transition type (using clone to not modify original data)
    clone <- CTR[which(CTR$model == evmod),]
    # isolate rates by transition type
    NC <- which(!is.na(str_extract(clone$rate, ".00..01.")))
    DC <- which(!is.na(str_extract(clone$rate, ".11..01.")))
    CN <- which(!is.na(str_extract(clone$rate, ".01..00.")))
    CD <- which(!is.na(str_extract(clone$rate, ".01..11.")))
    
    # make a summary table for each model
    OP <- as.data.frame(matrix(data=NA, nrow=4, ncol=ncol(avtran)))
    colnames(OP) <- colnames(CTR)[1:4]
    OP[,c(1,4)] <- clone[1,c(1,4)]
    OP[,2] <- c("N->C","C->N","C->D","D->C")
    
    # calculate mean of relevant rate estimates
    OP[,3] <- c(mean(clone$value[NC]), mean(clone$value[CN]), 
                mean(clone$value[CD]), mean(clone$value[DC]))
    
    # write in output table
    avtran[(i*4-3):(i*4),] <- OP
}


# calculate median estimates for plot labels
trn_rate <- c("N->C","C->N","C->D","D->C")
trn_medians <- as.data.frame(trn_rate)
for(i in unique(tbr2$rtype)){
    trn_medians$median[which(trn_medians$trn_rate == i)] <- median(log(tbr2$value[which(tbr2$rtype == i)]))
}

# round to zero estimates that are below detectable limit
#tbr5 <- tbr2
#tbr5[which(tbr5$value < 10^-6),3] <- 10^-6 # 2*10^-6 < 1/(total tree depth * total number of lineages) == 1/(167.662*2400) == 2.48*10^-6
# rates below this threshold are expected to be undetectable given current amount of data


# plot transition rate estimates
tbr2$rtype <- factor(tbr2$rtype, levels=c("N->C","C->N","C->D","D->C"))
pdf(file="Figure 5 - Transition rates estimates_25JUL.pdf", width = 12, height = 8)
ggplot(tbr2, aes(x = rtype, y = log(value+10^-10), fill = rtype)) + # added 10^-10 to the rate estimates to deal with zeros (72 rows)
    geom_flat_violin(position = position_nudge(x = .15, y = 0), alpha = .8) +
    geom_point(colour = "grey20", position = position_jitter(width = .02), shape = 16, 
               size = 1.2, alpha = 0.4) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.7) +
    geom_hline(aes(yintercept = log(2*10^-6), colour = '#5A4A6F'), linetype = 2) + # this marks the effective rate of 0, because 2*10^-6 transitions per lineage per million years is 1 transition per 5000 lineages per 100 million years, which is unlikely to be detected when I examine 2400 tips over 167.662myr (most lineages are MUCH younger than that)
    theme_light(base_size = 26) +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1)) +
    annotate("rect", xmin = 0.45, xmax = 4.6, ymin = -24, ymax = log(2*10^-6), alpha=0.5, fill="grey") +
    scale_fill_manual(values = cols <- c("#AA22BB","#EE4444","gold1","#33AAFF")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250","#EBB261", "#9D5A6C","#5A4A6F","#E47250","#EBB261", "#9D5A6C")) +
    labs(y = expression(Transition~rate~~~ln(lineage^-1~Myr^-1)), 
         x = "Direction of transition (from -> to)") +
    geom_label(label = round(trn_medians$median,2), data = trn_medians,
               x = c(1:4)-0.2, y = trn_medians$median,
               label.size = 0.5, label.padding = unit(0.3, "lines"), label.r = unit(0.1, "lines"),
               color = "black", fill = c("#AA22BB","#EE4444","gold1","#33AAFF"), alpha=0.3) +
    guides(fill = "none", color = "none") #+
    # Setting the limits of the y axis
    #scale_y_continuous(limits = c(0, .03)) 
dev.off()


## plot version with estimates averaged over hidden states
# calculate median estimates for plot labels
trn_rate <- c("N->C","C->N","C->D","D->C")
trn_medians <- as.data.frame(trn_rate)
for(i in unique(avtran$rate)){
    trn_medians$median[which(trn_medians$trn_rate == i)] <- median(log(avtran$value[which(avtran$rate == i)]))
}


avtran$rate <- factor(avtran$rate, levels=c("N->C","C->N","C->D","D->C"))
pdf(file="Figure 5 - Transition rates estimate means_25JUL.pdf", width = 12, height = 8)
ggplot(avtran, aes(x = rate, y = log(value), fill = rate)) +
    geom_flat_violin(position = position_nudge(x = .15, y = 0), alpha = .8) +
    geom_point(colour = "grey20", position = position_jitter(width = .02), shape = 16, 
               size = 1.2, alpha = 0.4) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.7) +
    theme_light(base_size = 26) +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1)) +
    scale_fill_manual(values = cols <- c("#AA22BB","#EE4444","gold1","#33AAFF")) +
    #scale_colour_manual(values = c("#5A4A6F", "#E47250","#EBB261", "#9D5A6C","#5A4A6F","#E47250","#EBB261", "#9D5A6C")) +
    labs(y = expression(Transition~rate~~~~ln(lineage^-1~Myr^-1)), 
         x = "Direction of transition (from -> to)") +
    geom_label(label = round(trn_medians$median,2), data = trn_medians,
               x = c(1:4)-0.2, y = trn_medians$median,
               label.size = 0.5, label.padding = unit(0.3, "lines"), label.r = unit(0.1, "lines"),
               color = "black", fill = c("#AA22BB","#EE4444","gold1","#33AAFF"), alpha=0.3) +
    guides(fill = "none", color = "none")
dev.off()


## 2.4 Omitting 50 worst performing trees -------------------------------------
# TL;DR - omitting 50 worst performing trees doesn't visibly change the results 

# order trees by best performing model (BIC)
srtmod <- mod_avg[order(mod_avg$BIC),]
slct <- srtmod[1:78,] # this eliminates the 50 worst performing trees (some have more than one model) 

# the below are trees that produce more than one supported model (in the full set and slct75, respectively)
srtmod[which(srtmod$tree %in% c("3120","6760","6375","4420","2555","2417","3865")),c(1:11)]
slct[which(slct$tree %in% c("3120","6760","6375","4420")), c(1:11)]

### 2.4.1 Evolutionary rates --------------------------------------------------
evol1 <- slct[,-c(1:11)]
evol1 <- evol1[which(is.na(str_extract(colnames(evol1), "q")))] # all rates that are NOT transition rates
evol1$model <- str_replace(rownames(evol1), ".RDS", "") # remove file suffix from model names

# find and then remove model/s based on MCC tree
for(i in 1:nrow(evol1)){if(strsplit(evol1$model[i], "_")[[1]][2] == "MCCtree"){drop_mcc <- i}}
evol1 <- evol1[-drop_mcc,]

# melt into long form
evo_1 <- pivot_longer(evol1, cols=1:ncol(evol1)-1, names_to = "rate", values_to = "value") %>%
    drop_na()

{type <- modtype <- rtype <- AP <- state <- character()
    
    type[which(!is.na(str_extract(evo_1$rate, "s")))] <- "spec"
    type[which(!is.na(str_extract(evo_1$rate, "mu")))] <- "ext"
    type[which(!is.na(str_extract(evo_1$rate, "lambda")))] <- "net.div"
    
    for (i in 1:nrow(evo_1)){
        modtype[i] <- str_extract(str_split(evo_1$model[i], "_")[[1]][1], "[0-9]")
        rtype[i] <- str_sub(evo_1$rate[i], start = 1L, end = -2L)}
    
    AP[which(str_sub(evo_1$rate, start = -3L, end = -2L) == "00")] <- "Noct"
    AP[which(str_sub(evo_1$rate, start = -3L, end = -2L) == "01")] <- "Cath"
    AP[which(str_sub(evo_1$rate, start = -3L, end = -2L) == "11")] <- "Diur"
    
    state <- str_sub(evo_1$rate, start = -1L, end = -1L)
    }

evo_1 <- cbind(evo_1, type, modtype, rtype, AP, state)
rm(type, modtype, rtype, AP, state)  

## removing transitions between states that do not exist in the model
tri_state <- which(str_sub(evo_1$model, start = -10L, end = -10L) == 3) # 3-state models
quad_state <- which(str_sub(evo_1$model, start = -10L, end = -10L) == 4)# 4-state models
hi_states <- which(evo_1$state %in% c("D","E","F","G","H"))   # transitions outside states A-C
vhi_states <- which(evo_1$state %in% c("E","F","G","H"))  # transitions outside states A-D
remove <- c(intersect(tri_state, hi_states), intersect(quad_state, vhi_states))
evo_1 <- evo_1[-remove,]
rm(tri_state, quad_state, hi_states, vhi_states, remove)

evo_1$rtype <- factor(evo_1$rtype, levels = c("s00","s01","s11","mu00","mu01","mu11","lambda00","lambda01","lambda11"))
ggplot(evo_1, aes(x = rtype, y = value)) +
    geom_flat_violin(aes(fill = AP), position = position_nudge(x = .15, y = 0), alpha = .8) +
    geom_point(aes(fill = AP), position = position_jitter(width = .05), shape = 21, size = 1.3, alpha = 0.3) + #
    geom_boxplot(aes(fill = AP), width = .15, outlier.shape = NA, alpha = 0.6) +
    theme_light() +
    theme(axis.text.x = element_text(vjust = 1, hjust = 0.5)) + #angle = 80, 
    scale_fill_manual(values = cols <- c("#33DD00","#FFCC00","#2233FF"))+#,"grey")) +
    # the below line should be used for coloring boxplots by AP with the argument 'colour='
    #scale_colour_manual(values = c("green3","gold2","dodgerblue3"))+#,"blue")) +
    guides(fill = "none", color = "none") +
    geom_rect(aes(xmin = 7 - 0.3, xmax = 9 + 0.55, ymin = 0 - 0.05, ymax = max(value, na.rm=TRUE) + 0.1),
              color = "red", fill = "transparent", size = 1.2) 


### 2.4.2 Transition rates ----------------------------------------------------
tran1 <- slct[,-c(1:11)]
tran1 <- tran1[which(!is.na(str_extract(colnames(tran1), "q")))] # transition rates only
tran1$model <- str_replace(rownames(tran1), ".RDS", "") # remove file suffix from model names

# find and then remove model/s based on MCC tree
for(i in 1:nrow(tran1)){if(strsplit(tran1$model[i], "_")[[1]][2] == "MCCtree"){drop_mcc <- i}}
tran1 <- tran1[-drop_mcc,]

# melt rates into long form using the 'model' column as denominator
Tr1 <- pivot_longer(tran1, cols=1:ncol(tran1)-1, names_to = "rate", values_to = "value") %>%
    drop_na()

# removing transitions between states that do not exist in the model
tri_state <- which(str_sub(Tr1$model, start = -10L, end = -10L) == 3) # 3-state models
quad_state <- which(str_sub(Tr1$model, start = -10L, end = -10L) == 4)# 4-state models
hi_states <- which(!is.na(str_extract(Tr1$rate, "[D-H]")))   # transitions outside states A-C
vhi_states <- which(!is.na(str_extract(Tr1$rate, "[E-H]")))  # transitions outside states A-D
remove <- c(intersect(tri_state, hi_states), intersect(quad_state, vhi_states))
Tr1 <- Tr1[-remove,]

modtype <- rtype <- rclass <- character()
for (i in 1:nrow(Tr1)) {
    prefix <- str_split(Tr1$model[i], "Mu")[[1]][1]
    num <- str_extract(str_split(Tr1$model[i], "_")[[1]][1], "[0-9]") # get # of states
    modtype[i] <- paste0(prefix,num)
    rtype[i] <- gsub("[A-H]_", "->", str_sub(Tr1$rate[i], start = 2L, end = -2L)) # transition type
    if (rtype[i] %in% c("00->00","01->01","11->11")) {
        rclass[i] <- "hidden change"
    } else {
        rclass[i] <- "character change"
    }
}
Tr1 <- cbind(Tr1, modtype, rtype, rclass)
Tr1$rtype <- sapply(Tr1$rtype, gsub, pattern = "00", replacement = "N")
Tr1$rtype <- sapply(Tr1$rtype, gsub, pattern = "01", replacement = "C") 
Tr1$rtype <- sapply(Tr1$rtype, gsub, pattern = "11", replacement = "D")
rm(modtype, rtype, rclass, tri_state, quad_state, hi_states, vhi_states, remove)

trns1 <- Tr1[which(Tr1$rclass == "character change"),] 
trns2 <- trns1
trns2[which(trns2$value < 2*10^-6),3] <- 2*10^-6

trns1$rtype <- factor(trns1$rtype, levels=c("N->C","C->N","C->D","D->C"))
ggplot(trns1, aes(x = rtype, y = log(value+0.0000000001), fill = rtype)) + # added 10^-10 to the rate estimates to deal with zeros (72 rows)
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(colour = "grey20", position = position_jitter(width = .05), shape = 21, size = 1.3, alpha = 0.5) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    geom_hline(aes(yintercept = log(10^-6), colour = 'red'), linetype = 2) + # this marks the effective rate of 0, because 10^-6 transitions per lineage per million years is 1 transition per 10000 lineages per 100 million years, which is unlikely to be detected when I examine <2500 lineages over 160 million years
    theme_light() +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1)) +
    annotate("rect", xmin = 0.45, xmax = 4.6, ymin = -24, ymax = log(2*10^-6), alpha=0.6, fill="grey") +
    scale_fill_manual(values = cols <- c("#AA22BB","#EE4444","gold1","#44DDFF")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    labs(y = "log(transition rate)", x = "Transition type (from -> to)") +
    #annotate("text", x=2.5, y=-16, label= "undetectable\nrange", size = 4.5) + 
    guides(fill = "none", color = "none") #+
# Setting the limits of the y axis
#scale_y_continuous(limits = c(0, .03)) 


# 3. Follow up analyses -------------------------------------------------------

library("hisse")

setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/Analyses/ResultsCluster")


## 3.1 Reconstruction based on fitted models ----------------------------------

# repeating steps from MuHiSSE analysis ("4-ClusterCode.R")

states <- data.frame(tree$tip.label, tree$tip.state, tree$tip.state)
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
}

mod <- readRDS(file="MuHiSSE2_tree6404.RDS")


# 4. [ open ] -----------------------------------------------------------------



# 5. Processing binary models -------------------------------------------------

## 5.1 Any daytime activity [N-CD] --------------------------------------------
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
pdf(file="Model support ordered_MCC marked N-CD.pdf", width = 10, height = 8)
#df1 <- score[score$states != 1,]
df2 <- subset(score, score$tree =="MCC")    # get only MCC models separately
#par(mfrow=c(2,1))
ggplot(score, aes(x = reorder(tree, AICc, FUN = max), y = AICc, colour = as.factor(states))) +
    geom_point(alpha= .8, size = 2.5) + 
    geom_point(data = df2, aes(tree, AICc), colour = 'red', shape = 1, size = 3) +
    scale_y_reverse() + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("tree variant") +
    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)

ggplot(score, aes(x = reorder(tree, BIC, FUN = max), y = BIC, colour = as.factor(states))) +
    geom_point(alpha= .8, size = 2.5) + 
    geom_point(data = df2, aes(tree, BIC), colour = 'red', shape = 1, size = 3) + 
    scale_y_reverse() +
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("tree variant") +
    scale_color_viridis_d(option = "inferno", begin = 0, end = .9)
dev.off()


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
pdf(file = "Model support aggregated N-CD.pdf", width = 8.5, height = 6)
ggplot(data = vars, aes(y = AICc, x = as.factor(states), fill = as.factor(states))) +
        geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
        geom_point(aes(y = AICc, x = as.factor(states), color = as.factor(states)), 
                   position = position_jitter(width = .1), size = 3, alpha = 0.6) +
        geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
        guides(fill = FALSE, color = FALSE) + # don't show a legend
        #labs(color = "States", fill = "States") + # in case a legend is needed
        xlab("No. of states") +
        facet_wrap(~type,) +
        expand_limits(x = 5.25) +
        scale_fill_viridis_d(option = "viridis", begin = 0, end = 1, alpha = 0.8, direction = -1) +
        scale_color_viridis_d(option = "viridis", begin = 0, end = 1, alpha = 1, direction = -1) +
        theme_minimal(base_size = 12)  

ggplot(data = vars, aes(y = BIC, x = as.factor(states), fill = as.factor(states))) +
       geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
       geom_point(aes(y = BIC, x = as.factor(states), color = as.factor(states)), 
                  position = position_jitter(width = .1), size = 3, alpha = 0.6) +
       geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
       guides(fill = FALSE, color = FALSE) + # don't show a legend
       #labs(color = "States", fill = "States") + # in case a legend is needed
       xlab("No. of states") +
       facet_wrap(~type,) +
       expand_limits(x = 5.25) +
       scale_fill_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 0.8, direction = -1) +
       scale_color_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 1, direction = -1) +
       theme_minimal(base_size = 12)
dev.off()

# plotting support by no. states with trendlines to show that the pattern is uniform accross trees
pdf(file = "Model support by state_trendlines N-CD.pdf", width = 7, height = 8.5)
ggplot(vars,aes(x=states, AICc, colour = reorder(tree, AICc, FUN = mean), group = tree)) +
    geom_point(alpha= .5, size = 2) +
    # use the below if need to plot MCC along with tree variants
    #geom_point(data = tert[tert$tree == "MCC",], aes(states, AICc), colour = 'grey10', shape = 4, size = 2) + 
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    theme_minimal() +
    guides(colour=FALSE) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .8, direction = -1)

ggplot(vars,aes(x=states, BIC, colour = reorder(tree, BIC, FUN = mean), group = tree)) +
    geom_point(alpha= .6, size = 2) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    theme_minimal() +
    guides(colour=FALSE) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .8, direction = -1)
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
pdf(file="model type performance by tree N-CD.pdf", width = 10, height = 8)
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


### 5.1.3 Rate estimates ------------------------------------------------------

*TODO*: FILL IN

## 5.2 Strict diurnality only [D-CN] ------------------------------------------
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
pdf(file="Model support ordered_MCC marked D-CN.pdf", width = 10, height = 8)
#df1 <- score[score$states != 1,]

df2 <- subset(score, score$tree =="MCC")    # get only MCC models separately
#par(mfrow=c(2,1))
ggplot(score, aes(x = reorder(tree, AICc, FUN = max), y = AICc, colour = as.factor(states))) +
    geom_point(alpha= .8, size = 2.5) + 
    geom_point(data = df2, aes(tree, AICc), colour = 'red', shape = 1, size = 3) +
    scale_y_reverse() + 
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("tree variant") +
    scale_color_viridis_d(option = "viridis", begin = 0, end = .9, direction = -1)

ggplot(score, aes(x = reorder(tree, BIC, FUN = max), y = BIC, colour = as.factor(states))) +
    geom_point(alpha= .8, size = 2.5) + 
    geom_point(data = df2, aes(tree, BIC), colour = 'red', shape = 1, size = 3) + 
    scale_y_reverse() +
    facet_wrap(~type) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    labs(color = "States") +
    xlab("tree variant") +
    scale_color_viridis_d(option = "inferno", begin = 0, end = .9)
dev.off()



# plotting support by no. states, aggregaetd over all variants
vars <- score[-which(score$tree == "MCC"),]
pdf(file = "Model support aggregated D-CN.pdf", width = 10, height = 8)
ggplot(data = vars, aes(y = AICc, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = AICc, x = as.factor(states), fill = as.factor(states)), 
               position = position_jitter(width = .1), shape = 21, size = 2.5, alpha = 0.6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    guides(fill = FALSE, color = FALSE) + # don't show a legend
    #labs(color = "States", fill = "States") + # in case a legend is needed
    xlab("No. of states") +
    facet_wrap(~type,) +
    expand_limits(x = 5.25) +
    scale_fill_viridis_d(option = "viridis", begin = 0, end = 1, alpha = 0.8, direction = -1) +
    scale_color_viridis_d(option = "viridis", begin = 0, end = 1, alpha = 1, direction = -1) +
    theme_minimal(base_size = 12)  

ggplot(data = vars, aes(y = BIC, x = as.factor(states), fill = as.factor(states))) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = BIC, x = as.factor(states), fill = as.factor(states)), 
               position = position_jitter(width = .1), shape = 21, size = 2.5, alpha = 0.6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.6) +
    guides(fill = FALSE, color = FALSE) + # don't show a legend
    #labs(color = "States", fill = "States") + # in case a legend is needed
    xlab("No. of states") +
    facet_wrap(~type,) +
    expand_limits(x = 5.25) +
    scale_fill_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 0.8, direction = -1) +
    scale_color_viridis_d(option = "inferno", begin = 0.1, end = 0.9, alpha = 1, direction = -1) +
    theme_minimal(base_size = 12)
dev.off()

# plotting support by no. states with trendlines to show that the pattern is uniform accross trees
pdf(file = "Model support by state_trendlines D-CN.pdf", width = 7, height = 8.5)
ggplot(vars,aes(x=states, AICc, colour = reorder(tree, AICc, FUN = mean), group = tree)) +
    geom_point(alpha= .8, size = 2) +
    # use the below if need to plot MCC along with tree variants
    #geom_point(data = tert[tert$tree == "MCC",], aes(states, AICc), colour = 'grey10', shape = 4, size = 2) + 
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    theme_minimal() +
    guides(colour=FALSE) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .8, direction = -1)

ggplot(vars,aes(x=states, BIC, colour = reorder(tree, BIC, FUN = mean), group = tree)) +
    geom_point(alpha= .8, size = 2) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust = 1)) +
    theme_minimal() +
    guides(colour=FALSE) +
    facet_wrap(~type,) + 
    scale_color_viridis_d(option = "turbo", alpha = .8, direction = -1)
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
pdf(file="model type performance by tree D-CN.pdf", width = 10, height = 8)
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


### 5.2.3 Rate estimates ------------------------------------------------------



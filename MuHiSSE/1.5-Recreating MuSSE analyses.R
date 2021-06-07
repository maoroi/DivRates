### Recreating the original MuSSE analyses with updated taxonomy and tree

require("ape")
require("diversitree")
require("phangorn")
require("phytools")
require("plotrix")
library("job")
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/New MuSSE")
set.seed(88)

## preparing data
act <- read.csv(file="ActivityData_MDD_v1_match.csv")   # 2438 species in the data
tree <- read.tree("tree6404.nex")                       # 2400 tips in the tree
tree <- phytools::force.ultrametric(tree)
dat <- act[which(act$Phylo_name %in% tree$tip.label),c(6,5)]  # remove unmatched data

states <- dat$AP
states[which(states=="Nocturnal")] <- 1
states[which(states=="Cathemeral")] <- 2
states[which(states=="Diurnal")] <- 3
states <- as.numeric(states)
names(states) <- dat$Phylo_name

freq <- c(0.3523062, 0.4836342, 0.4083013)
names(freq) <- c(1,2,3)

## MuSSE (3 state dataset)

job::job(simple_MuSSE = {
    lik <- make.musse(tree, states, 3, sampling.f=freq, strict=TRUE, control=list())
    p <- starting.point.musse(tree, 3)
    prior <- make.prior.exponential(1/(2 * (p[1] - p[4])))
    
    ## diversification equal for all states but transition rates follow ordered model
    lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1, mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q31 ~ 0)
    fit.base <- find.mle(lik.base, p[argnames(lik.base)])
    ## running an MCMC instead of ML
    samples.b <- mcmc(lik.base, coef(fit.base), nstep = 1000, w = 1, prior = prior, print.every = 50)
})
write.table(samples.b, file=paste('CIDrates_tree6404.csv', sep=''), 
            append = FALSE, quote = TRUE, sep=',', eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

## diversification unconstrained, ordered character-change model
job::job(full_MuSSE = {
    lik.free <- constrain(lik, q13 ~ 0, q31 ~ 0)        
    fit.free <- find.mle(lik.free, p[argnames(lik.free)])
    ## running an MCMC instead of ML
    samples.f <- mcmc(lik.free, coef(fit.free), nstep = 1000, w = 1, prior = prior, print.every = 50)
})
write.table(samples.f, file=paste('MCC',order[k],'_MuSSE_transition_diversification.csv', sep=''), 
            append = FALSE, quote = TRUE, sep=',', eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

# plot null model
pdf(file=paste('MCC',order[k],'_MuSSE_transitions_only.pdf', sep=''), height=6, width=8)
profiles.plot(samples.b[2:4], lwd = 1, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), 
              col.fill = alpha(c('dodgerblue3','#22dd11','darkgoldenrod2'), alpha=0.8), opacity = 0.2, n.br = 100)
dev.off()

# plot MuSSE model
pdf(file=paste('MCC',order[k],'_MuSSE_transition_diversification.pdf', sep=''), height=6, width=8)
profiles.plot(samples.f[2:4], lwd = 1, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), 
              col.fill = alpha(c('dodgerblue3','#22dd11','darkgoldenrod2'), alpha=0.8), opacity = 0.2, n.br = 100)
dev.off()

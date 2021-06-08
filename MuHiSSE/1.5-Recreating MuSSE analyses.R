### Recreating the original MuSSE analyses with updated taxonomy and tree

require("ape")
require("diversitree")
require("phangorn")
require("phytools")
require("plotrix")
require("scales")
library("job")
setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE/New MuSSE")
set.seed(88)

# function to taylor binomials in tree to v1 of the taxonomy
CorTax <- function(tree, data){
    tree$tip.label[which(tree$tip.label == "Equus_africanus")]  <- "Equus_asinus"   # asinus is domestic E. africanus
    tree$tip.label[which(tree$tip.label == "Pygeretmus_zhitkovi")]  <- "Pygeretmus_shitkovi"    # correct spelling
    tree$tip.label[which(tree$tip.label == "Aethomys_namaquensis")] <- "Micaelamys_namaquensis" # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Tamiops_macclellandii")]<- "Tamiops_mcclellandii"   # typo
    tree$tip.label[which(tree$tip.label == "Nesotragus_moschatus")] <- "Neotragus_moschatus"    # typo
    tree$tip.label[which(tree$tip.label == "Pseudalopex_culpaeus")] <- "Lycalopex_culpaeus"     # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Pseudalopex_griseus")]  <- "Lycalopex_griseus"      # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Pseudalopex_gymnocercus")]  <- "Lycalopex_gymnocercus"  # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Pseudalopex_gymnocercus")]  <- "Lycalopex_sechurae" # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Pseudalopex_vetulus")]  <- "Lycalopex_vetulus"      # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Lutra_maculicollis")]   <- "Hydrictis_maculicollis" # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Monodelphis_unistriatus")]  <- "Monodelphis_unistriata" # Latin grammar correction
    tree$tip.label[which(tree$tip.label == "Dactylonax_palpator")]  <- "Dactylopsila_palpator"  # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Elephantulus_revoilii")]<- "Elephantulus_revoili"   # spelling error
    tree$tip.label[which(tree$tip.label == "Zaglossus_bruijnii")]   <- "Zaglossus_bruijni"      # spelling error
    tree$tip.label[which(tree$tip.label == "Procolobus_badius")]    <- "Piliocolobus_badius"    # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Procolobus_kirkii")]    <- "Piliocolobus_kirkii"    # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Procolobus_preussi")]   <- "Piliocolobus_preussi"   # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Procolobus_pennantii")] <- "Piliocolobus_pennantii" # genus transfer v1
    tree$tip.label[which(tree$tip.label == "Procolobus_rufomitratus")]  <- "Piliocolobus_rufomitratus"   # genus transfer v1
    tree <- drop.tip(tree, which(!tree$tip.label %in% data$Phylo_name))
}

## preparing data
act <- read.csv(file="ActivityData_MDD_v1_match.csv")   # 2438 species in the data
tree <- read.tree("tree6404.nex")                       # 2400 tips in the tree
tree <- CorTax(tree, act)
tree <- force.ultrametric(tree)
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
lik <- make.musse(tree, states, 3, sampling.f=freq, strict=TRUE, control=list())
p <- starting.point.musse(tree, 3)
prior <- make.prior.exponential(1/(2 * (p[1] - p[4])))

## diversification equal for all states but transition rates follow ordered model
lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1, mu2 ~ mu1, mu3 ~ mu1, q13 ~ 0, q31 ~ 0)
fit.base <- find.mle(lik.base, p[argnames(lik.base)])
## running an MCMC instead of ML
samples.b <- mcmc(lik.base, coef(fit.base), nstep = 1000, w = 1, prior = prior, print.every = 50)
write.table(samples.b, file=paste('CIDrates_tree6404.csv', sep=''), 
            append = FALSE, quote = TRUE, sep=',', eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

## diversification unconstrained, ordered character-change model
lik.free <- constrain(lik, q13 ~ 0, q31 ~ 0)        
fit.free <- find.mle(lik.free, p[argnames(lik.free)])
## running an MCMC instead of ML
samples.f <- mcmc(lik.free, coef(fit.free), nstep = 1000, w = 1, prior = prior, print.every = 50)
write.table(samples.f, file=paste('MuSSErates_tree6404.csv', sep=''), 
            append = FALSE, quote = TRUE, sep=',', eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

# plot null model
pdf(file=paste('CID1_specext_tree6404.pdf', sep=''), height=6, width=8)
profiles.plot(samples.b[2:3], lwd = 1, col.line = c('purple','#22b211'), 
              col.fill = alpha(c('purple','#22dd11'), alpha=0.8), opacity = 0.2, n.br = 100)
dev.off()

# plot MuSSE model
pdf(file=paste('MuSSE_specrates_tree6404.pdf', sep=''), height=6, width=8)
profiles.plot(samples.f[2:4], lwd = 1, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), 
              col.fill = alpha(c('dodgerblue3','#22dd11','darkgoldenrod2'), alpha=0.8), opacity = 0.2, n.br = 60)
dev.off()

pdf(file=paste('MuSSE_extrates_tree6404.pdf', sep=''), height=6, width=8)
profiles.plot(samples.f[5:7], lwd = 1, col.line = c('dodgerblue4','#22b211','darkgoldenrod3'), 
              col.fill = alpha(c('dodgerblue3','#22dd11','darkgoldenrod2'), alpha=0.8), opacity = 0.2, n.br = 60)
dev.off()

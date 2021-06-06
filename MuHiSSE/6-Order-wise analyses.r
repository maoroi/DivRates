### Order-wise analyses (Artiodactyla, Carnivora, Primates, Rodentia, Terrestrial-only[no bats, whales, sirenia)

## consider also: all terrestrial mammals, Euarchontoglires, Laurasiatheria, Afrotheria

setwd("C:/Users/Roi Maor/Desktop/2nd Chapter/DivRates/MuHiSSE")
act <- read.csv(file="ActivityData_MDD_v1_match.csv")

# 1. Prune tree and trim dataset ----------------------------------------------
trees <- c("0028", "1166", "1845", "2966", "3865", "4496", "5221", "5333", "5455", "6404", "7198", "9128")
orders <-  c("Artiodactyla","Carnivora","Primates","Rodentia")
for (k in 1:length(trees)) {
    tfile <- paste0("treevar",trees[k],".nex")
    tree <- read.tree(tfile)
    for (ord in orders) {
        trim <-  getMRCA(tree, act$Phylo_name[which(act$Order == ord)])
        dat <- act[which(act$Order == ord),c(6,5)]

    }
}



# 2. calculate sampling fractions ---------------------------------------------


# 3. set up relevant models ---------------------------------------------------


# 4. generate input and shell files for cluster -------------------------------




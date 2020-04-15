# Updating Chapter 3 to include Faurby & Svenning (2015) phylogeny and Burgin et al. (2018) taxonomy
require("ape")
require("phangorn")
require("phytools")

## preparing data
read.csv("MamTax2018.csv") -> TAX
TAX[which(TAX$extinct. == 0),c(3,4,8,21,34,35)] -> tax       # living species and relevant columns only 

read.tree("Small_phylogeny_4125_species.nex") -> FStree
clades <- prop.part(FStree)
MCCtree <- maxCladeCred(FStree, part = clades)

read.csv('Appendix_1-ActivityData.csv') -> DATA
DATA[,c(1:4)] -> dat
dat[,3] <- gsub(' ', '_', dat[,3], ignore.case = FALSE)
dat <- cbind(dat[,c(1:3)], dat[,c(3:4)])
colnames(dat) <- c(colnames(dat[1:3]), "MSW3", colnames(dat[ncol(dat)]))
dat$Binomial <- "NA"
data <- dat

length(which(tax$SciName %in% data$Binomial)) == 2199
length(which(tax$SciName %in% MCCtree$tip.label)) == 3583
length(which(data$Binomial %in% MCCtree$tip.label)) = 2060

# retain unchanged binomials
unchanged <- which(data$MSW3 %in% tax$SciName)
data$Binomial[unchanged] <- data$MSW3[unchanged]

# match updated binomials (there are a few classes of changes)
updated <- which(data$MSW3 %in% tax$IfTransfer_oldSciName) # this wrongly includes "Leptonycteris_yerbabuenae" which is in MSW3
for (i in 1:length(updated)){
    data$Binomial[updated[i]] <- as.character(tax$SciName[which(tax$IfTransfer_oldSciName == data$MSW3[updated[i]])])
}
# this is a real bitch because they never mention the previous genus of some transfers, e.g. Bison_bison)
genus_transfer<- which(data$MSW3 %in% TAX$SciName[which(TAX$genusTransfersinceMSW3.==1)])
for (i in 1:length(genus_transfer)){
   # data$Binomial[genus_transfer[i]] <- TAX$SciName[which(TAX[,] == genus_transfer) %% nrow(TAX)]
}

which(TAX[,] == paste0(strsplit(as.character(TAX$TaxonomyNotes[90]), 'moved to ', fixed = TRUE)[[1]][2],"_bison") 


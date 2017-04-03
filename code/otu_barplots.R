# OTU bar plots to look at community membership

pack_used <- c('randomForest','ggplot2', 'pROC', 'knitr','dplyr','AUCRF', 'tidyr', 'caret', 'RColorBrewer', 'reshape2', 'wesanderson')
for (dep in pack_used){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), repos = 'http://cran.us.r-project.org', 
                     quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

meta <- 'data/raw/kws_metadata.tsv'
shared <- 'data/mothur/kws_final.an.shared'
tax <- 'data/mothur/kws_final.an.cons.taxonomy'
subsample <- 'data/mothur/kws_final.an.0.03.subsample.shared'

meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)
tax_file <- read.table(file='data/mothur/kws_final.an.cons.taxonomy', sep = '\t', header=T, row.names=1)

shared_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')

#make OTU abundance file
#remove numOTUs column 
#Create df with relative abundances

test <- subset(shared_file, select = -c(numOtus, label))
rel_abund <- 100*test/unique(apply(test, 1, sum)) #unique is bad


#Create vector of OTUs with median abundances >1%
OTUs_1 <- apply(rel_abund, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]
#get df of just top OTUs
rel_abund_top <- rel_abund[, OTUs_1]
rel_meta <- merge(meta_file, rel_abund_top, by.x='group', by.y="row.names")

source('code/Sum_OTU_by_Tax.R')
source('code/sum_shared.R')

#use this code to assign phyla to each OTU in the shared file 
shared_phyla <- get_tax_level_shared(subsample, tax, 2)
phyla_met <- merge(meta_file, shared_phyla, by.x='group', by.y='row.names')

#get median of all OTUs by location
phyla_loc <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), median)
phyla_upper <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), FUN= quantile, probs =0.75)
phyla_lower <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), FUN= quantile, probs =0.25)

#could do all three of these lines in one line with summarize
#could do dplyr way, start with melted DF 

#only get top 6 phyla - could order by top 6 phyla 
phyla_loc <- phyla_loc[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]
phyla_upper <- phyla_upper[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]
phyla_lower <- phyla_lower[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]

#get rel abundance - subsampled to 4231- define as variable (or calculate )
rownames(phyla_loc) <- phyla_loc$Group.1
phyla_loc <- phyla_loc[,-1]
phyla_abund <- 100*phyla_loc/4231

rownames(phyla_lower) <- phyla_lower$Group.1
phyla_lower <- phyla_lower[,-1]
phyla_lower <- 100*phyla_lower/4231

rownames(phyla_upper) <- phyla_upper$Group.1
phyla_upper <- phyla_upper[,-1]
phyla_upper <- 100*phyla_upper/4231

#gather AND PLOT OMG :D
#put rownames back in their own column - you can get around this w dplyr 
phyla_abund <- cbind(group=rownames(phyla_abund), phyla_abund)
rownames(phyla_abund) <- c()
phylanames <- colnames(phyla_abund[,1:7])
phylamelt <- melt(phyla_abund[, phylanames], id.vars=1)

phyla_abund_up <- cbind(group=rownames(phyla_upper), phyla_upper)
rownames(phyla_abund_up) <- c()
phylamelt_up <- melt(phyla_abund_up[, phylanames], id.vars=1)

phyla_abund_low <- cbind(group=rownames(phyla_lower), phyla_lower)
rownames(phyla_abund_low) <- c()
phylamelt_low <- melt(phyla_abund_low[, phylanames], id.vars=1)

#merge them riiite
names(phylamelt_up)[3] <- "upper"
phylamelt <- merge(phylamelt, phylamelt_up)

names(phylamelt_low)[3] <- "lower"
phylamelt <- merge(phylamelt, phylamelt_low)

#aaand heres the plot! IT WORKS
#set positions to order bars by R to L 
positions <- c("RB", "RS", "LB", "LS", "SS")
ggplot(phylamelt, aes(x=group, y=value, ymin=lower, ymax=upper, fill=variable)) + 
  geom_bar(position=position_dodge(), stat='identity') + 
  geom_errorbar(position=position_dodge(0.9), width=0.2) + theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title= element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=14)) +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("R Mucosa", "R Lumen", "L Mucosa", "L Lumen", "Stool")) +
  theme(axis.title.x=element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1)) + scale_fill_brewer(palette="Dark2", name="Phylum") +
  ylab("% Relative Abundance") 




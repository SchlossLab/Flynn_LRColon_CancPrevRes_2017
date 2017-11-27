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

#Create df with relative abundances
test <- subset(shared_file, select = -c(numOtus, label))
rel_abund <- 100*test/unique(apply(test, 1, sum)) 

#Create vector of OTUs with median abundances >1%
OTUs_1 <- apply(rel_abund, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]

#get df of just top OTUs
rel_abund_top <- rel_abund[, OTUs_1]
rel_meta <- merge(meta_file, rel_abund_top, by.x='group', by.y="row.names")

source('code/Sum_OTU_by_Tax.R')
source('code/sum_shared.R')

#this code makes the df for family level RA sorted by % abundance 
#add column of sites to rel_abund_top
locat <- meta_file$location
rel_abund_top <- cbind(site = locat, rel_abund_top)

#this is what i want! 
rel_fam <- sum_OTU_by_tax_level(2, rel_abund_top, tax_file)

#add a column of locations to the df
loc_names <- meta_file$location
rel_fam[(ncol(rel_fam)+1)] <- loc_names
colnames(rel_fam)[ncol(rel_fam)] <- "location"

fam_med <- aggregate(rel_fam[, 1:(ncol(rel_fam)-1)], list(rel_fam$location), median)


barplot(t(fam_med))

#use this code to assign phyla to each OTU in the shared file 
shared_phyla <- get_tax_level_shared(subsample, tax, 2)
phyla_met <- merge(meta_file, shared_phyla, by.x='group', by.y='row.names')

shared_fam <- get_tax_level_shared(subsample, tax, 5)
fam_met <- merge(meta_file, shared_fam, by.x='group', by.y='row.names')


#try to get the df organized to work as a boxplot - no median calculation 
phyla_test <- phyla_met[, c("location","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]

#get median of all OTUs by location
phyla_loc <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), median)
phyla_upper <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), FUN= quantile, probs =0.75)
phyla_lower <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), FUN= quantile, probs =0.25)

fam_loc <- aggregate(fam_met[,7:ncol(fam_met)], list(fam_met$location), median)

#reorganize to make next calculations easier 
rownames(fam_loc) <- fam_loc$Group.1
fam_loc <- fam_loc[,-1]
fam_abund <- 100*fam_loc/4231

#do this to remove 0s

fam_abund <- fam_abund[,colSums(fam_abund^2) !=0]


#could do all three of these lines in one line with summarize
#could do dplyr way, start with melted DF 

#only get top 6 phyla - could order by top 6 phyla 
phyla_loc <- phyla_loc[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]
phyla_upper <- phyla_upper[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]
phyla_lower <- phyla_lower[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]

#get rel abundance - subsampled to 4231- define as variable (or calculate )
subsampled_to <- 4231
rownames(phyla_loc) <- phyla_loc$Group.1
phyla_loc <- phyla_loc[,-1]
phyla_abund <- 100*phyla_loc/subsampled_to

rownames(phyla_lower) <- phyla_lower$Group.1
phyla_lower <- phyla_lower[,-1]
phyla_lower <- 100*phyla_lower/subsampled_to

rownames(phyla_upper) <- phyla_upper$Group.1
phyla_upper <- phyla_upper[,-1]
phyla_upper <- 100*phyla_upper/subsampled_to

RA <- function(x) 100*x/subsampled_to

phyla_RA <- data.frame(phyla_test[1], apply(phyla_test[2:ncol(phyla_test)],2, RA))
phylaRAnames <- colnames(phyla_RA[,1:7])
phyla_RAmelt <- melt(phyla_RA[, phylaRAnames], id.vars=1)


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

fam_abund2 <- cbind(group=rownames(fam_abund), fam_abund)
rownames(fam_abund2) <- c()
famnames <- colnames(fam_abund2)
fam_melt <- melt(fam_abund2[, famnames ], id.vars=1)


#this is a stacked bar chart, because there are like 33 groups. need to rank first then plot top 10 or so

ggplot(fam_melt, aes(x=variable, y=value, fill=variable)) +geom_bar(stat = 'identity') +facet_wrap(~group) +theme_bw()

#ok still need to sort these by top 10 OTUs 

fam_ordered <- fam_abund[order(fam_abund[,1:ncol(fam_abund)]),]

rownames(fam_abund2) <- fam_abund2$group
fam_abund2 <- fam_abund2[,-1]

decrease_sort <- function(x){
  sort(x, decreasing = TRUE)
}

fam_ordered <- t(apply(fam_abund2, 1, FUN = function(x) decrease_sort(x)))

#a ridiculously long way to get the top ten families and select them. 
#like a really ridiculously long way to go about that - come on kaitlin

LB_fam <- subset(fam_abund2, rownames(fam_abund2) == 'LB')
rownames(LB_fam) <- LB_fam$group
LB_fam <- LB_fam[,-1]
LB_fam <- t(apply(LB_fam, 1, FUN = function(x) decrease_sort(x)))

LB_melted <- melt(LB_fam)
LB_ten <- LB_melted[1:9,]
LB_six <- LB_melted[1:6,]

#use the top ten column to select on overall dataset then facet!
topten <- as.character(LB_ten$Var2)
topsix <- as.character(LB_six$Var2)

fam_ten <- fam_abund2[, topten]
fam_ten <- cbind(group=rownames(fam_ten), fam_ten)
rownames(fam_ten) <- c()
famtennames <- colnames(fam_ten)
fam_tenmelt <- melt(fam_ten[, famtennames ], id.vars=1)

#to get RA without median for boxplot 

fam_RA <- fam_met[, c("location",topsix)]
family_RA <- data.frame(fam_RA[1], apply(fam_RA[2:ncol(fam_RA)],2, RA))
famRAnames <- colnames(family_RA[,1:6])
fam_RAmelt <- melt(family_RA[, famRAnames], id.vars=1)

# boxplot for top family 
positions <- c("RB", "RS", "LB", "LS", "SS")
ggplot(fam_RAmelt, aes(x=location, y=value)) + geom_boxplot(aes(color=variable)) + 
  scale_color_discrete(guide=FALSE)+
  geom_boxplot(aes(fill=variable), outlier.shape=21, outlier.size=2.5) + theme_bw() + 
  theme(axis.text = element_text(size= 16), axis.title= element_text(size=18), legend.text=element_text(size=14), legend.title=element_text(size=16)) +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("R Mucosa", "R Lumen", "L Mucosa", "L Lumen", "Stool")) +
  theme(axis.title.x=element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1)) + scale_fill_brewer(palette="Dark2", name="Family") +
  ylab("% Relative Abundance") 

ggplot(fam_tenmelt, aes(x=variable, y=value, fill=variable)) +geom_bar(stat = 'identity') +facet_wrap(~group) +theme_bw()


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

#for boxplot version
positions <- c("RB", "RS", "LB", "LS", "SS")
ggplot(phyla_RAmelt, aes(x=location, y=value)) + geom_boxplot(aes(color=variable)) + 
  scale_color_discrete(guide=FALSE)+
  geom_boxplot(aes(fill=variable), outlier.shape=21, outlier.size=2.5) + theme_bw() + 
  theme(axis.text = element_text(size= 16), axis.title= element_text(size=18), legend.text=element_text(size=14), legend.title=element_text(size=16)) +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("R Mucosa", "R Lumen", "L Mucosa", "L Lumen", "Stool")) +
  theme(axis.title.x=element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1)) + scale_fill_brewer(palette="Dark2", name="Phylum") +
  ylab("% Relative Abundance") 

#to add log scale add +scale_y_continuous(trans = "log10")

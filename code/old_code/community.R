# OTU bar plots to look at community membership

pack_used <- c('randomForest','ggplot2', 'pROC', 'knitr','dplyr','AUCRF', 'tidyr', 'caret', 'RColorBrewer', 'reshape2', 'wesanderson', 'vegan')
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
rel_abund <- 100*test/unique(apply(test, 1, sum))

#Create vector of OTUs with median abundances >1%
OTUs_1 <- apply(rel_abund, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]
#get df of just top OTUs
rel_abund_top <- rel_abund[, OTUs_1]
rel_abund_top <- na.omit(rel_abund_top)
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

#only get top 6 phyla
phyla_loc <- phyla_loc[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]
phyla_upper <- phyla_upper[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]
phyla_lower <- phyla_lower[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]

#get rel abundance - subsampled to 4231
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
#put rownames back in their own column 
phyla_abund <- cbind(group=rownames(phyla_abund), phyla_abund)
rownames(phyla_abund) <- c()
phylanames <- colnames(phyla_abund[,1:7])
phylamelt <- melt(phyla_abund[, phylanames], id.vars=1)

phyla_abund_up <- cbind(group=rownames(phyla_upper), phyla_upper)
rownames(phyla_abund_up) <- c()
phylamelt_up <- melt(phyla_abund_up[, phylanames], id.vars=1)

#this isnt working but 
phyla_abund_low <- cbind(group=rownames(phyla_lower), phyla_lower)
rownames(phyla_abund_low) <- c()
phylamelt_low <- melt(phyla_abund_low[, phylanames], id.vars=1)

#merge them riiite
names(phylamelt_up)[3] <- "upper"
phylamelt <- merge(phylamelt, phylamelt_up)

names(phylamelt_low)[3] <- "lower"
phylamelt <- merge(phylamelt, phylamelt_low)

colors <- as.list(wes_palette("Darjeeling"))
colors <- colors + as.list(wes_palette("Darjeeling2"))

#aaand heres the plot! IT WORKS
ggplot(phylamelt, aes(x=group, y=value, ymin=lower, ymax=upper, fill=variable)) + geom_bar(position=position_dodge(), stat='identity') + 
  geom_errorbar(position=position_dodge(0.9), width=0.2) +theme_bw() + 
  scale_x_discrete(breaks=c("LB", "LS", "RB", "RS", "SS"), labels=c("L Mucosa", "L Lumen", "R Mucosa", "R Lumen", "Stool")) +theme(axis.title.x=element_blank(), panel.grid.major.x = element_blank(),
                                                                                                                                   panel.grid.minor.x = element_blank()) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1)) + scale_fill_brewer(palette="Dark2", name="Phylum") +ylab("% Relative Abundance")



#actually need to aggregate/medians first from shared file, then put into sumotubytaxa
source('code/tax_level.R')
gathered <- gather(rel_meta, key, abundance, contains("Otu"))

phyla_otu <- get_tax(5, df=tax_file)

gathered_otu <- merge(gathered, phyla_otu, by.x='key', by.y='row.names')

ggplot(gathered_otu, aes(x=location, y=abundance, fill=tax)) + stat_summary(fun.y='mean', geom='bar', position='dodge')

#+geom_errorbar(aes(ymin=value-low_melt, ymax=value+low_melt))


rel_location <- aggregate(rel_meta[, 7:ncol(rel_meta)], list(rel_meta$location), median)
rownames(rel_location) <- rel_location$Group.1
rel_location <- rel_location[,-1]
phylum <- sum_OTU_by_tax_level(5, rel_location, tax_file)

#gonna need some error bars
phy_upper <- aggregate(rel_meta[, 7:ncol(rel_meta)], list(rel_meta$location), FUN= quantile, probs =0.75)
phy_lower <- aggregate(rel_meta[, 7:ncol(rel_meta)], list(rel_meta$location), FUN= quantile, probs =0.25)
rownames(phy_lower) <- phy_lower$Group.1
rownames(phy_upper) <- phy_upper$Group.1
phy_lower <- phy_lower[,-1]
phy_upper <- phy_upper[,-1]
phy_upper_error <- sum_OTU_by_tax_level(5, phy_upper, tax_file)
phy_lower_error <- sum_OTU_by_tax_level(5, phy_lower, tax_file)

#do i need to melt before doing error bars? lets find out 
phy_lower_error <- cbind(group=rownames(phy_lower_error), phy_lower_error)
lower_names <- colnames(phy_lower_error[,1:8])
rownames(phy_lower_error) <- c()
low_melt <- melt(phy_lower_error[,lower_names], id.var=1)

phy_upper_error <- cbind(group=rownames(phy_upper_error), phy_upper_error)
upper_names <- colnames(phy_upper_error[,1:8])
rownames(phy_upper_error) <- c()
upper_melt <- melt(phy_upper_error[,upper_names], id.var=1)

#now try plotting this 
phylumnames <- cbind(group=rownames(phylum), phylum)
phylanames <- colnames(phylumnames[,1:8])
rownames(phylumnames) <- c()
med_phyla <- melt(phylumnames[, phylanames], id.var=1)

#merge
med_error <- merge(med_phyla, low_melt, by=c('variable', 'group'))
med_error <- merge(med_error, upper_melt, by=c('variable', 'group'))


low_melt$value <- as.numeric(low_melt$value)
upper_melt$value <- as.numeric(upper_melt$value)

#WOO 
ggplot(med_phyla, aes(x=group, y=value)) +geom_bar(aes(fill=variable), position='dodge', stat='identity') +geom_errorbar(aes(ymin=low_melt, ymax=upper_melt))

ggplot(med_error, aes(x=group, y=value.x)) +geom_bar(aes(fill=variable), position='dodge', stat='identity') +geom_errorbar(aes(ymin=value.y, ymax=value))


#add error bars with geom_error_bar. do i get them from summary or what 



#Ok lets forget the community stuff for right now. let's plot inv. simpson diversity by location 
#need to redo inverse simpson with new shared, missing some groups 

invsimp <- read.table("data/mothur/kws_final.an.groups.summary", sep ='\t', header = T)
simpmeta <- merge(invsimp, meta_file, by.x='group', by.y='group')

ggplot(simpmeta, aes(x=location, y=invsimpson)) + geom_jitter(col=as.factor(simpmeta$patient), size=2) + theme_bw() + ggtitle("Simpson diversity by sample location")

ggplot(simpmeta, aes(x=side, y=invsimpson)) + geom_jitter(col=simpmeta$patient) + theme_bw() + ggtitle("Simpson diversity by sample side")

ggplot(simpmeta, aes(x=site, y=invsimpson)) + geom_jitter(col=simpmeta$patient) + theme_bw() + ggtitle("Simpson diversity by sample site")


#thetayc distances

tyc <- read.table("data/mothur/kws_final.an.summary", sep = '\t', header = T, row.names=NULL)

#can i quantify the distances between samples from the same patient and samples from each side? 

alltyc <- read.table("data/process/allshared.summary", sep = '\t', header = T, row.names=NULL)
alltyc <- separate(alltyc, label, into= c('pt1', 'samp1'), sep="-", remove=F)
alltyc <- separate(alltyc, comparison, into= c('pt2', 'samp2'), sep="-", remove=F)
alltyc <- alltyc[-1]
alltyc <- alltyc[-7]

#ultimately want a plot of all points where pt1 == pt2 in one bar and all of the others in another column 
#unite and make column of 0/1 for matches? then can plot 1 and 0s 
#should i separate out lumen and mucosa ? sure or no not for now

alltyc["same_pt"] <- NA

for (i in 1:nrow(alltyc)){
  if (alltyc$pt1[i] == alltyc$pt2[i]){
    alltyc$same_pt[i] <- 1
  }
  else alltyc$same_pt[i] <- 0
}

alltyc[10] <- as.factor(alltyc[10])

ggplot(alltyc, aes(x=as.factor(same_pt), y=thetayc)) + geom_jitter(width=0.15, shape=21, size=3, fill='grey', col='black') + theme_bw()+
  theme(legend.position="none", axis.title.x=element_blank(), axis.text = element_text(size= 16), axis.title= element_text(size=18)) +
  scale_x_discrete(labels=c("interpersonal", "intrapersonal")) +
  theme(axis.title.x=element_blank()) +ylab("ThetaYC distance")+ 
  stat_summary(aes(x=as.factor(same_pt), y=thetayc), data = alltyc, fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.4, col = "black")


wilcox.test(thetayc ~ same_pt, data = alltyc)


#separate column values for comparisons
tyc <- separate(tyc, label, into= c('pt1', 'samp1'), sep="-", remove=F)
tyc <- separate(tyc, comparison, into= c('pt2', 'samp2'), sep="-", remove=F)
tyc <- tyc[-1]
tyc <- tyc[-7]

#subset data to only include patients comparisons to each other
tyc <- subset(tyc, pt1==pt2)
tyc <- unite_(tyc, "match", from=c('samp1', 'samp2'), sep="_", remove = F)

#now plot- this is a plot of thetayc distances in all pairings for all patients 
ggplot(tyc, aes(x=match, y=thetayc)) + geom_point(aes(col=match)) +facet_wrap(~pt1) +theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(tyc, aes(x=match, y=thetayc)) + geom_boxplot() +theme_bw()

lefttyc <- subset(tyc, samp1== 'LB' | samp1== 'LS')
lefttyc <- subset(lefttyc, samp2== 'LB' | samp2== 'LS' | samp2=='SS')

righttyc <- subset(tyc, samp1=='RB' | samp1 == 'RS')
righttyc <- subset(righttyc, samp2=='RB' | samp2 == 'RS' | samp2 == 'SS')

mucosatyc <- subset(tyc, samp1=='RB' | samp1 == 'LB')
mucosatyc <- subset(mucosatyc, samp2=='RB' | samp2 == 'LB')

stooltyc <- subset(tyc, match=='LB_RB' | match== 'LS_RS')

exittyc <- subset(tyc, samp2 == 'SS')
exitLtyc <- subset(exittyc, samp1 == 'LS' | samp1 == 'LB')
exitRtyc <- subset(exittyc, samp1 == 'LS' | samp1 == 'RS')


leftandrighttyc <- subset(tyc, match=='LB_LS' | match== 'RB_RS')

lvsr <- rbind(stooltyc, leftandrighttyc)

ggplot(lefttyc, aes(x=match, y=thetayc)) + geom_boxplot() +theme_bw() +ggtitle("Left")

ggplot(righttyc, aes(x=match, y=thetayc)) + geom_boxplot() +theme_bw() +ggtitle("Right")

ggplot(mucosatyc, aes(x=match, y=thetayc)) + geom_boxplot() +theme_bw()

ggplot(exittyc, aes(x=match, y=thetayc, fill = match)) +geom_boxplot() +theme_bw()+ scale_fill_manual(values=wes_palette("Darjeeling2")) +theme(legend.position="none")
exittyc$match <- as.factor(exittyc$match)
kruskal.test(thetayc ~ match, data = exittyc)

ggplot(stooltyc, aes(x=match, y=thetayc, fill = match)) + geom_boxplot() +theme_bw() + scale_fill_manual(values=c("blue3", "goldenrod1")) +theme(legend.position="none")
stooltyc <- stooltyc[-25,]
wilcox.test(thetayc ~ match, data=stooltyc, paired=T)

ggplot(leftandrighttyc, aes(x=match, y=thetayc, fill=match)) + geom_boxplot() +theme_bw() +scale_fill_manual(values=c("blue3", "goldenrod1")) +theme(legend.position="none")
leftandrighttyc <- leftandrighttyc[-25,]
wilcox.test(thetayc ~ match, data=leftandrighttyc, paired = T)

#thetayc plot for figure 3
#lvsr$match <- factor(lvsr$match, levels=lvsr$match[order(lvsr$samp1)])
tycpositions <- c("RB_RS", "LS_RS", "LB_RB", "LB_LS")
ggplot(lvsr, aes(x=match, y=thetayc)) + geom_point() + geom_jitter(width= 0.2) +theme_bw() + 
  theme(legend.position="none", axis.title.x=element_blank(), axis.text = element_text(size= 16), axis.title= element_text(size=18)) + 
  scale_x_discrete(limits = tycpositions, breaks = tycpositions,
                  labels=c("R Mucosa vs R Lumen", "L Lumen vs R Lumen", "L Mucosa vs R Mucosa", "L Mucosa vs L Lumen")) +
  ylab("ThetaYC distance") + stat_summary(aes(x=match, y=thetayc), data = lvsr, fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.4)

summary(aov(lvsr$thetayc~lvsr$match))
wilcox.test(thetayc ~ match, data=lvsr, paired =T)




exitpositions <- c("RB_SS", "RS_SS", "LB_SS", "LS_SS")
ggplot(exittyc, aes(x=match, y=thetayc)) + geom_point() + geom_jitter(width= 0.2) +theme_bw() +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text = element_text(size= 16), axis.title= element_text(size=18)) +
  scale_x_discrete(limits=exitpositions, breaks=exitpositions,
                  labels=c("R Mucosa vs Stool", "R Lumen vs Stool", "L Mucosa vs Stool", "L Lumen vs Stool")) +
  theme(axis.title.x=element_blank()) +ylab("ThetaYC distance")+ 
  stat_summary(aes(x=match, y=thetayc), data = exittyc, fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.4)

summary(aov(exittyc$thetayc ~ exittyc$match))

#cant i just do a fuckin box plot?
boxplot(tyc[,"thetayc"] ~ tyc[,"match"])


#get and plot median of each match, gonna need error bars too 
 
medians <- matrix(nrow=1, ncol=2)
colnames(medians) <- c("median", "iqr")
medians <- as.data.frame(medians)






# OTU bar plots to look at community membership
#for KWS stuff

install.packages("RColorBrewer")
library(RColorBrewer)
meta <- 'data/raw/kws_metadata.tsv'
shared <- 'data/mothur/kws_final.an.shared'
tax <- 'data/mothur/kws_final.an.cons.taxonomy'
subsample <- 'data/mothur/kws_final.an.0.03.subsample.shared'


meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)
tax_file <- read.table(file='data/mothur/kws_final.an.cons.taxonomy', sep = '\t', header=T, row.names=1)

shared_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')

#make OTU abundance file
#Create df with relative abundances
rel_abund <- 100*shared_file/unique(apply(shared_file, 1, sum))
#remove numOTUs column
rel_abund <- rel_abund[,-1] 
rel_abund <- rel_abund[,-1] 

#Create vector of OTUs with median abundances >1%
OTUs_1 <- apply(rel_abund, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]
#get df of just top OTUs
rel_abund_top <- rel_abund[, OTUs_1]
rel_abund_top <- na.omit(rel_abund_top)

source('code/Sum_OTU_by_Tax.R')
source('code/sum_shared.R')

#use this code to assign phyla to each OTU in the shared file 
shared_phyla <- get_tax_level_shared(subsample, tax, 2)
phyla_met <- merge(meta_file, shared_phyla, by.x='group', by.y='row.names')

#get median of all OTUs by location
phyla_loc <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), median)
phyla_upper <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), FUN= quantile, probs =0.75)
phyla_lower <- aggregate(phyla_met[, 7:ncol(phyla_met)], list(phyla_met$location), FUN= quantile, probs =0.25)

#remove phylA we dont want
phyla_loc <- phyla_loc[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]
phyla_upper <- phyla_upper[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]
phyla_lower <- phyla_lower[, c("Group.1","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]

#get rel abundance
rownames(phyla_loc) <- phyla_loc$Group.1
phyla_loc <- phyla_loc[,-1]
phyla_abund <- 100*phyla_loc/1900

rownames(phyla_lower) <- phyla_lower$Group.1
phyla_lower <- phyla_lower[,-1]
phyla_lower <- 100*phyla_lower/1900

rownames(phyla_upper) <- phyla_upper$Group.1
phyla_upper <- phyla_upper[,-1]
phyla_upper <- 100*phyla_upper/1900


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


family <- sum_OTU_by_tax_level(2, rel_abund_top, tax_file)

#actually need to aggregate/medians first from shared file, then put into sumotubytaxa

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

#now try plotting this 
phylumnames <- cbind(group=rownames(phylum), phylum)
phylanames <- colnames(phylumnames[,1:8])
rownames(phylumnames) <- c()
med_phyla <- melt(phylumnames[, phylanames], id.var=1)

#merge

#WOO FINALLY FUCKIN WORKS 
ggplot(med_phyla, aes(x=group, y=value)) +geom_bar(aes(fill=variable), position='dodge', stat='identity') +geom_errorbar(aes(ymin=low_melt, ymax=low_melt))

#add error bars with geom_error_bar. do i get them from summary or what 

#If you are trying to do a base version of this plot 
family_melt <- melt(family, id.var='group')
family_names <- cbind(group=rownames(family), family)
transposed_family <- t(family)
transposed_family <- transposed_family[-1,]

#ggplot(transposed_family, aes(x=colnames(transposed_family), y))

#colorbrewer
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

barplot(t(family), ylim=c(0,100), col=getPalette(nrow(family)))

#subset dataframe to get just sides 
rightside <- family_names[meta_file$side == 'right',]
rightside <- rightside[-1]
barplot(t(rightside), ylim=c(0,100), col=getPalette(nrow(family)), cex.axis = 0.8, cex.names = 0.8)

#right side is still too busy. lets just do one pt at a time

pt6 <- family_names[meta_file$patient == 6,]
pt6 <- pt6[-1]
barplot(t(pt6), ylim=c(0,100), col = getPalette(nrow(family)))


#Ok lets forget the community stuff for right now. let's plot inv. simpson diversity by location 
#need to redo inverse simpson with new shared, missing some groups 

invsimp <- read.table("data/mothur/kws_final.an.groups.summary", sep ='\t', header = T)
simpmeta <- merge(invsimp, meta_file, by.x='group', by.y='group')

ggplot(simpmeta, aes(x=location, y=invsimpson)) + geom_jitter(col=as.factor(simpmeta$patient), size=2) + theme_bw() + ggtitle("Simpson diversity by sample location")

ggplot(simpmeta, aes(x=side, y=invsimpson)) + geom_jitter(col=simpmeta$patient) + theme_bw() + ggtitle("Simpson diversity by sample side")

ggplot(simpmeta, aes(x=site, y=invsimpson)) + geom_jitter(col=simpmeta$patient) + theme_bw() + ggtitle("Simpson diversity by sample site")


#thetayc distances

tyc <- read.table("data/mothur/kws_final.an.summary", sep = '\t', header = T, row.names=NULL)

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



#cant i just do a fuckin box plot?
boxplot(tyc[,"thetayc"] ~ tyc[,"match"])


#get and plot median of each match, gonna need error bars too 
 
medians <- matrix(nrow=1, ncol=2)
colnames(medians) <- c("median", "iqr")
medians <- as.data.frame(medians)

#this was only working for medians. maybe try piping all data like nick did - this works but doesnt plot yet. 
for(i in (tyc$match)){
  compar <- subset(tyc, match==i)
  theta <- compar["thetayc"]
  med <- apply(theta, 2, median)
  iqr <- apply(theta, 2, IQR)
  medians[i, "median"] <- med
  medians[i, "iqr"] <- iqr
  }

medians <- as.data.frame(medians)
medians <- add_rownames(medians, var = "comparison")

ggplot(medians, aes(x=comparison, y=medians)) +geom_point()







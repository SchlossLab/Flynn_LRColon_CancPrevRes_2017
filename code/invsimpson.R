# Inverse simpson and theta yc plots 

pack_used <- c('ggplot2','dplyr', 'tidyr', 'RColorBrewer', 'reshape2', 'wesanderson')
for (dep in pack_used){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), repos = 'http://cran.us.r-project.org', 
                     quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)
tax_file <- read.table(file='data/mothur/kws_final.an.cons.taxonomy', sep = '\t', header=T, row.names=1)
invsimp <- read.table("data/mothur/kws_final.an.groups.summary", sep ='\t', header = T)
tyc <- read.table("data/mothur/kws_final.an.summary", sep = '\t', header = T, row.names=NULL)

simpmeta <- merge(invsimp, meta_file, by.x='group', by.y='group')

simpmeta$location <- factor(simpmeta$location, c("LB","RB", "LS", "RS", "SS"))
positions <- c("RB", "RS", "LB", "LS", "SS")
ggplot(simpmeta, aes(x=location, y=invsimpson, group =1)) +geom_point() +geom_jitter(width=0.2) +theme_bw() + ylab("Inverse Simpson Diversity") +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("R Mucosa", "R Lumen", "L Mucosa", "L Lumen", "Stool")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text = element_text(size= 16), axis.title= element_text(size=18)) +
  stat_summary(aes(x=location, y=invsimpson), data = simpmeta, fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.4)


#thetayc distances

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
ggplot(lvsr, aes(x=match, y=thetayc)) + geom_point() + geom_jitter(width= 0.2) +theme_bw() +
  scale_x_discrete(labels=c("L Mucosa vs L Lumen", "L Mucosa vs R Lumen", "L Lumen vs R Lumen", "R Mucosa vs R Lumen")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text = element_text(size= 16), axis.title= element_text(size=18)) +
  ylab("ThetaYC distance") + stat_summary(aes(x=match, y=thetayc), data = lvsr, fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.4)



ggplot(exittyc, aes(x=match, y=thetayc, color=match)) + geom_point() + geom_jitter(width= 0.5) +theme_bw() +scale_color_brewer(palette="Set1") +
  theme(legend.position="none") +scale_x_discrete(labels=c("L Mucosa vs Stool", "L Lumen vs Stool", "R Mucosa vs Stool", "R Lumen vs Stool")) +
  theme(axis.title.x=element_blank()) +ylab("ThetaYC distance")

#trying to get adonis to work for comparisons 

#for now just do paired wilcoxon with multiple comparisons 

pvalues <- c()

Atyc <- subset(tyc, match=='RB_RS' | match=='LS_RS')
pvalues <- c(pvalues, wilcox.test(thetayc~match, data=Atyc, paired=T)$p.value)

btyc <- subset(tyc, match=='RB_RS'| match=='LB_RB')
btyc <- btyc[-25,]
pvalues <- c(pvalues, wilcox.test(thetayc~match, data=btyc, paired=T)$p.value)

ctyc <- subset(tyc, match=='RB_RS'| match=='LB_LS')
ctyc <- ctyc[-25,]
pvalues <- c(pvalues, wilcox.test(thetayc~match, data=ctyc, paired=T)$p.value)

dtyc <- subset(tyc, match == 'LS_RS' | match == 'LB_RB')
dtyc <- dtyc[-25,]
pvalues <- c(pvalues, wilcox.test(thetayc~match, data=dtyc, paired=T)$p.value)

etyc <- subset(tyc, match == 'LS_RS' | match == 'LB_LS')
etyc <- etyc[-25,]
pvalues <- c(pvalues, wilcox.test(thetayc~match, data=etyc, paired=T)$p.value)

ftyc <- subset(tyc, match == 'LB_RB' | match == 'LB_LS')
pvalues <- c(pvalues, wilcox.test(thetayc~match, data=ftyc, paired=T)$p.value)

pvalues <- p.adjust(pvalues, method = "BH")

# now for exit comparisons

stoolpvalues <- c()

htyc <- subset(tyc, match=='RB_SS' | match=='RS_SS')
htyc <- htyc[-25,]
stoolpvalues <- c(stoolpvalues, wilcox.test(thetayc~match, data=htyc, paired=T)$p.value)

ityc <- subset(tyc, match=='RB_SS' | match=='LB_SS')
stoolpvalues <- c(stoolpvalues, wilcox.test(thetayc~match, data=ityc, paired=T)$p.value)

jtyc <- subset(tyc, match=='RB_SS' | match=='LS_SS')
stoolpvalues <- c(stoolpvalues, wilcox.test(thetayc~match, data=jtyc, paired=T)$p.value)

ktyc <- subset(tyc, match=='RS_SS' | match=='LB_SS')
ktyc <- ktyc[-25,]
stoolpvalues <- c(stoolpvalues, wilcox.test(thetayc~match, data=ktyc, paired=T)$p.value)

ltyc <- subset(tyc, match=='RS_SS' | match=='LS_SS')
ltyc <- ltyc[-25,]
stoolpvalues <- c(stoolpvalues, wilcox.test(thetayc~match, data=ltyc, paired=T)$p.value)

mtyc <- subset(tyc, match=='LB_SS' | match=='LS_SS')
stoolpvalues <- c(stoolpvalues, wilcox.test(thetayc~match, data=mtyc, paired=T)$p.value)

stoolpvalues <- p.adjust(stoolpvalues, method = "BH")

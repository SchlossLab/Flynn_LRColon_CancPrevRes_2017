###########################################################
# Builds Figure 3

###########################################################

pack_used <- c('ggplot2','dplyr', 'tidyr', 'RColorBrewer', 'reshape2', 'wesanderson', 'cowplot')
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
tyc <- read.table("data/mothur/kws_final.an.summary", sep = '\t', header = T, row.names=NULL)

#thetayc distances
#separate column values for comparisons
tyc <- separate(tyc, label, into= c('pt1', 'samp1'), sep="-", remove=F)
tyc <- separate(tyc, comparison, into= c('pt2', 'samp2'), sep="-", remove=F)
tyc <- subset(tyc, select = -c(row.names, X))

#subset data to only include patients comparisons to each other
tyc <- subset(tyc, pt1==pt2)
tyc <- unite_(tyc, "match", from=c('samp1', 'samp2'), sep="_", remove = F)

stooltyc <- subset(tyc, match=='LB_RB' | match== 'LS_RS')
leftandrighttyc <- subset(tyc, match=='LB_LS' | match== 'RB_RS')
lvsr <- rbind(stooltyc, leftandrighttyc)

exittyc <- subset(tyc, samp2 == 'SS')

#plots
tycpositions <- c("RB_RS", "LS_RS", "LB_RB", "LB_LS")
lvr_plot <- ggplot(lvsr, aes(x=match, y=thetayc)) + geom_boxplot(width=0.8) +theme_bw() + 
  theme(legend.position="none", axis.title.x=element_blank(), axis.text = element_text(size= 10), axis.title= element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  scale_x_discrete(limits = tycpositions, breaks = tycpositions,
                   labels=c("P Mucosa vs P Lumen", "D Lumen vs P Lumen", "D Mucosa vs P Mucosa", "D Mucosa vs D Lumen")) +
  ylab("ThetaYC distance") 



exitpositions <- c("RB_SS", "RS_SS", "LB_SS", "LS_SS")
exit_plot <- ggplot(exittyc, aes(x=match, y=thetayc)) + geom_boxplot(width=0.8) +theme_bw() +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text = element_text(size= 10), axis.title= element_text(size=12), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(limits=exitpositions, breaks=exitpositions,
                   labels=c("P Mucosa vs Stool", "P Lumen vs Stool", "D Mucosa vs Stool", "D Lumen vs Stool")) +
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

#####################################################################
# Intra / interpersonal comparison

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

inter_plot <- ggplot(alltyc, aes(x=as.factor(same_pt), y=thetayc)) + geom_boxplot(width=0.5) + theme_bw()+
  theme(legend.position="none", axis.title.x=element_blank(), axis.text = element_text(size= 10), axis.title= element_text(size=12)) +
  scale_x_discrete(labels=c("Interpersonal", "Intrapersonal")) +
  theme(axis.title.x=element_blank()) +ylab("ThetaYC distance")

inter_medians <- aggregate(thetayc ~ same_pt, alltyc, median)
  
wilcox.test(thetayc ~ same_pt, data = alltyc)

# build and export figure 
#export as PDF

plot_file <- '~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_3.pdf'
pdf(file=plot_file, width=12, height=11)
layout(matrix(c(1,
                2,
                3), 
              nrow=3, byrow = TRUE))

lvr_plot 
ggdraw(lvr_plot) + 
  draw_plot_label("A", size = 14)

exit_plot

inter_plot

dev.off()


#cowplot way- use this way!!


fig3 <- plot_grid(lvr_plot, exit_plot, inter_plot, labels = c("A", "B", "C"), ncol = 1, align = "v")  
save_plot('~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_3.pdf', fig3, ncol=1, nrow=3, base_width=8)



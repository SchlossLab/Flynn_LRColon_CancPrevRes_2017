##########################################
# Builds Figure 2: OTUs and Simpson diversity by location 
#
#
##########################################
library(reshape2)

meta <- 'data/raw/kws_metadata.tsv'
shared <- 'data/mothur/kws_final.an.shared'
tax <- 'data/mothur/kws_final.an.cons.taxonomy'
subsample <- 'data/mothur/kws_final.an.0.03.subsample.shared'

meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)
tax_file <- read.table(file='data/mothur/kws_final.an.cons.taxonomy', sep = '\t', header=T, row.names=1)
shared_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')

source('code/Sum_OTU_by_Tax.R')
source('code/sum_shared.R')

############################################
# OTU boxplot

#use this code to assign phyla to each OTU in the shared file 
shared_phyla <- get_tax_level_shared(subsample, tax, 2)
phyla_met <- merge(meta_file, shared_phyla, by.x='group', by.y='row.names')

#try to get the df organized to work as a boxplot - no median calculation 
phyla_test <- phyla_met[, c("location","Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")]

subsampled_to <- 4231
RA <- function(x) 100*x/subsampled_to

phyla_RA <- data.frame(phyla_test[1], apply(phyla_test[2:ncol(phyla_test)],2, RA))
phylaRAnames <- colnames(phyla_RA[,1:7])
phyla_RAmelt <- melt(phyla_RA[, phylaRAnames], id.vars=1)


#for boxplot version
positions <- c("RB", "RS", "LB", "LS", "SS")
phy_plot <- ggplot(phyla_RAmelt, aes(x=location, y=value)) + geom_boxplot(aes(color=variable)) + 
  scale_color_discrete(guide=FALSE)+
  geom_boxplot(aes(fill=variable), outlier.shape=21) + theme_bw() + 
  theme(axis.text = element_text(size= 10), axis.title= element_text(size=12), panel.border = element_blank()) + 
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("P Mucosa", "P Lumen", "D Mucosa", "D Lumen", "Stool")) +
  theme(axis.title.x=element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_brewer(palette="Dark2", name="Phylum") +
  ylab("% Relative Abundance") 


###########################################
#Simpson diversity plot
simps <- read.table(file='data/mothur/kws_final.an.groups.summary', header = T)
simpmeta <- merge(meta_file, simps)

simpmeta$location <- factor(simpmeta$location, c("LB","RB", "LS", "RS", "SS"))
positions <- c("RB", "RS", "LB", "LS", "SS")
simp_plot <- ggplot(simpmeta, aes(x=location, y=invsimpson, group =1)) +geom_point() +geom_jitter(width=0.2) +theme_bw() + ylab("Inverse Simpson Diversity") +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("P Mucosa", "P Lumen", "D Mucosa", "D Lumen", "Stool")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text = element_text(size= 10), axis.title= element_text(size=12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(aes(x=location, y=invsimpson), data = simpmeta, fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.4)



##### Save and export plots for publication

theme_set(theme_minimal())

fig2 <- plot_grid(
  plot_grid(
    phy_plot + theme(legend.position = "none")
    , simp_plot
    , ncol = 1
    , align = "hv"
    , labels=c("A", "B")
    , label_size = 16)
  , plot_grid(
    get_legend(phy_plot)
    , ggplot()
    , ncol =1)
  , rel_widths = c(7,3)
)

save_plot('~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_2.pdf', fig2, ncol=1, nrow=2, base_width=6, base_height = 3)


#looking for cancer OTUs in KWS samples

#first gotta find the OTU numbers for the cancer bugs

#fuso OTU numbers: Otu00179, Otu00184, Otu00314, Otu00472, Otu01440

meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)
tax_file <- read.table(file='data/mothur/kws_final.an.cons.taxonomy', sep = '\t', header=T, row.names=1)
subsampled_file <- read.table(file='data/mothur/kws_final.an.0.03.subsample.shared', sep = '\t', header=T, row.names=2)

shared_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')
subsampled_meta <- merge(meta_file, subsampled_file, by.x='group', by.y='row.names')


only_fuso <- shared_meta[, colnames(shared_meta) %in% c('group', 'location', 'Otu00179', 'Otu00184', 'Otu00314', 'Otu00472', 'Otu01440')]

only_subsampled_fuso <- subsampled_meta[, colnames(subsampled_meta) %in% c('group', 'location', 'Otu00179', 'Otu00184', 'Otu00314', 'Otu00472', 'Otu01440')]

only_fuso[,8] <- only_fuso[,3] + only_fuso[,4] + only_fuso[,5] + only_fuso[,6] + only_fuso[,7]

names(only_fuso)[8] <- "total_fuso"

ggplot(only_fuso, aes(x=location, y=total_fuso)) + geom_point() +geom_jitter()

ggplot(only_fuso, aes(x=location, y=Otu00179)) +geom_point() + geom_jitter()

ggplot(only_subsampled_fuso, aes(x=location, y=Otu00179)) +geom_point() + geom_jitter() +theme_bw()


fuso179 <- subsampled_meta[, colnames(subsampled_meta) %in% c("group", "patient", "location", "Otu00179")]

fuso179[,5] <- (fuso179[,4]/4321)*100
names(fuso179)[5] <- "Otu00179_relAbund"

positions <- c("RB", "RS", "LB", "LS", "SS")
fuso_plot <- ggplot(fuso179, aes(x=location, y=Otu00179_relAbund)) +geom_jitter(width=0.3) +theme_bw() +
  ylab("F. nucleatum Relative Abundance") +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("P Mucosa", "P Lumen", "D Mucosa", "D Lumen", "Stool")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text = element_text(size= 14), axis.title= element_text(size=14))

#IBD fuso OTU - F. varium 

fusoviv <- subsampled_meta[, colnames(subsampled_meta) %in% c("group", "patient", "location", "Otu00472")]
fusoviv[,5] <- (fusoviv[,4]/4321)*100
names(fusoviv)[5] <- "Otu00472_relAbund"

positions <- c("RB", "RS", "LB", "LS", "SS")
fusoviv_plot <- ggplot(fusoviv, aes(x=location, y=Otu00472_relAbund)) +geom_jitter(width=0.3) +theme_bw() +
  ylab("F. varium Relative Abundance") +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("P Mucosa", "P Lumen", "D Mucosa", "D Lumen", "Stool")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text = element_text(size= 14), axis.title= element_text(size=14))



#other oncogenic otus - Otu00152 is P. asaccharolytica. Otu00248 is p micra

p_asach <- shared_meta[, colnames(shared_meta) %in% c("Group", "patient", "location", "Otu00152")]

p_asach[,4] <- (p_asach[,3]/4321) *100
names(p_asach)[4] <- "Otu00152_abund"

p_152 <-subsampled_meta[, colnames(subsampled_meta) %in% c("Group", "patient", "location", "Otu00152")]

p_152[,4] <- (p_152[,3]/4321) *100
names(p_152)[4] <- "Otu152_abund"

ggplot(p_152, aes(x=location, y=Otu152_abund)) +geom_jitter() +theme_bw()

ggplot(p_asach, aes(x=location, y=Otu00152_abund)) +geom_jitter() +theme_bw()

positions <- c("RB", "RS", "LB", "LS", "SS")
porphy_plot <- ggplot(p_152, aes(x=location, y=Otu152_abund)) +geom_jitter(width=0.3) +theme_bw() +
  ylab("P. asacharolytica Relative Abundance") +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("P Mucosa", "P Lumen", "D Mucosa", "D Lumen", "Stool")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text = element_text(size= 14), axis.title= element_text(size=14))



p_micra <- shared_meta[, colnames(shared_meta) %in% c("Group", "patient", "location", "Otu00248")]

p_micra[,4] <- (p_micra[,3]/4321) *100
names(p_micra)[4] <- "Otu00248_abund"

ggplot(p_micra, aes(x=location, y=Otu00248_abund)) +geom_jitter() +theme_bw()

########################################Build figure 6 
#export as PDF

plot_file <- '~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_6.pdf'
pdf(file=plot_file, width=9, height=6)
layout(matrix(c(1,
                2), 
              nrow=2, byrow = TRUE))

fuso_plot

fusoviv_plot

porphy_plot

dev.off()


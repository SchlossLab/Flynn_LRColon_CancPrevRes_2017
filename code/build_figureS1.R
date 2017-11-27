#####################################################
# Builds Figure S1 - oncogenic OTUs 
#####################################################

library(cowplot)

meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)
tax_file <- read.table(file='data/mothur/kws_final.an.cons.taxonomy', sep = '\t', header=T, row.names=1)
subsampled_file <- read.table(file='data/mothur/kws_final.an.0.03.subsample.shared', sep = '\t', header=T, row.names=2)

shared_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')
subsampled_meta <- merge(meta_file, subsampled_file, by.x='group', by.y='row.names')

#Fusobacterium OTUs
#fuso OTU numbers: Otu00179, Otu00184, Otu00314, Otu00472, Otu01440

# Fusobacterium nucleatum plot - OTU00179
fuso179 <- subsampled_meta[, colnames(subsampled_meta) %in% c("group", "patient", "location", "Otu00179")]

fuso179[,5] <- (fuso179[,4]/4321)*100
names(fuso179)[5] <- "Otu00179_relAbund"

fuso_y_title <- expression(paste(italic("F. nucleatum"), " Rel. Abund."))
positions <- c("RB", "RS", "LB", "LS", "SS")
fuso_plot <- ggplot(fuso179, aes(x=location, y=Otu00179_relAbund)) +geom_jitter(width=0.3) +theme_bw() +
  ylab(fuso_y_title) +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("P Muc", "P Lum", "D Muc", "D Lum", "Feces")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.title.y=element_text(size=8),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_y_log10()


#IBD fuso OTU - F. varium OTU00472
fusoviv <- subsampled_meta[, colnames(subsampled_meta) %in% c("group", "patient", "location", "Otu00472")]
fusoviv[,5] <- (fusoviv[,4]/4321)*100
names(fusoviv)[5] <- "Otu00472_relAbund"

fusoviv_y_title <- expression(paste(italic("F. varium"), " Rel. Abund."))

positions <- c("RB", "RS", "LB", "LS", "SS")
fusoviv_plot <- ggplot(fusoviv, aes(x=location, y=Otu00472_relAbund)) +geom_jitter(width=0.3) +theme_bw() +
  ylab(fusoviv_y_title) +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("P Muc", "P Lum", "D Muc", "D Lum", "Feces")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.title.y=element_text(size=8), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_y_log10()


#Other oncogenic OTUs
#Otu00152 is P. asaccharolytica. 

p_152 <-subsampled_meta[, colnames(subsampled_meta) %in% c("Group", "patient", "location", "Otu00152")]

p_152[,4] <- (p_152[,3]/4321) *100
names(p_152)[4] <- "Otu152_abund"

porph_y_title <- expression(paste(italic("P. asacharolytica"), " Rel. Abund."))

positions <- c("RB", "RS", "LB", "LS", "SS")
porphy_plot <- ggplot(p_152, aes(x=location, y=Otu152_abund)) +geom_jitter(width=0.3) +theme_bw() +
  ylab(porph_y_title) +
  scale_x_discrete(limits = positions, breaks=positions, 
                   labels=c("P Muc", "P Lum", "D Muc", "D Lum", "Feces")) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.title.y=element_text(size=8),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_y_log10()


#Otu00248 is P.  micra
p_micra <- shared_meta[, colnames(shared_meta) %in% c("Group", "patient", "location", "Otu00248")]
p_micra[,4] <- (p_micra[,3]/4321) *100
names(p_micra)[4] <- "Otu00248_abund"

ggplot(p_micra, aes(x=location, y=Otu00248_abund)) +geom_jitter() +theme_bw() +scale_y_log10()

########################################
#Build figure S1 and export as PDF

figS1 <- plot_grid(fuso_plot, fusoviv_plot, porphy_plot, labels = c("A", "B", "C"), ncol = 1, align = "v")  
save_plot('~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_S1.pdf', figS1, ncol=1, nrow=3, base_width = 3, base_height=2)



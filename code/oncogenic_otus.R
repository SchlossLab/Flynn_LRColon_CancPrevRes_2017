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

ggplot(fuso179, aes(x=location, y=Otu00179_relAbund)) +geom_jitter(width=0.3) +theme_bw()

#other oncogenic otus - Otu00152 is P. asaccharolytica. Otu00248 is p micra

p_asach <- shared_meta[, colnames(shared_meta) %in% c("Group", "patient", "location", "Otu00152")]

p_asach[,4] <- (p_asach[,3]/4321) *100
names(p_asach)[4] <- "Otu00152_abund"

p_152 <-subsampled_meta[, colnames(shared_meta) %in% c("Group", "patient", "location", "Otu00152")]

p_152[,4] <- (p_152[,3]/4321) *100
names(p_152)[4] <- "Otu152_abund"

ggplot(p_152, aes(x=location, y=Otu152_abund)) +geom_jitter() +theme_bw()

ggplot(p_asach, aes(x=location, y=Otu00152_abund)) +geom_jitter() +theme_bw()

p_micra <- shared_meta[, colnames(shared_meta) %in% c("Group", "patient", "location", "Otu00248")]

p_micra[,4] <- (p_micra[,3]/4321) *100
names(p_micra)[4] <- "Otu00248_abund"

ggplot(p_micra, aes(x=location, y=Otu00248_abund)) +geom_jitter() +theme_bw()


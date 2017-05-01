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

ggplot(only_subsampled_fuso, aes(x=location, y=Otu00179)) +geom_point() + geom_jitter()


#ok probably now should get seq and blast each otu to figure out which is nucleatum 

#do this with subsampled file maybe too? 
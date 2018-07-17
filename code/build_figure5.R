#########################################################################
# Build figure 5 and 6 - abundance charts of feature selected OTUs from RF models
# Kaitlin Flynn, Schloss lab, updated 12-20-17

# PUT TAX FILE AND SHARED META UP IN HURR
source('code/random_functions.R')
source('code/tax_level.R')

library(AUCRF)
library(dplyr)
library(pROC)

#load in all of the files, get rel abund of >1% 

meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)
tax_file <- read.table(file='data/mothur/kws_final.an.cons.taxonomy', sep = '\t', header=T, row.names=1)

subs_file <- read.table(file='data/mothur/kws_final.an.0.03.subsample.shared', sep = '\t', header = T, row.names=2)

#make OTU abundance file
#Create df with relative abundances
shared_file <- subset(shared_file, select = -c(numOtus, label))
shared_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')
rel_abund <- 100*shared_file/unique(apply(shared_file, 1, sum))

#do rel abund calcs for subsampled
subs_file <- subset(subs_file, select = -c(numOtus, label))
subs_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')
subs_abund <- 100*subs_file/unique(apply(subs_file, 1, sum))


#Create vector of OTUs with median abundances >1%
OTUs_1 <- apply(rel_abund, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]

OTUs_sub <- apply(subs_abund, 2, max) > 1
OTU_list_sub <- colnames(subs_abund)[OTUs_sub]

#get df of just top OTUs
rel_abund_top <- rel_abund[, OTUs_1]
rel_meta <- merge(meta_file, rel_abund_top, by.x='group', by.y="row.names")

subs_abund_top <- subs_abund[, OTUs_sub]
subs_meta <- merge(meta_file, subs_abund_top, by.x='group', by.y="row.names")

seed <- 1
n_trees <- 2001

## run these with subs meta and use as input for rel abund plots

aucrf_data_left_bs <- auc_loc(subs_meta, "LB", "LS")
aucrf_data_LRbowel <- auc_loc(subs_meta, "LB", "RB")
aucrf_data_right_bs <- auc_loc(subs_meta, "RB", "RS")
aucrf_data_LRlumen <- auc_loc(subs_meta, "LS", "RS")

######################################################## build and export figure 
#export as PDF

plot_file <- '~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_5.pdf'
pdf(file=plot_file, width=6, height=7)
layout(matrix(c(1,
                2), 
              nrow=2, byrow = TRUE))

#try exporting as jpg for poster


jpeg('~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_5.jpg', width=6, height=7, units='in', res=300)
layout(matrix(c(1,
                2), 
              nrow=2, byrow = TRUE))

#RB vs RS

right_otu_feat <- colnames(aucrf_data_right_bs[2:6])
otu_taxa_right <- get_tax(1, right_otu_feat, tax_file)
otu_taxa_right <- separate(otu_taxa_right, tax_label, into = c("OTU", "otu_num"), sep = "\\(")
formatted4 <- lapply(1:nrow(otu_taxa_right), function(i) bquote(paste(italic(.(otu_taxa_right$OTU[i])), "(", .(otu_taxa_right$otu_num[i]), sep=" ")))
#Abundance stripchart or most predictive otus
rs_abunds <- shared_meta[shared_meta$location=='RS', right_otu_feat]/10000 + 1e-4
rb_abunds <- shared_meta[shared_meta$location=='RB', right_otu_feat]/10000 + 1e-4

par(mar=c(5, 11, 1, 1))
plot(1, type="n", ylim=c(0,length(right_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in right_otu_feat){
  stripchart(at=index-0.35, jitter(rs_abunds[,i], amount=1e-5), pch=21, bg="white", method="jitter", jitter=0.2, add=T, lwd=0.5)
  stripchart(at=index+0.35, jitter(rb_abunds[,i], amount=1e-5), pch=21, bg="gray29", method="jitter", jitter=0.2, add=T, lwd=0.5)
  segments(median(rs_abunds[,i]),index-0.8,median(rs_abunds[,i]),index, lwd=3)
  segments(median(rb_abunds[,i]),index+0.8,median(rb_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=do.call(expression,formatted4), las=1, line=-0.5, tick=F, cex.axis=0.9)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("P Muc", "P Lum"), pch=c(21, 21), pt.bg=c("gray29","white"), cex=0.7)

mtext('A', side=2, line=7.5, las=1, adj=2, padj=-4.5, cex=2, font=2)

#just LB vs LS 
left_otu_feat <- colnames(aucrf_data_left_bs[2:6])
otu_taxa_left <- get_tax(1, left_otu_feat, tax_file)
otu_taxa_left <- separate(otu_taxa_left, tax_label, into = c("OTU", "otu_num"), sep = "\\(")
formatted3 <- lapply(1:nrow(otu_taxa_left), function(i) bquote(paste(italic(.(otu_taxa_left$OTU[i])), "(", .(otu_taxa_left$otu_num[i]), sep=" ")))
#Abundance stripchart or most predictive otus
ls_abunds <- shared_meta[shared_meta$location=='LS', left_otu_feat]/10000 + 1e-4
lb_abunds <- shared_meta[shared_meta$location=='LB', left_otu_feat]/10000 + 1e-4

par(mar=c(5, 11, 1, 1))
plot(1, type="n", ylim=c(0,length(left_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in left_otu_feat){
  stripchart(at=index-0.35, jitter(ls_abunds[,i], amount=1e-5), pch=21, bg="white", method="jitter", jitter=0.2, add=T, lwd=0.5)
  stripchart(at=index+0.35, jitter(lb_abunds[,i], amount=1e-5), pch=21, bg="gray29", method="jitter", jitter=0.2, add=T, lwd=0.5)
  segments(median(ls_abunds[,i]),index-0.8,median(ls_abunds[,i]),index, lwd=3)
  segments(median(lb_abunds[,i]),index+0.8,median(lb_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted3), las=1, line=-0.5, tick=F, cex.axis=0.9)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("D Muc", "D Lum"), pch=c(21, 21), pt.bg=c("gray29","white"), cex=0.7)

mtext('B', side=2, line=7.5, las=1, adj=2, padj=-4.5, cex=2, font=2)

dev.off()


######### Build figure 6

####################################### Figure 6 - updated 
#export as PDF

plot_file <- '~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_6.pdf'
pdf(file=plot_file, width=6, height=7)
layout(matrix(c(1,
                2), 
              nrow=2, byrow = TRUE))

#for exporting as Jpeg for poster
jpeg('~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_6.jpg', width=6, height=7, units='in', res=300)
layout(matrix(c(1,
                2), 
              nrow=2, byrow = TRUE))


#Lb vs Rb
LRbowel_otu_feat <- colnames(aucrf_data_LRbowel[2:6])
otu_taxa_LRbowel <- get_tax(1, LRbowel_otu_feat, tax_file)
otu_taxa_LRbowel <- separate(otu_taxa_LRbowel, tax_label, into = c("OTU", "otu_num"), sep = "\\(")
formatted1 <- lapply(1:nrow(otu_taxa_LRbowel), function(i) bquote(paste(italic(.(otu_taxa_LRbowel$OTU[i])), "(", .(otu_taxa_LRbowel$otu_num[i]), sep=" ")))
#Abundance stripchart or most predictive otus 
lb_abunds <- shared_meta[shared_meta$location=='LB', LRbowel_otu_feat]/10000 + 1e-4
rblb_abunds <- shared_meta[shared_meta$location=='RB', LRbowel_otu_feat]/10000 + 1e-4

par(mar=c(5, 11, 1, 1))
plot(1, type="n", ylim=c(0,length(LRbowel_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in LRbowel_otu_feat){
  stripchart(at=index-0.35, jitter(lb_abunds[,i], amount=1e-5), pch=21, bg="gray29", method="jitter", jitter=0.2, add=T, lwd=0.5)
  stripchart(at=index+0.35, jitter(rblb_abunds[,i], amount=1e-5), pch=21, bg="white", method="jitter", jitter=0.2, add=T, lwd=0.5)
  segments(median(lb_abunds[,i]),index-0.8,median(lb_abunds[,i]),index, lwd=3)
  segments(median(rblb_abunds[,i]),index+0.8,median(rblb_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=do.call(expression,formatted1), las=1, line=-0.5, tick=F, cex.axis=0.9)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("D Muc", "P Muc"), pch=c(21, 21), pt.bg=c("gray29","white"), cex=0.7)

mtext('A', side=2, line=7.5, las=1, adj=2, padj=-4.5, cex=2, font=2)


#LS vs RS
LRlumen_otu_feat <- colnames(aucrf_data_LRlumen[2:5])
otu_taxa_LRlumen <- get_tax(1, LRlumen_otu_feat, tax_file)
otu_taxa_LRlumen <- separate(otu_taxa_LRlumen, tax_label, into = c("OTU", "otu_num"), sep = "\\(")
formatted <- lapply(1:nrow(otu_taxa_LRlumen), function(i) bquote(paste(italic(.(otu_taxa_LRlumen$OTU[i])), "(", .(otu_taxa_LRlumen$otu_num[i]), sep=" ")))
#Abundance stripchart or most predictive otus 
lsrs_abunds <- shared_meta[shared_meta$location=='LS', LRlumen_otu_feat]/10000 + 1e-4
rsls_abunds <- shared_meta[shared_meta$location=='RS', LRlumen_otu_feat]/10000 + 1e-4

par(mar=c(5, 11, 1, 1))
plot(1, type="n", ylim=c(0,length(LRlumen_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in LRlumen_otu_feat){
  stripchart(at=index-0.35, jitter(lsrs_abunds[,i], amount=1e-5), pch=21, bg="gray29", method="jitter", jitter=0.2, add=T, lwd=0.5)
  stripchart(at=index+0.35, jitter(rsls_abunds[,i], amount=1e-5), pch=21, bg="white", method="jitter", jitter=0.2, add=T, lwd=0.5)
  segments(median(lsrs_abunds[,i]),index-0.8,median(lsrs_abunds[,i]),index, lwd=3)
  segments(median(rsls_abunds[,i]),index+0.8,median(rsls_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=0.9)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("D Lum", "P Lum"), pch=c(21, 21), pt.bg=c("gray29","white"), cex=0.7)

mtext('B', side=2, line=7.5, las=1, adj=2, padj=-4.5, cex=2, font=2)

dev.off()



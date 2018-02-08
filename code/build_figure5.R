#########################################################################
# Build figure 5 - abundance charts of feature selected OTUs from RF models
# Kaitlin Flynn, Schloss lab, updated 12-20-17

# PUT TAX FILE AND SHARED META UP IN HURR
source('code/random_functions.R')
source('code/tax_level.R')
#from optimized model code that is run on the cluster, import file with OTUs and abundances 

r_top_otus <- read.table(file='data/process/r_imp_otus.tsv')

#####Relative abundance plots#####
#get top 5 OTUs and plot relative abundance

r_top_names <- colnames(r_top_otus[2:6])
r_top_labels <- get_tax(1, r_top_names, tax_file)

rs_abunds <- shared_meta[shared_meta$location=='RS', r_top_names]/10000 + 1e-4
rb_abunds <- shared_meta[shared_meta$location=='RB', r_top_names]/10000 + 1e-4

### this isn't pulling the correct top OTUs- probably because its only saved from one iteration of the model in the loop 

par(mar=c(5, 15, 1, 1))
plot(1, type="n", ylim=c(0,length(r_top_otus)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n", cex.lab=1.5)
index <- 1
for(i in r_top_otus){
  stripchart(at=index-0.35, jitter(rs_abunds[,i], amount=1e-5), pch=21, bg="purple", method="jitter", jitter=0.2, add=T, cex=1.2, lwd=0.5)
  stripchart(at=index+0.35, jitter(rb_abunds[,i], amount=1e-5), pch=21, bg="orange", method="jitter", jitter=0.2, add=T, cex=1.2, lwd=0.5)
  segments(median(rs_abunds[,i]),index-0.7,median(rs_abunds[,i]),index, lwd=3)
  segments(median(rb_abunds[,i]),index+0.7,median(rb_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=r_top_labels$tax_label, las=1, line=-0.5, tick=F, cex.axis=1.2)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"), cex.axis=1.2)
legend('topright', legend=c("Right mucosa", "Right lumen"), pch=c(21, 21), pt.bg=c("orange","purple"), cex=1.2)

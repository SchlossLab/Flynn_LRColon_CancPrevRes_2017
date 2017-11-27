#########################################################################
# Build figure 4 - Random Forest model building and ROCs
# Kaitlin Flynn, Schloss lab, updated 11-27-17

library(randomForest)
library(AUCRF)

#load in all of the files, get rel abund of >1% and/or >10%

meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)
tax_file <- read.table(file='data/mothur/kws_final.an.cons.taxonomy', sep = '\t', header=T, row.names=1)

#make OTU abundance file
#Create df with relative abundances
shared_file <- subset(shared_file, select = -c(numOtus, label))
shared_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')

rel_abund <- 100*shared_file/unique(apply(shared_file, 1, sum))

#Create vector of OTUs with median abundances >1%
OTUs_1 <- apply(rel_abund, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]
#get df of just top OTUs
rel_abund_top <- rel_abund[, OTUs_1]
rel_meta <- merge(meta_file, rel_abund_top, by.x='group', by.y="row.names")

print("loaded and organized data")

seed <- 1
n_trees <- 2001

source('code/random_functions.R')
source('code/tax_level.R')

# Modeling approach for leave-one-out
# iteratively, do the following:
#1) select 19 samples (patients) from the 20, hold out the other one
#2) train the model using AUCRF and cross validation.
#3) from that model, will have an AUCRF object that can be used to predict/test on the held out sample
#4) record the n of variables used in each model, the AUC of the training model and the mtry. Store in a list the value (0/1) that was predicted for the sample from the model
#5) do that 20x
#6) aggregate all of the models from the left-out-list to get the AUC table to plot a curve and get an AUC value for the whole model
# do that for all of the models we have

# Example of the code used for the distal models
# code needs to be run on a high performance computer cluster - takes too much memory and time

################ Run the following on HPCC / flux
testsub <- subset(rel_meta, location %in% c("LB", "LS"))
testsub$location <- factor(testsub$location)
levels(testsub$location) <- c(1:length(levels(testsub$location))-1)

#create empty list to store in
test_results <- data.frame(Patient = character(), zero = double(), one = double())
held_out_results <- data.frame()
for(p in unique(testsub$patient)){
  print(p)
  test_set <- subset(testsub, testsub$patient != p)
  held_out <- subset(testsub, testsub$patient == p)
  rf_testset <- AUCRF(location~., data=select(test_set, location, contains("Otu")), ntree=n_trees, pdel=0.05, ranking="MDA")
  aucrf_test <- AUCRFcv(rf_testset, nCV=10, M=20)
  test_held_out <- predict(aucrf_test$RFopt, held_out, type='prob')[,2]
  held_out_results <- as.data.frame(test_held_out)
  held_out_results$Patient <- p
  colnames(held_out_results) <- c("zero", "one", "Patient")
  test_results <- rbind(test_results, held_out_results)
}

left_optimized <- write.table(test_results, file = 'data/process/left_optimized.tsv', sep = '\t')

#################

#do for all models!

left_optimized_results <- read.table(file = 'data/process/left_optimized.tsv', sep = '\t')
right_optimized_results <- read.table(file = 'data/process/right_optimized.tsv', sep = '\t')

bowel_optimized_results <- read.table(file = 'data/process/bowel_optimized.tsv', sep = '\t')
stool_optimized_results <- read.table(file = 'data/process/stool_optimized.tsv', sep = '\t')

left_roc <- roc(testsub$location ~ left_optimized_results$one)
right_roc <- roc(testsubR$location ~ right_optimized_results$one)

muc_roc <- roc(testsubM$location ~ bowel_optimized_results$one)
stool_roc <- roc(testsubS$location ~ stool_optimized_results$one)

#########################################################################
# Save and export Figure 4

plot_file <- '~/Documents/Flynn_LRColon_XXXX_2017/submission/figure_4.pdf'
pdf(file=plot_file, width=4, height=6)
layout(matrix(c(1,
                2), 
              nrow=2, byrow = TRUE))

#plots for figures 
#Left and right plot
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(left_roc, col='red', add=T, lty=1)
plot(right_roc, col = 'blue', add=T, lty=1)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottomright', legend=c(
  sprintf('D Lum vs D Muc, AUC = 0.908'),
  sprintf('P Lum vs P Muc, AUC = 0.764')
),lty=1, lwd = 2, cex=0.7, col=c('red', 'blue'), bty='n')

mtext('A', side=2, line=2, las=1, adj=1.5, padj=-5, cex=1.5, font=2)

#stool and mucosa plot 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(muc_roc, col='darkgreen', add=T, lty=1)
plot(stool_roc, col = 'purple', add=T, lty=1)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottomright', legend=c(
  sprintf('D Muc vs P Muc, AUC = 0.850'),
  sprintf('D Lum vs P Lum, AUC = 0.580')
), lty=1, lwd=2, cex=0.7, col=c('darkgreen', 'purple'), bty='n')

mtext('B', side=2, line=2, las=1, adj=1.5, padj=-5, cex=1.5, font=2)

dev.off()



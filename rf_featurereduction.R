####optimizing Random Forest models by standardizing input OTUs

#Kaitlin Flynn, Schloss lab, updated 8 28 17
#load packages
pack_used <- c('randomForest','ggplot2', 'pROC', 'knitr','dplyr','AUCRF', 'tidyr', 'caret')
for (dep in pack_used){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), repos = 'http://cran.us.r-project.org', 
                     quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

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

seed <- 1
n_trees <- 2001

source('code/random_functions.R')
source('code/tax_level.R')


#####AUCRF#####################################################################################
###########Left mucosa (bowel) vs Right mucosa (bowel)
#testing AUCRFcv approach on unfiltered data
bowel_test <- subset(rel_meta, location %in% c("LB", "RB"))
bowel_test$location <- factor(bowel_test$location)
#change levels of variable of interest to 0/1
levels(bowel_test$location) <- c(1:length(levels(bowel_test$location))-1)
# create RF model
set.seed(seed)
rf_bowel_test <- AUCRF(location ~ ., data = select(bowel_test, location, contains("Otu")),
                       ntree = n_trees, pdel = 0.05, ranking = 'MDA')

#use aucrf object as input for AUCRFcv - cross validation command

aucrf_cv_bowel_unfiltered <- AUCRFcv(rf_bowel_test, nCV=10, M=100)

#use output of that in Optimal Set to get set of predictive OTUs and their probability of inclusion
optimal_bowel_unfiltered <- OptimalSet(aucrf_cv_bowel_unfiltered)

#then have file with OTUs and probabilities, select top 10, subset shared file for input

bowel_top10 <- as.vector(optimal_bowel_unfiltered$Name[1:10])

#subset shared to just be these columns then run model again 

bowel_top10_shared <- subset(rel_meta, select=colnames(rel_meta) %in% c('location',bowel_top10))
auc_bowel_top10 <- auc_loc(bowel_top10_shared, "LB", "RB")

######10fold cross v with reduced feature input
#10 fold cross validation for L vs R mucosa

iters <- 100
cv10f_aucs_muc10 <- c()
cv10f_all_resp_muc10 <- c()
cv10f_all_pred_muc10 <- c()
for(j in 1:iters){
  set.seed(j)
  sampling_muc <- sample(1:nrow(auc_bowel_top10),nrow(auc_bowel_top10),replace=F)
  cv10f_probs_muc10 <- rep(NA,39)
  for(i in seq(1,36,4)){
    train_muc <- auc_bowel_top10[sampling_muc[-(i:(i+3))],]
    test_muc <- auc_bowel_top10[sampling_muc[i:(i+3)],]
    set.seed(seed)
    temp_model_muc <- AUCRF(location~., data=train_muc, pdel=0.99, ntree=500)
    cv10f_probs_muc10[sampling_muc[i:(i+3)]] <- predict(temp_model_muc$RFopt, test_muc, type='prob')[,2]
  }
  cv10f_roc_muc10 <- roc(auc_bowel_top10$location~cv10f_probs_muc10)
  cv10f_all_pred_muc10 <- c(cv10f_all_pred_muc10, cv10f_probs_muc10)
  cv10f_all_resp_muc10 <- c(cv10f_all_resp_muc10, auc_bowel_top10$location)
  cv10f_aucs_muc10[j] <- cv10f_roc_muc10$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc_muc10 <- roc(cv10f_all_resp_muc10~cv10f_all_pred_muc10)


###########Left lumen (stool) vs Right lumen (stool)

lumen_test <- subset(rel_meta, location %in% c("LS", "RS"))
lumen_test$location <- factor(lumen_test$location)
levels(lumen_test$location) <- c(1:length(levels(lumen_test$location))-1)
set.seed(seed)
rf_lumen_test <- AUCRF(location ~ ., data = select(lumen_test, location, contains("Otu")),
                       ntree = n_trees, pdel = 0.05, ranking = 'MDA')

aucrf_cv_lumen_unfiltered <- AUCRFcv(rf_bowel_test, nCV=10, M=100)
optimal_lumen_unfiltered <- OptimalSet(aucrf_cv_lumen_unfiltered)

lumen_top10 <- as.vector(optimal_lumen_unfiltered$Name[1:10])
lumen_top10_shared <- subset(rel_meta, select=colnames(rel_meta) %in% c('location',lumen_top10))
auc_lumen_top10 <- auc_loc(lumen_top10_shared, "LS", "RS")


#10 fold cross validation for L vs R lumen
iters <- 100
cv10f_aucs_lum10 <- c()
cv10f_all_resp_lum10 <- c()
cv10f_all_pred_lum10 <- c()
for(j in 1:iters){
  set.seed(j)
  sampling_lum <- sample(1:nrow(auc_lumen_top10),nrow(auc_lumen_top10),replace=F)
  cv10f_probs_lum10 <- rep(NA,39)
  for(i in seq(1,36,4)){
    train_lum <- auc_lumen_top10[sampling_lum[-(i:(i+3))],]
    test_lum <- auc_lumen_top10[sampling_lum[i:(i+3)],]
    set.seed(seed)
    temp_model_lum <- AUCRF(location~., data=train_lum, pdel=0.99, ntree=500)
    cv10f_probs_lum10[sampling_lum[i:(i+3)]] <- predict(temp_model_lum$RFopt, test_lum, type='prob')[,2]
  }
  cv10f_roc_lum10 <- roc(auc_lumen_top10$location~cv10f_probs_lum10)
  cv10f_all_pred_lum10 <- c(cv10f_all_pred_lum10, cv10f_probs_lum10)
  cv10f_all_resp_lum10 <- c(cv10f_all_resp_lum10, auc_lumen_top10$location)
  cv10f_aucs_lum10[j] <- cv10f_roc_lum10$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc_lum10 <- roc(cv10f_all_resp_lum10~cv10f_all_pred_lum10)

########Now doing feature reduction to top 6 OTUs for the L vs L and R vs R models

#testing AUCRFcv approach on unfiltered data
left_test <- subset(rel_meta, location %in% c("LB", "LS"))
left_test$location <- factor(left_test$location)
#change levels of variable of interest to 0/1
levels(left_test$location) <- c(1:length(levels(left_test$location))-1)
# create RF model
set.seed(seed)
rf_left_test <- AUCRF(location ~ ., data = select(left_test, location, contains("Otu")),
                      ntree = n_trees, pdel = 0.05, ranking = 'MDA')

aucrf_cv_left_unfiltered <- AUCRFcv(rf_left_test, nCV=10, M=100)
optimal_left_unfiltered <- OptimalSet(aucrf_cv_left_unfiltered)

left_top6 <- as.vector(optimal_left_unfiltered$Name[1:6])
left_top6_shared <- subset(rel_meta, select=colnames(rel_meta) %in% c('location',left_top6))
auc_left_top6 <- auc_loc(left_top6_shared, "LB", "LS")


#10fold CV for L lumen vs L mucosa
iters <- 100
cv10f_aucs_left6 <- c()
cv10f_all_resp_left_bs6 <- c()
cv10f_all_pred_left_bs6 <- c()
for(j in 1:iters){
  set.seed(j)
  sampling <- sample(1:nrow(auc_left_top6),nrow(auc_left_top6),replace=F)
  cv10f_probs_left6 <- rep(NA,39)
  for(i in seq(1,36,4)){
    train_left_bs <- auc_left_top6[sampling[-(i:(i+3))],]
    test_left_bs <- auc_left_top6[sampling[i:(i+3)],]
    set.seed(seed)
    temp_model_left_bs <- AUCRF(location~., data=train_left_bs, pdel=0.99, ntree=500)
    cv10f_probs_left6[sampling[i:(i+3)]] <- predict(temp_model_left_bs$RFopt, test_left_bs, type='prob')[,2]
  }
  cv10f_roc_left_bs6 <- roc(auc_left_top6$location~cv10f_probs_left6)
  cv10f_all_pred_left_bs6 <- c(cv10f_all_pred_left_bs6, cv10f_probs_left6)
  cv10f_all_resp_left_bs6 <- c(cv10f_all_resp_left_bs6, auc_left_top6$location)
  cv10f_aucs_left6[j] <- cv10f_roc_left_bs6$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc_left_bs6 <- roc(cv10f_all_resp_left_bs6~cv10f_all_pred_left_bs6)

############R vs R with limited input
right_test <- subset(rel_meta, location %in% c("RB", "RS"))
right_test$location <- factor(right_test$location)
#change levels of variable of interest to 0/1
levels(right_test$location) <- c(1:length(levels(right_test$location))-1)
# create RF model
set.seed(seed)
rf_right_test <- AUCRF(location ~ ., data = select(right_test, location, contains("Otu")),
                      ntree = n_trees, pdel = 0.05, ranking = 'MDA')

aucrf_cv_right_unfiltered <- AUCRFcv(rf_right_test, nCV=10, M=100)
optimal_right_unfiltered <- OptimalSet(aucrf_cv_right_unfiltered)

right_top6 <- as.vector(optimal_right_unfiltered$Name[1:6])
right_top6_shared <- subset(rel_meta, select=colnames(rel_meta) %in% c('location',right_top6))
auc_right_top6 <- auc_loc(right_top6_shared, "RB", "RS")

#10fold CV for R lumen vs R mucosa
iters <- 100
cv10f_aucs_right6 <- c()
cv10f_all_resp_right_bs6 <- c()
cv10f_all_pred_right_bs6 <- c()
for(j in 1:iters){
  set.seed(j)
  sampling <- sample(1:nrow(auc_right_top6),nrow(auc_right_top6),replace=F)
  cv10f_probs_right6 <- rep(NA,39)
  for(i in seq(1,36,4)){
    train_right_bs <- auc_right_top6[sampling[-(i:(i+3))],]
    test_right_bs <- auc_right_top6[sampling[i:(i+3)],]
    set.seed(seed)
    temp_model_right_bs <- AUCRF(location~., data=train_right_bs, pdel=0.99, ntree=500)
    cv10f_probs_right6[sampling[i:(i+3)]] <- predict(temp_model_right_bs$RFopt, test_right_bs, type='prob')[,2]
  }
  cv10f_roc_right_bs6 <- roc(auc_right_top6$location~cv10f_probs_right6)
  cv10f_all_pred_right_bs6 <- c(cv10f_all_pred_right_bs6, cv10f_probs_right6)
  cv10f_all_resp_right_bs6 <- c(cv10f_all_resp_right_bs6, auc_right_top6$location)
  cv10f_aucs_right6[j] <- cv10f_roc_right_bs6$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc_right_bs6 <- roc(cv10f_all_resp_right_bs6~cv10f_all_pred_right_bs6)


##############Plots ##########

######left and right lumen and mucosa w limited input plot 

#####10fold plot left vs right mucosa and lumen plot 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='', cex.axis=1.5)
plot(cv10f_roc_muc10,col = 'green4', lwd=3, add=T, lty=1) #r vs l mucosa cross validation
plot(cv10f_roc_lum10, col = 'orange', lwd=3, add=T, lty=1) #r vs l lumen cross validation
mtext(side=2, text="True Positive (Sensitivity)", line=2.5, cex=1.2)
mtext(side=1, text="True Negative (Specificity)", line=2.5, cex=1.2)
legend('bottom', legend=c(sprintf('D mucosa vs P mucosa 10-fold CV, AUC = 0.912'),
                          sprintf('D lumen vs P lumen 10-fold CV, AUC = 0.6243'),
                          sprintf('Mucosa model vs Lumen model', roc.test(cv10f_all_resp_lum10~cv10f_all_pred_lum10)$p.value)
                          # sprintf('OOB vs 10-fold CV: p=%.2g', roc.test(otu_euth_roc,cv10f_roc)$p.value)
),lty=c(1, 1), lwd=2, col=c('green4', 'orange'), bty='n', cex=1.2)


#Lumen vs mucosa plot 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='', cex.axis=1.5)
plot(cv10f_roc_right_bs6, col='blue', lwd=3, add=T, lty=1)
#plot(cv10f_roc, col = 'purple', lwd=3, add=T, lty=1)
plot(cv10f_roc_left_bs6, col = 'red', lwd=3, add=T, lty=1)
mtext(side=2, text="True Positive (Sensitivity)", line=2.5, cex=1.5)
mtext(side=1, text="True Negative (Specificity)", line=2.5, cex=1.5)
legend('bottom', legend=c(#sprintf('Lumen vs Mucosa, 10-fold CV, AUC = 0.925'),
  sprintf('D Lumen vs D Mucosa, 10-fold CV, AUC ='),
  sprintf('P Lumen vs P Mucosa, 10-fold CV, AUC = ')
  #sprintf('OOB vs Leave-1-out: p=%.2g', roc.test(otu_euth_roc,LOO_roc)$p.value),
  #sprintf('OOB vs 10-fold CV: p=%.2g', roc.test(otu_euth_roc,cv10f_roc)$p.value)
),lty=c(1, 1, 1), lwd=3, col=c('red', 'blue'), bty='n', cex=1.2)







#####P value for comparing AUC/models#######

pval <- roc.test(cv10f_roc_muc10, cv10f_roc_lum10)





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

aucrf_cv_bowel_unfiltered <- AUCRFcv(rf_bowel_test, nCV=10, M=n_trees)

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



##############Plots ##########

######left and right lumen and mucosa




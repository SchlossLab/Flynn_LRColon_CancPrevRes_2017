# random forest and AUCRF of left vs right, stool vs mucosa


#load packages- code swiped from nick
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
shared_file <- shared_file[,-1]
shared_file <- shared_file[,-1]
rel_abund <- 100*shared_file/unique(apply(shared_file, 1, sum))

#Create vector of OTUs with median abundances >1%
OTUs_1 <- apply(rel_abund, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]
#get df of just top OTUs
rel_abund_top <- rel_abund[, OTUs_1]
rel_abund_top <- na.omit(rel_abund_top)

rel_meta <- merge(meta_file, rel_abund_top, by.x='group', by.y="row.names")

seed <- 1
n_trees <- 2001

#then can run random forest and just select the variables we want to predict on- but might be better to subset first...

left_bs <- subset(rel_meta, location %in% c("LB", "LS"))
left_bs$location <- factor(left_bs$location)

rf_left <- randomForest(location ~ ., data = select(left_bs, location, contains("Otu")), importance=T, ntree=n_trees)
#OOB estimate of  error rate: 10.26%
#Confusion matrix:
 # LB LS class.error
#LB 16  3   0.1578947
#LS  1 19   0.0500000

#so for left, these OTUS are higher in LB than in stool, also (see below) higher in LB than RB 
left_importance <- sort(importance(rf_left)[,1], decreasing = T)

right_bs <- subset(rel_meta, location %in% c("RB", "RS"))
right_bs$location <- factor(right_bs$location)

rf_right <- randomForest(location ~ ., data = select(right_bs, location, contains("Otu")), importance=T, ntree=n_trees)
#OOB estimate of  error rate: 53.85%
#Confusion matrix:
 # RB RS class.error
#RB 11  9   0.4500000
#RS 12  7   0.6315789

LR_bowel <- subset(rel_meta, location %in% c("LB", "RB"))
LR_bowel$location <- factor(LR_bowel$location)

rf_bowel <- randomForest(location ~ ., data = select(LR_bowel, location, contains("Otu")), importance=T, ntree=n_trees)

#OOB estimate of  error rate: 25.64%
#Confusion matrix:
 # LB RB class.error
#LB 13  6   0.3157895
#RB  4 16   0.2000000

#sort importance to show which OTUs are higher in Left than right bowel
bowel_importance <- sort(importance(rf_bowel)[,1], decreasing = T)


LR_lumen <- subset(rel_meta, location %in% c("LS", "RS"))
LR_lumen$location <- factor(LR_lumen$location)
rf_lumen <- randomForest(location ~ ., data = select(LR_lumen, location, contains("Otu")), importance=T, ntree=n_trees)

#OOB estimate of  error rate: 69.23%
#Confusion matrix:
 # LS RS class.error
#LS  7 13   0.6500000
#RS 14  5   0.7368421

#to get what is higher in LS vs RS 
lumen_importance <- sort(importance(rf_lumen)[,1], decreasing=T)

#ok now move on to looking at aucrf and roc curves 

#predicted
otu_euth_probs <- predict(rf_early_aucrf$RFopt, type = 'prob')

#compare real with predicted with ROC
otu_euth_roc <- roc(Predict_early_euth_df$Euth_Early ~ otu_euth_probs[ , 2])

#draw blank plot, plot each line over it individually 

#also run AUCRF
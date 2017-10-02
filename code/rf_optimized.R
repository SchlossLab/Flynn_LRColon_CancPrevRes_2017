### optimizing RF models by doing leave-one-out, per reviewers suggestions on marc's papers

#Kaitlin Flynn, Schloss lab, updated 9-23-17


#6) aggregate all of the models from the left-out-list to get the AUC table to plot a curve and get an AUC value for the whole model
# do that for all of the models we have

#Load packages, files, make OTU table 
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

#iteratively, do the following:
#1) select 19 samples (patients) from the 20, hold out the other one
#2) train the model using AUCRF and cross validation.
#3) from that model, will have an AUCRF object that i can use to predict/test on the held out sample
#4) record the # of variables used in each model, the AUC of the training model and the mtry. Store in a list the value (0/1) that was predicted for the sample from the model
#5) do that 20x

#testing 

testsub <- subset(rel_meta, location %in% c("LB", "LS"))
testsub$location <- factor(testsub$location)
levels(testsub$location) <- c(1:length(levels(testsub$location))-1)
p <- 6
for(p in testsub$patient){
  test_set <- subset(testsub, testsub$patient != p)
  held_out <- subset(testsub, testsub$patient == p)
  rf_testset <- AUCRF(location~., data=select(test_set, location, contains("Otu")), ntree=n_trees, pdel=0.05, ranking="MDA")
  aucrf_test <- AUCRFcv(rf_testset, nCV=10, M=20)
  test_held_out <- predict(aucrf_test$RFopt, held_out, type='prob')
  #store output in something here 
  }

aucrf_cv_left_bs <- AUCRFcv(rf_aucrf, nCV=10, M=20)

#wait ok lets try it using base randomForest package

rf_only  <- randomForest(location ~ ., data = select(test_set, location, contains("Otu")), importance = T, ntree=n_trees)
#add cross validation in here before testing held out set?
held_out_fixed <- held_out[,-3]
held_prediction <- predict(rf_only, held_out)

#so wait, isn't this already what I'm doing? 

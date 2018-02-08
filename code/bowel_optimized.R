### optimizing RF models by doing leave-one-out, per reviewers suggestions on marc's papers

#Kaitlin Flynn, Schloss lab, updated 11-27-17

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


print("loaded and organized data")

seed <- 1
n_trees <- 2001

source('code/random_functions.R')
source('code/tax_level.R')

testsub <- subset(subs_meta, location %in% c("LB", "RB"))
testsub$location <- factor(testsub$location)
levels(testsub$location) <- c(1:length(levels(testsub$location))-1)

print("about to start model loop")
#create empty list to store in
test_results <- data.frame(Patient = character(), zero = double(), one = double())
held_out_results <- data.frame()
for(p in unique(testsub$patient)){
  print(p)
  test_set <- subset(testsub, testsub$patient != p)
  held_out <- subset(testsub, testsub$patient == p)
  rf_testset <- AUCRF(location~., data=select(test_set, location, contains("Otu")), ntree=n_trees, pdel=0.05, ranking="MDA")
  aucrf_test <- AUCRFcv(rf_testset, nCV=10, M=20)
  test_held_out <- predict(aucrf_test$RFopt, held_out, type='prob')
  held_out_results <- as.data.frame(test_held_out)
  held_out_results$Patient <- p
  colnames(held_out_results) <- c("zero", "one", "Patient")
  test_results <- rbind(test_results, held_out_results)
}

bowel_optimized <- write.table(test_results, file = 'data/process/bowel_optimized.tsv', sep = '\t')

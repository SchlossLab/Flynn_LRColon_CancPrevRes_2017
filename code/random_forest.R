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

#all mucosa vs all lumen

all_lum <- subset(rel_meta, site %in% c("stool", "mucosa"))
all_lum$site <- factor(all_lum$site)
rf_all <- randomForest(site ~ ., data=select(all_lum, site, contains("Otu")), importance = T, ntree=n_trees)

# exit vs all lumen
exit_lum <- subset(rel_meta, site %in% c("stool", "exit"))
exit_lum$site <- factor(exit_lum$site)
rf_exitlum <- randomForest(site ~ ., data=select(exit_lum, site, contains("Otu")), importance=T, ntree=n_trees)

# exit vs all mucosa
exit_muc <- subset(rel_meta, site %in% c("mucosa", "exit"))
exit_muc$site <- factor(exit_muc$site)
rf_exitmuc <- randomForest(site ~ ., data=select(exit_muc, site, contains("Otu")), importance=T, ntree=n_trees)

#exit vs L lumen
exit_Llum <- subset(rel_meta, location %in% c("LS", "SS"))
exit_Llum$location <- factor(exit_Llum$location)
rf_exitLlum <- randomForest(location ~ ., data=select(exit_Llum, location, contains("Otu")), importance=T, ntree=n_trees)

#exit vs R lumen
exit_Rlum <- subset(rel_meta, location %in% c("RS", "SS"))
exit_Rlum$location <- factor(exit_Rlum$location)
rf_exitRlum <- randomForest(location ~ ., data=select(exit_Rlum, location, contains("Otu")), importance=T, ntree=n_trees)



#now do aucrf model 
# remove cage/inocula OTUs from prediction dataframes
#1 is stool, 0 is mucosa
classification_labels <- levels(left_bs$location) # save level labels as a vector for reference
levels(left_bs$location) <- c(1:length(levels(left_bs$location))-1) # convert levels to numeric based on number of levels

classification_labelsLR <- levels(LR_bowel$location) # save level labels as a vector for reference
levels(LR_bowel$location) <- c(1:length(levels(LR_bowel$location))-1) # convert levels to numeric based on number of levels

classification_labelsR <- levels(right_bs$location) # save level labels as a vector for reference
levels(right_bs$location) <- c(1:length(levels(right_bs$location))-1) # convert levels to numeric based on number of levels

classification_labelslumen <- levels(LR_lumen$location) # save level labels as a vector for reference
levels(LR_lumen$location) <- c(1:length(levels(LR_lumen$location))-1) # convert levels to numeric based on number of levels

classification_labelsall <- levels(all_lum$site)
levels(all_lum$site) <- c(1:length(levels(all_lum$site)) - 1)

classification_labelsexitlum <- levels(exit_lum$site)
levels(exit_lum$site) <- c(1:length(levels(exit_lum$site))-1)

classification_labelsexitmuc <- levels(exit_muc$site)
levels(exit_muc$site) <- c(1:length(levels(exit_muc$site))-1)

classification_labelsexitLlum <- levels(exit_Llum$location) # save level labels as a vector for reference
levels(exit_Llum$location) <- c(1:length(levels(exit_Llum$location))-1) 

classification_labelsexitRlum <- levels(exit_Rlum$location) # save level labels as a vector for reference
levels(exit_Rlum$location) <- c(1:length(levels(exit_Rlum$location))-1) 

# create RF model
set.seed(seed)
rf_left_aucrf <- AUCRF(location ~ ., data = select(left_bs, location, contains("Otu")),
                        ntree = n_trees, pdel = 0.05, ranking = 'MDA')
set.seed(seed)
rf_LRbowel_aucrf <- AUCRF(location ~ ., data = select(LR_bowel, location, contains("Otu")),
                       ntree = n_trees, pdel = 0.05, ranking = 'MDA')
set.seed(seed)
rf_right_aucrf <- AUCRF(location ~ ., data = select(right_bs, location, contains("Otu")),
                          ntree = n_trees, pdel = 0.05, ranking = 'MDA')
set.seed(seed)
rf_LRlumen_aucrf <- AUCRF(location ~ ., data = select(LR_lumen, location, contains("Otu")),
                        ntree = n_trees, pdel = 0.05, ranking = 'MDA')
set.seed(seed)
rf_all_aucrf <- AUCRF(site ~ ., data = select(all_lum, site, contains("Otu")),
                          ntree = n_trees, pdel = 0.05, ranking = 'MDA')
set.seed(seed)
rf_exitlum_aucrf <- AUCRF(site ~ ., data = select(exit_lum, site, contains("Otu")),
                      ntree = n_trees, pdel = 0.05, ranking = 'MDA')
set.seed(seed)
rf_exitLlum_aucrf <- AUCRF(location ~ ., data = select(exit_Llum, location, contains("Otu")),
                          ntree = n_trees, pdel = 0.05, ranking = 'MDA')
set.seed(seed)
rf_exitRlum_aucrf <- AUCRF(location ~ ., data = select(exit_Rlum, location, contains("Otu")),
                           ntree = n_trees, pdel = 0.05, ranking = 'MDA')


#this is just left stool vs mucosa  
otu_left_probs <- predict(rf_left_aucrf$RFopt, type = 'prob')

otu_LRbowel_probs <- predict(rf_LRbowel_aucrf$RFopt, type ='prob')

otu_right_probs <- predict(rf_right_aucrf$RFopt, type = 'prob')

otu_LRlumen_probs <- predict(rf_LRlumen_aucrf$RFopt, type = 'prob')

otu_all_probs <- predict(rf_all_aucrf$RFopt, type = 'prob')

otu_exitlum_probs <- predict(rf_exitlum_aucrf$RFopt, type = 'prob')

otu_exitmuc_probs <- predict(rf_exitmuc_aucrf$RFopt, type = 'prob')

otu_exitLlum_probs <- predict(rf_exitLlum_aucrf$RFopt, type ='prob')

otu_exitRlum_probs <- predict(rf_exitRlum_aucrf$RFopt, type ='prob')

all_left_probs <- data.frame(obs = left_bs$location,
                             pred = otu_left_probs[,2])

all_LRbowel_probs <- data.frame(obs = LR_bowel$location,
                             pred = otu_LRbowel_probs[,2])

all_right_probs <- data.frame(obs = right_bs$location,
                                pred = otu_right_probs[,2])

all_LRlumen_probs <- data.frame(obs = LR_lumen$location,
                              pred = otu_LRlumen_probs[,2])

all_all_probs <- data.frame(obs = all_lum$site,
                                pred = otu_all_probs[,2])

all_exitlum_probs <- data.frame(obs = exit_lum$site,
                            pred = otu_exitlum_probs[,2])

all_exitmuc_probs <- data.frame(obs = exit_muc$site,
                                pred = otu_exitmuc_probs[,2])

all_exitLlum_probs <- data.frame(obs = exit_Llum$location,
                                pred = otu_exitLlum_probs[,2])

all_exitRlum_probs <- data.frame(obs = exit_Rlum$location,
                                 pred = otu_exitRlum_probs[,2])

#compare real with predicted with ROC
otu_left_roc <- roc(left_bs$location ~ otu_left_probs[ , 2])
left_otu_feat <- rf_left_aucrf$Xopt

otu_LRbowel_roc <- roc(LR_bowel$location ~ otu_LRbowel_probs[ , 2])
LRbowel_otu_feat <- rf_LRbowel_aucrf$Xopt

otu_right_roc <- roc(right_bs$location ~ otu_right_probs[ , 2])
right_otu_feat <- rf_right_aucrf$Xopt


otu_LRlumen_roc <- roc(LR_lumen$location ~ otu_LRlumen_probs[ , 2])
LRlumen_otu_feat <- rf_LRlumen_aucrf$Xopt

otu_all_roc <- roc(all_lum$site ~ otu_all_probs[,2])
all_otu_feat <- rf_all_aucrf$Xopt

otu_exitlum_roc <- roc(exit_lum$site ~ otu_exitlum_probs[,2])
exitlum_feat <- rf_exitlum_aucrf$Xopt

otu_exitmuc_roc <- roc(exit_muc$site ~ otu_exitmuc_probs[,2])
exitmuc_feat <- rf_exitmuc_aucrf$Xopt

otu_exitLlum_roc <- roc(exit_Llum$location ~ otu_exitLlum_probs[,2])
exitLlum_feat <- rf_exitLlum_aucrf$Xopt

otu_exitRlum_roc <- roc(exit_Rlum$location ~ otu_exitRlum_probs[,2])
exitRlum_feat <- rf_exitRlum_aucrf$Xopt

#really should make that a function 

#generate entire figure just of exit comparisons ? 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(otu_exitlum_roc, col='red', lwd=2, add=T, lty=1) #all lumen vs exit 
plot(otu_exitmuc_roc, col='blue', lwd=2, add=T, lty=1) #all mucosa vs exit
plot(otu_exitLlum_roc, col='green4', lwd=2, add=T, lty=1) #left lumen vs exit 
plot(otu_exitRlum_roc, col='purple', lwd=2, add=T, lty=1) #right lumen vs left lumen 
#plot(otu_all_roc, col = 'pink', lwd=2, add =T, lty=1)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottom', legend=c(sprintf('All lumen vs exit, AUC = 0.882'),
                          sprintf('All mucosa vs exit, AUC = 0.991'),
                          sprintf('L lumen vs exit, AUC = 0.802'),
                          sprintf('R lumen vs exit, AUC = 0.934')
                          # sprintf('all lumen vs all mucosa, AUC = 0.922')#,
                          #                               sprintf('OOB vs Leave-1-out: p=%.2g', roc.test(otu_euth_roc,LOO_roc)$p.value),
                          #                               sprintf('OOB vs 10-fold CV: p=%.2g', roc.test(otu_euth_roc,cv10f_roc)$p.value)
), lty=1, lwd=2, col=c('red','blue', 'green4', 'purple'), bty='n')




#draw blank plot, plot each line over it individually 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(otu_left_roc, col='red', lwd=2, add=T, lty=1) #left stool vs mucosa
plot(otu_LRbowel_roc, col='blue', lwd=2, add=T, lty=1) #left mucosa vs right mucosa
plot(otu_right_roc, col='green4', lwd=2, add=T, lty=1) #right stool vs mucosa
plot(otu_LRlumen_roc, col='purple', lwd=2, add=T, lty=1) #right lumen vs left lumen 
plot(otu_all_roc, col = 'pink', lwd=2, add =T, lty=1)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottom', legend=c(sprintf('L lumen vs L mucosa, AUC = 0.984', otu_left_roc$auc),
                               sprintf('L mucosa vs R mucosa, AUC = 0.926',otu_LRbowel_roc$auc),
                               sprintf('R lumen vs R mucosa, AUC = 0.860',otu_right_roc$auc),
                               sprintf('R lumen vs L lumen, AUC = 0.773', otu_LRlumen_roc$auc),
                                sprintf('all lumen vs all mucosa, AUC = 0.922')#,
                               #                               sprintf('OOB vs Leave-1-out: p=%.2g', roc.test(otu_euth_roc,LOO_roc)$p.value),
                               #                               sprintf('OOB vs 10-fold CV: p=%.2g', roc.test(otu_euth_roc,cv10f_roc)$p.value)
),lty=1, lwd=2, col=c('red','blue','green4', 'purple', 'pink'), bty='n')


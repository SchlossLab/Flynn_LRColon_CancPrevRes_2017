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
shared_meta <- merge(meta_file, shared_file, by.x='group', by.y='row.names')

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

#gender test, all M vs all F 
#but i really want to look at m mucosa vs f mucosa

gender_mucosa <- subset(rel_meta, site == 'mucosa')
gender_mucosa$gender <- factor(gender_mucosa$gender)
rf_gender_muc <- randomForest(gender~., data=select(gender_mucosa, gender, contains("Otu")), importance = T, ntree=n_trees)

gender_lumen <- subset(rel_meta, site == 'stool')
gender_lumen$gender <- factor(gender_lumen$gender)
rf_gender_lum <- randomForest(gender~., data=select(gender_lumen, gender, contains("Otu")), importance = T, ntree=n_trees)

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

classification_labelsgendermuc <- levels(gender_mucosa$gender) # save level labels as a vector for reference
levels(gender_mucosa$gender) <- c(1:length(levels(gender_mucosa$gender))-1)

classification_labelsgenderlum <- levels(gender_lumen$gender) # save level labels as a vector for reference
levels(gender_lumen$gender) <- c(1:length(levels(gender_lumen$gender))-1)


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
set.seed(seed)
rf_gendermuc_aucrf <- AUCRF(gender~., data = select(gender_mucosa, gender, contains("Otu")), ntree=n_trees, pdel=0.05, ranking = 'MDA')

set.seed(seed)
rf_genderlum_aucrf <- AUCRF(gender~., data = select(gender_lumen, gender, contains("Otu")), ntree=n_trees, pdel=0.05, ranking = 'MDA')



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

otu_gendermuc_probs <- predict(rf_gendermuc_aucrf$RFopt, type = 'prob')

otu_genderlum_probs <- predict(rf_genderlum_aucrf$RFopt, type = 'prob')


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

gender_muc_probs <- data.frame(obs=gender_mucosa$gender, pred=otu_gendermuc_probs[,2])

gender_lum_probs <- data.frame(obs=gender_lumen$gender, pred=otu_genderlum_probs[,2])



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

otu_gendermuc_roc <- roc(gender_mucosa$gender ~ otu_gendermuc_probs[,2])
gendermuc_feat <- rf_gendermuc_aucrf$Xopt

otu_genderlum_roc <- roc(gender_lumen$gender ~ otu_genderlum_probs[,2])
genderlum_feat <- rf_genderlum_aucrf$Xopt

#really should make that a function 

aucrf_data_allum <- all_lum[,c('site',all_otu_feat)]
aucrf_data_LRbowel <- LR_bowel[, c('location', LRbowel_otu_feat)]
aucrf_data_LRlumen <- LR_lumen[, c('location', LRlumen_otu_feat)]


#10 fold cross validation for all lumen vs mucosa 
iters <- 100
cv10f_aucs <- c()
cv10f_all_resp <- c()
cv10f_all_pred <- c()
for(j in 1:iters){
  set.seed(j)
  sampling <- sample(1:nrow(aucrf_data_allum),nrow(aucrf_data_allum),replace=F)
  cv10f_probs <- rep(NA,78)
  for(i in seq(1,77,7)){
    train <- aucrf_data_allum[sampling[-(i:(i+6))],]
    test <- aucrf_data_allum[sampling[i:(i+6)],]
    set.seed(seed)
    temp_model <- AUCRF(site~., data=train, pdel=0.99, ntree=500)
    cv10f_probs[sampling[i:(i+6)]] <- predict(temp_model$RFopt, test, type='prob')[,2]
  }
  cv10f_roc <- roc(aucrf_data_allum$site~cv10f_probs)
  cv10f_all_pred <- c(cv10f_all_pred, cv10f_probs)
  cv10f_all_resp <- c(cv10f_all_resp, aucrf_data_allum$site)
  cv10f_aucs[j] <- cv10f_roc$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc <- roc(cv10f_all_resp~cv10f_all_pred)

#10 fold cross validation for L vs R mucosa

iters <- 100
cv10f_aucs_muc <- c()
cv10f_all_resp_muc <- c()
cv10f_all_pred_muc <- c()
for(j in 1:iters){
  set.seed(j)
  sampling_muc <- sample(1:nrow(aucrf_data_LRbowel),nrow(aucrf_data_LRbowel),replace=F)
  cv10f_probs_muc <- rep(NA,39)
  for(i in seq(1,36,4)){
    train_muc <- aucrf_data_LRbowel[sampling_muc[-(i:(i+3))],]
    test_muc <- aucrf_data_LRbowel[sampling_muc[i:(i+3)],]
    set.seed(seed)
    temp_model_muc <- AUCRF(location~., data=train_muc, pdel=0.99, ntree=500)
    cv10f_probs_muc[sampling_muc[i:(i+3)]] <- predict(temp_model_muc$RFopt, test_muc, type='prob')[,2]
  }
  cv10f_roc_muc <- roc(aucrf_data_LRbowel$location~cv10f_probs_muc)
  cv10f_all_pred_muc <- c(cv10f_all_pred_muc, cv10f_probs_muc)
  cv10f_all_resp_muc <- c(cv10f_all_resp_muc, aucrf_data_LRbowel$location)
  cv10f_aucs_muc[j] <- cv10f_roc_muc$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc_muc <- roc(cv10f_all_resp_muc~cv10f_all_pred_muc)


#10 fold cross validation for L vs R lumen

iters <- 100
cv10f_aucs_lum <- c()
cv10f_all_resp_lum <- c()
cv10f_all_pred_lum <- c()
for(j in 1:iters){
  set.seed(j)
  sampling_lum <- sample(1:nrow(aucrf_data_LRlumen),nrow(aucrf_data_LRlumen),replace=F)
  cv10f_probs_lum <- rep(NA,39)
  for(i in seq(1,36,4)){
    train_lum <- aucrf_data_LRlumen[sampling_lum[-(i:(i+3))],]
    test_lum <- aucrf_data_LRlumen[sampling_lum[i:(i+3)],]
    set.seed(seed)
    temp_model_lum <- AUCRF(location~., data=train_lum, pdel=0.99, ntree=500)
    cv10f_probs_lum[sampling_lum[i:(i+3)]] <- predict(temp_model_lum$RFopt, test_lum, type='prob')[,2]
  }
  cv10f_roc_lum <- roc(aucrf_data_LRlumen$location~cv10f_probs_lum)
  cv10f_all_pred_lum <- c(cv10f_all_pred_lum, cv10f_probs_lum)
  cv10f_all_resp_lum <- c(cv10f_all_resp_lum, aucrf_data_LRlumen$location)
  cv10f_aucs_lum[j] <- cv10f_roc_lum$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc_lum <- roc(cv10f_all_resp_lum~cv10f_all_pred_lum)




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




#Lumen vs mucosa plot 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(otu_left_roc, col='red', lwd=2, add=T, lty=1) #left stool vs mucosa
#plot(otu_LRbowel_roc, col='blue', lwd=2, add=T, lty=1) #left mucosa vs right mucosa
plot(otu_right_roc, col='blue', lwd=2, add=T, lty=1) #right stool vs mucosa
#plot(otu_LRlumen_roc, col='purple', lwd=2, add=T, lty=1) #right lumen vs left lumen 
plot(otu_all_roc, col = 'purple', lwd=2, add =T, lty=1)
plot(cv10f_roc, col = 'purple', lwd=2, add=T, lty=2)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottom', legend=c(sprintf('L lumen vs L mucosa, AUC = 0.984', otu_left_roc$auc),
                               #sprintf('L mucosa vs R mucosa, AUC = 0.926',otu_LRbowel_roc$auc),
                               sprintf('R lumen vs R mucosa, AUC = 0.860',otu_right_roc$auc),
                               #sprintf('R lumen vs L lumen, AUC = 0.773', otu_LRlumen_roc$auc),
                                sprintf('all lumen vs all mucosa, AUC = 0.922'),
                                sprintf('all vs all 10-fold CV, AUC = 0.925')#,
                               #                               sprintf('OOB vs Leave-1-out: p=%.2g', roc.test(otu_euth_roc,LOO_roc)$p.value),
                               #                               sprintf('OOB vs 10-fold CV: p=%.2g', roc.test(otu_euth_roc,cv10f_roc)$p.value)
),lty=c(1, 1, 1, 2), lwd=2, col=c('red','blue', 'purple', 'purple'), bty='n')


#left vs right mucosa and lumen plot 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
#plot(otu_left_roc, col='red', lwd=2, add=T, lty=1) #left stool vs mucosa
plot(otu_LRbowel_roc, col='green4', lwd=2, add=T, lty=1) #left mucosa vs right mucosa
#plot(otu_right_roc, col='blue', lwd=2, add=T, lty=1) #right stool vs mucosa
plot(otu_LRlumen_roc, col='orange', lwd=2, add=T, lty=1) #right lumen vs left lumen 
plot(cv10f_roc_muc,col = 'green4', lwd=2, add=T, lty=2) #r vs l mucosa cross validation
plot(cv10f_roc_lum, col = 'orange', lwd=2, add=T, lty=2) #r vs l lumen cross validation
plot(otu_gendermuc_roc, col = 'pink', lwd=2, add=T, lty=1) #gender mucosa test
#plot(otu_all_roc, col = 'purple', lwd=2, add =T, lty=1)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottom', legend=c(#sprintf('L lumen vs L mucosa, AUC = 0.984', otu_left_roc$auc),
                          sprintf('L mucosa vs R mucosa, AUC = 0.926',otu_LRbowel_roc$auc),
                          #sprintf('R lumen vs R mucosa, AUC = 0.860',otu_right_roc$auc),
                          sprintf('L lumen vs R lumen, AUC = 0.773', otu_LRlumen_roc$auc),
                          sprintf('L mucosa vs R mucosa 10-fold CV, AUC = 0.912'),
                          sprintf('L lumen vs R lumen 10-fold CV, AUC = 0.7551')
                          #sprintf('all lumen vs all mucosa, AUC = 0.922')#,
                          #                               sprintf('OOB vs Leave-1-out: p=%.2g', roc.test(otu_euth_roc,LOO_roc)$p.value),
                          #                               sprintf('OOB vs 10-fold CV: p=%.2g', roc.test(otu_euth_roc,cv10f_roc)$p.value)
),lty=c(1, 1, 2, 2), lwd=2, col=c('green4', 'orange', 'green4', 'orange'), bty='n')


#predict gender of mucosa?! new plot 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(otu_gendermuc_roc, col = 'pink', lwd=2, add=T, lty=1) #gender mucosa test
plot(otu_genderlum_roc, col = 'blue', lwd=2, add=T, lty=1)
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
legend('bottom', legend=c(#sprintf('L lumen vs L mucosa, AUC = 0.984', otu_left_roc$auc),
  sprintf('L mucosa vs R mucosa, AUC = 0.926',otu_LRbowel_roc$auc)
  #                               sprintf('OOB vs Leave-1-out: p=%.2g', roc.test(otu_euth_roc,LOO_roc)$p.value),
),lty=c(1, 1, 2, 2), lwd=2, col=c('pink', 'blue'), bty='n')






#importance plots for OTUs

#do this for all combinations??

tax_function <- 'code/tax_level.R'
source(tax_function)

#add nicks looping back in. store all plots in a list. then call from the list. but will need a list of rf_names first 
n_features <- 10
#just going to start with rf_left
importance_sorted_rfleft <- sort(importance(rf_left)[,1], decreasing = T)
top_important_OTU_rfleft <- data.frame(head(importance_sorted_rfleft, n_features))
colnames(top_important_OTU_rfleft) <- 'Importance'
top_important_OTU_rfleft$OTU <- rownames(top_important_OTU_rfleft)
otu_taxa_rfleft <- get_tax(1, top_important_OTU_rfleft$OTU, tax_file)

#importance_plot_day_rfleft <- 
  ggplot(data = top_important_OTU_rfleft, aes(x = factor(OTU), y = Importance)) + 
  geom_point() + scale_x_discrete(limits = rev(top_important_OTU_rfleft$OTU),
                                  labels = rev(paste(otu_taxa_rfleft[,1],' (',
                                                     rownames(otu_taxa_rfleft),')',
                                                     sep=''))) +
  labs(x= '', y = '% Increase in MSE') + theme_bw() + coord_flip() + ggtitle('LB vs LS')

  
#RB vs RS
importance_sorted_rfright <- sort(importance(rf_right)[,1], decreasing = T)
top_important_OTU_rfright <- data.frame(head(importance_sorted_rfright, n_features))
colnames(top_important_OTU_rfright) <- 'Importance'
top_important_OTU_rfright$OTU <- rownames(top_important_OTU_rfright)
otu_taxa_rfright <- get_tax(1, top_important_OTU_rfright$OTU, tax_file)
  
  #importance_plot_day_rfleft <- 
ggplot(data = top_important_OTU_rfright, aes(x = factor(OTU), y = Importance)) + 
    geom_point() + scale_x_discrete(limits = rev(top_important_OTU_rfright$OTU),
                                    labels = rev(paste(otu_taxa_rfright[,1],' (',
                                                       rownames(otu_taxa_rfright),')',
                                                       sep=''))) +
    labs(x= '', y = '% Increase in MSE') + theme_bw() + coord_flip() + ggtitle('RB vs RS')
  
#all lumen vs mucosa 

importance_sorted_rfall <- sort(importance(rf_all)[,1], decreasing = T)
top_important_OTU_rfall <- data.frame(head(importance_sorted_rfall, n_features))
colnames(top_important_OTU_rfall) <- 'Importance'
top_important_OTU_rfall$OTU <- rownames(top_important_OTU_rfall)
otu_taxa_rfall <- get_tax(1, top_important_OTU_rfall$OTU, tax_file)

#importance_plot_day_rfleft <- 
ggplot(data = top_important_OTU_rfall, aes(x = factor(OTU), y = Importance)) + 
  geom_point() + scale_x_discrete(limits = rev(top_important_OTU_rfall$OTU),
                                  labels = rev(paste(otu_taxa_rfall[,1],' (',
                                                     rownames(otu_taxa_rfall),')',
                                                     sep=''))) +
  labs(x= '', y = '% Increase in MSE') + theme_bw() + coord_flip() + ggtitle('All mucosa vs lumen')


#L mucosa vs R mucosa 
importance_sorted_rfbowel <- sort(importance(rf_bowel)[,1], decreasing = T)
top_important_OTU_rfbowel <- data.frame(head(importance_sorted_rfbowel, n_features))
colnames(top_important_OTU_rfbowel) <- 'Importance'
top_important_OTU_rfbowel$OTU <- rownames(top_important_OTU_rfbowel)
otu_taxa_rfbowel <- get_tax(1, top_important_OTU_rfbowel$OTU, tax_file)

#importance_plot_day_rfleft <- 
ggplot(data = top_important_OTU_rfbowel, aes(x = factor(OTU), y = Importance)) + 
  geom_point() + scale_x_discrete(limits = rev(top_important_OTU_rfbowel$OTU),
                                  labels = rev(paste(otu_taxa_rfbowel[,1],' (',
                                                     rownames(otu_taxa_rfbowel),')',
                                                     sep=''))) +
  labs(x= '', y = '% Increase in MSE') + theme_bw() + coord_flip() + ggtitle('L mucosa vs R mucosa')


#L lumen vs R lumen 
importance_sorted_rflumen <- sort(importance(rf_lumen)[,1], decreasing = T)
top_important_OTU_rflumen <- data.frame(head(importance_sorted_rflumen, n_features))
colnames(top_important_OTU_rflumen) <- 'Importance'
top_important_OTU_rflumen$OTU <- rownames(top_important_OTU_rflumen)
otu_taxa_rflumen <- get_tax(1, top_important_OTU_rflumen$OTU, tax_file)

#importance_plot_day_rfleft <- 
ggplot(data = top_important_OTU_rflumen, aes(x = factor(OTU), y = Importance)) + 
  geom_point() + scale_x_discrete(limits = rev(top_important_OTU_rflumen$OTU),
                                  labels = rev(paste(otu_taxa_rflumen[,1],' (',
                                                     rownames(otu_taxa_rflumen),')',
                                                     sep=''))) +
  labs(x= '', y = '% Increase in MSE') + theme_bw() + coord_flip() + ggtitle('L lumen vs R lumen')



#make all of this shit a function

#also make RA plots of the choice OTUs 

#all_otu_feat holds the important OTUs for all lumen vs all mucosa 

all_otu_feat <- rev(all_otu_feat[1:5])
otu_taxa_all <- get_tax(1, all_otu_feat, tax_file)
#Abundance stripchart or most predictive otus - Niel's code 
lumen_abunds <- shared_meta[shared_meta$site=='stool', all_otu_feat]/10000 + 1e-4
mucosa_abunds <- shared_meta[shared_meta$site=='mucosa', all_otu_feat]/10000 + 1e-4

par(mar=c(4, 9, 1, 1))
plot(1, type="n", ylim=c(0,length(all_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in all_otu_feat){
  stripchart(at=index-0.35, jitter(lumen_abunds[,i], amount=1e-5), pch=21, bg="royalblue1", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(mucosa_abunds[,i], amount=1e-5), pch=21, bg="orange", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(mean(lumen_abunds[,i]),index-0.7,mean(lumen_abunds[,i]),index, lwd=3)
  segments(mean(mucosa_abunds[,i]),index+0.7,mean(mucosa_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=otu_taxa_all$tax_label, las=1, line=-0.5, tick=F, cex.axis=0.8)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("mucosa", "lumen"), pch=c(21, 21), pt.bg=c("orange","royalblue1"), cex=0.7)
mtext('C', at=1e-7, font=2, side=3)





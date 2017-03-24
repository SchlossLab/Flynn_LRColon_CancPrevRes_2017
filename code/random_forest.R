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

#build randomForest model for each location comparison using randomize_loc function 
rf_left <- randomize_loc(rel_meta, "LB", "LS") #OOB 10.26%
rf_right <- randomize_loc(rel_meta, "RB", "RS") #OOB 53%
rf_bowel <- randomize_loc(rel_meta, "LB", "RB") #OOB 25.64%
rf_lumen <- randomize_loc(rel_meta, "LS", "RS") #OOB 69.23%
rf_exitRlum <- randomize_loc(rel_meta, "RS", "SS")
rf_exitLlum <- randomize_loc(rel_meta, "LS", "SS")

#and for each site
rf_all <- randomize_site(rel_meta, "mucosa", "stool")
rf_exitlum <- randomize_site(rel_meta, "stool", "exit")
rf_exitmuc <- randomize_site(rel_meta, "mucosa", "exit")



#for left, these OTUS are higher in LM than in LL, also (see below) higher in LB than RB 
left_importance <- sort(importance(rf_left)[,1], decreasing = T)

#sort importance to show which OTUs are higher in Left than right bowel
bowel_importance <- sort(importance(rf_bowel)[,1], decreasing = T)

#to get what is higher in LS vs RS 
lumen_importance <- sort(importance(rf_lumen)[,1], decreasing=T)


# create RF model with AUCRF outputs top OTUs
aucrf_data_left_bs <- auc_loc(rel_meta, "LB", "LS")
aucrf_data_LRbowel <- auc_loc(rel_meta, "LB", "RB")
aucrf_data_right_bs <- auc_loc(rel_meta, "RB", "RS")
aucrf_data_LRlumen <- auc_loc(rel_meta, "LS", "RS")
aucrf_data_allum <- auc_site(rel_meta, "mucosa", "stool")

#need to fix specific function for these to output just the rf object 
rf_exitlum_aucrf <- auc_site(rel_meta, "stool", "exit") # not working 
rf_exitLlum_aucrf <- auc_loc(rel_meta, "LB", "SS")
rf_exitRlum_aucrf <- auc_loc(rel_meta, "RB", "SS")



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

#10fold CV for L lumen vs L mucosa
iters <- 100
cv10f_aucs <- c()
cv10f_all_resp_left_bs <- c()
cv10f_all_pred_left_bs <- c()
for(j in 1:iters){
  set.seed(j)
  sampling <- sample(1:nrow(aucrf_data_left_bs),nrow(aucrf_data_left_bs),replace=F)
  cv10f_probs <- rep(NA,39)
  for(i in seq(1,36,4)){
    train_left_bs <- aucrf_data_left_bs[sampling[-(i:(i+3))],]
    test_left_bs <- aucrf_data_left_bs[sampling[i:(i+3)],]
    set.seed(seed)
    temp_model_left_bs <- AUCRF(location~., data=train_left_bs, pdel=0.99, ntree=500)
    cv10f_probs[sampling[i:(i+3)]] <- predict(temp_model_left_bs$RFopt, test_left_bs, type='prob')[,2]
  }
  cv10f_roc_left_bs <- roc(aucrf_data_left_bs$location~cv10f_probs)
  cv10f_all_pred_left_bs <- c(cv10f_all_pred_left_bs, cv10f_probs)
  cv10f_all_resp_left_bs <- c(cv10f_all_resp_left_bs, aucrf_data_left_bs$location)
  cv10f_aucs[j] <- cv10f_roc_left_bs$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc_left_bs <- roc(cv10f_all_resp_left_bs~cv10f_all_pred_left_bs)

#10fold CV for R lumen vs R mucosa
iters <- 100
cv10f_aucs <- c()
cv10f_all_resp_right_bs <- c()
cv10f_all_pred_right_bs <- c()
for(j in 1:iters){
  set.seed(j)
  sampling <- sample(1:nrow(aucrf_data_right_bs),nrow(aucrf_data_right_bs),replace=F)
  cv10f_probs <- rep(NA,39)
  for(i in seq(1,36,4)){
    train_right_bs <- aucrf_data_right_bs[sampling[-(i:(i+3))],]
    test_right_bs <- aucrf_data_right_bs[sampling[i:(i+3)],]
    set.seed(seed)
    temp_model_right_bs <- AUCRF(location~., data=train_right_bs, pdel=0.99, ntree=500)
    cv10f_probs[sampling[i:(i+3)]] <- predict(temp_model_right_bs$RFopt, test_right_bs, type='prob')[,2]
  }
  cv10f_roc_right_bs <- roc(aucrf_data_right_bs$location~cv10f_probs)
  cv10f_all_pred_right_bs <- c(cv10f_all_pred_right_bs, cv10f_probs)
  cv10f_all_resp_right_bs <- c(cv10f_all_resp_right_bs, aucrf_data_right_bs$location)
  cv10f_aucs[j] <- cv10f_roc_right_bs$auc #stores aucs for all iterations, can use to calc IQR
}
cv10f_roc_right_bs <- roc(cv10f_all_resp_right_bs~cv10f_all_pred_right_bs)

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



#generate entire figure just of exit comparisons 
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
plot(cv10f_roc_right_bs, col='blue', lwd=2, add=T, lty=1)
plot(cv10f_roc, col = 'purple', lwd=2, add=T, lty=1)
plot(cv10f_roc_left_bs, col = 'red', lwd=2, add=T, lty=1)
mtext(side=2, text="True Positive (Sensitivity)", line=2.5)
mtext(side=1, text="True Negative (Specificity)", line=2.5)
legend('bottom', legend=c(sprintf('Lumen vs Mucosa, 10-fold CV, AUC = 0.925'),
                          sprintf('L Lumen vs L Mucosa, 10-fold CV, AUC =0.980'),
                          sprintf('R Lumen vs R Mucosa, 10-fold CV, AUC = 0.797')
                               #sprintf('OOB vs Leave-1-out: p=%.2g', roc.test(otu_euth_roc,LOO_roc)$p.value),
                               #sprintf('OOB vs 10-fold CV: p=%.2g', roc.test(otu_euth_roc,cv10f_roc)$p.value)
),lty=c(1, 1, 1), lwd=2, col=c('purple','red', 'blue'), bty='n')


#left vs right mucosa and lumen plot 
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
plot(cv10f_roc_muc,col = 'green4', lwd=2, add=T, lty=1) #r vs l mucosa cross validation
plot(cv10f_roc_lum, col = 'orange', lwd=2, add=T, lty=1) #r vs l lumen cross validation
mtext(side=2, text="True Positive (Sensitivity)", line=2.5)
mtext(side=1, text="True Negative (Specificity)", line=2.5)
legend('bottom', legend=c(sprintf('L mucosa vs R mucosa 10-fold CV, AUC = 0.912'),
                          sprintf('L lumen vs R lumen 10-fold CV, AUC = 0.7551')
                          # sprintf('OOB vs Leave-1-out: p=%.2g', roc.test(otu_euth_roc,LOO_roc)$p.value),
                          # sprintf('OOB vs 10-fold CV: p=%.2g', roc.test(otu_euth_roc,cv10f_roc)$p.value)
),lty=c(1, 1), lwd=2, col=c('green4', 'orange'), bty='n')


#importance plots for OTUs
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

ggplot(data = top_important_OTU_rflumen, aes(x = factor(OTU), y = Importance)) + 
  geom_point() + scale_x_discrete(limits = rev(top_important_OTU_rflumen$OTU),
                                  labels = rev(paste(otu_taxa_rflumen[,1],' (',
                                                     rownames(otu_taxa_rflumen),')',
                                                     sep=''))) +
  labs(x= '', y = '% Increase in MSE') + theme_bw() + coord_flip() + ggtitle('L lumen vs R lumen')


#all_otu_feat holds the important OTUs for all lumen vs all mucosa 

all_otu_feat <- rev(all_otu_feat[1:5])
otu_taxa_all <- get_tax(1, all_otu_feat, tax_file)
#Abundance stripchart or most predictive otus 
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


#just LB vs LS 

left_otu_feat <- rev(left_otu_feat[1:5])
otu_taxa_left <- get_tax(1, left_otu_feat, tax_file)
#Abundance stripchart or most predictive otus
ls_abunds <- shared_meta[shared_meta$location=='LS', left_otu_feat]/10000 + 1e-4
lb_abunds <- shared_meta[shared_meta$location=='LB', left_otu_feat]/10000 + 1e-4

par(mar=c(4, 9, 1, 1))
plot(1, type="n", ylim=c(0,length(left_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in left_otu_feat){
  stripchart(at=index-0.35, jitter(ls_abunds[,i], amount=1e-5), pch=21, bg="lightblue", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(lb_abunds[,i], amount=1e-5), pch=21, bg="yellow", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(mean(ls_abunds[,i]),index-0.7,mean(ls_abunds[,i]),index, lwd=3)
  segments(mean(lb_abunds[,i]),index+0.7,mean(lb_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=otu_taxa_left$tax_label, las=1, line=-0.5, tick=F, cex.axis=0.8)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("left mucosa", "left lumen"), pch=c(21, 21), pt.bg=c("yellow","lightblue"), cex=1)

#RB vs RS

right_otu_feat <- rev(right_otu_feat[1:5])
otu_taxa_right <- get_tax(1, right_otu_feat, tax_file)
#Abundance stripchart or most predictive otus
rs_abunds <- shared_meta[shared_meta$location=='RS', right_otu_feat]/10000 + 1e-4
rb_abunds <- shared_meta[shared_meta$location=='RB', right_otu_feat]/10000 + 1e-4

par(mar=c(4, 10, 1, 1))
plot(1, type="n", ylim=c(0,length(right_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in right_otu_feat){
  stripchart(at=index-0.35, jitter(rs_abunds[,i], amount=1e-5), pch=21, bg="purple", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(rb_abunds[,i], amount=1e-5), pch=21, bg="orange", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(mean(rs_abunds[,i]),index-0.7,mean(rs_abunds[,i]),index, lwd=3)
  segments(mean(rb_abunds[,i]),index+0.7,mean(rb_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=otu_taxa_right$tax_label, las=1, line=-0.5, tick=F, cex.axis=0.8)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("right mucosa", "right lumen"), pch=c(21, 21), pt.bg=c("orange","purple"), cex=1)


#Lb vs Rb
LRbowel_otu_feat <- rev(LRbowel_otu_feat[1:5])
otu_taxa_LRbowel <- get_tax(1, LRbowel_otu_feat, tax_file)
#Abundance stripchart or most predictive otus 
lb_abunds <- shared_meta[shared_meta$location=='LB', LRbowel_otu_feat]/10000 + 1e-4
rblb_abunds <- shared_meta[shared_meta$location=='RB', LRbowel_otu_feat]/10000 + 1e-4

par(mar=c(4, 10, 1, 1))
plot(1, type="n", ylim=c(0,length(LRbowel_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in LRbowel_otu_feat){
  stripchart(at=index-0.35, jitter(lb_abunds[,i], amount=1e-5), pch=21, bg="darkgreen", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(rblb_abunds[,i], amount=1e-5), pch=21, bg="pink", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(mean(lb_abunds[,i]),index-0.7,mean(lb_abunds[,i]),index, lwd=3)
  segments(mean(rblb_abunds[,i]),index+0.7,mean(rblb_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=otu_taxa_LRbowel$tax_label, las=1, line=-0.5, tick=F, cex.axis=0.8)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("left mucosa", "right mucosa"), pch=c(21, 21), pt.bg=c("darkgreen","pink"), cex=1)


#LS vs RS
LRlumen_otu_feat <- rev(LRlumen_otu_feat[1:5])
otu_taxa_LRlumen <- get_tax(1, LRlumen_otu_feat, tax_file)
#Abundance stripchart or most predictive otus 
lsrs_abunds <- shared_meta[shared_meta$location=='LS', LRlumen_otu_feat]/10000 + 1e-4
rsls_abunds <- shared_meta[shared_meta$location=='RS', LRlumen_otu_feat]/10000 + 1e-4

par(mar=c(4, 10, 1, 1))
plot(1, type="n", ylim=c(0,length(LRlumen_otu_feat)*2), xlim=c(1e-4,3), log="x", ylab="", xlab="Relative Abundance (%)", xaxt="n", yaxt="n")
index <- 1
for(i in LRlumen_otu_feat){
  stripchart(at=index-0.35, jitter(lsrs_abunds[,i], amount=1e-5), pch=21, bg="brown", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(rsls_abunds[,i], amount=1e-5), pch=21, bg="magenta", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(mean(lsrs_abunds[,i]),index-0.7,mean(lsrs_abunds[,i]),index, lwd=3)
  segments(mean(rsls_abunds[,i]),index+0.7,mean(rsls_abunds[,i]),index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=otu_taxa_LRlumen$tax_label, las=1, line=-0.5, tick=F, cex.axis=0.8)
axis(1, at=c(1e-4, 1e-3, 1e-2, 1e-1, 1), label=c("0", "0.1", "1", "10", "100"))
legend('topright', legend=c("left lumen", "right lumen"), pch=c(21, 21), pt.bg=c("brown","magenta"), cex=1)


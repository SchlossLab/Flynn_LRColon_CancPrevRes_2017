#Random Forest functions
#Kaitlin Flynn, Schloss lab, 2016-present

#randomForest to compare two locations
randomize_loc <- function(yourdata, samp1, samp2){
  subsetted <- subset(yourdata, location %in% c(samp1, samp2))
  subsetted$location <- factor(subsetted$location)
  rf_out <- randomForest(location ~ ., data = select(subsetted, location, contains("Otu")), importance = T, ntree=n_trees)
}

#randomForest to compare 2 sites 
randomize_site <- function(yourdata, samp1, samp2){
  subsetted <- subset(yourdata, site %in% c(samp1, samp2))
  subsetted$site <- factor(subsetted$site)
  rf_out <- randomForest(site ~ ., data = select(subsetted, site, contains("Otu")), importance = T, ntree=n_trees)
}

#generate AUCRF model for 2 locations and return top OTUs
auc_loc <- function(yourdata, samp1, samp2){
  subsetted <- subset(yourdata, location %in% c(samp1, samp2))
  subsetted$location <- factor(subsetted$location)
  #change levels of variable of interest to 0/1
  levels(subsetted$location) <- c(1:length(levels(subsetted$location))-1)
  # create RF model
  set.seed(seed)
  rf_aucrf <- AUCRF(location ~ ., data = select(subsetted, location, contains("Otu")),
                    ntree = n_trees, pdel = 0.05, ranking = 'MDA')
  otu_probs <- predict(rf_aucrf$RFopt, type = 'prob') #pull out predictive OTUs
  all_probs <- data.frame(obs = subsetted$location, pred = otu_probs[,2])
  #compare real with predicted with ROC
  otu_roc <- roc(subsetted$location ~ otu_probs[,2])
  otu_feat <- rf_aucrf$Xopt
  aucrf_data <- subsetted[, c('location', otu_feat)]
  return(aucrf_data)
  #also should store this data somewhere. maybe as unique values or within a list?
}

#generate AUCRF model for 2 sites and return top OTUs
auc_site <- function(yourdata, samp1, samp2){
  subsetted <- subset(yourdata, site %in% c(samp1, samp2))
  subsetted$site <- factor(subsetted$site)
  #change levels of variable of interest to 0/1
  levels(subsetted$site) <- c(1:length(levels(subsetted$site))-1)
  # create RF model
  set.seed(seed)
  rf_aucrf <- AUCRF(site ~ ., data = select(subsetted, site, contains("Otu")),
                    ntree = n_trees, pdel = 0.05, ranking = 'MDA')
  otu_probs <- predict(rf_aucrf$RFopt, type = 'prob') #pull out predictive OTUs
  all_probs <- data.frame(obs = subsetted$site, pred = otu_probs[,2])
  #compare real with predicted with ROC
  otu_roc <- roc(subsetted$location ~ otu_probs[,2])
  otu_feat <- rf_aucrf$Xopt
  aucrf_data <- subsetted[, c('site', otu_feat)]
  return(aucrf_data)
}


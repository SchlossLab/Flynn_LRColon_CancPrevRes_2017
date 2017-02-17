##Created by Pat Schloss and Marc Sze
#Set of functions to get the taxonomy of an OTU and insert into shared file. 
get_tax_substring <- function(tax, tax_level){
  substring <- unlist(strsplit(tax, ";"))[tax_level]
  paste(substring, collapse='.')
}
get_tax_name <- function(tax_file, tax_level){
  
  
  tax_data <- read.table(file=tax_file, header=T, stringsAsFactors=F)
  taxonomy <- tax_data$Taxonomy
  taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
  taxonomy <- gsub('"', '', taxonomy)
  
  tax_substring <- sapply(taxonomy, get_tax_substring, tax_level)
  
  names(tax_substring) <- tax_data$OTU
  
  tax_substring
}
get_tax_level_shared <- function(shared_file, tax_file, tax_level){
  
  shared_otus <- read.table(file=shared_file, header=T, stringsAsFactors=F, row.names=2)[,-c(1,2)]
  is_present <- apply(shared_otus, 2, sum) > 0
  shared <- shared_otus[,is_present]
  
  taxonomy <- get_tax_name(tax_file, tax_level)
  taxonomy <- taxonomy[colnames(shared)]
  unique_taxa <- levels(as.factor(taxonomy))
  
  shared_tax_level <- NULL
  
  for(ut in unique_taxa){
    otus <- names(taxonomy[taxonomy %in% ut])
    sub_shared <- shared_otus[,colnames(shared_otus) %in% otus]
    
    if(is.null(dim(sub_shared))){
      shared_tax_level <- cbind(shared_tax_level, sub_shared)
    } else {
      tax_level_count <- apply(sub_shared, 1, sum)
      shared_tax_level <- cbind(shared_tax_level, tax_level_count)
    }
  }
  colnames(shared_tax_level) <- unique_taxa
  rownames(shared_tax_level) <- rownames(shared)
  return(shared_tax_level)
}
##Created by Amanda Elmore
##use along with get_tax_level_shared to make a spread and regular relative abundance file with metadata merged
#Usage:  get_tax_rabund(shared_file, tax_file, metadata, subsample=2400, tax_level, select_taxa, type)
#shared_file = subsampled shared file
#tax_level (2=phylum)
#select_taxa = vector of taxa to pick and others go into "other". If you want everything make select='all'.
#
get_tax_rabund <- function(shared_file, tax_file, metadata_file, subsample=2400, tax_level, select_taxa=c("Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria"), type='gathered'){
  
  metadata<- read.table(file=metadata_file, header=T, stringsAsFactors = FALSE)
  
  #Create gathered phylum level with metadata to create graphs
  shared <- as.data.frame(get_tax_level_shared(shared_file, tax_file, tax_level))
  rabund <- shared/subsample  #abundide by the number of sequences that I subsampled to to get percent of each phyla
  #remove others except for c("Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria","Fusobacteria")
  
  if (select_taxa!='all'){
    rabund <- subset(rabund, select=select_taxa)
    rabund$other <- 1-rowSums(rabund)
  }
  rabund$id <- rownames(rabund)
  #create label
  label <-  if(tax_level==6) "genus" else if(tax_level==5) "family" else if(tax_level==4) "order" else if(tax_level==3) "class" else if(tax_level==2) "phylum"
  #create gathered version of shared file
  gathered <- gather(rabund, key, "abundance", -(length(rabund[1,])))
  colnames(gathered)[colnames(gathered)=="key"] <- label
  
  #merge with metadata
  gathered <- merge(metadata, gathered, by="id")  #gathered version of shared file
  rabund <- merge(metadata, rabund, by="id")  #shared file with phylum level abundances
  if(type=='gathered') return(gathered) else return(rabund)
}
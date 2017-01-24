#get taxonomy - df with columns for each level - 
# level - 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum), 6 (kingdom)
## - convert taxonomy text list to df, remove percentages
## - subset df by desired tax level, then replace any unclassified with next level up
get_tax <- function(tax_level=1,row_list=rownames(tax_file),df=tax_file){
     if (tax_level %in% c(1:5)){
          otu_list <- as.character(row_list)
          otu_taxonomy <- df[rownames(df) %in% otu_list,]
          otu_taxonomy <- otu_taxonomy[as.numeric(factor(otu_list)), ]
          taxonomy <- strsplit(as.character(otu_taxonomy$Taxonomy),';',fixed=TRUE)
          taxonomy <- lapply(taxonomy, `[`, seq_len(6))
          taxonomy <- data.frame(do.call('rbind', taxonomy))
          taxonomy <- sapply(taxonomy,gsub,pattern="\\(\\d*\\)",replacement="")
          level <- 7-tax_level
          if(length(row_list) < 2){
            tax_out <- as.character(taxonomy[level])
            for (i in level:2){
              next_level <- i-1
              tax_out[is.na(tax_out)] <- as.character(taxonomy[next_level])
            }
          }else{
           tax_out <- as.character(taxonomy[ ,level])
            for (i in level:2){
              next_level <- i-1
              tax_out[is.na(tax_out)] <- 
               as.character(taxonomy[is.na(tax_out),next_level])
            }
          }
          tax_out <- gsub('_unclassified', '', tax_out)
          label_out <- paste0(tax_out, ' (', 
                              gsub('tu0*', 'TU ', rownames(otu_taxonomy)),')')
          return(data.frame(tax= tax_out, tax_label = label_out, row.names=rownames(otu_taxonomy),
                            stringsAsFactors = FALSE))
     } else {print(
          'Error: Level must be 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum)')
     }
}
####
# function to output a dataframe with OTU rel abund summed by tax group per sample
#
# input file with tax level class and file with samples OTU rel abund
####

sum_OTU_by_tax_level <- function(TAX_LEVEL, OTU_DF, tax_df){
     source('code/tax_level.R')
     TAX_DF <-get_tax(TAX_LEVEL,names(OTU_DF), tax_df)
     tax_levels <- as.character(unique(TAX_DF$tax))
     OUTPUT_DF <- data.frame(rep(0,length(rownames(OTU_DF))), row.names=rownames(OTU_DF))
     for (i in 1:length(tax_levels)){
          OTU_by_level <- rownames(TAX_DF)[TAX_DF$tax %in% tax_levels[i]]
          if (length(OTU_by_level)>1){
               level_column <- apply(OTU_DF[,names(OTU_DF)[names(OTU_DF) %in% OTU_by_level]],1,sum)
          } else {
               level_column <- OTU_DF[,names(OTU_DF)[names(OTU_DF) %in% OTU_by_level]]
          }     
          OUTPUT_DF[,i] <- level_column
          colnames(OUTPUT_DF)[i] <- tax_levels[i]
     }
     return(OUTPUT_DF)
}

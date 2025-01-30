#Custom functions used in RREMseq manuscript analysis

#calculate CpG coverage using data exported from Seqmonk
UniqueCpG_Seqmonk = function(methylation_call, SampleStart = NA, cov = NA)
{
  library(dplyr)
  library(data.table)
  library(tidyverse)
  
  Sample_fil = as.data.frame(methylation_call$Start)
  EndCol = ncol(methylation_call)
  SampleData = methylation_call[,SampleStart:EndCol]
  Sample_fil = cbind(Sample_fil, SampleData)
  all_sample = as.list(SampleStart:EndCol)
  SampleEnd = ncol(Sample_fil)
  SampleName = colnames(Sample_fil[2:SampleEnd])
  
  uniqueCpG = lapply(all_sample, function(comp){
    totalCpG =  length(na.omit(methylation_call[,comp]))
    methylation_call_fil = methylation_call[!is.na(methylation_call[,comp]), ]
    sample_position = as.data.frame(methylation_call_fil$Start)
    Sample_CpGdup = CountCpGduplicate(sample_position[,1]) 
    dupCpG_num = length(Sample_CpGdup$duplicatedCount)
    uniqueCpG = totalCpG - dupCpG_num
    return(uniqueCpG)
  })
  
  names(uniqueCpG) <- SampleName
  
  uniqueCpG.df <- uniqueCpG %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Sample")
  
  colnames(uniqueCpG.df) = c("Sample", "UniqueCpGs")
  
  #write.csv(uniqueCpG.df, file = "./UniqueCpGcov.csv", row.names = F)
  uniqueCpG.df$cov = cov
  output = uniqueCpG.df
  return(output)
}




#calculate average CpG calls using data exported from Seqmonk
get_AverageCalls <- function(CpGicall_count, ) {
  CpGicall_count_fil = CpGicall_count[apply(CpGicall_count, 1, 
                                            function(x) all(x != 0)), ]
  CpGicall_count_means <- CpGicall_count_fil %>% 
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    t() %>%
    as.data.frame()
  return(CpGicall_count_means)
}

#Calculate genomic elements coverage from Seqmonk data
FeatureCount <- function(df,TotalFeature,Cov) {
  counts <- colSums(!is.na(df))
  result <- data.frame(Sample = names(counts), Non_NA_Count = counts)
  result <- extract_before_underscore(result, "Sample", "Group")
  result$percentage <- result$Non_NA_Count/TotalFeature
  result$Cov <- Cov
  return(result)
}

#extract methylation coverage within TregSE

FilTregSE = function(CovFile, TregSE,ExportName)
{
  library(dplyr)
  library(data.table)
  library(tidyverse)
  
  
  cov <- CovFile
  
  probe <- paste(cov[,1], cov[,2], sep = ":")
  probe <- paste(probe, cov[,2], sep = "-")
  probe <- as.data.frame(probe)
  probe <- mutate_all(probe, funs(gsub("chr", "Chr", .)))
  cov.location <- cbind(probe,cov)
  
  TregSE.cov <- merge(cov.location,
                      TregSE,
                      by.x = "probe",
                      by.y = "Probe", all = F)
  TregSE.cov.export <- TregSE.cov[,2:7]
  
  write.table(TregSE.cov.export, file = ExportName, sep = "\t", row.names = F, col.names = F,quote = F)
  return(TregSE.cov.export)
}


#calculate average genomic elements methylation status 
ma <- function(x, n = 50){stats::filter(x, rep(1 / n, n), sides = 2)} #sliding window smoothing

get_FeatureBeta_mean <- function(FeatureMethyl, position){
  FeatureMethyl_smooth = apply(FeatureMethyl, 2, ma)
  FeatureMethyl_smooth = na.omit(FeatureMethyl_smooth)
  FeatureMethyl_mean = apply(FeatureMethyl_smooth, 1, mean)
  FeatureMethyl_mean = na.omit(FeatureMethyl_mean)
  return(FeatureMethyl_mean)
}

#calculate the range of genomic elements methylation status 
get_FeatureBeta_range <- function(FeatureMethyl, position){
  FeatureMethyl_smooth = apply(FeatureMethyl, 2, ma)
  FeatureMethyl_smooth = na.omit(FeatureMethyl_smooth)
  FeatureMethyl_max <- apply(FeatureMethyl_smooth, 1, max)
  FeatureMethyl_min <- apply(FeatureMethyl_smooth, 1, min)
  FeatureMethyl_range <- cbind(FeatureMethyl_min,FeatureMethyl_max)
  colnames(FeatureMethyl_range) <- c("min", "max")
  FeatureMethyl_range <- as.data.frame(FeatureMethyl_range)
  return(FeatureMethyl_range)
}



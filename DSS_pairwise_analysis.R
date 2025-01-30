library(stringr)
library(genomation)
library(DSS)
library(bsseq)
library(data.table)

#Exmaple analysis script for DSS pair-wise analysis of Tconv vs Treg cells

#import cov files
files <- list.files(path=".", pattern = "*.gz") # get file names
samples <- str_match(files,"(\\w+?)_S") # get dataframe of sample names
files <- as.list(files)
names <- as.data.frame(samples[,2]) # sample names
file.vector <- unlist(files)

#create dataset for DSS
BSobj <- read.bismark(files = file.vector,
                      colData = names,
                      loci = NULL,
                      rmZeroCov = FALSE,
                      strandCollapse = FALSE,
                      nThread = 1,
                      verbose = TRUE)

design <- as.data.frame(as.matrix(names))
design[,2] <- str_sub(design[,1],1,nchar(design[,1])-1)
colnames(design) <- c("Sample", "Group")
rownames(design) <- NULL
design[,2] <- as.factor(design[,2])


sampleNames(BSobj) <- c("P1_Tconv","P1_Treg","P2_Tconv","P2_Treg","P3_Tconv","P3_Treg","P4_Tconv","P4_Treg","P5_Tconv","P5_Treg", 
                        "P6_Tconv","P6_Treg","P7_Tconv","P7_Treg")



SCRIPT_Tconv_vsTreg <- DMLtest(BSobj, 
                                group1 = as.character(subset(design, Group %in% "Tconv")[,1]), 
                                group2 = as.character(subset(design, Group %in% "Treg")[,1]), 
                                smoothing = T, equal.disp = FALSE)

SCRIPT_Tconv_vsTreg.fdr05 <- subset(SCRIPT_Tconv_vsTreg, fdr <= 0.05)




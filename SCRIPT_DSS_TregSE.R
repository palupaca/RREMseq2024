library(stringr)
library(genomation)
library(DSS)
library(bsseq)
library(data.table)


#load DSS results

load("./SCRIPT_DSS_pairwise_all.RData")

DMLtest.trimmed <- SCRIPT_Tconv_vsTreg[,1:4]
DMLtest.trans <- trunc(DMLtest.trimmed[,3:4] * 100)
DMLtest.trans <- cbind(DMLtest.trimmed[,1:2], DMLtest.trans[,1:2])
colnames(DMLtest.trans)[3:4] <- c("Tconv", "Treg")

#export for seqmonk 
write.table(DMLtest.trans[,c(1,2,3)], file = "./!DSS_values_Tconv.txt", sep = "\t", row.names = F, quote = F)
write.table(DMLtest.trans[,c(1,2,4)], file = "./!DSS_values_Treg.txt", sep = "\t", row.names = F, quote = F)


#reshape DSS results
probe <- paste(SCRIPT_Tconv_vsTreg[,1], SCRIPT_Tconv_vsTreg[,2], sep = ":")
probe <- paste(probe, SCRIPT_Tconv_vsTreg[,2], sep = "-")
probe <- as.data.frame(probe)
probe <- paste0("Chr",probe[,1])

SCRIPT_Tconv_vsTreg.locations <- cbind(probe,SCRIPT_Tconv_vsTreg)

# get Treg-SE locations for CDF plot
Treg.SE.locations <- fread("./Treg_SE_CpGs.txt", select = c(1)) # from SeqMonk
Treg.SE.locations <- as.data.frame(Treg.SE.locations)

#filter for DSS values within TregSE
SCRIPT_Tconv_vsTreg.TregSE.DSS <- merge(Treg.SE.locations,
                                      SCRIPT_Tconv_vsTreg.locations,
                                      by.x = "Probe",
                                      by.y = 1, all = F)

par(pty="s", cex.lab = 1.5, cex.axis = 1.1, las = 1)

#plot CDF
plot(ecdf(SCRIPT_Tconv_vsTreg.TregSE.DSS[,4]),
     xlab = "beta",
     ylab = "Cumulative proportion",
     main = "Cumulative distribution function \nTreg-SE",
     col = "blue", lwd = 3, verticals = T, do.points = F)
lines(ecdf(SCRIPT_Tconv_vsTreg.TregSE.DSS[,5]), col = "red", lwd = 3, verticals = T, do.points = F)
legend("bottomright",
       legend = c("Tconv", "Treg"),
       col = c("blue", "red"),
       pch = 15)

#plot bar plots
plot(boxplot(Beta~CellType, data = SCRIPT_DSS_plot, range=0.0, horizontal=FALSE, varwidth=TRUE, notch=FALSE,
             boxwex = 0.6, col = c("blue", "red"), sep = ":", lex.order = TRUE)) 


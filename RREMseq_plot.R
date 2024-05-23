##master plotting script for the RREM-seq paper 

library(ggsci)
library(ggplot2)
library(ggfortify)
library(stringr)
library(pheatmap)
library(dplyr)
library(data.table)
library(tidyverse)
library(svglite)

#General setup for plotting

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  title = element_text(size = 20),
  legend.text = element_text(size = 18)
)

options(
  ggplot2.discrete.colour = ggsci::scale_colour_nejm,
  ggplot2.discrete.fill = ggsci::scale_fill_nejm
)




#plot global CpG coverage
CpGcov <- as.data.frame(fread("./uniqueCpG.sum.csv"))
CpGcov <- as.data.frame(fread("./uniqueCpGs_comp.csv"))

#reorder data
CpGcov <- CpGcov %>%
  mutate(Group = fct_relevel(Group, 
                            "WGBS", "WGEM-seq", "RRBS", 
                            "RREM-seq"))



dat <- summarise(group_by(CpGcov, Group),
                 mean = mean(UniqueCpGs),
                 sd = sd(UniqueCpGs))
p1 <- ggplot() + 
  geom_bar(data = dat,
           aes(y = mean, x = Group, fill = Group), stat="identity", width=0.75) + 
  geom_errorbar(data = dat,
                aes(y = mean, x = Group,
                    ymin= mean - sd,
                    ymax= mean + sd), stat="identity", width=0.4) + 
  geom_dotplot(data =  CpGcov, aes(y= UniqueCpGs, x = Group, fill = Group),binaxis='y', stackdir='center')  +
  ylim(c(0,25))+
  labs(y= "Unique CpGs/M") +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill='transparent'),
    #legend.box.background = element_rect(fill='transparent')
  )+
  theme(axis.line = element_line(color="black", size = 0.5))

p1+My_Theme+theme(legend.position="none")

#plot library conc. 

library.conc <- as.data.frame(fread("./libConc.csv"))

library.conc <- library.conc %>%
  mutate(Group = fct_relevel(Group, 
                             "WGBS_25ng", "WGEM-seq_25ng", "RRBS_25ng", 
                             "RREM-seq_25ng", "RREM-seq_2ng"))

library.conc <- library.conc %>%
  mutate(Methods = fct_relevel(Methods, 
                             "WGBS", "WGEM-seq", "RRBS", 
                             "RREM-seq"))
dat <- summarise(group_by(library.conc, Group),
                 mean = mean(Conc.),
                 sd = sd(Conc.))
dat <- as.data.frame(dat)

dat$Methods <- c("WGBS", "WGEM-seq", "RRBS", "RREM-seq", "RREM-seq", "RREM-seq")

dat <- dat %>%
  mutate(Methods = fct_relevel(Methods, 
                               "WGBS", "WGEM-seq", "RRBS", 
                               "RREM-seq"))

p2 <- ggplot() + 
  geom_bar(data = dat,
           aes(y = mean, x = Group, fill = Methods), stat="identity", width=0.75) + 
  geom_errorbar(data = dat,
                aes(y = mean, x = Group,
                    ymin= mean - sd,
                    ymax= mean + sd), stat="identity", width=0.4) + 
  geom_dotplot(data =  library.conc, aes(y= Conc., x = Group, fill = Methods),binaxis='y', stackdir='center')  +
  ylim(c(0,15))+
  labs(y= "LibraryConc(ng/ul)", x="") +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill='transparent'),
    #legend.box.background = element_rect(fill='transparent')
  )+
  theme(axis.line = element_line(color="black", size = 1))

libconc.all <- p2+My_Theme+theme(legend.position="none")+My_Theme+scale_fill_nejm()




#plot controls
phageControl <- as.data.frame(fread("./phageControl.csv"))

phageControl <- phageControl %>%
  mutate(GroupName = fct_relevel(GroupName, 
                                 "WGBS_25ng", "WGEM-seq_25ng", "RRBS_25ng", 
                                 "RREM-seq_25ng","RREM-seq_2ng"))
phageControl$CytosineConversionRate

dat_phage <- summarise(group_by(phageControl, GroupName),
                 mean = mean(CytosineConversionRate),
                 sd = sd(CytosineConversionRate))

dat_phage$Group <- c("WGBS", "WGEM-seq", "RRBS", "RREM-seq", "RREM-seq")

dat_phage <- dat_phage %>%
  mutate(Group = fct_relevel(Group, 
                             "WGBS", "WGEM-seq", "RRBS", "RREM-seq"))

p.phage <- ggplot() + 
  geom_bar(data = dat_phage,
           aes(y = mean, x = GroupName, fill = Group), stat="identity", width=0.75) + 
  geom_errorbar(data = dat_phage,
                aes(y = mean, x = GroupName,
                    ymin= mean - sd,
                    ymax= mean + sd), stat="identity", width=0.4) + 
  geom_dotplot(data =  phageControl, aes(y= CytosineConversionRate, x = GroupName, fill = Group),binaxis='y', stackdir='center')  +
  ylim(c(0,100))+
  labs(y= "CytosineConversionRate/%", x="") +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill='transparent'),
    #legend.box.background = element_rect(fill='transparent')
  )+
  theme(axis.line = element_line(color="black", size = 1))

p.phage+My_Theme+theme(legend.position="none")+My_Theme+scale_fill_nejm()



#plot CpG coverage
CpGcov <- CpGcov %>%
  mutate(GroupName = fct_relevel(GroupName, 
                                 "WGBS_25ng", "WGEM-seq_25ng", "RRBS_25ng", 
                                 "RREM-seq_25ng","RREM-seq_2ng"))

dat.CpGcov <- summarise(group_by(CpGcov, GroupName),
                        mean = mean(UniqueCpGs),
                        sd = sd(UniqueCpGs))

dat.CpGcov$Group <- c("WGBS", "WGEM-seq", "RRBS", "RREM-seq", "RREM-seq")

dat.CpGcov <- dat.CpGcov %>%
  mutate(Group = fct_relevel(Group, 
                             "WGBS", "WGEM-seq", "RRBS", "RREM-seq"))

p.CpGcovs <- ggplot() + 
  geom_bar(data = dat.CpGcov,
           aes(y = mean, x = GroupName, fill = Group), stat="identity", width=0.75) + 
  geom_errorbar(data = dat.CpGcov,
                aes(y = mean, x = GroupName,
                    ymin= mean - sd,
                    ymax= mean + sd), stat="identity", width=0.4) + 
  geom_dotplot(data =  CpGcov, aes(y= UniqueCpGs, x = GroupName, fill = Group),binaxis='y', stackdir='center')  +
  ylim(c(0,22))+
  labs(y= "Unique CpGs(M)", x="") +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill='transparent'),
    #legend.box.background = element_rect(fill='transparent')
  )+
  theme(axis.line = element_line(color="black", size = 1))

p.CpGcovs+My_Theme+theme(legend.position="none")+My_Theme+scale_fill_nejm()


CpGiref <- as.data.frame(fread("./cpgIslandExt.txt")) #from UCSC
CpGi <- as.data.frame(fread("./SCRIPT_run2n3_CpGiprobes.txt"))
CpGi <- na.omit(CpGi$P4_Tconv)                     
length(CpGi) / 22564
22564 -length(CpGi)

promoters <- read.table("./promoters.txt", quote = "", sep = "\t", header = T)
promoters <- na.omit(promoters$all)                     
length(promoters) / 63631
63631 -length(promoters)


#code by Sara Di Rienzi, 2024


library("dplyr")
library("ggpubr")
library("reshape2")
library("stringr")
library("data.table")


#Makes Fig2 or FigS5A 
#Fig S3C is similarly made from this plots

#start with genes desired from Supplemental Table 4,and subset the columns:
hormones<-read.delim(file="Fig2.txt",header=TRUE,sep = "\t")  #or load FigS5A.txt for Supplemental Fig 5A
GeTMMULDMindexes<-str_which(names(hormones), "GeTMM.ULDM")
GeTMMU6475indexes<-str_which(names(hormones), "GeTMM.U6475")
GeTMMILDMindexes<-str_which(names(hormones), "GeTMM.ILDM")
GeTMMI6475indexes<-str_which(names(hormones),"GeTMM.I6475") 
GeTMMI17938indexes<-str_which(names(hormones), "GeTMM.I1793")

hormones$GeTMMULDM6475<-rowMeans(hormones[,c(GeTMMULDMindexes, GeTMMU6475indexes)])
hormones$GeTMMILDM6475<-rowMeans(hormones[,c(GeTMMILDMindexes, GeTMMI6475indexes)])
hormones$GeTMMILDM17938<-rowMeans(hormones[,c(GeTMMILDMindexes, GeTMMI6475indexes)])
hormones$GeTMMI6475I17938<-rowMeans(hormones[,c(GeTMMI6475indexes,GeTMMI17938indexes)])
hormones$GeTMMI6475U6475<-rowMeans(hormones[c(GeTMMI6475indexes,GeTMMU6475indexes)])
hormones$GeTMMILDMULDM<-rowMeans(hormones[c(GeTMMILDMindexes,GeTMMULDMindexes)])

hormones$GeTMMaverage<-rowMeans(hormones[,c(GeTMMULDMindexes,GeTMMU6475indexes,GeTMMILDMindexes,GeTMMI6475indexes,GeTMMI17938indexes)])

U6475foldchange<-str_which(names(hormones),"U6475.ULDM4log2FoldChange")
I6475foldchange<-str_which(names(hormones),"I6475.ILDM4log2FoldChange")
I17938foldchange<-str_which(names(hormones),"I17938.ILDM4log2FoldChange")
I6475I17938foldchange<-str_which(names(hormones),"I6475.I17938log2FoldChange")
I6475U6475foldchange<-str_which(names(hormones),"I6475.U6475log2FoldChange")
ILDMULDMfoldchange<-str_which(names(hormones),"ILDM4.ULDM4log2FoldChange")

hormoneshort<-cbind(hormones[,c(1:8,U6475foldchange,I6475foldchange,I17938foldchange,I6475I17938foldchange,I6475U6475foldchange,ILDMULDMfoldchange)])
names(hormoneshort)[c(9:14)]<-c("U6475-ULDM4","I6475-ILDM4","I17938-ILDM4","I6475-I17938","I6475-U6475","ILDM4-ULDM4")

#need to collaspse on all of the fold changes
melted<-reshape2::melt(hormoneshort,id.vars=names(hormoneshort)[c(1:8)], variable.name="Comparison", value.name="FoldChange")
melted$Comparison = factor(melted$Comparison,c("U6475-ULDM4","I6475-ILDM4","I17938-ILDM4","I6475-I17938","ILDM4-ULDM4","I6475-U6475" ))

#pull out geTMM data to use to color the bars
GeTMMdata<-cbind(hormones[,c(1,116:122)])
names(GeTMMdata)[c(2:7)]<-c("U6475-ULDM4","I6475-ILDM4","I17938-ILDM4","I6475-I17938",
                            "I6475-U6475","ILDM4-ULDM4")
melted2<-reshape2::melt(GeTMMdata,id.vars=names(GeTMMdata)[c(1,8)], variable.name="Comparison", value.name="GeTMM")

merged1<-merge(melted,melted2,by=c("EnsemblID","Comparison"),all.x=TRUE,all.y=FALSE)

#pull out the DEG part to use as alpha overlay
hormoneshort2<-cbind(hormones[,c(1,9:14)])
names(hormoneshort2)[c(2:7)]<-c("U6475-ULDM4","I6475-ILDM4","I17938-ILDM4","I6475-I17938","I6475-U6475","ILDM4-ULDM4")
melted3<-reshape2::melt(hormoneshort2,id.vars=names(hormoneshort2)[1], variable.name="Comparison", value.name="DEG")
melted3$DEG<-abs(as.numeric(melted3$DEG))

merged2<-merge(merged1,melted3,by=c("EnsemblID","Comparison"),all.x=TRUE,all.y=FALSE)

merged2$InflGroup<-paste(merged2$Cluster,merged2$Gene,sep=" ")

#fix the order
merged2$InflGroup<-factor(merged2$InflGroup,rev(levels(as.factor((merged2$InflGroup)))))


merged2<-as.data.table(merged2)

plot<-ggplot(merged2,aes(InflGroup,as.numeric(FoldChange),fill=as.numeric(log10(GeTMM)),alpha=DEG)) +
  geom_bar(stat="identity",show.legend = TRUE,colour="black")+
  facet_wrap(~Comparison,ncol=6)+
  coord_flip()+
  scale_fill_gradient(low="blue", high="red")+
  ylab("Log 2 Fold change")+
  xlab("Hormone genes")+
  theme_bw() +
  theme(
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=16),
    axis.title.x= element_text(size=16),
    axis.text.x= element_text(size=12, colour = "black"),
    axis.text.y= element_text(size=12, colour = "black"), 
    legend.position="none", #turnon to see the legend
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size=16,face="bold")
  )

plot

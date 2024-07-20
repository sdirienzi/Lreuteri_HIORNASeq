#code by Sara Di Rienzi, 2024

#to make plots to show DEG function enrichments

library("dplyr")
library("ggpubr")
library("viridis")
library("gage")
library("pathview")
library("gageData")
library("reshape2")




DEGfunctions<-read.delim(file="Fig1C_DEGfunctions.txt",header=TRUE,sep="\t")


DEGfunctionsmelt<-reshape2::melt(DEGfunctions,id.vars=names(DEGfunctions)[c(1,5:7)],variable.name="Comparison",value.name = "Count")
names(DEGfunctionsmelt)[1]="Group"

DEGfunctionsmelt$Group = factor(DEGfunctionsmelt$Group,rev(c("Immune response","Metal/stress response",
                                                         "Signaling","Membrane component",
                                                         "Cell fate/growth","Hormone secretion",
                                                         "Mucus","Nutrient metabolism/response")))

colours=rev(c("orange","gold","black","royalblue","slateblue1","#F0027F","tan","seagreen1"))
colored<-setNames(colours,levels(DEGfunctionsmelt$Group))

ggplot(data=DEGfunctionsmelt,aes(x=Comparison,y=Count,fill=Group))+
  geom_bar(stat="identity",position = "stack")+
  labs(y="Pathway Count",x="Treatment") +
  scale_fill_manual(name = "Group",values=colored)+
  theme(plot.title = element_text(size = 16,hjust = 0.5)) + 
  theme_classic() +
  theme(
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x = element_text(size=20, colour="black"),
    axis.title.y= element_text(size=18, colour="black"),
    axis.text.x= element_text(size = 20, colour="black",angle=45,vjust=0.7),
    axis.text.y= element_text(size=20, colour="black"),
    legend.text=element_text(size=18))

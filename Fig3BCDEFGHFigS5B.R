#code by Sara Di Rienzi, 2024


require(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(data.table)
library(stringr)
library(ggsci)
library(DescTools)


data<-read.delim(file="Fig3BCDEFFigS5B.txt",sep="\t",header=TRUE)
data$IndTreat<-paste(data$ShortI,data$Treatment,sep="")
data$IndTreat<-factor(data$IndTreat,c("ULDM4","U6475","ILDM4","I6475","I17938"))

melted<-reshape2::melt(data,id.vars=names(data)[c(1,2,3,11)],variable.name="Hormone",value.name="Concentration")
meltable<-as.data.table(melted)


total<-ggplot(melted,aes(x=IndTreat, y=as.numeric(Concentration),fill=as.factor(IndTreat))) +  
  geom_boxplot(show.legend = FALSE)+
  geom_point(position=position_jitterdodge(), size=2,pch=21,show.legend = FALSE)+
  facet_wrap(~Hormone,strip.position="top",scales = "free_y",ncol = 4) + 
  xlab("Treatment") +
  ylab(bquote('pg/ml')) +
  scale_fill_npg()+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=14),
    axis.title.x= element_text(size=14),
    axis.text.x= element_text(size=12, colour = "black", angle=45,vjust=0.6),
    plot.title = element_text(size=14,face="bold"),
    axis.text.y= element_text(size=12, colour = "black"),
    panel.border = element_blank(),
    strip.text = element_text(size=14,face="bold"),
    strip.background = element_blank()
   )

total


#do stats on the data

#first do ANOVA and then Dunnett's

hormonelist<-levels(as.factor(meltable$Hormone))
pvaluelist<-c()
#i<-1
for (i in 1:length(1:7)) {
hormonename<-hormonelist[i]
aovdata<-meltable[Induction == "Induced" & Hormone == hormonename]
oneway<-aov(Concentration ~ Treatment, data=aovdata)
summaryaov<-summary(oneway)
pval<-summaryaov[[1]]$`Pr(>F)`[1]
pvaluelist<-c(pvaluelist,pval)
}
names(pvaluelist)<-hormonelist
#Induced:
# Amylin    C.Peptide      Ghrelin    GIP.total        MCP.1           PP          PYY 
# 4.788223e-04 2.760552e-02 1.076050e-01 3.454448e-02 9.130631e-05 2.605713e-01 1.160966e-05 

pvaluelistuninuced<-c()
for (i in 1:length(1:7)) {
  hormonename<-hormonelist[i]
  aovdata<-meltable[Induction == "Uninduced" & Hormone == hormonename]
  oneway<-aov(Concentration ~ Treatment, data=aovdata)
  summaryaov<-summary(oneway)
  pval<-summaryaov[[1]]$`Pr(>F)`[1]
  pvaluelistuninuced<-c(pvaluelistuninuced,pval)
}
names(pvaluelistuninuced)<-hormonelist
pvaluelistuninuced
# Amylin C.Peptide   Ghrelin GIP.total     MCP.1        PP       PYY 
# NaN       NaN 0.3818823 0.3739010 0.3822538 0.3016092       NaN 

#make subplot of Amylin, Ghrelin, GIP, PP, PYY
#make subplot of #Cpep and MCP


#do Dunnett's test with LDM4 on induced as ref level for Amylin, Cpep, GIP, MCP, PYY

library("DescTools")
DunnettTest(Concentration~as.factor(Treatment), data=meltable[Induction=="Induced"&Hormone=="PYY"], control="LDM4")

#Amylin
#                diff   lwr.ci   upr.ci    pval    
# 17938-LDM4 18.51667 11.43036 25.60298 0.00053 ***
# 6475-LDM4  17.44333 10.35702 24.52964 0.00073 ***

#GIP
# diff    lwr.ci   upr.ci   pval    
# 17938-LDM4 537.7733  90.81323 984.7334 0.0241 *  
# 6475-LDM4  367.4233 -79.53677 814.3834 0.0969 .  

#PYY
# diff   lwr.ci    upr.ci    pval    
# 17938-LDM4 94.52 76.76827 112.27173 9.1e-06 ***
# 6475-LDM4  74.67 56.91827  92.42173 3.6e-05 ***

DunnettTest(Concentration~as.factor(Treatment), data=meltable[Induction=="Induced"&Hormone=="MCP.1"], control="LDM4")
# diff    lwr.ci    upr.ci    pval    
# 17938-LDM4 -373.6000 -484.4818 -262.7182 0.00013 ***
# 6475-LDM4  -382.1633 -493.0452 -271.2815 0.00011 ***

DunnettTest(Concentration~as.factor(Treatment), data=meltable[Induction=="Induced"&Hormone=="C.Peptide"], control="LDM4")
# diff    lwr.ci   upr.ci   pval    
# 17938-LDM4 28.26333 -3.138077 59.66474 0.0722 .  
# 6475-LDM4  39.59000  8.188590 70.99141 0.0198 * 
#   



#Sertonin plot
data<-read.delim(file="Fig3GH.txt",sep="\t",header=TRUE)
data$TreatmentName<-factor(data$TreatmentName,c("ULDM4","U6475","ILDM4","I6475","I17938"))


shapevalues<-c(21,22)
names(shapevalues)<-c("1","2")


serotonin<-ggplot(data,aes(x=TreatmentName, y=as.numeric(Serotonin),fill=as.factor(TreatmentName))) +  
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+ 
  geom_point(aes(shape=as.factor(Experiment)),position=position_jitterdodge(), size=2,show.legend = FALSE)+
  facet_wrap(~Side,strip.position="top",scales = "fixed",ncol = 2) + 
  scale_shape_manual(values=shapevalues)+
  xlab("Treatment") +
  ylab(bquote('ng/ml')) +
  scale_y_continuous(expand=expansion(mult= c(0.03,.35)))+
  scale_fill_npg()+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=14),
    axis.title.x= element_text(size=14),
    axis.text.x= element_text(size=12, colour = "black", angle=45,vjust=0.6),
    plot.title = element_text(size=14,face="bold"),
    axis.text.y= element_text(size=12, colour = "black"),
    panel.border = element_blank(),
    strip.text = element_text(size=14,face="bold"),
    strip.background = element_blank(),
    legend.position="none"
  )

serotonin
datatable<-as.data.table(data)
DunnettTest(Serotonin~as.factor(Treatment), data=datatable[Induction=="Induced" &Side=="Basolateral"], control="Media")
#serotonin, apical
# diff    lwr.ci   upr.ci   pval    
# 17938-Media 193.0713 78.995168 307.1474 0.0017 ** 
#   6475-Media  122.6796  8.603501 236.7558 0.0350 *  

#serotonin, basolateral 
# $Media
# diff   lwr.ci   upr.ci    pval    
# 17938-Media 584.0407 426.8387 741.2428 3.5e-07 ***
#   6475-Media  461.7389 304.5369 618.9409 6.3e-06 ***


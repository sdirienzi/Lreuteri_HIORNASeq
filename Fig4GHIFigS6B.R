#code by Sara Di Rienzi, 2024

require(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(data.table)
library(stringr)
library(ggsci)
library(lme4)
library(emmeans)


data<-as.data.table(read.delim(file="Fig4FGHIFigS6B.txt",sep="\t",header=TRUE))

data$Treatment<-factor(data$Treatment,c("LDM4","6475"))
data$Hormone<-factor(data$Hormone, c("AVP", "LHB", "Adipolin", "KISS1"))


shapevalues<-c(23,24,25)
names(shapevalues)<-c("LG5","LG6","LG7")
hormoneplot<-ggplot(data,aes(x=Treatment, y=as.numeric(Concentration),fill=as.factor(Treatment))) +  
  geom_boxplot(show.legend = FALSE,color="black",alpha=0.4,outlier.shape = NA)+
  geom_point(aes(shape = LifeGift),position=position_jitterdodge(), size=2,show.legend = FALSE)+
  facet_wrap(~Hormone,strip.position="top",scales = "free_y",ncol = 3) + 
  scale_shape_manual(values=shapevalues)+
  xlab("Treatment") +
  ylab(bquote('pg/ml')) +
  scale_y_continuous(expand=expansion(mult= c(0.03,.35)))+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF","#E64B35FF", "#4DBBD5FF"))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=12),
    axis.title.x= element_text(size=12),
    axis.text.x= element_text(size=10, colour = "black", angle=45,vjust=0.6),
    plot.title = element_text(size=12,face="bold"),
    axis.text.y= element_text(size=10, colour = "black"),
    panel.border = element_blank(),
    strip.text = element_text(size=12,face="bold"),
    panel.spacing=unit(2,"lines"),
    strip.background = element_blank()
  )

hormoneplot
#switch for the hormone being tested
model<-lmer(Concentration~Treatment+(1|LifeGift), data=subset(data,data$Hormone=="KISS1"), control=lmerControl(optimizer = "bobyqa"))
nullmodel<-lmer(Concentration ~ (1|LifeGift),data=subset(data,data$Hormone=="KISS1"), control=lmerControl(optimizer = "bobyqa"))
anova(model, nullmodel)

performance::r2(model)
Modelemmeans<-emmeans(model,pairwise~Treatment)
Modelemmeans



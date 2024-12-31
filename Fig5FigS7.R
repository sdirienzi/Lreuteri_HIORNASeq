#code by Sara Di Rienzi, 2024

library(ggplot2)
library(dplyr) # easier data wrangling 
library(viridis) # colour blind friendly palette, works in B&W also
library(ggExtra) 
library(tidyr)
library(reshape2)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(limma)
library(gridExtra)
library(cluster)    # clustering algorithms
library(factoextra)
library(purrr)
library(Biostrings)
library(vegan)
library(phyloseq)
library(gridExtra)
library(data.table)
library(stringr)
library(ggpubr)


#Fig 5A, PCoA using imputed data

metabolitedata<-read.delim(file="Fig5FigS7MediaBatchnormImputedData.txt",header = TRUE,sep="\t", na.strings = "NA",check.names = FALSE)

metadata<-read.delim(file="Fig5FigS7metadata.txt",header=TRUE,sep="\t")

otumatrix<-data.matrix(cbind(metabolitedata[,2:258]))
rownames(otumatrix)<-metabolitedata[,1]
otumatrix<-otumatrix[complete.cases(otumatrix), ]
sampletable<-otu_table(otumatrix,taxa_are_rows=FALSE)

metadata<-sample_data(metadata)
sample_names(metadata)<-metadata$Sample

phyloseq<-merge_phyloseq(sampletable,metadata)
bray_pcoa = ordinate(phyloseq, "PCoA", "bray")
jaccard_pcoa = ordinate(phyloseq, "PCoA", "jaccard")
ordplot<-plot_ordination(phyloseq ,bray_pcoa, justDF=TRUE )
ordplotJac<-plot_ordination(phyloseq ,jaccard_pcoa, justDF=TRUE )
plot_ordination(phyloseq ,bray_pcoa) #axis 1 = 55.3 axis 2 =33.5
plot_ordination(phyloseq ,jaccard_pcoa) #axis 1 = 38.9, #axis. 2 = 30.7
ordtable<-as.data.table(ordplot) 
ordtable_Jac<-as.data.table(ordplotJac)

pdf(file="MDS_bray.pdf",width=4, height=3)
ggplot(ordtable, aes(x=Axis.1,y=Axis.2,fill=Treatment))+
  geom_point(size= 2, pch=21, colour="black") +
  stat_ellipse()+
  xlab("Axis 1 (55.3%)")+
  ylab("Axis 2 (33.5%)")+
  scale_fill_manual(values=c("gold","orchid1","purple2")) +
  theme_classic() + theme(
    axis.line.x = element_line(colour = "black", size = 0.5), 
    axis.line.y = element_line(colour = "black", size = 0.5), 
    axis.text.x= element_text(size=8, colour = "black"),
    axis.text.y= element_text(size=8, colour = "black"),
    axis.title = element_text(size=10, colour="black"),
      legend.position="none"
  )
dev.off()

dm<-as.data.frame((sampletable))
metadata<-read.delim(file="Fig5FigS7metadata.txt",header=TRUE,sep="\t")


adonis2(dm ~ Treatment + Replicate, data = metadata,method="bray")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dm ~ Treatment + Replicate, data = metadata, method = "bray")
# Df SumOfSqs      R2       F Pr(>F)    
# Treatment  2  1.24402 0.82212 51.8812  0.001 ***
#   Replicate  1  0.10131 0.06695  8.4501  0.002 ** 
#   Residual  14  0.16785 0.11092                   
# Total     17  1.51318 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#supplemental Fig 5B, PCoA using imputed data but without malate and fumarate (chem IDs 409 and 330)

metabolitedatanomalnofum<-read.delim(file="Fig5FigS7MediaBatchnormImputedData_NoMalateOrFumarate.txt",header = TRUE,sep="\t", na.strings = "NA",check.names = FALSE)
metadata<-read.delim(file="Fig5FigS7metadata.txt",header=TRUE,sep="\t")

otumatrix2<-data.matrix(cbind(metabolitedatanomalnofum[,2:256]))
rownames(otumatrix2)<-metabolitedatanomalnofum[,1]
otumatrix2<-otumatrix2[complete.cases(otumatrix2), ]
sampletable<-otu_table(otumatrix2,taxa_are_rows=FALSE)
metadata<-sample_data(metadata)
sample_names(metadata)<-metadata$Sample

phyloseq<-merge_phyloseq(sampletable,metadata)
bray_pcoa = ordinate(phyloseq, "PCoA", "bray")
ordplot<-plot_ordination(phyloseq ,bray_pcoa, justDF=TRUE )
plot_ordination(phyloseq ,bray_pcoa) #axis 1 = 60.2 axis 2 =24.9
ordtable<-as.data.table(ordplot) 

pdf(file="MDS_bray_nomalateorfumarate.pdf",width=4, height=3)
ggplot(ordtable, aes(x=Axis.1,y=Axis.2,fill=Treatment))+
  geom_point(size= 2, pch=21, colour="black") +
  stat_ellipse()+
  xlab("Axis 1 (60.2%)")+
  ylab("Axis 2 (24.9%)")+
  scale_fill_manual(values=c("gold","orchid1","purple2")) +
  theme_classic() + theme(
    axis.line.x = element_line(colour = "black", size = 0.5), 
    axis.line.y = element_line(colour = "black", size = 0.5), 
    axis.text.x= element_text(size=8, colour = "black"), 
    axis.text.y= element_text(size=8, colour = "black"),
    axis.title = element_text(size=10, colour="black"),
    legend.position="none"
  )
dev.off()

dm<-as.data.frame((sampletable))
metadata<-read.delim(file="Fig5FigS7metadata.txt",header=TRUE,sep="\t")


adonis2(dm ~ Treatment + Replicate, data = metadata,method="bray")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dm ~ Treatment + Replicate, data = metadata, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# Treatment  2  0.74400 0.83291 86.242  0.001 ***
#   Replicate  1  0.08886 0.09948 20.601  0.001 ***
#   Residual  14  0.06039 0.06761                  
# Total     17  0.89325 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Supplemental Fig 7A, B
#correlation matrix
#build the metabolitematrix as for the heatmap but need to use imputed data
metabolitedata<-read.delim(file="Fig5FigS7MediaBatchnormImputedData.txt",header = TRUE,sep="\t", na.strings = "NA",check.names = FALSE)

metabolitematrix<-as.matrix(metabolitedata[,2:258])

#samplenames
rownames(metabolitematrix)<-metabolitedata[,1]

cormatrix<-cor(t(metabolitematrix), method="spearman") #transposed allow you to compare the strain. the other way the metabolites

cormatrix[lower.tri(cormatrix, diag=TRUE)]<-NA

correlationmatrix<-cormatrix

LR6475subindexes<-str_which(rownames(correlationmatrix),"6475")
LR17938subindexes<-str_which(rownames(correlationmatrix),"17938")
LDM4subindexes<-str_which(rownames(correlationmatrix),"LDM4")

LR6475subvalues<-unlist(as.list(correlationmatrix[LR6475subindexes,LR6475subindexes]))
LR17938subvalues<-unlist(as.list(correlationmatrix[LR17938subindexes,LR17938subindexes]))
LDM4subvalues<-unlist(as.list(correlationmatrix[LDM4subindexes,LDM4subindexes]))
LDM4LR6475subvalues<-unlist(as.list(correlationmatrix[LDM4subindexes, LR6475subindexes]))
LDM4LR17938subvalues<-unlist(as.list(correlationmatrix[LDM4subindexes,LR17938subindexes]))
LR6475LR17938subvalues<-unlist(as.list(correlationmatrix[LR6475subindexes,LR17938subindexes]))

corredata<-data.frame("LDM4-LDM4" = LDM4subvalues, "LR6475-LR6475" = LR6475subvalues,
                      "LR17938-17938" = LR17938subvalues, "LDM4-6475" = LDM4LR6475subvalues, 
                      "LDM4-17938" = LDM4LR17938subvalues, "LR6475-LR17938" =LR6475LR17938subvalues )


melt1<-reshape2::melt(corredata, variable.name="Comparison", value.name="Correlation")

boxplotcorrelations<-ggboxplot(melt1,x="Comparison",y="Correlation",
                               add="jitter", linewidth=0.5, add.params = list(size = 0.5))+
  theme(
    axis.line.x = element_line(colour = "black", size = 0.5), 
    axis.line.y = element_line(colour = "black", size = 0.5), 
    axis.text.x= element_text(size=8, colour = "black", angle=45,vjust=0.7),
    axis.text.y= element_text(size=8, colour = "black"),
    axis.title = element_text(size=10, colour="black"),
    legend.position="none"
  )


pdf(file="metabolitecormatix.pdf",width=4,height=3)
boxplotcorrelations
dev.off()


library("rlist")
my_comparisons <- list()
names(corredata)
for(i in 1: (length(names(corredata))-1) ){
  for(j in 2:(length(names(corredata)) )){
    if (i !=j & j > i) {
      print(c("i", i, "j", j))
      a<-names(corredata)[i]
      b<-names(corredata)[j]
      compx<-c(a,b)
      my_comparisons<-list.append(my_comparisons,compx)
    }
  }
}

statdata<-compare_means(Correlation~Comparison, melt1, method = "wilcox", p.adjust.method = "holm")   #using holm to reduce false positives
write.table(statdata, file="stats.txt",sep="\t")

mycomps<-list()
indexes<-which(statdata$p.adj<0.05)
for ( i in 1:length(indexes) ) {
  index<-indexes[i]
  comp1<-statdata$group1[index]
  comp2<-statdata$group2[index]
  compy<-c(comp1,comp2)
  mycomps<-list.append(mycomps,compy)
}


#here are the means of the comparisons
melt2<-na.omit(melt1)
aggregate(melt2[,2], list(melt2$Comparison), mean, na.action = na.omit)
#use the stats file to create the next two files manually. 

#this file shows the mean differences of the comparisons; make this from the above calculations
cormapmeandiff<-as.matrix(read.delim(file="FigS7meandifferences.txt",sep="\t",header=TRUE,row.names = 1))

#this file is a matrix of which comparisons had significant mean differences; make this from the above calculations
cormapstats<-as.matrix(read.delim(file="FigS7correlationplotstats.txt",sep="\t",header=TRUE,row.names=1))

col_fun= circlize::colorRamp2(c(0, -1,-1.3, -6,-16), c("navyblue","cornflowerblue","white", "lightcoral","maroon")) #redblue
col_fun= circlize::colorRamp2(breaks=c(0, -1.4,-2,-4), c("navyblue","white","lightcoral","maroon")) #redblue


circleheat<-Heatmap(cormapstats, name = "p.adj", rect_gp = gpar(type = "none"),
                    col = col_fun,
                    show_row_names = TRUE, show_column_names = TRUE, 
                    row_dend_reorder = FALSE, column_dend_reorder = FALSE,
                    row_names_gp = gpar(fontsize = 10),
                    column_names_gp = gpar(fontsize = 10),
                    cluster_rows =FALSE, cluster_columns = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.circle(x = x, y = y, r = abs(cormapmeandiff[i, j]) * min(unit.c(width, height)*.5), 
                                  gp = gpar(fill = col_fun(log10(cormapstats[i, j])), col = NA))
                      grid.text(cormapstats[i, j], x, y, gp = gpar(fontsize = 9))
                    })
pdf(file="metabolitemap.pdf", width=5, height=4)
circleheat
dev.off()



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Supplemental Figure 7C

metabolitesummarizeddata<-read.delim(file="Fig5FigS7FoldChangeData.txt",header = TRUE,sep="\t", na.strings = "NA",check.names = FALSE)

names(metabolitesummarizeddata)
length(which(metabolitesummarizeddata$LR6475_LDM4_qval <0.05 & metabolitesummarizeddata$LR6475_LDM4_pval<0.05)) #185
length(which(metabolitesummarizeddata$LR17938_LDM4_qval <0.05 & metabolitesummarizeddata$LR17938_LDM4_pval <0.05)) #159

LR6475sig<-which(metabolitesummarizeddata$LR6475_LDM4_qval <0.05 & metabolitesummarizeddata$LR6475_LDM4_pval<0.05)
LR17938sig<-which(metabolitesummarizeddata$LR17938_LDM4_qval <0.05 & metabolitesummarizeddata$LR17938_LDM4_pval <0.05)

#length(LR17938sig) <- max(length(LR6475sig), length(LR17938sig))
notsharedsig6475<-setdiff(LR6475sig,LR17938sig)
write.table(metabolitesummarizeddata[notsharedsig6475,],file="6475onlysig.txt", sep = "\t",row.names=FALSE)
notsharedsig17938<-setdiff(LR17938sig,LR6475sig)
write.table(metabolitesummarizeddata[notsharedsig17938,],file="17938onlysig.txt", sep = "\t",row.names=FALSE)

nonoverlapping<-c(notsharedsig6475,notsharedsig17938)
overlapping<-union(LR6475sig,LR17938sig)

overlappedids<-metabolitesummarizeddata[overlapping,]

write.table(overlappedids, file="overlappingmetabolites.txt",sep="\t",row.names=FALSE)


#do wide to long
metabolitesummarizeddata<-read.delim(file="Fig5FigS7FoldChangeData.txt",header = TRUE,sep="\t", na.strings = "NA",check.names = FALSE)

metabolitesummarizedlong<-reshape2::melt(metabolitesummarizeddata, id.vars=names(metabolitesummarizeddata)[13:17], 
                               measure.vars=names(metabolitesummarizeddata)[4:5],
                               variable.name="Strain", value.name= "ControlEnrichment")



pdf("metaboliteloadabs.pdf",width=2, height=2.8)
ggplot(metabolitesummarizedlong,  aes(x=Strain, y=abs(log(ControlEnrichment)),fill=Strain)) +
  geom_violin(linewidth = 0.5, color = "black")+
  ylab("Absolute log of enrichment")+
  xlab("Comparison")+
  scale_x_discrete(labels=c("LR6475-LDM4", "LR17938-LDM4"))+
  scale_fill_manual(values=c("purple2","orchid1"))+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.x= element_text(size=10),
    axis.title.y= element_text(size=10),
    axis.text.x= element_text(size=8, colour = "black",angle=45,vjust=0.5),
    axis.text.y= element_text(size=8, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    legend.position="none", plot.title=element_blank())
dev.off()

wilcox.test(abs(log(ControlEnrichment))~Strain,data=metabolitesummarizedlong)
# data:  abs(log(ControlEnrichment)) by Strain
# W = 39051, p-value = 0.0003447
# alternative hypothesis: true location shift is not equal to 0

Lreu6475<-subset(metabolitesummarizedlong,metabolitesummarizedlong$Strain=="LR6475_LDM4")
mean(abs(log(Lreu6475$ControlEnrichment))) #0.9932056

Lreu17938<-subset(metabolitesummarizedlong,metabolitesummarizedlong$Strain=="LR17938_LDM4")
mean(abs(log(Lreu17938$ControlEnrichment))) #0.8438429




#Fig 5C

#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html#plot-the-heamtap-list
#https://www.nature.com/articles/nbt1205-1499


#read in metabolite table:

metabolitedatanotimp<-read.delim(file="Fig5MediaBatchnormData.txt",header = TRUE,sep="\t", na.strings = "NA",check.names = FALSE)
metadata<-read.delim(file="Fig5FigS7metadata.txt",header=TRUE,sep="\t")
filteredfunctiondata<-read.delim(file="Fig5HeatMapMetabolites.txt",header = TRUE,sep="\t", na.strings = "NA",check.names = FALSE ) #selected metabolites sig compared to control and between Lreu strains
metabolitemetadata<-read.delim(file="Fig5metabolitemetadata.txt",header = TRUE,sep="\t", na.strings = "NA",check.names = FALSE )
metabolitemetadata$CHEM_ID<-as.character(metabolitemetadata$CHEM_ID)
filteredfunctiondata<-subset(filteredfunctiondata,filteredfunctiondata$Plot=="Y")

#subset the metabolitedata to just those in filteredfunctiondata
keeplist<-sort(filteredfunctiondata$ChemicalID)
keeppositions<-which(colnames(metabolitedatanotimp[,c(1:258)]) %in% as.character(keeplist),arr.ind = T)


#make a matrix of the data to plot, currently sorted by chemcial ID
metabolitematrix<-as.matrix(cbind(metabolitedatanotimp[,keeppositions]))

#samplenames
rownames(metabolitematrix)<-metabolitedatanotimp[,1]


#change the colnames to the order of the groupings
filteredfunctiondatatoplot<-subset(filteredfunctiondata,filteredfunctiondata$ChemicalID %in% colnames(metabolitematrix))
#get the plot order
PlotOrderNames<-filteredfunctiondatatoplot$PlotOrder[order(filteredfunctiondatatoplot$ChemicalID)]
colnames(metabolitematrix)<-PlotOrderNames
#make a matrix of the data to plot, currently sorted by chemcial ID

#now sort by Plot order
sorted_mat <- metabolitematrix[, order(as.numeric(colnames(metabolitematrix)))]

#then convert to biochem name
sortednamesBioChem <- filteredfunctiondatatoplot$BiochemicalName[order(filteredfunctiondatatoplot$SubPathway)]
colnames(sorted_mat)<-sortednamesBioChem


#transpose, center and scale the data
myX<- t(scale(sorted_mat,center=FALSE,scale=TRUE)) #scale works on the columns; here the columns are metabolites but the heatmap will have the samples as columns
#so to compared across metabolites, do scale, then transpose


#Add sample annotations

table(metadata$Treatment)
table(metadata$Replicate)

metause<- data.frame(Treatment = metadata$Treatment,Replicate=metadata$Replicate)


treatment.cols<- c("purple2","gold","orchid1")
replicate.cols<-c("gray","black")

treatment.cols.assigned<- setNames(treatment.cols, c("Lreu6475", "LDM4","Lreu17938"))
replicate.cols.assigned<- setNames(replicate.cols,unique(as.character(metause$Replicate)))
myha<- HeatmapAnnotation(df = metause, col = list(Treatment = treatment.cols.assigned,
                                                  Replicate = replicate.cols.assigned))

genelistmeta<-data.frame(Group = filteredfunctiondatatoplot$SubPathway)
length(unique(filteredfunctiondatatoplot$SubPathway))

#to color the annotations
palettelength<-length(unique(filteredfunctiondatatoplot$SubPathway))
palette <- c(brewer.pal("Paired",n=12), brewer.pal("Paired",n=palettelength-12))
group.cols<-c(palette )
sortedgroup <- filteredfunctiondatatoplot$SubPathway[order(filteredfunctiondatatoplot$PlotOrder)]
group.cols.assigned<-setNames(group.cols,c(unique(sortedgroup)))
geneha<-rowAnnotation(df=genelistmeta,col=list(Group = group.cols.assigned))

heatmap<-Heatmap(myX, name = "MediaBatchnormData", 
                    show_row_names = TRUE, show_column_names = TRUE, 
                    column_dend_reorder = TRUE, row_dend_reorder = FALSE,
                   cluster_rows = FALSE,
                    clustering_distance_columns = "euclidean",
                    top_annotation = myha,
                  left_annotation = geneha,
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 10),
                 show_heatmap_legend = TRUE,
                 width = unit(5, "cm"), height = unit(18, "cm")
                 )
pdf(file="heatmap_small.pdf",width=12, height=9)
draw(heatmap)
dev.off()







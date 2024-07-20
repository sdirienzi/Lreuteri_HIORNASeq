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


#correlationmatrix work is down further; the first part makes a heatmap based on all DEGs not just those annotated

#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html#plot-the-heamtap-list
#https://www.nature.com/articles/nbt1205-1499

#get list of genes
genelist<-read.delim(file="FigS4genes.txt",header = TRUE,sep="\t", na.strings = "NA")
allDEGS<-read.delim(file="FigS4DEGs.txt",header=TRUE,sep="\t",na.strings="NA")
allDEGSIds<-data.frame("EnsemblID" =allDEGS$EnsemblID)
suptable2<-read.delim(file = "FigS4supplementaltable.txt",header = TRUE,sep = "\t")
metadata<-read.delim(file="FigS4metadata.txt",header=TRUE,sep="\t")
fcdata<-read.delim(file="FigS4completeDESeq2output.txt",header=TRUE,sep="\t")
fcs22<-data.frame("Gene" = fcdata$Gene, "U6475" = fcdata$U6475.ULDM4log2FoldChange, "I6475" = fcdata$I6475.ILDM4log2FoldChange,
                  "I17938" = fcdata$I17938.ILDM4log2FoldChange, "I647517938" = fcdata$I6475.I17938log2FoldChange ,
                  "I6475U6475" = fcdata$I6475.U6475log2FoldChange, "IU" = fcdata$ILDM4.ULDM4log2FoldChange)



#get expression values for each sample
allgenes<-merge(allDEGSIds,genelist,by.x="EnsemblID",by.y="EnsemblID",all=TRUE)
heatmapdata<-merge(allgenes,suptable2,by.x=c("EnsemblID"),by.y=c("EnsemblID"),all.x=TRUE,all.y=FALSE)
heatmapdat<-merge(heatmapdata,fcs22,by.x="EnsemblID",by.y="Gene",all.x=TRUE,all.y=FALSE)
heatmapdata<-heatmapdat

heatmapdatamelt<-reshape2::melt(heatmapdata,id.vars=names(heatmapdata[c(1:14,53:61) ]),measure.vars =names(heatmapdata[15:44]), variable.name = "Sample", value.name ="Expression" )

genematrix<-cbind(heatmapdata[,c(15:44)])


rownames(genematrix)<-heatmapdata$EnsemblID#[c(1:30)]

mybatch <- c("A","A","A","B","B","B","A","A","A","B","B","B","A","A","A","B","B","B",
           "A","A","A","B","B","B","A","A","A","B","B","B")




mydesign<-model.matrix(~TreatmentInduction,metadata)

batchedgenematrix<-removeBatchEffect(genematrix,batch = mybatch,design = mydesign)

mymatrix<-batchedgenematrix  #genematrix is the uncorrected matrix

#for the heatmap plot, rescale the expression values
myX<- t(scale(t(mymatrix),center=TRUE,scale=TRUE))

#get the number of clusters
kmean_withinss <- function( k) {
  cluster <- kmeans(myX, k)
  return (cluster$tot.withinss)
}

max_k<-20
wss<-sapply(2:max_k,kmean_withinss)
elbow <-data.frame(2:max_k, wss)
ggplot(elbow, aes(x = X2.max_k, y = wss)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 20, by = 1))


# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(myX, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(myX))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")


set.seed(123)
gap_stat <- clusGap(myX, FUN = kmeans, nstart = 25,
                    K.max = 15, B = 100,iter.max=30)
# Print the result
print(gap_stat, method = "firstmax")

fviz_gap_stat(gap_stat)




metause<- data.frame(Induction = metadata$Induction, Treatment = metadata$Treatment,Replicate=metadata$Replicate)

genelistmeta<-data.frame(Group = heatmapdata$MiddleGroup, #[c(1:30)]
                         U6475=heatmapdata$U6475.ULDM, #[c(1:30)]
                         I6475=heatmapdata$I6475.ILDM, #[c(1:30)]
                         I17938=heatmapdata$I1793.ILDM,
                         I6475I17938 = heatmapdata$I6475.I1793,
                         I6475U6475 = heatmapdata$I6475.U6475,
                         IU = heatmapdata$ILDM.ULDM) #[c(1:30)]


genelistmeta2<-data.frame(Group = heatmapdata$MiddleGroup, #[c(1:30)]
                         U6475=heatmapdata$U6475, #[c(1:30)]
                         I6475=heatmapdata$I6475, #[c(1:30)]
                         I17938=heatmapdata$I17938,
                         I6475I17938 = heatmapdata$I647517938,
                         I6475U6475 = heatmapdata$I6475U6475,
                         IU = heatmapdata$IU) #[c(1:30)]

table(metadata$Induction)
table(metadata$Treatment)
table(metadata$Replicate)
table(genelistmeta$Group)
table(genelistmeta$I6475)

induction.cols<-c("blue", "red")
treatment.cols<- c("gold","purple2","orchid1")
replicate.cols<-c("gray","black")
group.cols<-c("orange","slateblue1","seagreen1","pink","black",
              "blue","#F0027F","steelblue1","gold","tan", "aquamarine", "gray47","khaki1",
              "chocolate1","deepskyblue")
U6475.cols<-c("white","maroon","navyblue")
I6475.cols<-c( "navyblue","white","maroon")
I17938.cols<-c("white", "navyblue","maroon")
I6475I17938.cols<-c("white", "maroon","navyblue")
I6475U6475.cols<-c("white", "maroon","navyblue")
IU.cols<-c("white", "maroon","navyblue")


induction.cols.assigned<- setNames(induction.cols, unique(as.character(metadata$Induction)))
treatment.cols.assigned<- setNames(treatment.cols, unique(as.character(metadata$Treatment)))
replicate.cols.assigned<- setNames(replicate.cols,unique(as.character(metadata$Replicate)))
myha<- HeatmapAnnotation(df = metause, col = list(Induction = induction.cols.assigned, 
                                                  Treatment = treatment.cols.assigned,
                                                  Replicate = replicate.cols.assigned))

group.cols.assigned<-setNames(group.cols,c(unique(as.character(genelistmeta$Group))[2:11], 
                                           "Cadherin","Metal response","Molecular transduction",
                                           "Immune response","Stress response"))
U6475.cols.assigned<-setNames(U6475.cols,unique(as.character(genelistmeta$U6475)))
I6475.cols.assigned<-setNames(I6475.cols,unique(as.character(genelistmeta$I6475)))
I17938.cols.assigned<-setNames(I17938.cols,unique(as.character(genelistmeta$I17938)))
I6475I17938.cols.assigned<-setNames(I6475I17938.cols,c(0,1,-1))
I6475U6475.cols.assigned<-setNames(I6475U6475.cols,c(0,1,-1))
IU.cols.assigned<-setNames(IU.cols, c(0,1,-1))

geneha<-rowAnnotation(df=genelistmeta,col=list(Group = group.cols.assigned, 
                                                U6475 = U6475.cols.assigned,
                                               I6475 =I6475.cols.assigned,
                                               I17938 = I17938.cols.assigned,
                                               I6475I17938 = I6475I17938.cols.assigned,
                                               I6475U6475 = I6475U6475.cols.assigned,
                                               IU = IU.cols.assigned))

bigheatmap<-Heatmap(myX, name = "rlog RNAseq\ncounts scaled\n & centered", 
        show_row_names = FALSE, show_column_names = TRUE, 
        row_dend_reorder = TRUE, column_dend_reorder = TRUE, 
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        top_annotation = myha,
        left_annotation = geneha,
        row_km=6,column_km = 5)

bigheatimage<-draw(bigheatmap)

bigheatimage

geneclusterorder<-row_order(bigheatimage)
Genes<-data.frame("EnsemblID" = heatmapdata$EnsemblID, "Gene" = heatmapdata$GeneSymbol, "GeneDescription"=heatmapdata$GeneDescription, "MiddleGroup"=heatmapdata$MiddleGroup,
                  "FinalGroup"=heatmapdata$FinalGroup,"Subtype"=heatmapdata$Subtype,"Cluster" = NA)

Genes$Cluster[geneclusterorder$`1`]<-"1"
Genes$Cluster[geneclusterorder$`2`]<-"2"
Genes$Cluster[geneclusterorder$`3`]<-"3"
Genes$Cluster[geneclusterorder$`4`]<-"4"
Genes$Cluster[geneclusterorder$`5`]<-"5"
Genes$Cluster[geneclusterorder$`6`]<-"6"



kmeans(myX,6)



#do correlation among samples:
#Start the same as for the big heatmap:
#if you use t(myX), will get cor map for genes, not samples
#correlation is pearson 
#adjust rowannotations:


#repeat with all genes, minus the lowly expressed ones:
#need rlog values of all genes

rlogallgenes<-read.delim(file="FigS3ABconvertedrlogcounts.txt",header=TRUE,sep="\t")
#remove the lower quartile of expressed genes?
maxGene <- apply(rlogallgenes,1,max)
# remove bottom 25% lowly expressed genes, which inflate the PPC
rlogtopgenes <- rlogallgenes[which(maxGene > quantile(maxGene)[1] ) ,] 


mybatch <- c("A","A","A","B","B","B","A","A","A","B","B","B","A","A","A","B","B","B",
             "A","A","A","B","B","B","A","A","A","B","B","B")


batchedgenematrixrlog<-removeBatchEffect(rlogtopgenes,batch = mybatch,design = mydesign)
myXall<- t(scale(t(batchedgenematrixrlog),center=TRUE,scale=TRUE))

allgenecormatrix<-cor((myXall))
allgenecormatrix[lower.tri(allgenecormatrix, diag=TRUE)]<-NA
correlationmatrix<-allgenecormatrix

ULDMsubindexes<-str_which(rownames(correlationmatrix), "ULDM")
U6475subindexes<-str_which(rownames(correlationmatrix),"U6475")
ILDMsubindexes<-str_which(rownames(correlationmatrix),"ILDM")
I6475subindexes<-str_which(rownames(correlationmatrix),"I6475")
I17938subindexes<-str_which(rownames(correlationmatrix),"I1793")

ULDMsubvalues<-unlist(as.list(correlationmatrix[ULDMsubindexes,ULDMsubindexes]))
U6475subvalues<-unlist(as.list(correlationmatrix[U6475subindexes,U6475subindexes]))
ILDMsubvalues<-unlist(as.list(correlationmatrix[ILDMsubindexes,ILDMsubindexes]))
I6475subvalues<-unlist(as.list(correlationmatrix[I6475subindexes,I6475subindexes]))
I17938subvalues<-unlist(as.list(correlationmatrix[I17938subindexes,I17938subindexes]))
ULDMU6475subvalues<-unlist(as.list(correlationmatrix[ULDMsubindexes,U6475subindexes]))
ILDMI6475subvalues<-unlist(as.list(correlationmatrix[ILDMsubindexes,I6475subindexes]))
ILDMI17938subvalues<-unlist(as.list(correlationmatrix[ILDMsubindexes,I17938subindexes]))
I6475I17938subvalues<-unlist(as.list(correlationmatrix[I6475subindexes,I17938subindexes]))


corredata<-data.frame("ULDM4-ULDM4" = ULDMsubvalues,"U6475-U6475" = U6475subvalues,
                      "ILDM4-ILDM4" = ILDMsubvalues,"I6475-I6475" = I6475subvalues,
                      "I17938-I17938" =I17938subvalues,
                      "ULDM-U6475" = ULDMU6475subvalues, "ILDM-I6475" = ILDMI6475subvalues,
                      "ILDM-I17938" = ILDMI17938subvalues, "I6475-I17938" = I6475I17938subvalues)
#write.table(corredata,file="correlationvalues.txt",sep="\t")
melt1<-reshape2::melt(corredata, variable.name="Comparison", value.name="Correlation")

boxplotcorrelations<-ggboxplot(melt1,x="Comparison",y="Correlation",
                               add="jitter")
boxplotcorrelations

library("ggpubr")
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

statdata<-compare_means(Correlation~Comparison, melt1, method = "t.test", p.adjust.method = "holm")   #using holm to reduce false positives
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

#this file shows the mean differences of the comparisons
cormapmeandiff<-as.matrix(read.delim(file="FigS3Bmeandifferences.txt",sep="\t",header=TRUE,row.names = 1))

#this file is a matrix of which comparisons had significant mean differences
cormapstats<-as.matrix(read.delim(file="FigS3Bcorrelationplotstats.txt",sep="\t",header=TRUE,row.names=1))


#change color scheme
#add mean values
col_fun= circlize::colorRamp2(c(0, -1,-1.3, -6,-16), c("navyblue","cornflowerblue","white", "lightcoral","maroon")) #redblue
col_fun= circlize::colorRamp2(breaks=c(0, -1.4,-2,-4), c("navyblue","white","lightcoral","maroon")) #redblue


circleheat<-Heatmap(cormapstats, name = "p.adj", rect_gp = gpar(type = "none"),
        col = col_fun,
        show_row_names = TRUE, show_column_names = TRUE, 
        row_dend_reorder = FALSE, column_dend_reorder = FALSE, 
        cluster_rows =FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.circle(x = x, y = y, r = abs(cormapmeandiff[i, j]) * min(unit.c(width, height)), 
                      gp = gpar(fill = col_fun(log10(cormapstats[i, j])), col = NA))
          grid.text(format.pval(cormapstats[i, j],digits = max(1, getOption("digits") - 3)), x, y, gp = gpar(fontsize = 10))
        })

circleheat

#only need the upper right triangle of the plot
ggarrange(boxplotcorrelations,circleheat)

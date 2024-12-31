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


#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html#plot-the-heamtap-list
#https://www.nature.com/articles/nbt1205-1499


genelist<-read.delim(file="FigS4genes.txt",header = TRUE,sep="\t", na.strings = "NA")
suptable2<-read.delim(file = "FigS4supplementaltable.txt",header = TRUE,sep = "\t")
metadata<-read.delim(file="FigS4metadata.txt",header=TRUE,sep="\t")
fcdata<-read.delim(file="FigS4completeDESeq2output.txt",header=TRUE,sep="\t")
fcs22<-data.frame("Gene" = fcdata$Gene, "U6475" = fcdata$U6475.ULDM4log2FoldChange, "I6475" = fcdata$I6475.ILDM4log2FoldChange,
                  "I17938" = fcdata$I17938.ILDM4log2FoldChange, "I647517938" = fcdata$I6475.I17938log2FoldChange ,
                  "I6475U6475" = fcdata$I6475.U6475log2FoldChange, "IU" = fcdata$ILDM4.ULDM4log2FoldChange)


#get expression values for each sample
heatmapdata<-merge(genelist,suptable2,by.x="EnsemblID",by.y="EnsemblID",x.all=TRUE,y.all=FALSE)
heatmapdat<-merge(heatmapdata,fcs22,by.x="EnsemblID",by.y="Gene",x.all=TRUE,y.all=FALSE)
heatmapdata<-heatmapdat


heatmapdatamelt<-reshape2::melt(heatmapdata,id.vars=names(heatmapdata[c(1:14,56:61) ]),measure.vars =names(heatmapdata[15:44]), variable.name = "Sample", value.name ="Expression" )

genematrix<-cbind(heatmapdata[,c(15:44)])

rownames(genematrix)<-heatmapdata$Gene#[c(1:30)]

mybatch <- c("A","A","A","B","B","B","A","A","A","B","B","B","A","A","A","B","B","B",
           "A","A","A","B","B","B","A","A","A","B","B","B")




mydesign<-model.matrix(~TreatmentInduction,metadata)

batchedgenematrix<-removeBatchEffect(genematrix,batch = mybatch,design = mydesign)

mymatrix<-batchedgenematrix  #genematrix is the uncorrected matrix

#for the heatmap plot, rescale the expression values
myX<- t(scale(t(mymatrix),center=TRUE,scale=TRUE)) #scale works on the columns; here the columns are samples and the heatmap will have the samples as columns
#so to compared across metabolites, transpose, do scale, then transpose back

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

U6475.cols<-c("navyblue","white","maroon")
I6475.cols<-c("white", "maroon","navyblue")
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

group.cols.assigned<-setNames(group.cols,c(unique(as.character(genelistmeta$Group)), 
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
        row_km=8,column_km = 5)

bigheatimage<-draw(bigheatmap) #Supplemental Fig 4



geneclusterorder<-row_order(bigheatimage)
Genes<-data.frame("EnsemblID" = heatmapdata$EnsemblID, "Gene" = heatmapdata$GeneSymbol, "GeneDescription"=heatmapdata$GeneDescription, "MiddleGroup"=heatmapdata$MiddleGroup,
                  "FinalGroup"=heatmapdata$FinalGroup,"Subtype"=heatmapdata$Subtype,"Cluster" = NA)

Genes$Cluster[geneclusterorder$`1`]<-"1"
Genes$Cluster[geneclusterorder$`2`]<-"2"
Genes$Cluster[geneclusterorder$`3`]<-"3"
Genes$Cluster[geneclusterorder$`4`]<-"4"
Genes$Cluster[geneclusterorder$`5`]<-"5"
Genes$Cluster[geneclusterorder$`6`]<-"6"
Genes$Cluster[geneclusterorder$`7`]<-"7"
Genes$Cluster[geneclusterorder$`8`]<-"8"



kmeans(myX,8)


#if you want to plot smaller heatmaps by cluster
#break into separate plots based on cluster
clusters<-levels(as.factor(heatmapdata$Cluster))
counter<-c(seq(1,length(clusters),1))

#then break into subgroups with gene names printed
#sub on Final Group,
subtype.cols<-c("brown","black","lavender","royalblue3","palegreen1","yellow3","antiquewhite","burlywood")
subtype.cols.assigned<-setNames(subtype.cols,unique(as.character(heatmapdata$Subtype)))
U6475.cols<-c("navyblue","white","maroon")
I6475.cols<-c("white", "maroon","navyblue")
I17938.cols<-c("white", "navyblue","maroon")
induction.cols.assigned<- setNames(induction.cols, unique(as.character(metadata$Induction)))
treatment.cols.assigned<- setNames(treatment.cols, unique(as.character(metadata$Treatment)))
replicate.cols.assigned<- setNames(replicate.cols,unique(as.character(metadata$Replicate)))
myha<- HeatmapAnnotation(df = metause, col = list(Induction = induction.cols.assigned, 
                                                  Treatment = treatment.cols.assigned,
                                                  Replicate = replicate.cols.assigned))

U6475.cols.assigned<-setNames(U6475.cols,unique(as.character(genelistmeta$U6475)))
I6475.cols.assigned<-setNames(I6475.cols,unique(as.character(genelistmeta$I6475)))
I17938.cols.assigned<-setNames(I17938.cols,unique(as.character(genelistmeta$I17938)))


i<-2
for(i in 1:(length(counter))) {
index<-which(heatmapdata$Cluster==clusters[i])
heatmapdatameltsub<-reshape2::melt(heatmapdata[index,],id.vars=names(heatmapdata[1:14]),measure.vars =names(heatmapdata[15:44]), variable.name = "Sample", value.name ="Expression" )
genematrixsub<-cbind(heatmapdata[index,c(15:44)])
rownames(genematrixsub)<-heatmapdata$Gene[index]

batchedgenematrixsub<-removeBatchEffect(genematrixsub,batch = mybatch,design = mydesign)

myXsub<-batchedgenematrixsub

myXsub<- t(scale(t(genematrixsub),center=TRUE,scale=TRUE))


genelistmetsub<-data.frame(Group = heatmapdata$FinalGroup[index], #[c(1:30)]
                           Subtype = heatmapdata$Subtype[index],
                         U6475=heatmapdata$U6475.ULDM[index], #[c(1:30)]
                         I6475=heatmapdata$I6475.ILDM[index], #[c(1:30)]
                         I17938=heatmapdata$I1793.ILDM[index],
                         I6475I17938 = heatmapdata$I6475.I1793[index],
                         I6475U6475 = heatmapdata$I6475.U6475[index],
                         IU = heatmapdata$ILDM.ULDM[index]) #[c(1:30)]

genelistmetasub<-droplevels(genelistmetsub)

induction.cols<-c("blue", "red")
treatment.cols<- c("gold","purple2","orchid1")
replicate.cols<-c("gray","black")


genehasub<-rowAnnotation(df=genelistmetasub,col=list(Group = group.cols.assigned, 
                                               Subtype =subtype.cols.assigned,
                                               U6475 = U6475.cols.assigned,
                                               I6475 =I6475.cols.assigned,
                                               I17938 = I17938.cols.assigned,
                                               I6475I17938 = I6475I17938.cols.assigned,
                                               I6475U6475 = I6475U6475.cols.assigned,
                                               IU = IU.cols.assigned))

groupname<-clusters[i]
groupname<-gsub("/","",groupname)

myplot<-Heatmap(myXsub, name = "rlog RNAseq\ncounts scaled", 
          show_row_names = TRUE, show_column_names = TRUE, 
          row_dend_reorder = TRUE, column_dend_reorder = TRUE, 
          clustering_distance_rows = "pearson",
          clustering_distance_columns = "euclidean",
          clustering_method_rows = "complete",
          clustering_method_columns = "complete",
          top_annotation = myha,
          left_annotation = genehasub)


pdf(file=paste("batchcorrectedheatmap8_cluster",groupname,".pdf",sep=""),width=12, height=8)
print(myplot)
dev.off()
}






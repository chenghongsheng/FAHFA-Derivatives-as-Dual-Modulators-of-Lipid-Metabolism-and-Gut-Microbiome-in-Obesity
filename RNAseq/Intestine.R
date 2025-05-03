setwd("path/to/your/working/directory")

library(DESeq2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(pheatmap)
library(ViSEAGO)
library(EnhancedVolcano)
library(stringr)
library(ggfortify)
library(ggforce)
library(GSVA)
library(rstatix)

##### load data ######

metadata.sub <- read.csv("./metadata.csv", row.names=1)
metadata.sub$Group<-factor(metadata$Group,levels=c(
  'Control', 'LIDPAD', 'Sema','12SAHSA','12TAHSA','12SAASA'  ,'12TAASA', '12HDTZSA'))


count.sub <- read.table("./Fcount.txt", row.names=1, header = T)
  
#sanity check

all(rownames(metadata.sub) == colnames(count.sub))

#### set color ######

color<-c('black','red',brewer.pal(6,'Set2'))
pie(rep(1,8),col=color)

##### DESeq-first trial #####
#conducts normalization
dds <- DESeqDataSetFromMatrix(countData = count.sub,
                                    colData = metadata.sub,
                                    design = ~ Group)
dds<-DESeq(dds)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

vsd

pcadat<-plotPCA(vsd,intgroup="Group",ntop=3000,returnData=T)

all(rownames(pcadat)==rownames(metadata.sub)) #sanity check
pcadat$Group<-metadata.sub$Group

percentVar.vsd<-round(100*attr(pcadat,"percentVar"))

pcadat$name<-str_split_i(pcadat$name,'_',1)
#pcadat$name<-str_split_i(pcadat$name,'\\.',2)

pcadat<-subset(pcadat,!pcadat$Group %in% c('Control','Sema'))

color2<-color[c(2,4:8)]

pca<-ggplot(pcadat, aes(PC1, PC2, fill=Group)) +
  geom_point(size=5,pch=21,stroke=0.5)+ 
  #geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd[2],"% variance")) + 
  #stat_ellipse()+
  geom_mark_ellipse(aes(fill = Group,color = Group))+
  theme_light()+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = color2)+
  scale_color_manual(values = color2)+
  coord_fixed()

dev.new()
pdf("pca.FAHFA.pdf")
pca
dev.off()

#write.csv(pcadat,'pcadat.csv')

###### DEGs Liver ####
Mmu.dataset<-useDataset('mmusculus_gene_ensembl',mart=useMart("ensembl"))
Genemap<-getBM(attributes = c('ensembl_gene_id','external_gene_name',"gene_biotype"), 
               filters='ensembl_gene_id',
               values=rownames(count),mart=Mmu.dataset)
genesymbols <- tapply(Genemap$external_gene_name, 
                      Genemap$ensembl_gene_id, paste, collapse="; ")

# LIDPAD vs control
res.LPvCtrl <- as.data.frame(results(dds, contrast=c('Group','LIDPAD','Control'),alpha=0.05)) 
res.LPvCtrl<-res.LPvCtrl[order(res.LPvCtrl$padj),]
res.LPvCtrl$symbol<-genesymbols[rownames(res.LPvCtrl)]
sig.LPvCtrl<-res.LPvCtrl[which(res.LPvCtrl$pvalue<0.01 &
                                       (abs(res.LPvCtrl$log2FoldChange)>1)),]

write.csv(sig.LPvCtrl,"sig.LPvCtrl.csv")
write.csv(res.LPvCtrl,"res.LPvCtrl.csv")

# Sema vs LIDPAD
res.semavLP <- as.data.frame(results(dds, contrast=c('Group','Sema','LIDPAD'),alpha=0.05)) 
res.semavLP<-res.semavLP[order(res.semavLP$padj),]
res.semavLP$symbol<-genesymbols[rownames(res.semavLP)]
sig.semavLP<-res.semavLP[which(res.semavLP$pvalue<0.01 &
                                 (abs(res.semavLP$log2FoldChange)>1)),]

write.csv(sig.semavLP,"sig.semavLP.csv")
write.csv(res.semavLP,"res.semavLP.csv")


# 12SAHSA vs LIDPAD
res.12SAHSAvLP <- as.data.frame(results(dds, contrast=c('Group','12SAHSA','LIDPAD'),alpha=0.05)) 
res.12SAHSAvLP<-res.12SAHSAvLP[order(res.12SAHSAvLP$padj),]
res.12SAHSAvLP$symbol<-genesymbols[rownames(res.12SAHSAvLP)]
sig.12SAHSAvLP<-res.12SAHSAvLP[which(res.12SAHSAvLP$pvalue<0.01 &
                                 (abs(res.12SAHSAvLP$log2FoldChange)>1)),]

write.csv(sig.12SAHSAvLP,"sig.12SAHSAvLP.csv")
write.csv(res.12SAHSAvLP,"res.12SAHSAvLP.csv")


# 12TAHSA vs LIDPAD
res.12TAHSAvLP <- as.data.frame(results(dds, contrast=c('Group','12TAHSA','LIDPAD'),alpha=0.05)) 
res.12TAHSAvLP<-res.12TAHSAvLP[order(res.12TAHSAvLP$padj),]
res.12TAHSAvLP$symbol<-genesymbols[rownames(res.12TAHSAvLP)]
sig.12TAHSAvLP<-res.12TAHSAvLP[which(res.12TAHSAvLP$pvalue<0.01 &
                                       (abs(res.12TAHSAvLP$log2FoldChange)>1)),]

write.csv(sig.12TAHSAvLP,"sig.12TAHSAvLP.csv")
write.csv(res.12TAHSAvLP,"res.12TAHSAvLP.csv")


# 12SAASA vs LIDPAD
res.12SAASAvLP <- as.data.frame(results(dds, contrast=c('Group','12SAASA','LIDPAD'),alpha=0.05)) 
res.12SAASAvLP<-res.12SAASAvLP[order(res.12SAASAvLP$padj),]
res.12SAASAvLP$symbol<-genesymbols[rownames(res.12SAASAvLP)]
sig.12SAASAvLP<-res.12SAASAvLP[which(res.12SAASAvLP$pvalue<0.01 &
                                       (abs(res.12SAASAvLP$log2FoldChange)>1)),]

write.csv(sig.12SAASAvLP,"sig.12SAASAvLP.csv")
write.csv(res.12SAASAvLP,"res.12SAASAvLP.csv")


# 12TAASA vs LIDPAD
res.12TAASAvLP <- as.data.frame(results(dds, contrast=c('Group','12TAASA','LIDPAD'),alpha=0.05)) 
res.12TAASAvLP<-res.12TAASAvLP[order(res.12TAASAvLP$padj),]
res.12TAASAvLP$symbol<-genesymbols[rownames(res.12TAASAvLP)]
sig.12TAASAvLP<-res.12TAASAvLP[which(res.12TAASAvLP$pvalue<0.01 &
                                       (abs(res.12TAASAvLP$log2FoldChange)>1)),]

write.csv(sig.12TAASAvLP,"sig.12TAASAvLP.csv")
write.csv(res.12TAASAvLP,"res.12TAASAvLP.csv")


# 12HDTZSA vs LIDPAD
res.12HDTZSAvLP <- as.data.frame(results(dds, contrast=c('Group','12HDTZSA','LIDPAD'),alpha=0.05)) 
res.12HDTZSAvLP<-res.12HDTZSAvLP[order(res.12HDTZSAvLP$padj),]
res.12HDTZSAvLP$symbol<-genesymbols[rownames(res.12HDTZSAvLP)]
sig.12HDTZSAvLP<-res.12HDTZSAvLP[which(res.12HDTZSAvLP$pvalue<0.01 &
                                       (abs(res.12HDTZSAvLP$log2FoldChange)>1)),]

write.csv(sig.12HDTZSAvLP,"sig.12HDTZSAvLP.csv")
write.csv(res.12HDTZSAvLP,"res.12HDTZSAvLP.csv")


#plot heatmap ####
DEG.all<-c(rownames(sig.LPvCtrl),
                 rownames(sig.semavLP),
           rownames(sig.12SAHSAvLP),
           rownames(sig.12TAHSAvLP),
           rownames(sig.12SAASAvLP),
           rownames(sig.12TAASAvLP),
           rownames(sig.12HDTZSAvLP)
                 )
DEG.all<-DEG.all[!duplicated(DEG.all)]

hm_mat<-assay(vsd)
hm_mat<-subset(hm_mat,rownames(hm_mat) %in% DEG.all)

coldat_hm<-metadata.sub[order(metadata.sub$Group),,drop=F]

hm_mat<-hm_mat[,rownames(coldat_hm)]
hm_mat<-as.data.frame(hm_mat)

hm_mat<-as.data.frame(t(hm_mat))
hm_mat$group<-coldat_hm$Group

hm_mat_long <- hm_mat %>%
  tidyr::gather('gene','normalized.count',1:1144) %>%
  group_by(group,gene)%>%
  summarise(average=mean(normalized.count)) %>%
  tidyr::spread(gene,average)

hm_mat_long<-as.data.frame(hm_mat_long)
rownames(hm_mat_long)<-hm_mat_long$group
hm_mat_long$group<-NULL
hm_mat_long<-as.data.frame(t(hm_mat_long))


coldat_hm.used<-as.data.frame(colnames(hm_mat_long))
coldat_hm.used$Group<-coldat_hm.used$`colnames(hm_mat_long)`
rownames(coldat_hm.used)<-coldat_hm.used$`colnames(hm_mat_long)`
coldat_hm.used$`colnames(hm_mat_long)`<-NULL

hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
break_hm = seq(-2, 2,length.out=100)

my_color_annotation<-list(Group= c('Control'='black','LIDPAD'='red',
                                   'Sema'='#66C2A5',
                                   '12SAHSA'='#FC8D62',
                                   '12TAHSA'='#8DA0CB',
                                   '12SAASA'='#E78AC3',
                                   '12TAASA'='#A6D854',
                                   '12HDTZSA'='#FFD92F'   
                                   ))



pheatmap(hm_mat_long,scale="row",border_color = NA,color = hm_color,
                     show_rownames = F,show_colnames = F,
                     cluster_rows = T,cluster_cols =F,
                     annotation_col = coldat_hm.used,
                     breaks = break_hm,
                     annotation_colors = my_color_annotation,
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     angle_col = 45,
        # cutree_rows = 5,
        # gaps_col = 10,
                     fontsize_col = 5) 



#plot heatmap- down DEGs in TAASA ####

down.DEG.TAASA <- sig.12TAASAvLP %>%
  subset(log2FoldChange<0) %>%
  rownames(.)

write.csv(down.DEG.TAASA,'down.DEG.TAASA.csv')

down.DEG.SAHSA <- res.12SAHSAvLP %>%
  subset(log2FoldChange<0 & pvalue<0.05) %>%
  rownames(.)


write.csv(down.DEG.SAHSA,'down.DEG.SAHSA.csv')

## selected genes
sel.lipid.gene <- read.csv("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/FAHFA/Animal experiment/RNAseq/analysis/hm_mat_long.sub.csv", row.names=1)
lipid_genes <- read.csv("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/FAHFA/Animal experiment/RNAseq/analysis/lipid_genes.csv")



hm_mat_long.sub<-hm_mat_long %>%
  mutate(symbol = genesymbols[rownames(.)]) %>%
  subset(symbol %in% rownames(sel.lipid.gene))

rownames(hm_mat_long.sub)<-hm_mat_long.sub$symbol
hm_mat_long.sub$symbol <- NULL


rowdat<-lipid_genes %>%
  subset(Gene %in% rownames(sel.lipid.gene)) %>%
  mutate(pos=rep('Yes', nrow(.))) %>%
  tidyr::spread(Function, pos)



rownames(rowdat) <- rowdat$Gene
rowdat$Gene<-NULL
rowdat<-rowdat[,c(7,6,9,2,5,4,8,3,1)]

#write.csv(rowdat,'rowdat.csv')
#rowdat <- read.csv("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/FAHFA/Animal experiment/RNAseq/analysis/rowdat.csv", row.names=1)

#write.csv(hm_mat_long %>%
#            mutate(symbol = genesymbols[rownames(.)]) %>%
#            subset(symbol %in% lipid_genes$Gene),'hm_mat_long.sub.csv')

hm_mat_long.sub<-hm_mat_long.sub[rownames(rowdat),]

p_heatmap_lipid<-pheatmap(hm_mat_long.sub,scale="row",border_color = 'white',color = hm_color,
         show_rownames = T,show_colnames = F,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm.used,
         annotation_row = rowdat,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
        clustering_distance_cols = "correlation",
        # angle_col = 0,
        cellwidth = 10,cellheight = 10,
         # cutree_rows = 5,
         # gaps_col = 10,
         fontsize_col = 1,annotation_names_row=T) 


dev.new()
pdf("p_heatmap_lipid.pdf",height = 10)
p_heatmap_lipid
dev.off()

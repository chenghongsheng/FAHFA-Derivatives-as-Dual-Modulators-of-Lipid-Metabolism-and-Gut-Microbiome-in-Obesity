setwd("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/FAHFA/Animal experiment/16s/16S_B2")

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(vegan)
library(pairwiseAdonis)
library(ALDEx2)
library(ggpubr)
library(rstatix)
library(emmeans)
library(stringr)
library(pheatmap)
library(RColorBrewer)

library(phyloseq)
library(Biostrings)
library(dplyr)

#load sample from batch 1 and batch 2 and combine ####

#create metadata####

metadata <- read.delim("./metadata.txt", row.names=1)
metadata<-sample_data(metadata)

#taxa
taxa.b1 <- read.csv("./taxa.b1.csv", row.names=1)
taxa.b1$seq<-rownames(taxa.b1)

taxa.b2 <- read.csv("./taxa.b2.csv", row.names=1)
taxa.b2$seq<-rownames(taxa.b2)

taxa.merge<-dplyr::full_join(taxa.b1,taxa.b2,by='seq')
write.csv(taxa.merge,'taxa.merge.csv')

taxa.merge <- read.csv("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/FAHFA/Animal experiment/16s/16S_B2/taxa.merge.fix.csv", row.names=1)
taxa.merge[taxa.merge=='']<-NA

#seqtab 
seqtab.nochim.b1<-readRDS("seqtab.nochim.b1.rds") %>%
  t(.) %>%
  data.frame(.) %>%
  mutate(seq=rownames(.))

seqtab.nochim.b2<-readRDS("seqtab.nochim.b2.rds") %>%
  t(.) %>%
  data.frame(.) %>%
  mutate(seq=rownames(.))


seqtab.nochim.merged<-dplyr::full_join(seqtab.nochim.b1,seqtab.nochim.b2,by='seq')
seqtab.nochim.merged[is.na(seqtab.nochim.merged)]<-0
rownames(seqtab.nochim.merged) <- seqtab.nochim.merged$seq
seqtab.nochim.merged$seq<-NULL

length(intersect(rownames(seqtab.nochim.merged),rownames(taxa.merge)))
seqtab.nochim.merged<-seqtab.nochim.merged[rownames(taxa.merge),]
colnames(seqtab.nochim.merged)<-str_replace(colnames(seqtab.nochim.merged),pattern='\\.','-')

seqtab.nochim.merged<-as.data.frame(t(seqtab.nochim.merged))
rownames(seqtab.nochim.merged)<-str_replace(rownames(seqtab.nochim.merged),pattern='\\.','-')


################# set unique ID ########
giveUniqueID<-function(x){
  x<-x %>%
    mutate(ID = case_when(Genus=='na'~paste0('Family_',Family),
                          TRUE ~ Genus))
  x<-x %>%
    mutate(ID = case_when(ID=='Family_na'~paste0('Order_',Order),
                          TRUE ~ ID))
  x<-x %>%
    mutate(ID = case_when(ID=='Order_na'~paste0('Class_',Class),
                          TRUE ~ ID))
  
  x<-x %>%
    mutate(ID = case_when(ID=='Class_na'~paste0('Phylum_',Phylum),
                          TRUE ~ ID))
  
  x<-x %>%
    mutate(ID = case_when(ID=='Phylum_na'~paste0('Kingdom_',Kingdom),
                          TRUE ~ ID))
}

phyloseq_counts_to_df <- function(ps) {
  df <- as.data.frame(t(otu_table(ps)))
  df <- cbind("#OTU ID" = rownames(df), df)  # Cheeky way to add a comment to first line
  
  # Note: we use OTU ID above since at the moment this function is meanly used
  # to pipe files into picrust, which needs the 'OTU ID'.
  return(df)
}

write_counts_to_file <- function(ps, filename) {
  df <- phyloseq_counts_to_df(ps)
  
  write.table(df,
              file=filename,
              sep = "\t",
              row.names = FALSE)
}

write_seqs_to_file <- function(ps, filename) {
  df <- as.data.frame(refseq(ps))
  
  unlink(filename)  # make sure we delete before we concatenate below
  
  names <- rownames(df)
  for (i in 1:nrow(df)) {
    cat(paste0(">", names[i]), file=filename, sep="\n", append=TRUE)
    cat(df[i,1], file=filename, sep="\n", append=TRUE)
  }
  message(paste0("Wrote ", nrow(df), " ASV sequences to file ", filename))
}

# create ps data
ps <- phyloseq(otu_table(seqtab.nochim.merged, taxa_are_rows=FALSE), 
               tax_table(as.matrix(taxa.merge)))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, metadata,dna)

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps_percent <- transform_sample_counts(ps, function(x) 100*x/sum(x))

saveRDS(ps_percent,file="./ps.merge_percent.RDS")
saveRDS(ps,file="./ps.merge.RDS")


#######load data############

ps.genus<-tax_glom(ps,taxrank = 'Genus',NArm = F)

ps.genus_percent <- transform_sample_counts(ps.genus, function(x) 100*x/sum(x))

saveRDS(ps.genus_percent,file="./ps.merge.genus_percent.RDS")
saveRDS(ps.genus,file="./ps.merge.genus.RDS")

bar_Genus <- plot_bar(ps.genus, fill = "Genus")
relabundance<-as.data.frame(bar_Genus$data)
relabundance<-subset(relabundance,relabundance$Kingdom=="Bacteria" & !relabundance$Abundance==0)

relabundance<-relabundance %>% #collapse genus
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  summarise(percent=sum(Abundance))

relabundance<-spread(relabundance,'Sample','percent')

relabundance[,c(1:6)][is.na(relabundance[,c(1:6)])]<-'na'

relabundance<-giveUniqueID(relabundance)
write.csv(relabundance,'relabundance.merge.csv')


coldat<-as.data.frame(bar_Genus$data) %>%
  select(Sample,Group)
coldat<-coldat[!duplicated(coldat$Sample),]


coldat$Group<-factor(coldat$Group,levels=c('Control', 'LIDPAD', 'Sema','12SAHSA','12TAHSA','12SAASA'  ,'12-TAASA', '12-HDTZSA'))
coldat<-coldat[order(coldat$Group),]


relabundance<-gather(relabundance,"Sample",'percent',7:30)
relabundance<-full_join(relabundance,coldat,by='Sample')

n <- length(levels(factor(relabundance$ID)))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(888)
col<-sample(col_vector, n,replace = T)
pie(rep(1,n), col=col)

relabundance<-arrange(relabundance,rev(Group),desc(percent))

df.mean <-relabundance %>%
  select(ID,percent) %>%
  group_by(ID)%>%
  mutate(percent=replace_na(percent,0)) %>%
  summarise(mean=mean(percent)) %>%
  arrange(desc(mean))

relabundance$ID<-factor(relabundance$ID,levels =df.mean$ID[!duplicated(df.mean$ID)])


Genus_bar_FMT_used <-ggplot(relabundance,aes(fill=ID,y=percent,x=Sample))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity',show.legend=T,width = 0.9)+
  scale_fill_manual(values=col)+
  guides(fill=guide_legend(reverse=T))+
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "None",
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  xlab("")+
  ylab("Relative Abundance (%)")

dev.new()
pdf('Genus_bar_FMT_used.pdf',width = 20)
Genus_bar_FMT_used
dev.off()

#Top20
relabundance_long<-relabundance %>%
  select(Genus,Family,Order,Class,Phylum,Kingdom,Sample,percent,ID) %>%
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  spread('Sample','percent')

colnames(relabundance_long)

relabundance_long[,c(8:31)][is.na(relabundance_long[,c(8:31)])]<-0
relabundance_long$sum<-rowSums(as.matrix(relabundance_long[,c(8:31)]))
relabundance_long<-relabundance_long %>%
  arrange(desc(sum))

Top20 = relabundance_long$ID[1:20]

relabundance.top20<-subset(relabundance,relabundance$ID  %in% Top20)

Genus_bar_FMT.top20 <-ggplot(relabundance.top20,aes(fill=ID,y=percent,x=Sample))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity',show.legend=T)+
  scale_fill_manual(values=col)+
  guides(fill=guide_legend(reverse=T))+
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "None",
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  xlab("")+
  ylab("Relative Abundance (%)")

dev.new()
pdf('Genus_bar_FMT.top20.pdf',width = 15)
Genus_bar_FMT.top20
dev.off()

### set color #####

color<-c('black','red',brewer.pal(6,'Set2'))
pie(rep(1,8),col=color)


#alpha diversity############
#alpha-diversity

alpha<-plot_richness(ps.genus) #do not use Chao1 index from DADA2
diversity_index<-as.data.frame(alpha$data)

diversity_index$Group <-factor(diversity_index$Group,
                               levels=c('Control', 'LIDPAD', 'Sema','12SAHSA','12TAHSA','12SAASA'  ,'12-TAASA', '12-HDTZSA'))

diversity_index<-diversity_index %>%
  select(-c(se)) %>%
  spread(variable,value)

write.csv(diversity_index,file="diversity_index.csv")

Shannon <- ggplot(diversity_index,aes(x=Group,y=Shannon,color=Group))+
  geom_point(size=3)+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  stat_summary(fun = mean,geom = "crossbar", width = 0.5,color='black')+
  #scale_y_log10()+
  scale_color_manual(values=color)+
  ggtitle("bar is mean")+
  xlab("")+
  ylab("Shannon")

dev.new()
pdf('Shannon.pdf')
Shannon
dev.off()



InvSimpson <- ggplot(diversity_index,aes(x=Group,y=InvSimpson,color=Group))+
  geom_point(size=3)+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  stat_summary(fun = mean,geom = "crossbar", width = 0.5,color='black')+
  scale_y_log10()+
  #scale_color_manual(values=color)+
  ggtitle("bar is mean")+
  xlab("")+
  ylab("InvSimpson")

dev.new()
pdf('InvSimpson.pdf')
InvSimpson
dev.off()

#ANOVA 
Shannon_aov<-anova_test(diversity_index,dv=Shannon, between = Group)
Shannon_aov #sig
#       Effect DFn DFd     F     p p<.05   ges
#1      Group   4  10 4.392 0.026     * 0.637

pwc_Shannon<-pairwise_t_test(diversity_index,Shannon ~ Group)
pwc_Shannon


InvSimpson_aov<-anova_test(diversity_index,dv=InvSimpson, between = Group)
InvSimpson_aov #non-sig

pwc_InvSimpson<-pairwise_t_test(diversity_index,InvSimpson ~ Group)
pwc_InvSimpson


########### beta diversity ################
# Plot PCoA of Sample using bray-curtis dissimilarity 
# need to collapse to genus level 

#ps.genus_percent <- readRDS("./ps.genus_percent.RDS")

ps.genus_percent@sam_data$Group<-factor(ps.genus_percent@sam_data$Group,
                                        levels=c('Control', 'LIDPAD', 'Sema','12SAHSA','12TAHSA','12SAASA'  ,'12-TAASA', '12-HDTZSA'))

ps.genus.ord <-ordinate(ps.genus_percent, "PCoA", "bray")

PCoA <- plot_ordination(ps.genus_percent, ps.genus.ord, 
                        type="samples",
                        color='Group',
                        title="Treatment effect")

all(rownames(PCoA$data)==rownames(ps.genus_percent@sam_data))
PCoA$data$Group<-ps.genus_percent@sam_data$Group


PCoA  <- PCoA  + 
  geom_point(aes(color=Group),size=5)+
  scale_color_manual(values=color)+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1))+
  ggforce::geom_mark_ellipse(aes(fill = Group))+
  scale_fill_manual(values=color)+
  #facet_wrap(~Treatment,1)+
  scale_y_continuous(limits = c(-0.3,0.65))+
  xlab("PCoA 1 (37.7%)")+
  ylab("PCoA 2 (21.5%)")

dev.new()
pdf('PCoA.beta.diversity.pdf')
PCoA 
dev.off()


#Permanova using BC dissimilarity
metadata <- as(sample_data(ps.genus_percent), "data.frame")

dist<-phyloseq::distance(ps.genus_percent, method="bray")

permanova_ps <- adonis2(dist~ Group , 
                        #strata=metadata$Sema_ID,
                        data=metadata, 
                        permutations=9999)
permanova_ps
#adonis2(formula = dist ~ Group, data = metadata, permutations = 9999)
#Df SumOfSqs      R2      F Pr(>F)    
#Group     7   2.1522 0.62211 3.7628  1e-04 ***
#  Residual 16   1.3074 0.37789                  
#Total    23   3.4596 1.00000                    

pairwise.adonis(dist, metadata$Group,p.adjust.m ="BH") #NS


# differential analysis ########
ps.genus<-readRDS('./ps.merge.genus.rds')
ps.genus_percent<-readRDS('./ps.merge.genus_percent.rds')

bar_Genus <- plot_bar(ps.genus, fill = "Genus")


coldat<-as.data.frame(bar_Genus$data) %>%
  select(Sample,Group)
coldat<-coldat[!duplicated(coldat$Sample),]
rownames(coldat)<-coldat$Sample

count<-as.data.frame(bar_Genus$data)
count<-subset(count,count$Kingdom=="Bacteria" & !count$Abundance==0)

count<-count %>% #collapse genus
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  summarise(count=sum(Abundance))

count<-spread(count,'Sample','count')

count[,c(1:6)][is.na(count[,c(1:6)])]<-'na'

count<-giveUniqueID(count)

count<-as.data.frame(count[,c(coldat$Sample,'ID')])
count[is.na(count)]<-0
rownames(count)<-count$ID
count$ID<-NULL


#interpret aldex2 results#
# rab.all is median log2 relative abundance
# diff.btw is log2 fold difference
# diff.win is log2 variance
# effect is the ratio of the two
# we.ep - Expected P value of Welch’s t test
# we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t test
# wi.ep - Expected P value of Wilcoxon rank test
# wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test


## LP vs ctrl ####
coldata_LPvCtrl <- subset(coldat, coldat$Group %in% c("Control","LIDPAD"))
coldata_LPvCtrl$Group<-factor(coldata_LPvCtrl$Group,levels=c('Control','LIDPAD'))
count_LPvCtrl<- count[,rownames(coldata_LPvCtrl)]
aldex_LPvCtrl <- aldex(count_LPvCtrl, coldata_LPvCtrl$Group, mc.samples=1000, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LPvCtrl<-subset(aldex_LPvCtrl,abs(aldex_LPvCtrl$effect)>1 & 
                          abs(aldex_LPvCtrl$diff.btw)>1)

write.csv(aldex_LPvCtrl,'aldex_LPvCtrl.csv' )
write.csv(sig_LPvCtrl,'sig_LPvCtrl.csv' )


## Sema vs LP ####
coldata_SemavLP <- subset(coldat, coldat$Group %in% c("Sema","LIDPAD"))
coldata_SemavLP$Group<-factor(coldata_SemavLP$Group,levels=c('LIDPAD','Sema'))
count_SemavLP<- count[,rownames(coldata_SemavLP)]
aldex_SemavLP <- aldex(count_SemavLP, coldata_SemavLP$Group, mc.samples=1000, test="t", effect=TRUE,
                       include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_SemavLP<-subset(aldex_SemavLP,abs(aldex_SemavLP$effect)>1 & 
                      abs(aldex_SemavLP$diff.btw)>1)

write.csv(aldex_SemavLP,'aldex_SemavLP.csv' )
write.csv(sig_SemavLP,'sig_SemavLP.csv' )


## 12TAASA vs LP ####
coldata_12TAASAvLP <- subset(coldat, coldat$Group %in% c("12-TAASA","LIDPAD"))
coldata_12TAASAvLP$Group<-factor(coldata_12TAASAvLP$Group,levels=c('LIDPAD','12-TAASA'))
count_12TAASAvLP<- count[,rownames(coldata_12TAASAvLP)]
aldex_12TAASAvLP <- aldex(count_12TAASAvLP, coldata_12TAASAvLP$Group, mc.samples=1000, test="t", effect=TRUE,
                       include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_12TAASAvLP<-subset(aldex_12TAASAvLP,abs(aldex_12TAASAvLP$effect)>1 & 
                      abs(aldex_12TAASAvLP$diff.btw)>1)

write.csv(aldex_12TAASAvLP,'aldex_12TAASAvLP.csv' )
write.csv(sig_12TAASAvLP,'sig_12TAASAvLP.csv' )

## 12HDTZSA vs LP ####
coldata_12HDTZSAvLP <- subset(coldat, coldat$Group %in% c("12-HDTZSA","LIDPAD"))
coldata_12HDTZSAvLP$Group<-factor(coldata_12HDTZSAvLP$Group,levels=c('LIDPAD','12-HDTZSA'))
count_12HDTZSAvLP<- count[,rownames(coldata_12HDTZSAvLP)]
aldex_12HDTZSAvLP <- aldex(count_12HDTZSAvLP, coldata_12HDTZSAvLP$Group, mc.samples=1000, test="t", effect=TRUE,
                       include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_12HDTZSAvLP<-subset(aldex_12HDTZSAvLP,abs(aldex_12HDTZSAvLP$effect)>1 & 
                      abs(aldex_12HDTZSAvLP$diff.btw)>1)

write.csv(aldex_12HDTZSAvLP,'aldex_12HDTZSAvLP.csv' )
write.csv(sig_12HDTZSAvLP,'sig_12HDTZSAvLP.csv' )



## 12SAASA vs LP ####
coldata_12SAASAvLP <- subset(coldat, coldat$Group %in% c("12SAASA","LIDPAD"))
coldata_12SAASAvLP$Group<-factor(coldata_12SAASAvLP$Group,levels=c('LIDPAD','12SAASA'))
count_12SAASAvLP<- count[,rownames(coldata_12SAASAvLP)]
aldex_12SAASAvLP <- aldex(count_12SAASAvLP, coldata_12SAASAvLP$Group, mc.samples=1000, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_12SAASAvLP<-subset(aldex_12SAASAvLP,abs(aldex_12SAASAvLP$effect)>1 & 
                          abs(aldex_12SAASAvLP$diff.btw)>1)

write.csv(aldex_12SAASAvLP,'aldex_12SAASAvLP.csv' )
write.csv(sig_12SAASAvLP,'sig_12SAASAvLP.csv' )



## 12SAHSA vs LP ####
coldata_12SAHSAvLP <- subset(coldat, coldat$Group %in% c("12SAHSA","LIDPAD"))
coldata_12SAHSAvLP$Group<-factor(coldata_12SAHSAvLP$Group,levels=c('LIDPAD','12SAHSA'))
count_12SAHSAvLP<- count[,rownames(coldata_12SAHSAvLP)]
aldex_12SAHSAvLP <- aldex(count_12SAHSAvLP, coldata_12SAHSAvLP$Group, mc.samples=1000, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_12SAHSAvLP<-subset(aldex_12SAHSAvLP,abs(aldex_12SAHSAvLP$effect)>1 & 
                         abs(aldex_12SAHSAvLP$diff.btw)>1)

write.csv(aldex_12SAHSAvLP,'aldex_12SAHSAvLP.csv' )
write.csv(sig_12SAHSAvLP,'sig_12SAHSAvLP.csv' )


## 12TAHSA vs LP ####
coldata_12TAHSAvLP <- subset(coldat, coldat$Group %in% c("12TAHSA","LIDPAD"))
coldata_12TAHSAvLP$Group<-factor(coldata_12TAHSAvLP$Group,levels=c('LIDPAD','12TAHSA'))
count_12TAHSAvLP<- count[,rownames(coldata_12TAHSAvLP)]
aldex_12TAHSAvLP <- aldex(count_12TAHSAvLP, coldata_12TAHSAvLP$Group, mc.samples=1000, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_12TAHSAvLP<-subset(aldex_12TAHSAvLP,abs(aldex_12TAHSAvLP$effect)>1 & 
                         abs(aldex_12TAHSAvLP$diff.btw)>1)

write.csv(aldex_12TAHSAvLP,'aldex_12TAHSAvLP.csv' )
write.csv(sig_12TAHSAvLP,'sig_12TAHSAvLP.csv' )


# heatmap using pheatmap #####


Dgenus <- c(rownames(sig_LPvCtrl),
            rownames(sig_SemavLP),
            rownames(sig_12TAASAvLP),
            rownames(sig_12HDTZSAvLP),
            rownames(sig_12SAASAvLP),
            rownames(sig_12SAHSAvLP),
            rownames(sig_12TAHSAvLP)
                        )

Dgenus <- Dgenus[!duplicated(Dgenus)]

bar_Genus <- plot_bar(ps.genus_percent, fill = "Genus")

relabundance<-as.data.frame(bar_Genus$data)
relabundance<-subset(relabundance,relabundance$Kingdom=="Bacteria" & !relabundance$Abundance==0)

relabundance<-relabundance %>% #collapse genus
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  summarise(percent=sum(Abundance))

relabundance<-spread(relabundance,'Sample','percent')

relabundance[,c(1:6)][is.na(relabundance[,c(1:6)])]<-'na'

relabundance<-giveUniqueID(relabundance)

relabundance<-as.data.frame(relabundance[,c(coldat$Sample,'ID')])
relabundance[is.na(relabundance)]<-0
rownames(relabundance)<-relabundance$ID
relabundance$ID<-NULL


rel.Dgenus <-subset(relabundance,rownames(relabundance) %in% Dgenus)

coldat.hm<-coldat[,'Group',drop=F]
coldat.hm$Group<-factor(coldat.hm$Group,levels=c(
  'Control', 'LIDPAD', 'Sema',
  '12SAHSA', '12TAHSA',
   '12SAASA', '12-TAASA','12-HDTZSA' 
))
coldat.hm<-coldat.hm[order(coldat.hm$Group),,drop=F]

rel.Dgenus<-rel.Dgenus[,rownames(coldat.hm)]

hm_color<- colorRampPalette(rev(brewer.pal(11,'RdBu')))(100)

breaks=seq(-1.5,1.5,length.out=100)

my_color_annotation<-list(Group= c('Control'="black",
                                     "LIDPAD"="red",
                                     'Sema'="#66C2A5",
                                   '12SAHSA'='#FC8D62',
                                   '12TAHSA'='#8DA0CB',
                                   '12SAASA'='#E78AC3',
                                   '12-TAASA'='#A6D854',
                                   '12-HDTZSA'='#FFD92F'
                                        ))


Dgenus_hm<-pheatmap(rel.Dgenus,scale="row",border_color = NA,color = hm_color,
                    show_rownames = T,show_colnames = F,
                    cluster_rows = T,cluster_cols = F,
                    #annotation_row = rowdat,
                    breaks = breaks,
                    annotation_col = coldat.hm,
                    annotation_colors = my_color_annotation,
                    clustering_distance_rows = "maximum",
                    cellwidth = 10,cellheight = 10,
                    angle_col = 45,gaps_col = c(9,15))

dev.new()
pdf('Dgenus_hm.pdf', width=10,height=10)
Dgenus_hm
dev.off()

# export for  picrust #####

ps.clean <- subset_taxa(ps.genus, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("mitochondria"))
ps.clean


ps.clean.trim <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= 2}, prune=TRUE)
ps.clean.trim

sample_data(ps.clean.trim)

dir.create("picrust2", showWarnings = FALSE)
write_counts_to_file(ps.clean.trim, filename = "picrust2/raw_counts.tsv")
write_seqs_to_file(ps.clean.trim, filename = "picrust2/seqs.fna")

# plot SCFA ####

bar_Genus <- plot_bar(ps.genus_percent, fill = "Genus")

coldat<-as.data.frame(bar_Genus$data) %>%
  select(Sample,Group)
coldat<-coldat[!duplicated(coldat$Sample),]
coldat$Group<-factor(coldat$Group,levels=c('Control', 'LIDPAD', 'Sema','12SAHSA','12TAHSA','12SAASA'  ,'12-TAASA', '12-HDTZSA'))
coldat<-coldat[order(coldat$Group),]


relabundance<-as.data.frame(bar_Genus$data)
relabundance<-subset(relabundance,relabundance$Kingdom=="Bacteria" & !relabundance$Abundance==0)

relabundance<-relabundance %>% #collapse genus
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  summarise(percent=sum(Abundance))

relabundance<-spread(relabundance,'Sample','percent')

relabundance[,c(1:6)][is.na(relabundance[,c(1:6)])]<-'na'

relabundance<-giveUniqueID(relabundance)

relabundance<-as.data.frame(relabundance[,c(coldat$Sample,'ID')])
relabundance[is.na(relabundance)]<-0
rownames(relabundance)<-relabundance$ID
relabundance$ID<-NULL


## blautia abundance #####
blautia.abundance <-relabundance %>%
  subset(rownames(.)=='Blautia') %>%
  tidyr::gather(Sample,Percent,1:24) %>%
  left_join(.,coldat,by='Sample')

blautia<-ggplot(blautia.abundance,aes(x=Group,y=Percent))+
  stat_summary(fun=mean,geom='bar',fill='white',color='black',linewidth=1)+
  stat_summary(fun.data=mean_se,geom='errorbar',width=0.5,linewidth=1)+
  geom_dotplot(stackdir = 'center',binaxis = 'y',fill='black')+
  theme_light()

dev.new()
pdf("blautia.abundance.pdf")
blautia
dev.off()


pwt.blautia<-pairwise_t_test(blautia.abundance, Percent ~ Group)


## roseburia abundance #####
roseburia.abundance <-relabundance %>%
  subset(rownames(.)=='Roseburia') %>%
  tidyr::gather(Sample,Percent,1:24) %>%
  left_join(.,coldat,by='Sample')

roseburia<-ggplot(roseburia.abundance,aes(x=Group,y=Percent))+
  stat_summary(fun=mean,geom='bar',fill='white',color='black',linewidth=1)+
  stat_summary(fun.data=mean_se,geom='errorbar',width=0.5,linewidth=1)+
  geom_dotplot(stackdir = 'center',binaxis = 'y',fill='black')+
  theme_light()

dev.new()
pdf("roseburia.abundance.pdf")
roseburia
dev.off()


pwt.roseburia<-pairwise_t_test(roseburia.abundance, Percent ~ Group)

## NK4A136 abundance #####
NK4A136.abundance <-relabundance %>%
  subset(rownames(.)=='Lachnospiraceae NK4A136 group') %>%
  tidyr::gather(Sample,Percent,1:24) %>%
  left_join(.,coldat,by='Sample')

NK4A136<-ggplot(NK4A136.abundance,aes(x=Group,y=Percent))+
  stat_summary(fun=mean,geom='bar',fill='white',color='black',linewidth=1)+
  stat_summary(fun.data=mean_se,geom='errorbar',width=0.5,linewidth=1)+
  geom_dotplot(stackdir = 'center',binaxis = 'y',fill='black')+
  theme_light()

dev.new()
pdf("NK4A136.abundance.pdf")
NK4A136
dev.off()

## multiple genera #####

SCFA.abundance<-relabundance %>%
  subset(rownames(.) %in% c('Blautia','Roseburia','Lachnospiraceae NK4A136 group')) %>%
  mutate(Genus = rownames(.)) %>%
  tidyr::gather(Sample,Percent,1:24) %>%
  left_join(.,coldat,by='Sample')

p_SCFA<-ggplot(SCFA.abundance,aes(x=Group,y=Percent))+
  stat_summary(fun=mean,geom='bar',fill='white',color='black',linewidth=1)+
  stat_summary(fun.data=mean_se,geom='errorbar',width=0.5,linewidth=1)+
  geom_dotplot(stackdir = 'center',binaxis = 'y',fill='black')+
  theme_light()+
  facet_grid(Genus~.,scales='free')

dev.new()
pdf("SCFA.abundance.pdf")
p_SCFA
dev.off()


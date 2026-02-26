library(Seurat)
library(BayesSpace)
library(cowplot)
library(SingleCellExperiment)
library(scales)
library(dplyr)
library(reshape2)
library(harmony)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggsci)
library(cowplot)
library(dittoSeq)
library(ggpubr)
library(msigdbr)
library(viridis)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(ggrepel)
library(viridis)
library(gghalves)
library(ggunchained)
library(scales)
library(magrittr)

for (i in 1:9) {
  # i=4
  print(paste("now procession", names(dir)[i]), sep="")
  Autism[[i]] <- NormalizeData(Autism[[i]])
  Autism[[i]] <- FindVariableFeatures(Autism[[i]], selection.method = "vst")
  Autism[[i]] <- ScaleData(Autism[[i]], verbose = FALSE)
  Autism[[i]] <- FindVariableFeatures(Autism[[i]], selection.method = "vst",nfeatures = 2000)
  Autism[[i]] <- RunPCA(Autism[[i]], npcs = 15, verbose = FALSE)
  Autism[[i]] <- RunUMAP(Autism[[i]], reduction = "pca", dims = 1:15)
  Autism[[i]] <- FindNeighbors(Autism[[i]], reduction = "pca", dims = 1:15)
  Autism[[i]] <- FindClusters(Autism[[i]], resolution = 0.1,verbose = T)
  
  pdf(paste("./figures/Seurat/nCount_", names(dir)[i],"_spatial.pdf",sep=""),width = 10,height = 10)
  p<-SpatialFeaturePlot(object = Autism[[i]], 
                        images ='slice1',
                        features = "nCount_Spatial",image.alpha = 0,pt.size.factor = 1.25,
                        crop = F)+scale_fill_continuous(type  = "viridis")+NoLegend()
  print(p)
  dev.off()
  
  pdf(paste("./figures/Seurat/nFeature_", names(dir)[i],"_spatial.pdf",sep=""),width = 10,height = 10)
  p<-SpatialFeaturePlot(object = Autism[[i]], 
                        images ='slice1',
                        features = "nFeature_Spatial",image.alpha = 0,pt.size.factor = 1.25,
                        crop = F)+scale_fill_continuous(type = "viridis")+NoLegend()
  print(p)
  dev.off()
}

ASD.Bys <- list()
ASD.BysE <- list()

set.seed(123)
for (i in for (i in 1:9)) {
  set.seed(123)
  bys <- Autism[[i]]
  bys <- DietSeurat(bys, graphs = "pca")
  bys <- as.SingleCellExperiment(bys)
  colData(bys) <- cbind(colData(bys), Autism[[i]]@images$slice1@coordinates)
  bys <- spatialPreprocess(bys, platform="Visium", n.PCs=15, n.HVGs=2000, log.normalize=T)
  ASD.Bys[[i]] <- spatialCluster(bys, q=24, platform="Visium", d=15, init.method="mclust", model="t", gamma=2,
                                 nrep=1000, burn.in=100,save.chain=TRUE)
  
  BysCP <- clusterPlot(ASD.Bys[[i]], label = "spatial.cluster", is.enhanced = F, platform = "Visium", palette = cols)
  pdf(paste("./figures/Bayes/ASD_ClusterPlot_Bayes_", names(dir)[i], ".pdf", sep =""), width = 7.5, height = 7.5)
  print(BysCP)
  dev.off()
  print(paste("BayesCluster", names(dir)[i], "is over", sep =" "))
}

for (i in 1:9) {
  j = length(table(ASD.Bys[[i]]$spatial.cluster))
  ASD.BysE[[i]] <- spatialEnhance(ASD.Bys[[i]], q=j, platform="Visium", d=15,
                                 model="t", gamma=2,
                                 jitter_prior=0.3, jitter_scale=3.5,
                                 nrep=1000, burn.in=100,
                                 save.chain=T)
  BysEh <- clusterPlot(ASD.BysE[[i]],label = "spatial.cluster" , palette = cols) + ggtitle("enhanced")
  pdf(paste("./figures/Bayes/BayesEnhance/ASD_ClusterPlot_BayesEnhance_", names(dir)[i],
            ".pdf", sep =""), width = 7.5, height = 7.5)
  print(BysEh)
  dev.off()
  print(paste("BayesEnhance", names(dir)[i], "is over", sep =" "))
}


newname <- c('IsoCTX', 'OLF', 'ENT&OLF','CTXsp',
                      'HIP','STR','TH','HY',
                      'MB', 'Fiber tracts', 'VC')
names(newname) <- c('Isocortex', 'Olfactory areas', 'Entorhinal & Olfactory areas', 'Cortical subplate',
  'Hippocampal region', 'Striatum', 'Thalamus', 'Hypothalamus', 'Midbrain',
  'Fiber tracts', 'ventricular systems')
ASD.brain <- RenameIdents(ASD.brain, newname)
ASD.brain$regionSimp <- ASD.brain@active.ident
ASD.brain$regionSimp <- factor(ASD.brain$regionSimp, levels = rev())


for (i in 1:3){
  # Region SpatialDimPlot
  STDimP <- SpatialDimPlot(Autism[[i]], group.by = "RegionSub", cols = ann_cols_sub, pt.size = 1.25, image.alpha = 0, crop = F)
  pdf(paste("./figures/Region/ASD_RegionSpatialDimPlot_sub", names(dir)[i], ".pdf", sep =""), width = 5, height = 5)
  print(STDimP)
  dev.off()

}




ASD.brain$Region <- factor(ASD.brain$Region, levels = c('Isocortex', 'Olfactory areas', 'Entorhinal & Olfactory areas', 'Cortical subplate',
                                                        'Hippocampal region', 'Striatum', 'Thalamus', 'Hypothalamus', 'Midbrain', 
                                                        'Fiber tracts', 'ventricular systems'))




ASD.brain$RegionSub <- factor(ASD.brain$RegionSub, levels = c('IsoCTX_L1', 'IsoCTX_L2/3/4', 'IsoCTX_L5', 'IsoCTX_L6', 
                                                              'OLF', 'ENT&OLF', 'OLF_AON', 'OLF_PIR',
                                                              'CTXsp', 'HIP_DG', 'HIP_CA1', 'HIP_CA2/3',
                                                              'STR', 'STR_CP', 'STR_sAMY',
                                                              'TH', 'TH_DORpm', 'TH_DORsm', 'HY',
                                                              'MBmot','MB_VTA', 'Fiber tracts', 'VC'))
ASD.brain$RegionSub_group <- paste(ASD.brain$RegionSub, ASD.brain$group, sep = '@')
ASD.brain$Region_group <- paste(ASD.brain$regionSimp, ASD.brain$group, sep = '@')



for (i in 1:9) {
  bys <- ASD.Bys[[i]]
  bys$region <- Autism[[i]]$RegionSub
  p <- clusterPlot(bys, label = "region", is.enhanced = F, platform = "Visium", palette = byscol, color = NA)
  pdf(paste("./figures/Region/ASD_RegionSpatialDimPlot_bys_n", names(dir)[i], ".pdf", sep =""), width = 5, height = 5)
  print(p)
  dev.off()
}

#fiber
pdf('./figures/Region/Spatialfeatureplots_markers_Mbp.pdf', height = 10, width = 10)
p<-SpatialFeaturePlot(object = Autism[[i]], images ='slice1',features = "Mbp",image.alpha = 0,
                   pt.size.factor = 1.25,crop = F)+NoLegend()+scale_fill_viridis()
print(p)
dev.off()
#hip
pdf('./figures/Region/Spatialfeatureplots_markers_Hpca.pdf', height = 10, width = 10)
p<-SpatialFeaturePlot(object = Autism[[i]], images ='slice1',features = "Hpca",image.alpha = 0,
                      pt.size.factor = 1.25,crop = F)+NoLegend()+scale_fill_viridis()
print(p)
dev.off()
#ctx
pdf('./figures/Region/Spatialfeatureplots_markers_Stx1a.pdf', height = 10, width = 10)
p<-SpatialFeaturePlot(object = Autism[[i]], images ='slice1',features = "Stx1a",image.alpha = 0,
                      pt.size.factor = 1.25,crop = F)+NoLegend()+scale_fill_viridis()
print(p)
dev.off()
#TH
pdf('./figures/Region/Spatialfeatureplots_markers_Prkcd.pdf', height = 10, width = 10)
p<-SpatialFeaturePlot(object = Autism[[i]], images ='slice1',features = "Prkcd",image.alpha = 0,
                      pt.size.factor = 1.25,crop = F)+scale_fill_viridis()
print(p)
dev.off()
#HY
pdf('./figures/Region/Spatialfeatureplots_markers_Gpx3.pdf', height = 10, width = 10)
p<-SpatialFeaturePlot(object = Autism[[i]], images ='slice1',features = "Gpx3",image.alpha = 0,
                      pt.size.factor = 1.25,crop = F)+scale_fill_viridis()
print(p)
dev.off()
#OLF
pdf('./figures/Region/Spatialfeatureplots_markers_Lmo3.pdf', height = 10, width = 10)
p<-SpatialFeaturePlot(object = Autism[[i]], images ='slice1',features = "Lmo3",image.alpha = 0,
                      pt.size.factor = 1.25,crop = F)+scale_fill_viridis()
print(p)
dev.off()

ASD.brain <- merge(Autism[[1]], y=c(Autism[[2]], Autism[[3]], Autism[[4]], Autism[[5]], Autism[[6]], Autism[[7]],
                                    Autism[[8]], Autism[[9]]), project = "ASD_ST")
#ASD.brain <- SCTransform(ASD.brain, assay = 'Spatial')
ASD.brain <- NormalizeData(ASD.brain, assay = "Spatial")
ASD.brain <- ScaleData(ASD.brain, verbose = FALSE, assay = "Spatial")
ASD.brain <- FindVariableFeatures(ASD.brain, assay = "Spatial")
DefaultAssay(ASD.brain) <- 'Spatial'
ASD.brain <- RunPCA(ASD.brain)
ASD.brain <- RunHarmony(ASD.brain,group.by.vars = "Sample")
ElbowPlot(ASD.brain,reduction = "harmony")

ASD.brain <- RunUMAP(ASD.brain,reduction = "harmony",dims = 1:30)
ASD.brain <- FindNeighbors(ASD.brain, reduction = "harmony", dims = 1:30)
ASD.brain <- FindClusters(ASD.brain, 
                          resolution =0.1,
                          method ='igraph',  
                          verbose = T)
save(ASD.brain, file = './object/ASD.brain_renameIso.RData')

DimPlot(ASD.brain, group.by = "Sample", pt.size = 1)
DimPlot(ASD.brain, group.by = "Region", pt.size = 1, cols = byscol)

#UMAP plot of all spots
pdf('./figures/Dimplot_all_subregion_new.pdf', height = 6, width = 8)
p<-DimPlot(ASD.brain, group.by = "RegionSub", pt.size = 1, cols = byscol) +ggtitle(' ')
print(p)
dev.off()

#vlnplot markers
# ASD.brain$RegionSub <- factor(ASD.brain$RegionSub, levels = sort(unique(ASD.brain$RegionSub)))
# Idents(ASD.brain) <- ASD.brain$RegionSub
p<-VlnPlot(ASD.brain, features = c('Pak7','Myl4','Epop','Cobl','Stx1a',#IsoCTX
                                   'Mbp', #FB
                                   'Adora2a',#CP
                                   'Lmo3','Slc30a3','Nptx1','Prkcd',#TH
                                   'Hpca', #HIP
                                   'Ddc', 'Slc6a3'#MB
                                   ), group.by = 'RegionSub',stack = T, sort = F, 
           flip = T, split.by = 'RegionSub')+NoLegend() + scale_fill_manual(values = byscol)
pdf('./figures/VlnPlot_markers_new.pdf', height = 5, width = 7.5)
print(p)
dev.off()

regionmarkers <- FindAllMarkers(ASD.brain, only.pos = T)
regionmarkers10 <- regionmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

p<-DotPlot(ASD.brain, features = c('Epop','Hs3st2', 'Stx1a', #IsoCTX 'nrgn'
                                'Lmo3', 'Ptprt','Frzb',#OLF/CTX 'plxnd1'
                                'Sst', 'Hpcal1', 'Nptxr', #CTXsp 'lypd1'
                                'Hpca', 'Cabp7', 'Wipf3', #HIP
                                'Ppp1r1b', 'Gpr88', 'Tmem158', #STR  'Adora2a'
                                'Syt9','Synpo2', 'Prkcd',#TH 'Ramp3', 
                                'Baiap3', 'Dlk1',  'Gpx3',#HY 'Hcrt',
                                'Slc6a11', 'Gata3', 'Slc6a3', #MB 'Ddc'
                                'Mbp', 'Tmem88b', 'Cldn11', #Fiber 'S1pr5',
                                'Ttr', 'Enpp2', 'Slc16a8' #VC 'Folr1',
), group.by = 'regionSimp')+coord_flip()+
  scale_color_gradientn(colours = c("#330066","#336699","#66CC66","#FFCC33"))+
  theme(#panel.grid = element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 0.9,vjust=0.9))
p
pdf('./figures/Seurat/DotPlot_markers.pdf', height = 6.4, width = 8)
print(p)
dev.off()


ASD.brain@meta.data$RegionAll <- ASD.brain@meta.data$RegionSub
ASD.brain@meta.data[which(ASD.brain@meta.data$RegionAll == 'OLF_AON'),]$RegionAll <- 'OLF'
ASD.brain@meta.data[which(ASD.brain@meta.data$RegionAll == 'OLF_PIR'),]$RegionAll <- 'OLF'
ASD.brain@meta.data[which(ASD.brain@meta.data$RegionAll == 'STR_CP'),]$RegionAll <- 'STR'
ASD.brain@meta.data[which(ASD.brain@meta.data$RegionAll == 'STR_sAMY'),]$RegionAll <- 'STR'
ASD.brain@meta.data[which(ASD.brain@meta.data$RegionAll == 'TH_DORpm'),]$RegionAll <- 'TH'
ASD.brain@meta.data[which(ASD.brain@meta.data$RegionAll == 'TH_DORsm'),]$RegionAll <- 'TH'

p<-dittoBarPlot(ASD.brain, 'RegionAll', group.by = 'Sample', color.panel = byscol,
                x.reorder = c(1,7,4,2,8,5,3,9,6))+xlab('')+
  theme(axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 10))
pdf('./figures/ST/Region/dittoBarPlot.pdf', height = 5, width = 6)
print(p)
dev.off()

#####DEG#####
#AN
Allreg <- sort(unique(ASD.brain$RegionAll))

Idents(ASD.brain) <- ASD.brain$RegionAll
Markers_AN <- list()
for (i in 1:length(Allreg)) {
  markers <- FindMarkers(object = ASD.brain, ident.1 = 'Autism', ident.2 = 'Normal', min.pct = 0.1,
                         logfc.threshold = 0.25, group.by = 'group',subset.ident = Allreg[i])
  markers <- markers[which(markers$p_val < 0.05),]
  Markers_AN[[i]] <- markers
}
    
names(Markers_AN) <- Allreg

Markers_AN_merge <- data.frame(matrix(ncol = 8, nrow = 1))
colnames(Markers_AN_merge) <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","gene","RegionSub","Region" )
for (i in 1:length(Markers_AN)) {
  markers <- Markers_AN[[i]]
  markers$gene <- row.names(markers)
  markers$RegionSub <- names(Markers_AN[i])
  markers$Region <- strsplit(names(Markers_AN[i]),'_')[[1]][1]
  row.names(markers) <- NULL
  Markers_AN_merge <- rbind(Markers_AN_merge,markers)
}
Markers_AN_merge <- Markers_AN_merge[-1,]

Markers_AN_merge[which(Markers_AN_merge$RegionSub == 'MBmot'),]$Region <- 'MB'

names(table(Markers_AN_merge$Region))

Markers_AN_merge <- Markers_AN_merge[which(!Markers_AN_merge$gene %in% c('Hbb-bs', 'Hbb-bt', 'Hba-a2')),]
write.csv(Markers_AN_merge, file = './results/DEG/Markers_AN_merge.csv')

freq <- as.data.frame(table(Markers_AN_merge[which(abs(Markers_AN_merge$avg_log2FC)>1),]$gene))
freq_genes <- top_n(freq, 20, wt = Freq)$Var1
table(Markers_AN_merge$Region)

#####DEG plots#####
topgene_all <- Markers_AN_merge[1,]
topgene_all <- topgene_all[-1,]
for (i in 1:length(Allreg)) {
  topgenes <- filter(Markers_AN_merge,RegionSub==Allreg[i]) %>% distinct(gene,.keep_all = T) %>% top_n(5,abs(avg_log2FC))
  topgene_all <- rbind(topgene_all, topgenes)
}

topgene_all$generegion <- paste(topgene_all$gene, topgene_all$RegionSub, sep = '-')
Markers_AN_merge$generegion <- paste(Markers_AN_merge$gene, Markers_AN_merge$RegionSub, sep = '-')
Markers_AN_merge$size <- case_when(!(Markers_AN_merge$generegion %in% topgene_all$generegion)~ 1,
                                   Markers_AN_merge$generegion %in% topgene_all$generegion ~ 2)

Markers_AN_merge_volo <- filter(Markers_AN_merge, size == 1)
# volocano
# genes <- c('Ttr', 'Cdkn1a','Igfbp7')#, 'B2m', 'Hbb-bs'
# genes <- c('Ttr', 'Dgkz', 'Cdkn1a','Igfbp7')

p<-ggplot()+
  geom_jitter(data = Markers_AN_merge_volo, # 绘制所有数据点
              aes(x = RegionSub, y = avg_log2FC, color = RegionSub),
              size = 0.85,
              width =0.4, show.legend = F) +
  scale_color_manual(values = byscol)+
  geom_jitter(data = topgene_all, # 绘制top10数据点
              aes(x = RegionSub, y = avg_log2FC, color = RegionSub),
              size = 1,
              width =0.35) +
  scale_color_manual(values = byscol)+
  geom_tile(data = topgene_all, # 绘制中心分组标记图
            aes(x = RegionSub,
                y = 0,
                fill = RegionSub),
            height=0.4,
            # width = 2.5,
            color = "black",
            alpha = 1,
            show.legend = F)+
  scale_fill_manual(values = byscol) +
  geom_text_repel(data = topgene_all,  # 这里的filter很关键，筛选你想要标记的基因
                  aes(x = RegionSub, y = avg_log2FC, label = gene),
                  size = 4, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 5),
                  color = 'black',
                  force = 1.2,
                  arrow = arrow(length = unit(0.008, "npc"),
                                type = "open", ends = "last")) +
  labs(x="Regions", y="Average logFC") +
  geom_text(data=Markers_AN_merge, # 绘制中心分组标记图文本注释
            aes(x=RegionSub, y=0, label=RegionSub), size = 4, color ="white") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, color = 'black'),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 13))
pdf('./figures/ST/DEG/mergevolocanoAN_new.pdf', width = 15, height = 7.5)
print(p)
dev.off()

#####DEG barplot num#####
Markers_AN_merge$group <- 'stable'
Markers_AN_merge[which(Markers_AN_merge$avg_log2FC>0),]$group <- 'up'
Markers_AN_merge[which(Markers_AN_merge$avg_log2FC<0),]$group <- 'down'
Markers_AN_bar <- Markers_AN_merge %>%group_by(Region) %>% 
  summarise(up = sum(group == 'up', na.rm = T),
            down = -sum(group == 'down', na.rm = T))
# Markers_AN_bar$region <- gsub('_.*$', '', Markers_AN_bar$RegionSub)
# Markers_AN_bar[which(Markers_AN_bar$region == 'MBmot'),]$region <- 'MB'
mainreg <- unique(Markers_AN_merge$Region)

p<-ggplot() + 
  geom_col(data = Markers_AN_bar, aes(Region, up, fill = mainreg, width = 0.7),
           position = "dodge") +
  geom_col(data = Markers_AN_bar, aes(Region, down, fill = mainreg, width = 0.7),
           position = "dodge") +
  scale_fill_manual(values = maincol)+# coord_flip() + 
  geom_hline(yintercept = 0, col="black", size=0.6) +
  geom_text(data = Markers_AN_bar, 
            aes(Region, up, label = up),
            position = position_dodge(1), hjust = 0.5,vjust = -0.5, size = 5)+
  geom_text(data = Markers_AN_bar, 
            aes(Region, down, label = down),
            position = position_dodge(1), hjust = 0.5,vjust = 1.5, size = 5)+
  theme_classic()+ #Class R Graph style
  theme(panel.grid = element_blank(), legend.title=element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.line.x=element_line(linetype=1,color="black",size=0.5),
        axis.line.y=element_line(linetype=1,color="black",size=0.5))+ylim(-2100,1000)

pdf("./figures/DEG/Barplot_AN_nums.pdf", width = 11, height = 6)
print(p)
dev.off()

Markers_AC_merge$group <- 'stable'
Markers_AC_merge[which(Markers_AC_merge$avg_log2FC>0),]$group <- 'up'
Markers_AC_merge[which(Markers_AC_merge$avg_log2FC<0),]$group <- 'down'
Markers_AC_bar <- Markers_AC_merge %>%group_by(Region) %>% 
  summarise(up = sum(group == 'up', na.rm = T),
            down = -sum(group == 'down', na.rm = T))

p<-ggplot() + 
  geom_col(data = Markers_AC_bar, aes(Region, up, fill = 'up', width = 0.5),
           position = "dodge") +
  geom_col(data = Markers_AC_bar, aes(Region, down, fill = 'down', width = 0.5),
           position = "dodge") + coord_flip() + geom_hline(yintercept = 0, col="black", size=0.75) +
  geom_text(data = Markers_AC_bar, 
            aes(Region, up,fill = 'up', label = up),
            position = position_dodge(1), hjust = 0,size = 4)+
  geom_text(data = Markers_AC_bar, 
            aes(Region, down,fill = 'down', label = down),
            position = position_dodge(1), hjust = 1,size = 4)+
  scale_fill_manual(values = c("#3182bd","#de2d26"))+
  theme_classic()+ #Class R Graph style
  theme(panel.grid = element_blank(), legend.title=element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.line.x=element_line(linetype=1,color="black",size=0.5),
        axis.line.y=element_line(linetype=1,color="black",size=0.5))+ylim(-1200,1000)

pdf("./figures/DEG/Barplot_AC_nums.pdf", width = 8, height = 4)
print(p)
dev.off()

#####DEG logFC heatmap#####
Markers_AN_FC_hm <- Markers_AN_merge[, c(2, 6, 7)]
Markers_AN_FC_hm <- dcast(Markers_AN_FC_hm, gene~RegionSub, value.var = 'avg_log2FC')
# Markers_AN_FC_hm[is.na(Markers_AN_FC_hm)] <- 0
row.names(Markers_AN_FC_hm) <- Markers_AN_FC_hm$gene
Markers_AN_FC_hm <- Markers_AN_FC_hm[,-1]


Markers_AN_FC_hm_count <- Markers_AN_FC_hm %>%
  mutate(across(everything(), ~case_when(
    . > 0 ~ 'up',
    . < 0 ~ 'down',
    . == 0 ~ 'stable',  # 确保这里是逻辑表达式
    TRUE ~ NA_character_  # 处理NA值，保留为NA
  )))
Markers_AN_FC_hm_count <- as.data.frame(t(Markers_AN_FC_hm_count))
countlist <- lapply(Markers_AN_FC_hm_count, table)
countdf <- do.call(rbind, lapply(countlist, data.frame))
countdf$gene <- row.names(countdf)
countdf$gene <- gsub('\\..*', '', countdf$gene)
countdf_d <- dcast(countdf, gene~Var1, value.var = 'Freq')
countdf_d[is.na(countdf_d)] <- 0

 
countdf_d_s <- countdf_d[which(countdf_d$stable < 2),]

countdf_d_s$gene
write.table(countdf_d_s, file = './results/DEG/logFC/AN_countdf_d_s1.txt', row.names = F)
Markers_AN_FC_hm_data <- Markers_AN_FC_hm[which(row.names(Markers_AN_FC_hm) %in% countdf_d_s$gene),]
Markers_AN_FC_hm_data <- Markers_AN_FC_hm_data[!row.names(Markers_AN_FC_hm_data) %in% c('Hbb-bs', 'Hba-a2'), ]
p<-pheatmap(Markers_AN_FC_hm_data,
         color = viridis(60,  begin = 0.2, end = 1),
         main = 'Autism vs Normal',
         na_col = 'grey',
         cluster_rows = T,
         cluster_cols = F)

pdf("./figures/ST/DEG/LogFC_HM_AN_.pdf", width = 6, height = 6)
print(p)
dev.off()

#####Upset#####
Markers_AN_upset <- dcast(Markers_AN_merge,gene~RegionSub,value.var = 'p_val')
Markers_AN_upset[,2:length(Markers_AN_upset)] <- ifelse(Markers_AN_upset[,2:length(Markers_AN_upset)] < 0.05, 1, 0)
Markers_AN_upset[is.na(Markers_AN_upset)] <- 0

library(UpSetR)
p <- upset(Markers_AN_upset, nsets = 17, nintersects = 300, order.by = 'freq',
           sets.x.label="Total DEGs per region", mainbar.y.label="Unique and overlapping DEGs sets",
           point.size = 4, line.size = 1.5, text.scale = c(2,2,1.5,1.5,2,2),
           #main.bar.color = c(pal_npg()(4)),
           sets.bar.color = pal_igv()(17)
)

save(Markers_AN, Markers_AC, file = './object/DEG/MarkersAN&AC.RData')
save(Markers_AN_merge, Markers_AC_merge, file = './object/DEG/MarkersAN&AC_merge.RData')

#####enrich#####

genelist <- unique(Markers_AN_merge[which(Markers_AN_merge$avg_log2FC>0),]$gene)
genelist <- unique(Markers_AN_merge$gene)
genelist <- unique(Markers_AN_merge[which(abs(Markers_AN_merge$avg_log2FC)>1),]$gene)

#write.csv(genelist, file = './results/DEG/enrich.GO_AN_genelist_logFC1.csv', row.names = F)

Region <- unique(Markers_AN_merge$Region)
enrich.go_df_mergeAN <- enrich.go_r_df[1,]
enrich.go_df_mergeAN <- enrich.go_df_mergeAN[-1,]
for (i in 1:length(Region)) {
  print(paste(Region[i] ,'begin', sep = ' '))
  region_now <- Region[i]
  genelist <- unique(Markers_AN_merge[which(Markers_AN_merge$Region == region_now),]$gene)
  enrich.go_r <- enrichGO(gene = genelist,
                        OrgDb = 'org.Mm.eg.db',
                        keyType = 'SYMBOL',
                        ont = 'ALL',
                        pAdjustMethod = 'fdr',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
  
  enrich.go_r_df <- as.data.frame(enrich.go_r)
  enrich.go_r_df$region <- Region[i]
  enrich.go_df_mergeAN <- rbind(enrich.go_df_mergeAN, enrich.go_r_df)
  print(paste(Region[i] ,'over', sep = ' '))
}
write.csv(enrich.go_df_mergeAN, file = './results/DEG/enrich.go_df_mergeAN.csv', row.names = F)

enrich.go <- enrichGO(gene = genelist,
                      OrgDb = 'org.Mm.eg.db',
                      keyType = 'SYMBOL',
                      ont = 'ALL',
                      pAdjustMethod = 'fdr',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

enrich.go.df <- as.data.frame(enrich.go)
write.csv(enrich.go.df, file = './results/DEG/enrich.GO_AN_ALL.csv')

index <- which(row.names(enrich.go@result) %in% row.names(top_n(enrich.go@result, 10, -p.adjust)))
des <- enrich.go@result$Description[index]

save(enrich.go, enrich.go.down, file = './object/AN_enrich.RData')
save(enrich.go, enrich.go.down, file = './object/AC_enrich.RData')


enrich.go.df_s <- top_n(enrich.go.df, 10, -p.adjust)
enrich.go.df_s$labelx <- rep(0.2, nrow(enrich.go.df_s))
enrich.go.df_s$labely <- seq(nrow(enrich.go.df_s), 1)
enrich.go.df_s$`-log10qvalue` <- -log10(enrich.go.df_s$qvalue)

p<-ggplot(data = enrich.go.df_s, aes(Count, reorder(Description,Count))) +
  geom_col(aes(x=Count,y=reorder(Description,Count),fill=`-log10qvalue`),width = 0.8) +
  geom_bar(stat="identity",alpha=0.1,width = 0.8) + 
  geom_text(aes(x=labelx,y=labely,label = Description),size=5, hjust =0)+
  theme_classic()+
  scale_fill_viridis(alpha = 0.9, begin = 0.7, end = 0.9, direction = -1, option = 'B')+
  #scale_fill_distiller(palette="Oranges", direction = 1)+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 16),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.line.y = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 16))+
  xlab("Count")+ylab('GO terms')+
  scale_x_continuous(expand = c(0.01,0))
pdf('./figures/DEG/enrich_GO_barplot.pdf', width = 6, height = 4)
print(p)
dev.off()

FINAL_GO <- read_csv("results/DEG/enrich/metascape_AN_ALL/Enrichment_GO/_FINAL_GO.csv")
enrich.go.df_s <- top_n(FINAL_GO, 15, -LogP)
enrich.go.df_s <- enrich.go.df_s[order(enrich.go.df_s$LogP),]
enrich.go.df_s$labelx <- rep(0.1, nrow(enrich.go.df_s))
enrich.go.df_s$labely <- seq(nrow(enrich.go.df_s), 1)
enrich.go.df_s$`-LogP` <- -enrich.go.df_s$LogP

p<-ggplot(data = enrich.go.df_s, aes(`-LogP`, reorder(Description,`-LogP`))) +
  geom_col(aes(x=`-LogP`,y=reorder(Description,`-LogP`),fill=`Z-score`),width = 0.8) +
  geom_bar(stat="identity",alpha=0.1,width = 0.8) + 
  geom_text(aes(x=labelx,y=labely,label = Description),size=5.5, hjust =0)+
  theme_classic()+
  scale_fill_viridis(alpha = 0.9, begin = 0.7, end = 0.9, direction = -1, option = 'B')+
  #scale_fill_distiller(palette="Oranges", direction = 1)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.5),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5),
        axis.text.x = element_text(colour = 'black', size = 14),
        axis.ticks.x = element_line(colour = 'black', linewidth = 0.5),
        axis.title.x = element_text(colour = 'black', size = 18),
        axis.title.y = element_text(colour = 'black', size = 18))+
  xlab("-LogP")+ylab('GO terms')+
  scale_x_continuous(expand = c(0.01,0))
pdf('./figures/DEG/Enrich/enrich_GO_barplot_metascope.pdf', width = 9, height = 6)
print(p)
dev.off()

pdf('./figures/DEG/enrich.dotplot_AN.pdf', width = 5, height = 5)
pdf('./figures/DEG/enrich.dotplot_AC.pdf', width = 6, height = 5)

p <- dotplot(enrich.go,showCategory = des, orderBy = 'p.adjust') +
  scale_color_gradientn(colours = c("#2166AC", "#B2182B"))
print(p)
dev.off()

R.utils::setOption('clusterProfiler.download.method', 'auto')
ENTREZID <- bitr(unique(Markers_AN_merge$gene), fromType = 'SYMBOL', toType = "ENTREZID", OrgDb = 'org.Mm.eg.db')[2]

enrich.kegg <- enrichKEGG(gene = ENTREZID$ENTREZID,  
                          keyType = 'kegg',  
                          organism = 'mmu',  #hsa 
                          pAdjustMethod = 'fdr',  #指定 p 值校正方法
                          pvalueCutoff = 0.05,  #
                          qvalueCutoff = 0.2)
enrich.kegg_df <- as.data.frame(enrich.kegg)
write.csv(enrich.kegg_df, file = './enrich.KEGG_AN.csv')

save(ENTREZID, file = './KEGG_ENTREZID.RData')

#GSEA


Markers_AN$OLF$avg_log2FC
rank <- Markers_AN$OLF$avg_log2FC[order(Markers_AN$OLF$avg_log2FC, decreasing = T)]
names(rank) <- row.names(Markers_AN$OLF)[order(Markers_all$OLF$avg_log2FC, decreasing = T)]

Markers_AC$STR$avg_log2FC
rank <- Markers_AC$STR$avg_log2FC[order(Markers_AC$STR$avg_log2FC, decreasing = T)]
names(rank) <- row.names(Markers_AC$STR)[order(Markers_AC$STR$avg_log2FC, decreasing = T)]
#####GSEA_KEGG#####
geneset_KEGG = msigdbr(species = "Mus musculus",#Homo sapiens
                       category = "C2", 
                       subcategory = "KEGG"
) %>% dplyr::select(gs_name,gene_symbol)
gsea_geneset_KEGG <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)

Agsea_res_KEGG <- fgsea(pathways = gsea_geneset_KEGG, 
                        stats = rank,
                        minSize=5,
                        maxSize=500,
                        nperm=1000
)
agea_rf_KEGG<-as.data.frame(Agsea_res_KEGG)
agea_rf_KEGG$leadingEdge<-as.character(agea_rf_KEGG$leadingEdge)


geneset_GOBP = msigdbr(species = "Mus musculus",#Homo sapiens
                       category = "C5", 
                       subcategory = "BP"
) %>% dplyr::select(gs_name,gene_symbol)
gsea_geneset_GOBP <- geneset_GOBP %>% split(x = .$gene_symbol, f = .$gs_name)

#####GSEA_GOBP#####
GSEA_res_GOBP <- fgsea(pathways = gsea_geneset_GOBP, 
                        stats = rank,
                        minSize=5,
                        maxSize=500,
                        nperm=1000
)
GSEA_res_GOBP<-as.data.frame(GSEA_res_GOBP)
GSEA_res_GOBP$leadingEdge<-as.character(GSEA_res_GOBP$leadingEdge)
write.csv(GSEA_res_GOBP, file = './results/DEG/enrich/GSEA_res_GOBP_AN.csv')

#####GSEA plot#####
GO_secleted<-c('GOBP_RESPONSE_TO_TEMPERATURE_STIMULUS', 
               'GOBP_POSITIVE_REGULATION_OF_DEFENSE_RESPONSE', 
               'GOBP_RESPONSE_TO_STEROID_HORMONE', 
               'GOBP_NCRNA_PROCESSING', 
               'GOBP_RESPONSE_TO_BMP')


a<-data.frame(NA,NA,NA)
colnames(a)<-c("Rank","Enrichment Score","pathway")
a<-na.omit(a)
for (i in 1:length(GO_secleted)) {
  p<-plotEnrichment(gsea_geneset_GOBP[[GO_secleted[i]]],
                    rank,ticksSize = 1)
  pdata<-data.frame(p$data)
  pdata$pathway<-paste(GO_secleted[i],' p=',signif(GSEA_res_GOBP[grepl(paste0(GO_secleted[i],"$"),GSEA_res_GOBP$pathway),"pval"],2)," NES=",round(GSEA_res_GOBP[grep(GO_secleted[i],GSEA_res_GOBP$pathway),"NES"],2)) 
  colnames(pdata)<-c("Rank","Enrichment Score","pathway")
  #write.table(pdata,paste("GSEA/",GO_secleted[i],'.txt',collapse = '',sep = ""),sep = "\t",quote = F,row.names = F)
  a<-rbind(a,pdata)
  rm(pdata)
}

{
  df<-a
  # plot
  up <- ggplot(df,aes(x =Rank,y = `Enrichment Score`,color = pathway)) +
    geom_hline(yintercept = 0,color = 'grey',size = 1) +#
    geom_line(show.legend = T,size = 1) +#
    scale_color_manual(values = col1) +
    theme_bw(base_size = 18) +#ylim(-1,1)+
    theme(panel.grid = element_blank(),
          axis.ticks.length = unit(0.25,'cm'),
          panel.border = element_rect(size = 2),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),axis.text = element_text(size = 15,colour = "black"),
          axis.title.y = element_text(size = 15,color = "black"),
          axis.ticks.x = element_blank(),legend.position="right") +
    ylab('Enrichment score') 
  
  # line
  mid <- ggplot(df,aes(x = Rank,y = pathway)) +
    geom_vline(aes(xintercept = Rank,
                   color = pathway),
               show.legend = F,size = 1) +
    scale_color_manual(values = col1) +
    theme_classic(base_size = 18) +
    facet_wrap(~pathway,ncol = 1) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.ticks.length = unit(0.25,'cm'),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(color = 'white'),
          axis.title.y = element_text(color = 'white'))
  # heatmap
  result<-matrix(0,length(rank),2)
  result<-as.data.frame(result)
  colnames(result)<-c('Rank','y')
  result$Rank<-1:length(rank)
  bot <- ggplot(result,aes(x = Rank,y = 1,fill = Rank)) +
    geom_col(width = 1,show.legend = F) +
    scale_fill_gradient2(low = 'red',mid = 'white',high = 'blue',
                         midpoint = max(result$Rank)/2) +
    scale_x_continuous(breaks = seq(0,max(result$Rank),5000),
                       labels = seq(0,max(result$Rank),5000)/1000) +
    theme_classic(base_size = 18) +
    coord_cartesian(expand = 0) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(size = 10,color = "black"),
          axis.title.x = element_text(size = 15,color = "black"),
          axis.ticks.length = unit(0.25,'cm')) +
    xlab('Ranked in ordered dataset\n(x1000)') +
    labs(caption = c("High", "Low")) +
    theme(plot.caption = element_text(hjust=c(0, 1),color = c('red','blue')))

  p<-up + mid + bot + plot_layout(nrow = 3,heights = c(1,0.5,0.2))
  rm(up, mid, bot)
}

pdf('./figures/DEG/GSEA_HIP_AN_DEG.pdf', width = 13, height = 6)
print(p)
dev.off()

GSEA_BARPLOT <- function(x){
  result_df <- x %>% top_n(n = 4, wt = NES)
  result_df$`-log10(pval)` <- -log10(result_df$pval)
  
  result_df$pathway <- gsub("_", " ", result_df$pathway)
  result_df$pathway <- tolower(result_df$pathway)
  result_df$pathway <- str_to_title(result_df$pathway) 
  result_df$pathway <- substr(result_df$pathway, 10, nchar(result_df$pathway))
  result_df <- result_df[order(result_df$NES),]
  result_df$pathway <- factor(result_df$pathway, levels = result_df$pathway)
  
  p<-ggplot(result_df, aes(x = pathway, y = NES, fill = `-log10(pval)`)) +
    geom_col(width = 0.7) +
    coord_flip() +  # 将坐标轴翻转，使 pathway 名字在 y 轴
    scale_fill_gradient(low = "#6BAED6", high = "#084594", name = "-log10(pval)") +  # 颜色映射
    labs(x = "", y = "NES") +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 14, angle = 0, colour = 'black'),# 调整 y 轴标签
      axis.text.x = element_text(size = 12, angle = 0, colour = 'black'),
      axis.title.x = element_text(size = 14, angle = 0, colour = 'black')
    )
  return(p)
}

GSEA_df_RE_CTX_L1 <- GSEA_df_RE_CTX_L1[which(GSEA_df_RE_CTX_L1$pval<0.05),]
p <- GSEA_BARPLOT(GSEA_df_RE_CTX_L1)
pdf("./figures/ST/DEG/GSEA_RE_CTX_L1_AN_barplot.pdf", width = 8, height = 3)
print(p)
dev.off()

GSEA_df_RE_HIP_DG <- GSEA_df_RE_HIP_DG[which(GSEA_df_RE_HIP_DG$pval<0.05),]
p <- GSEA_BARPLOT(GSEA_df_RE_HIP_DG)
pdf("./figures/ST/DEG/GSEA_RE_HIP_DG_AN_barplot.pdf", width = 9, height = 2)
print(p)
dev.off()

#####GSEA_section_SASP######
GSEA_PLOT <- function(GO_secleted, GSEA_res, gsea_geneset, var, gsecol){
  genelist <- Markers_AN_merge[which(Markers_AN_merge$RegionSub == var),]
  rank <- genelist$avg_log2FC[order(genelist$avg_log2FC, decreasing = T)]
  names(rank) <- genelist$gene[order(genelist$avg_log2FC, decreasing = T)]
  
  a<-data.frame(NA,NA,NA)
  colnames(a)<-c("Rank","Enrichment Score","pathway")
  a<-na.omit(a)
  for (i in 1:length(GO_secleted)) {
    p<-plotEnrichment(gsea_geneset[[GO_secleted[i]]],
                      rank,ticksSize = 1)
    pdata<-data.frame(p$data)
    pdata$pathway<-paste(GO_secleted[i],' p=',signif(GSEA_res[grepl(paste0(GO_secleted[i],"$"),GSEA_res$pathway),"pval"],2),
                         " NES=",round(GSEA_res[grep(GO_secleted[i],GSEA_res$pathway),"NES"],2)) 
    colnames(pdata)<-c("Rank","Enrichment Score","pathway")
    #write.table(pdata,paste("GSEA/",GO_secleted[i],'.txt',collapse = '',sep = ""),sep = "\t",quote = F,row.names = F)
    a<-rbind(a,pdata)
    rm(pdata)
  }
  
  {
    df<-a
    # plot
    up <- ggplot(df,aes(x =Rank,y = `Enrichment Score`,color = pathway)) +
      geom_hline(yintercept = 0,color = 'grey',size = 1) +#
      geom_line(show.legend = T,size = 1) +#
      scale_color_manual(values = gsecol) +
      theme_bw(base_size = 18) +#ylim(-1,1)+
      theme(panel.grid = element_blank(),
            axis.ticks.length = unit(0.25,'cm'),
            panel.border = element_rect(size = 2),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),axis.text = element_text(size = 15,colour = "black"),
            axis.title.y = element_text(size = 15,color = "black"),
            axis.ticks.x = element_blank(),legend.position="right") +
      ylab('Enrichment score') 
    
    # line
    mid <- ggplot(df,aes(x = Rank,y = pathway)) +
      geom_vline(aes(xintercept = Rank,
                     color = pathway),
                 show.legend = F,size = 1) +
      scale_color_manual(values = gsecol) +
      theme_classic(base_size = 18) +
      facet_wrap(~pathway,ncol = 1) +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            axis.ticks.length = unit(0.25,'cm'),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.y = element_line(color = 'white'),
            axis.title.y = element_text(color = 'white'))
    # heatmap
    result<-matrix(0,length(rank),2)
    result<-as.data.frame(result)
    colnames(result)<-c('Rank','y')
    result$Rank<-1:length(rank)
    bot <- ggplot(result,aes(x = Rank,y = 1,fill = Rank)) +
      geom_col(width = 1,show.legend = F) +
      scale_fill_gradient2(low = '#C94741',mid = 'white',high = '#3783BB',
                           midpoint = max(result$Rank)/2) +
      scale_x_continuous(breaks = seq(0,max(result$Rank),5000),
                         labels = seq(0,max(result$Rank),5000)/1000) +
      theme_classic(base_size = 18) +
      coord_cartesian(expand = 0) +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_text(size = 10,color = "black"),
            axis.title.x = element_text(size = 15,color = "black"),
            axis.ticks.length = unit(0.25,'cm')) +
      xlab('Ranked in ordered dataset\n(x1000)') +
      labs(caption = c("High", "Low")) +
      theme(plot.caption = element_text(hjust=c(0, 1),color = c('#C94741','#3783BB')))
    
    p<-up + mid + bot + plot_layout(nrow = 3,heights = c(1,0.2,0.2))
    return(p)
  }
  
}


GO_secleted <- c('REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP',
                 'REACTOME_CELLULAR_SENESCENCE')



pdf("./figures/ST/DEG/GSEA_RE_CTX_L5_AN_SASP.pdf", width = 16, height = 6)
print(GSEA_PLOT(GO_secleted,GSEA_df_RE_CTX_L5,geneset_RE, 'CTX_L5', c('#C94741','#3783BB')))
dev.off()


pdf("./figures/ST/DEG/GSEA_RE_HIP_CA23_AN_SASP.pdf", width = 16, height = 6)
print(GSEA_PLOT(GO_secleted,GSEA_df_RE_HIP_CA23,geneset_RE, 'HIP_CA2/3', c('#C94741','#3783BB')))
dev.off()

pdf("./figures/ST/DEG/GSEA_RE_MBmot_AN_SASP.pdf", width = 16, height = 6)
print(GSEA_PLOT(GO_secleted,GSEA_df_RE_AN[['MBmot']], geneset_RE, 'MBmot', c('#C94741','#3783BB')))
dev.off()

pdf("./figures/ST/DEG/GSEA_RE_MB_VTA_AN_SASP.pdf", width = 16, height = 6)
print(GSEA_PLOT(GO_secleted,GSEA_df_RE_AN[['MB_VTA']], geneset_RE, 'MB_VTA', c('#C94741','#3783BB')))
dev.off()


igfbp7_exp <- data.frame(FetchData(object = ASD.brain, vars = c('Igfbp7', 'Region', 'RegionAll', 'group')))
igfbp7_exp_hm_df <- igfbp7_exp %>% group_by(RegionAll, group) %>% 
  summarise(mean_exp = mean(Igfbp7, na.rm = T))
igfbp7_exp_hm_df <- dcast(igfbp7_exp_hm_df,group~RegionAll,value.var = 'mean_exp')
row.names(igfbp7_exp_hm_df) <- igfbp7_exp_hm_df$group
igfbp7_exp_hm_df <- igfbp7_exp_hm_df[,-1]
igfbp7_exp_hm_df <- igfbp7_exp_hm_df[c(1,3,2),]
igfbp7_exp_hm_df <- t(igfbp7_exp_hm_df)

p<-pheatmap(igfbp7_exp_hm_df,
         show_rownames = T,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows= F,
         fontsize_row = 8,
         fontsize_col = 10,
         cols <- viridis(8),
         border_color = 'grey',
         legend = T,
         angle_col = '90',
         #legend_labels = 'regulon activity',
         clustering_method = 'single'
) 
pdf('./figures/igfbp7/igfbp7_heatmap_exp.pdf', height = 3.5, width = 3)
print(p)
dev.off()

DoHeatmap(ASD.brain, features = Gata2.genes, group.by = 'group', 
          group.colors = byscol, size = 3, angle = 90)+scale_fill_viridis()+NoLegend()



p<-VlnPlot(ASD.brain, features = "Igfbp7", pt.size = 0, group.by = 'RegionAll', split.by = 'group',
        cols = pal_npg()(3)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))#+
  # stat_compare_means(comparisons = list(c('Control', 'CCL4')))+
  # ylim(-0.1, 1.5)+scale_fill_lancet()
pdf('./figures/igfbp7/igfbp7_vlnplot_exp.pdf', height = 5, width = 15)
print(p)
dev.off()

p<-VlnPlot(ASD.brain, features = 'Igfbp7',group.by = 'group', sort = T, pt.size = 0, cols = pal_npg()(3))+
  stat_compare_means(comparisons = list(c('Autism', 'Normal'), c('Autism', 'Control')),label = 'p.signif')+
  ylim(-0.1, 7)
pdf('./figures/igfbp7/igfbp7_vlnplot_single.pdf', height = 5, width = 4.5)
print(p)
dev.off()


####Cdkn1a####
#####SpatialFeaturePlot####
p<-SpatialFeaturePlot(ASD.brain, features = 'Cdkn1a', image.alpha = 0, ncol = 3, 
                      pt.size.factor = 1.25, crop = F) & scale_fill_viridis()

pdf('./figures/SpatialFeaturePlot_Cdkn1a.pdf', height = 20, width = 10)
print(p)
dev.off()

#####heatmap#####
Cdkn1a_exp <- data.frame(FetchData(object = ASD.brain, vars = c('Cdkn1a', 'Region', 'RegionAll', 'group')))
Cdkn1a_exp_hm_df <- Cdkn1a_exp %>% group_by(RegionAll, group) %>% 
  summarise(mean_exp = mean(Cdkn1a, na.rm = T))
Cdkn1a_exp_hm_df <- dcast(Cdkn1a_exp_hm_df,group~RegionAll,value.var = 'mean_exp')
row.names(Cdkn1a_exp_hm_df) <- Cdkn1a_exp_hm_df$group
Cdkn1a_exp_hm_df <- Cdkn1a_exp_hm_df[,-1]
Cdkn1a_exp_hm_df <- Cdkn1a_exp_hm_df[c(1,3,2),]
Cdkn1a_exp_hm_df <- t(Cdkn1a_exp_hm_df)

p<-pheatmap(Cdkn1a_exp_hm_df,
            show_rownames = T,
            show_colnames = T,
            cluster_cols = F,
            cluster_rows= F,
            fontsize_row = 8,
            fontsize_col = 10,
            cols <- viridis(8),
            border_color = 'grey',
            legend = T,
            angle_col = '90',
            #legend_labels = 'regulon activity',
            clustering_method = 'single'
) 
pdf('./figures/Cdkn1a/Cdkn1a_heatmap_exp.pdf', height = 3.5, width = 3)
print(p)
dev.off()


p<-VlnPlot(ASD.brain, features = "Cdkn1a", pt.size = 0, group.by = 'RegionAll', split.by = 'group',
           cols = pal_npg()(3)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))#+
# stat_compare_means(comparisons = list(c('Control', 'CCL4')))+
# ylim(-0.1, 1.5)+scale_fill_lancet()
pdf('./figures/Cdkn1a/Cdkn1a_vlnplot_exp.pdf', height = 5, width = 15)
print(p)
dev.off()


####CellChat####
library(CellChat)
library(patchwork)
library(future)
set.seed(123)
showDatabaseCategory(CellChatDB.mouse)
future::plan("multisession", workers = 24) # do parallel

ASD.group <- list()
ASD.group[[1]] <- subset(ASD.brain, subset = group == 'Autism')
ASD.group[[2]] <- subset(ASD.brain, subset = group == 'Normal')
ASD.group[[3]] <- subset(ASD.brain, subset = group == 'Control')
save(ASD.group, file = './object/ASD.group.RData')

cellchat_list <- list()
for (i in 1:3) {
  sce <- ASD.group[[i]]
  cellchat <- createCellChat(object = sce, meta = sce@meta.data, group.by = "Region")
  cellchat@DB <- CellChatDB.mouse
  cellchat <- subsetData(cellchat)
  
  future::plan("multisession", workers = 24)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat_list[[i]] <- cellchat
}
rm(sce, cellchat)
names(cellchat_list) <- c('Autism', 'Normal', 'Control')
save(cellchat_list, file = './object/cellchat_list.RData')
rm(cellchat_list)

cellchat_db_mm <- CellChatDB.mouse
cellchat_db_mm <- as.data.frame(cellchat_db_mm$interaction)

object.list <- list()
object.list[[1]] <- cellchat_list[[1]]
object.list[[2]] <- cellchat_list[[2]]
object.list[[2]] <- cellchat_list[[3]]
names(object.list) <- c('Autism', 'Normal')
names(object.list) <- c('Autism', 'Control')

#####merge plots#####
#object.list <- list(Autism = cellchat_list[[1]], Normal = cellchat_list[[2]])
cellchatAll <- mergeCellChat(cellchat_list, add.names = names(cellchat_list),cell.prefix = TRUE)
cellchatAN <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
cellchatAC <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

save(cellchatAll, file = './object/cellchatAll_merge.RData')
save(cellchatAN, file = './object/cellchatAN_merge.RData')

gg1 <- compareInteractions(cellchatAll, show.legend = F, group = c(1,2,3),
                           color.use = c("#E64B35", "#4DBBD5","#00A087"))
gg2 <- compareInteractions(cellchatAll, show.legend = F, group = c(1,2,3), measure = "weight",
                           color.use = c("#E64B35", "#4DBBD5","#00A087"))
gg1 + gg2
pdf('./figures/cellchat/barplot.pdf', width = 8, height = 5)
print(gg1 + gg2)
dev.off()

gg1 <- compareInteractions(cellchatAN, show.legend = F, group = c(1,2),
                           color.use = c("#E64B35", "#4DBBD5"))
gg2 <- compareInteractions(cellchatAN, show.legend = F, group = c(1,2), measure = "weight",
                           color.use = c("#E64B35", "#4DBBD5"))
gg1 + gg2
pdf('./figures/cellchat/barplot_AN.pdf', width = 4, height = 3)
print(gg1 + gg2)
dev.off()

gg1 <- compareInteractions(cellchatAC, show.legend = F, group = c(1,2),
                           color.use = c("#E64B35", "#4DBBD5"))
gg2 <- compareInteractions(cellchatAC, show.legend = F, group = c(1,2), measure = "weight",
                           color.use = c("#E64B35", "#4DBBD5"))
gg1 + gg2
pdf('./figures/cellchat/barplot_AC.pdf', width = 4, height = 3)
print(gg1 + gg2)
dev.off()

groupSize1 <- as.numeric(table(cellchat_list[[1]]@idents))
groupSize2 <- as.numeric(table(cellchat_list[[2]]@idents))
groupSize3 <- as.numeric(table(cellchat_list[[3]]@idents))

#nets
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_circle(cellchat_list[[1]]@net$count, vertex.weight = groupSize1, weight.scale = T, 
                        label.edge= F, title.name = "Number of interactions in Autism")
gg2 <- netVisual_circle(cellchat_list[[2]]@net$count, vertex.weight = groupSize2, weight.scale = T, 
                        label.edge= F, title.name = "Number of interactions in Normal")
gg3 <- netVisual_circle(cellchat_list[[3]]@net$count, vertex.weight = groupSize3, weight.scale = T, 
                        label.edge= F, title.name = "Number of interactions in Control")
pdf('./figures/cellchat/netAN_n.pdf', width = 12, height = 6)
print(gg2)
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_circle(cellchat_list[[1]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                        label.edge= F, title.name = "Interaction weights/strength in Autism")
gg2 <- netVisual_circle(cellchat_list[[2]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                        label.edge= F, title.name = "Interaction weights/strength in Normal")
gg3 <- netVisual_circle(cellchat_list[[3]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                        label.edge= F, title.name = "Interaction weights/strength in Control")
pdf('./figures/cellchat/netAN_s.pdf', width = 12, height = 6)
print(gg2)
dev.off()



gg1 <- netVisual_heatmap(cellchat_list[[1]], measure = "weight")
gg2 <- netVisual_heatmap(cellchat_list[[2]], measure = "weight")
gg3 <- netVisual_heatmap(cellchat_list[[3]], measure = "weight")

gg1 + gg2 + gg3+ plot_layout(byrow = F)
pdf('./figures/cellchat/heatmapAN_s.pdf', width = 12, height = 4.5)
print(gg1 + gg2)
dev.off()

gg1 <- netVisual_heatmap(cellchat_list[[1]], measure = "count")
gg2 <- netVisual_heatmap(cellchat_list[[2]], measure = "count")
gg3 <- netVisual_heatmap(cellchat_list[[3]], measure = "count")

gg1 + gg2 + gg3+ plot_layout(byrow = F)
pdf('./figures/cellchat/heatmapAC_n.pdf', width = 12, height = 4.5)
print(gg1 + gg3)
dev.off()

i = 3
cellchat_list[[i]] <- netAnalysis_computeCentrality(cellchat_list[[i]])

pathway.union <- union(cellchat_list[[1]]@netP$pathways, cellchat_list[[2]]@netP$pathways)

gg1 <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "outgoing", 
                                         signaling = pathway.union, title = names(cellchat_list)[1], 
                                         width = 7.5, height = 10, color.use = maincol)
gg2 <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "outgoing", 
                                         signaling = pathway.union, title = names(cellchat_list)[2], 
                                         width = 7.5, height = 10, color.use = maincol)
gg1+gg2
pdf('./figures/cellchat/outAN_heatmap.pdf', width = 10, height = 7.5)
print(gg1 + gg2)
dev.off()



gg1 <- rankNet(cellchatAN, mode = "comparison", comparison = c(1,2), stacked = T, do.stat = TRUE,
               color.use = c("#E64B35", "#4DBBD5"), do.flip = T)#, measure = 'count'
gg2 <- rankNet(cellchatAN, mode = "comparison", comparison = c(1,2), stacked = F, do.stat = TRUE,
               color.use = c("#E64B35", "#4DBBD5"))#, measure = 'count'
gg1 + gg2
pdf('./figures/cellchat/compareAC.pdf', width = 10, height = 5)
print(gg1 + gg2)
dev.off()
#ggsave('./cellchat/plots/compare.png', width = 10, height = 5)


p <- netVisual_bubble(cellchatAll, sources.use = 'Isocortex', targets.use = c(1:11),  comparison = c(1, 2), 
                      angle.x = 45)
p <- netVisual_bubble(cellchatAN, sources.use = 'Isocortex', targets.use = c(1:10),  comparison = c(1, 2), 
                      angle.x = 90,color.text = c("#E64B35", "#4DBBD5"), font.size = 12)

p <- netVisual_bubble(cellchatAC, sources.use = 'Isocortex', targets.use = c(1:10),  comparison = c(1, 2), 
                      angle.x = 90,color.text = c("#E64B35", "#4DBBD5"), font.size = 12)

pdf('./figures/cellchat/dotplotAC.pdf', width = 10, height = 7.5)
print(p)
dev.off()

cellchat_list$Autism@netP$pathways
levels(cellchat_list$Autism@idents)
pathways.show <- c("APP")
netVisual_aggregate(cellchat_list$Autism, signaling = pathways.show,  
                    vertex.receiver = c(3,4,5,6,7,9),layout = "hierarchy")


netAnalysis_contribution(cellchat_list$Autism, signaling = 'APP')
netVisual_aggregate(cellchat_list$Autism, signaling = 'APP', layout = "circle")
netAnalysis_contribution(cellchat_list$Autism, signaling = 'COLLAGEN')
netVisual_aggregate(cellchat_list$Autism, signaling = 'COLLAGEN', layout = "circle")

netAnalysis_contribution(cellchat_list$Normal, signaling = 'NOTCH')
netVisual_aggregate(cellchat_list$Normal, signaling = 'NOTCH', layout = "circle")


####Scenic####
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(KernSmooth)

mat = ASD.brain@assays$Spatial@counts
mat = ASD.brain@assays$SCT@meta.features
mat = t(as.matrix(mat))
write.csv(mat, file = "./scenic/ASD_exp.csv")
rm(mat)

#return to R
loom <- open_loom('./scenic/CM_SCENIC.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)
close_loom(loom)
#
sub_regulonAUC <- regulonAUC[,match(colnames(ASD.brain),colnames(regulonAUC))]
dim(sub_regulonAUC)
identical(colnames(sub_regulonAUC), colnames(ASD.brain))

#groups
cellGroups <- data.frame(row.names = colnames(ASD.brain), 
                         celltype = ASD.brain$RegionAll)

cellGroups <- data.frame(row.names = colnames(ASD.brain), 
                         celltype = ASD.brain$group)


cellsPerGroup <- split(rownames(cellGroups), 
                       cellGroups[,'celltype'])

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)
#
# ASD.brain@meta.data$GroupRegionSub <- paste(ASD.brain@meta.data$group, 
#                                              ASD.brain@meta.data$RegionAll, sep = ".")

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

write.csv(regulonActivity_byGroup, file = './results/scenic/regulonActivity_byGroup.csv')

#####merge seurat#####
regulonsToPlot = c('Sox4(+)','Rreb1(+)','Klf16(+)','Cebpb(+)','Maz(+)','Gata2(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
ASD.brain@meta.data = cbind(ASD.brain@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))

# Vis
DotPlot(ASD.brain, features = unique(regulonsToPlot), group.by = 'group') #+ RotatedAxis()
p1 <- RidgePlot(ASD.brain, features = regulonsToPlot, cols = byscol, group.by = 'RegionAll',ncol = 3)
pdf('./figures/scenic/RidgePlot_scenic_region.pdf', width = 18, height = 10)
print(p1)
dev.off()

# RidgePlot(ASD.brain, features = regulonsToPlot, cols = maincol, group.by = 'Region',ncol = 3)
SpatialFeaturePlot(ASD.brain,features = regulonsToPlot, image.alpha = 0, crop = F) #& scale_fill_viridis()

p2 <- SpatialFeaturePlot(ASD.brain,features = 'Klf16(+)', image.alpha = 0, crop = F, 
                   images = c('slice1.1', 'slice1.4', 'slice1.7'), pt.size.factor = 1.2) & 
  scale_fill_gradient2(limits = c(0, 0.15), breaks = c(0, 0.05, 0.10, 0.15), low = "#440154", 
                       mid = "#21908C",high = "#FDE725", midpoint = 0.075)
pdf('./figures/scenic/SpatialFeaturePlot_scenic_Klf16.pdf', width = 15, height = 5)
print(p2)
dev.off()

p3 <- SpatialFeaturePlot(ASD.brain,features = 'Gata2(+)', image.alpha = 0, crop = F, 
                   images = c('slice1.1', 'slice1.4', 'slice1.7'), pt.size.factor = 1.2) & 
  scale_fill_gradient2(limits = c(0, 0.6), breaks = c(0, 0.2, 0.4, 0.6), low = "#440154", 
                       mid = "#21908C",high = "#FDE725", midpoint = 0.3)
pdf('./figures/scenic/SpatialFeaturePlot_scenic_Gata2.pdf', width = 18, height = 10)
print(p3)
dev.off()

p4 <- SpatialFeaturePlot(ASD.brain,features = 'Maz(+)', image.alpha = 0, crop = F, 
                         images = c('slice1.1', 'slice1.4', 'slice1.7'), pt.size.factor = 1.2) & 
  scale_fill_gradient2(limits = c(0, 0.3), breaks = c(0, 0.1, 0.2, 0.3), low = "#440154", 
                       mid = "#21908C",high = "#FDE725", midpoint = 0.15)
pdf('./figures/scenic/SpatialFeaturePlot_scenic_Maz.pdf', width = 18, height = 10)
print(p4)
dev.off()

regulonActivity_byGroup_Scaled <- t(regulonActivity_byGroup_Scaled)
hm <- pheatmap(regulonActivity_byGroup_Scaled,
               show_rownames = T,
               show_colnames = T,
               cluster_cols = T,
               cluster_rows= F,
               fontsize_row = 10,
               fontsize_col = 6,
               cols <- viridis(8),
               border_color = 'grey',
               legend = T,
               #legend_labels = 'regulon activity',
               clustering_method = 'single'
)
pdf('./figures/scenic/heatmap_scenic_group.pdf', width = 9, height = 4.5)
print(hm)
dev.off()

rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellGroups[colnames(sub_regulonAUC), 'celltype']) 
rss=na.omit(rss) 
#viridis(10, alpha = 1, begin = 0, end = 1, direction = 1, option = "A")
rssPlot <- plotRSS(rss) #, col.low = "#440154FF", col.mid = "#21908CFF", col.high = "#FDE725FF"
pdf('./figures/scenic/rssPlot.pdf', width = 4, height = 5)
print(rssPlot)
dev.off()

write.csv(rss, file = './results/scenic/scenic_rss.csv')

#gata2 heatmap
p<-DotPlot(ASD.brain, features = 'Gata2(+)', group.by = 'Region_group')
pdata <- p$data
strsplit(as.character(pdata$id), '@')[[2]][1]
pdata$regionsub <- as.vector(sapply(as.character(pdata$id), function(x) strsplit(x,'@')[[1]])[1,])
pdata$group <- as.vector(sapply(as.character(pdata$id), function(x) strsplit(x,'@')[[1]])[2,])
pdata <- pdata[,c(1,6,7)]
hdata <- dcast(pdata, group~regionsub,value.var = 'avg.exp')
row.names(hdata) <- hdata$group
hdata <- hdata[,-1]
pheatmap(hdata,
         show_rownames = T,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows= F,
         fontsize_row = 10,
         fontsize_col = 6,
         cols <- viridis(8),
         border_color = 'grey',
         legend = T)


list <- as.vector(sapply(as.character(pdata$id), function(x) strsplit(x,'@')[[1]])[1,])

# VlnPlot(ASD.brain, features = 'Gata2(+)', group.by = 'RegionSub_group')
#c('Rin3', 'Tsc22d3', 'Ntsr2', 'Igfbp7')
VlnPlot(ASD.brain, features = "Sox17", pt.size = 0.2, group.by = 'group',
        cols = pal_npg()(3)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  stat_compare_means(comparisons = list(c('Autism', 'Normal'), c('Autism', 'Control')),label = 'p.signif')+
  ylim(-0.1, 10)

#gata2 genes heatmap
p<-DotPlot(ASD.brain, features = Gata2.genes, group.by = 'Sample')
pdata <- p$data
#pdata <- pdata[,c(1,3,4)]
pdata$id <- paste(substr(pdata$id, 8,8), substr(pdata$id, 11,11), sep = '')
# list_group <- list(Autism = c('A1', 'A2', 'A3'),
#                    Normal = c('N1', 'N2', 'N2'),
#                    Control = c('C1', 'C2', 'C3'))
# 
# DotPlot(ASD.brain, features = 'Gata2(+)', cluster.idents = T, group.by = 'Sample')
pdata$group <- substr(pdata$id, 1,1)

intersect(Gata2.genes,unique(Markers_AN_merge$gene))

pdata$id <- factor(pdata$id, levels = c('A1', 'A2', 'A3', 'N1', 'N2', 'N3', 'C1', 'C2', 'C3'))
pdata$features.plot

gata_df_hm <- dcast(pdata,id~features.plot,value.var = 'avg.exp.scaled')
row.names(gata_df_hm) <- gata_df_hm$id
gata_df_hm <- gata_df_hm[,-1]

p<-pheatmap(t(gata_df_hm),
            show_rownames = T,
            show_colnames = T,
            cluster_cols = F,
            cluster_rows= F,
            fontsize_row = 8,
            fontsize_col = 10,
            cols <- viridis(64),
            border_color = 'grey',
            legend = T,
            angle_col = '90',
            #legend_labels = 'regulon activity',
            clustering_method = 'single'
) 



pdata$features.plot <- factor(pdata$features.plot, levels = rev(sort(as.character(unique(pdata$features.plot)))))

pdata$group <- factor(pdata$group, levels = c('A', 'N', 'C'))

p <- ggplot(pdata,aes(x=id,y=features.plot))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_bw()+labs(x = NULL, y = NULL)+
  theme(panel.grid = element_blank())+
  facet_wrap(~ group, scales = 'free_x') +
  #labs(x = "Sample", y = "Gene", fill = "Expression")+
  scale_color_gradientn(colours = c("#330066","#336699","#66CC66","#FFCC33"))

pdf('./figures/scenic/Dotplot_Gata2genes.pdf', width = 6.5, height = 4)
print(p)
dev.off()


library(igraph)
library(ggraph)
edges <- data.frame(
  from = c("Gata2"),
  to = Gata2.genes,
  weight = c(2.2150758277183638, 1.8327814225564607, 1.8064873153958, 2.227560265553162, 2.006536971497071, 
             3.286460451299887,3.6716947899891834, 35.25542525983002, 4.856995901887817) 
)

edges <- data.frame(
  from = c("Gata2"),
  to = Gata2.genes,
  weight = c(rep(1,9)) 
)

g <- graph_from_data_frame(d = edges, directed = FALSE)
V(g)$color <- c('red') 
V(g)$size <- 30 
V(g)$label.cex <- 1.5
V(g)$group <- c("top", rep("down", 9))  

plot(g, 
     vertex.label = V(g)$name,
     vertex.color = V(g)$color, 
     vertex.size = V(g)$size, 
     vertex.label.cex = V(g)$label.cex, 
     edge.width = E(g)$weight, 
     edge.color = "grey", 
     main = "Gene Interaction Network") 

p<-ggraph(g, layout = 'stress') + 
  geom_edge_link(aes(width = weight), color = "black") + 
  geom_node_point(aes(color = group), size = 5,) + 
  geom_node_text(aes(label = name), vjust = 1.5) + 
  scale_edge_width(range = c(0.2, 2)) + 
  theme_graph() + 
  ggtitle("Gata2 Target Genes")+
  scale_color_manual(values = c("top" = "#E64B35", "down" = "#4DBBD5"))

pdf('./figures/ST/scenic/Gata2_Target_Genes.pdf', width = 4.5, height = 4)
print(p)
dev.off()


####SPOTlight####

library(Matrix)
library(data.table)
library(Seurat)
library(dplyr)
library(SPOTlight.old)
library(igraph)
library(RColorBrewer)
library(ggplot2)

Idents(ASD.sce.qc) <- ASD.sce.qc$celltype_main
cluster_markers_all <- ASD.sce.qc@misc$cluster_markers_all <- Seurat::FindAllMarkers(object = ASD.sce.qc, 
                                                                                      assay = "SCT",
                                                                                      slot = "data",
                                                                                      verbose = TRUE, 
                                                                                      only.pos = TRUE, 
                                                                                      logfc.threshold = 1,
                                                                                      min.pct = 0.9)

spotlight_ls <- spotlight_deconvolution(
  se_sc = ASD.sce.qc,
  counts_spatial = ASD.brain@assays$Spatial@counts,
  clust_vr = "celltype_main", 
  cluster_markers = cluster_markers_all, 
  cl_n = 100, 
  hvg = 1000, 
  ntop = 50, 
  transf = "uv", 
  method = "nsNMF",
  min_cont = 0 
)

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
# decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(ASD.brain)
cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]



decon_df <- decon_mtrx %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

ASD.brain@meta.data <- ASD.brain@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")
save(decon_df, decon_mtrx, decon_mtrx_sub, file = './object/spotlight/spotlight_decon.RData')


SPOTlight.old::spatial_scatterpie(se_obj = ASD.brain,
                                  cell_types_all = cell_types_all,
                                  img_path = "./data/ST/A1_Visium/spatial/tissue_lowres_image.png",
                                  pie_scale = 0.4)


pdf('./sc/spotlight/A2_merge.pdf', width = 15, height = 7.5)
print(sp)
dev.off()

SPOTlight.old::spatial_scatterpie(se_obj = ASD.brain,
                                  cell_types_all = cell_types_all,
                                  img_path = "./data/ST/A1_2/spatial/tissue_lowres_image.png",
                                  cell_types_interest = 'Schwann.cells',
                                  pie_scale = 0.4)

sccol <- c("#7E6148FF","#4DBBD5FF","#5050FFFF","#3C5488FF","#B09C85FF",
           "#F39B7FFF","#E64B35FF","#00A087FF","#8491B4FF","#91D1C2FF")



sccol <- c("#ea5c6f","#f7905a","#e187cb","#fb948d","#e2b159","#ebed6f",
           "#7ee7bb","#64cccf","#a9dce6","#e4b7d6")


sp <- SPOTlight.old::spatial_scatterpie(se_obj = ASD.brain,
                                        cell_types_all = cell_types_all,
                                        img_path = "./figures/blank.png",
                                        slice = 'slice1.1',
                                        pie_scale = 0.4)+scale_fill_manual(values = sccol)

pdf('./figures/ST/spotlight/spotlight_A2.pdf', width = 7.5, height = 5)
print(sp)
dev.off()

sp <- SPOTlight.old::spatial_scatterpie(se_obj = ASD.brain,
                                        cell_types_all = cell_types_all,
                                        img_path = "./figures/blank.png",
                                        slice = 'slice1.4',
                                        pie_scale = 0.4)+scale_fill_manual(values = sccol)

pdf('./figures/ST/spotlight/spotlight_N2.pdf', width = 7.5, height = 5)
print(sp)
dev.off()


DotPlot(ASD.brain, features = cell_types_all, group.by = 'group')

spotmeta <- ASD.brain@meta.data[,c(18,19,21,38:47)]

spotmeta_GRS <- ASD.brain@meta.data[,c(21,38:47)]
spotmeta_GRS  <- aggregate(.~GroupRegionSub, mean, data=spotmeta_GRS)

spotmeta_Group <- ASD.brain@meta.data[,c(18,38:47)]
spotmeta_Group  <- aggregate(.~group, mean, data=spotmeta_Group)

spotmeta_region <- ASD.brain@meta.data[,c(19,38:47)]
spotmeta_region  <- aggregate(.~RegionAll, mean, data=spotmeta_region)

rowSums(spotmeta_region[,c(2:11)])

save(ASD.brain, file = './object/ASD.brian.spotlight.RData')

ASD.brain$GroupRegionSub

save(spotmeta, spotmeta_Group, spotmeta_GRS, spotmeta_region, file = './object/spotmeta.objects.RData')

spotmeta_region1 <- melt(spotmeta_region)

spotmeta_GRS1 <- melt(spotmeta_GRS)
xx <- as.data.frame(t(as.data.frame(strsplit(spotmeta_GRS1$GroupRegionSub, '[.]'))))

spotmeta_GRS1$region <- xx$V2
spotmeta_GRS1$group <- factor(xx$V1, levels = c('Autism', 'Normal', 'Control'))

spotmeta_GRS1$group <- gsub('Autism', 'Positive', spotmeta_GRS1$group)
spotmeta_GRS1$group <- gsub('Normal', 'Negative', spotmeta_GRS1$group)
spotmeta_GRS1$group <- factor(spotmeta_GRS1$group, levels = c('Positive', 'Negative', 'Control'))

p<-ggplot(data = spotmeta_GRS1, #[which(spotmeta_GRS1$group != 'Control'),]
          aes(x=group, y=value, fill=variable))+
  geom_bar(stat = "identity",position = "stack")+scale_fill_manual(values = sccol)+theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, colour = 'black'),
        axis.text.y = element_text(colour = 'black'),
        legend.position = 'none')+
  facet_grid(~region)

pdf('./figures/ST/spotlight/AllREG_spotlight_sn.pdf', width = 15, height = 5)
print(p)
dev.off()

spotmeta_Group1 <- melt(spotmeta_Group)
xx <- as.data.frame(t(as.data.frame(strsplit(spotmeta_Group1$GroupRegionSub, '[.]'))))

spotmeta_Group1$region <- xx$V2
spotmeta_Group1$group <- factor(xx$V1, levels = c('Autism', 'Normal', 'Control'))


p<-ggplot(data = spotmeta_Group1, #[which(spotmeta_GRS1$group != 'Control'),]
          aes(x=group, y=value, fill=variable))+labs(x = '', y = 'Fraction of cells')+
  geom_bar(stat = "identity",position = "stack")+scale_fill_manual(values = sccol)+theme_classic()+
  theme(axis.text.x = element_text(size = 12,color = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.y = element_text(size = 12))

pdf('./figures/ST/spotlight/AllREG_spotlight_group_sn.pdf', width = 4, height = 5)
print(p)
dev.off()

####hdWGCNA####
library(hdWGCNA)
library(proxy)
enableWGCNAThreads(nThreads = 8)
set.seed(123)

# image_df <- do.call(rbind, lapply(names(ASD.brain@images), function(x){
#   ASD.brain@images[[x]]@coordinates
# }))

modules_list <- list()
for (i in 4:9) {
  i = 3
  sce <- Autism[[i]]
  sce$imagerow <- sce@images$slice1@coordinates$imagerow
  sce$imagecol <- sce@images$slice1@coordinates$imagecol
  sce$row <- sce@images$slice1@coordinates$row
  sce$col <- sce@images$slice1@coordinates$col
  
  sce <- SetupForWGCNA(sce, gene_select = "fraction", fraction = 0.05, wgcna_name = "vis")
  
  sce <- MetaspotsByGroups(sce, group.by = c("Region"), ident.group = "Region", assay = 'Spatial', slot = 'counts')
  sce  <- NormalizeMetacells(sce)
  # m_obj <- GetMetacellObject(sce)
  
  sce  <- SetDatExpr(sce, assay = 'Spatial', group.by=NULL, group_name = NULL)
  sce <- TestSoftPowers(sce)
  plot_list <- PlotSoftPowers(sce)
  wrap_plots(plot_list, ncol=2)
  
  sce <- ConstructNetwork(sce, tom_name='test', overwrite_tom=TRUE, soft_power=9)
  PlotDendrogram(sce, main='Spatial hdWGCNA dendrogram')
  
  sce <- ModuleEigengenes(sce)
  sce <- ModuleConnectivity(sce)
  sce <- ResetModuleNames(sce, new_name = "SM")
  modules <- GetModules(sce) %>% subset(module != 'grey')
  modules_list[[i]] <- modules
  print(paste(i, 'is over', sep = ' '))
}

save(modules_list, file = './object/modules_list.RData')
save(plot_list, file = './object/WGCNA_plot_list.RData')

filtered_list <- lapply(modules_list, function(df) {
  subset(df, gene_name == "Igfbp7")[,c(1:3)]
})
combined_df <- do.call(rbind, filtered_list)

igfbp7_list <- list()
for (i in 1:9) {
  igfbp7_list[[i]] <- subset(modules_list[[i]], color == combined_df$color[i])$gene_name
}

save(igfbp7_list, file = './object/igfbp7_list.RData')
common_elements <- Reduce(intersect, igfbp7_list[1:3])

venn_list <- igfbp7_list[1:3]
names(venn_list) <- c('A1', 'A2', 'A3')
library(ggvenn)
p<-ggvenn(
  data = venn_list,         
  columns = NULL,           
  show_elements = F,        
  label_sep = "\n",         
  show_percentage = T,      
  digits = 1,               
  fill_color = pal_npg()(3), 
  fill_alpha = 0.7,         
  stroke_color = "white",   
  stroke_alpha = 0.5,       
  stroke_size = 0.5,        
  stroke_linetype = "blank", 
  set_name_color = "black", 
  set_name_size = 6,        
  text_color = "black",     
  text_size = 5             
)
pdf('./figures/ST/hdWGCNA/Venn_STWGCNA_igfbp7_intersects.pdf', width = 4, height = 5)
print(p)
dev.off()



SpatialFeaturePlot(Autism[[5]], features = 'Mapt', pt.size.factor = 3)
SpatialFeaturePlot(Autism[[3]], features = 'B2m', pt.size.factor = 3)

DotPlot(ASD.brain, features = common_elements, group.by = 'group')

# common_elements <- common_elements[-15]#delete igfbp7
WGCNA <- data.frame(wgcna = common_elements)
WGCNA <- as.list(WGCNA)
ASD.brain <- AddModuleScore(ASD.brain, features = WGCNA, ctrl = 100, name = 'WGCNA')
DotPlot(ASD.brain, features = 'WGCNA1', group.by = 'group')
VlnPlot(ASD.brain, features = 'WGCNA1', group.by = 'group', sort = T)
DotPlot(ASD.brain, features = WGCNA, group.by = 'group')
SpatialFeaturePlot(ASD.brain, features = 'WGCNA1', pt.size.factor = 3, ncol = 3)

WGCNA_int <- intersect(unique(Markers_AN_merge[which(Markers_AN_merge$avg_log2FC>0),]$gene), WGCNA$wgcna)
DotPlot(ASD.brain, features = WGCNA_int, group.by = 'group')

write.table(WGCNA$wgcna, file = './results/wgcna/WGCNA_Spatial_intersect.txt',sep = "\t",
            quote = F, col.names = F,row.names = F)

save(SFARI, SFARI_gene, SFARI_DEG, file = './object/sfari.RData')

save(igfbp7_list, SFARI_gene, SFARI_DEG, file = './object/sfari.RData')

library(ggrepel)
library(ggpubr)
tmp <- data.frame(igfbp7 = ASD.brain@assays$SCT@data %>% .['Igfbp7',], 
                  WGCNA = ASD.brain@assays$SCT@data %>% .['Mt1',])
# , group = ASD.brain$group
# tmp <- subset(tmp, group == 'Autism')
ggscatter(tmp ,x = 'igfbp7', y = 'WGCNA',
                add = "reg.line", conf.int = TRUE,color = col1[1],
                add.params = list(fill = col1[20]),
                ggtheme = theme_minimal()
)+ stat_cor(method = "spearman",
            color='black',
            p.accuracy = 0.001
)+theme_bw()+theme(panel.grid.major = element_blank() ,axis.text.x = element_text(size = 12,color = "black",angle = 0,hjust = 1,vjust = 1),
                   axis.text.y = element_text(size = 12,color = "black"),
                   axis.title = element_text(size = 15,color = 'black'),
                   legend.text = element_text(size = 15,color = "black"))                         

######STRING RESULTS######
String_GOBP <- read.table('./results/sc Igfbp7 exp/enrichment.Process.tsv', sep="\t")
String_REACTOME <- read.table('./results/sc Igfbp7 exp/enrichment.RCTM.tsv', sep="\t")
String_GOCC <- read.table('./results/sc Igfbp7 exp/enrichment.Component.tsv', sep="\t")
String_KEGG <- read.table('./results/sc Igfbp7 exp/enrichment.KEGG.tsv', sep="\t")

colnames(String_GOBP) <- c('ID','Description',	'count',	'background gene count','strength',
                           'signal', 'FDR', 'matching proteins in your network (IDs)',	
                           'matching proteins in your network (labels)')
String_GOBP <- String_GOBP[,c(2,3,6,7)]
String_GOBP$`-log10(FDR)` <- -log10(String_GOBP$FDR)
String_GOBP <- String_GOBP %>% top_n(n = 5, wt = signal)
String_GOBP <- String_GOBP[order(String_GOBP$signal),]
String_GOBP$Description <- factor(String_GOBP$Description, levels = String_GOBP$Description)
p<-ggplot(String_GOBP, aes(x = Description, y = signal, fill = `-log10(FDR)`)) +
  geom_col() +
  coord_flip() +  
  scale_fill_gradient(low = "#6BAED6", high = "#084594", name = "-log10(FDR)") + 
  labs(x = "", y = "signal") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 14, angle = 0, colour = 'black'),
    axis.text.x = element_text(size = 12, angle = 0, colour = 'black'),
    axis.title.x = element_text(size = 14, angle = 0, colour = 'black')
  )
pdf('./figures/ST/hdWGCNA/STRING_IGFBP7_GOBP.pdf', width = 11, height = 3)
print(p)
dev.off()


colnames(String_REACTOME) <- c('ID','Description',	'count',	'background gene count','strength',
                           'signal', 'FDR', 'matching proteins in your network (IDs)',	
                           'matching proteins in your network (labels)')
String_REACTOME <- String_REACTOME[,c(2,3,6,7)]
String_REACTOME$`-log10(FDR)` <- -log10(String_REACTOME$FDR)
String_REACTOME <- String_REACTOME %>% top_n(n = 5, wt = signal)
String_REACTOME <- String_REACTOME[order(String_REACTOME$signal),]
String_REACTOME$Description <- factor(String_REACTOME$Description, levels = String_REACTOME$Description)
p<-ggplot(String_REACTOME, aes(x = Description, y = signal, fill = `-log10(FDR)`)) +
  geom_col() +
  coord_flip() +  
  scale_fill_gradient(low = "#6BAED6", high = "#084594", name = "-log10(FDR)") +  
  labs(x = "", y = "signal") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 14, angle = 0, colour = 'black'),
    axis.text.x = element_text(size = 12, angle = 0, colour = 'black'),
    axis.title.x = element_text(size = 14, angle = 0, colour = 'black')
  )
pdf('./figures/ST/hdWGCNA/STRING_IGFBP7_REACT.pdf', width = 14, height = 2.5)
print(p)
dev.off()

colnames(String_GOCC) <- c('ID','Description',	'count',	'background gene count','strength',
                           'signal', 'FDR', 'matching proteins in your network (IDs)',	
                           'matching proteins in your network (labels)')
String_GOCC <- String_GOCC[,c(2,3,6,7)]
String_GOCC$`-log10(FDR)` <- -log10(String_GOCC$FDR)
String_GOCC <- String_GOCC %>% top_n(n = 5, wt = signal)
String_GOCC <- String_GOCC[order(String_GOCC$signal),]
String_GOCC$Description <- factor(String_GOCC$Description, levels = String_GOCC$Description)
p<-ggplot(String_GOCC, aes(x = Description, y = signal, fill = `-log10(FDR)`)) +
  geom_col() +
  coord_flip() +  
  scale_fill_gradient(low = "#6BAED6", high = "#084594", name = "-log10(FDR)") +  
  labs(x = "", y = "signal") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 14, angle = 0, colour = 'black'),
    axis.text.x = element_text(size = 12, angle = 0, colour = 'black'),
    axis.title.x = element_text(size = 14, angle = 0, colour = 'black')
  )
pdf('./figures/ST/hdWGCNA/STRING_IGFBP7_GOCC.pdf', width = 9, height = 3)
print(p)
dev.off()


colnames(String_KEGG) <- c('ID','Description',	'count',	'background gene count','strength',
                           'signal', 'FDR', 'matching proteins in your network (IDs)',	
                           'matching proteins in your network (labels)')
String_KEGG <- String_KEGG[,c(2,3,6,7)]
String_KEGG$`-log10(FDR)` <- -log10(String_KEGG$FDR)
String_KEGG <- String_KEGG %>% top_n(n = 5, wt = signal)
String_KEGG <- String_KEGG[order(String_KEGG$signal),]
String_KEGG$Description <- factor(String_KEGG$Description, levels = String_KEGG$Description)
p<-ggplot(String_KEGG, aes(x = Description, y = signal, fill = `-log10(FDR)`)) +
  geom_col() +
  coord_flip() +  
  scale_fill_gradient(low = "#6BAED6", high = "#084594", name = "-log10(FDR)") +  
  labs(x = "", y = "signal") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 14, angle = 0, colour = 'black'),
    axis.text.x = element_text(size = 12, angle = 0, colour = 'black'),
    axis.title.x = element_text(size = 14, angle = 0, colour = 'black')
  )
pdf('./figures/ST/hdWGCNA/STRING_IGFBP7_KEGG.pdf', width = 7.5, height = 2.5)
print(p)
dev.off()

####cellchat v2####
# devtools::install_github("sqjin/CellChat")
library(CellChat)


cellchat_list <- list()
sample <- unique(ASD.brain$Sample)
dirname <- c('A1_Visium', 'A1_2', 'A1_3', 'N1_Visium', 'N1_2', 'N1_3', 'C3_1', 'C3_2', 'C3_3')
slices <- c('slice1', 'slice1.1', 'slice1.2', 'slice1.3', 'slice1.4', 'slice1.5',
            'slice1.6', 'slice1.7', 'slice1.8')

for (i in 1:9) {
  Brain_ST <- subset(ASD.brain, Sample == sample[i])
  data.input = Seurat::GetAssayData(Brain_ST, slot = "data", assay = "SCT")
  meta = data.frame(labels = Brain_ST$RegionAll, #
                    row.names = names(Idents(Brain_ST))) # manually create a dataframe consisting of the cell labels
  unique(meta$labels)

  spatial.locs = Seurat::GetTissueCoordinates(Brain_ST[[slices[i]]], scale = NULL,
                                              cols = c("imagerow", "imagecol"))
  scale.factors = jsonlite::fromJSON(txt =
                                       file.path(paste('./data/ST/', dirname[i],'/spatial', sep = ''), 'scalefactors_json.json'))
  scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                       fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
  )
  cellchat <- createCellChat(object = data.input,
                             meta = meta,
                             group.by = "labels", #
                             datatype = "spatial", ###
                             coordinates = spatial.locs,
                             scale.factors = scale.factors)

  cellchat@DB <- CellChatDB.mouse

  cellchat <- subsetData(cellchat)
  # future::plan("multisession", workers = 24)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,
                                type = "truncatedMean", trim = 0.1,
                                distance.use = TRUE,
                                scale.distance = 0.01)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  cellchat_list[[i]] <- cellchat
  cellchat_list[[i]] <- netAnalysis_computeCentrality(cellchat_list[[i]])
  print(paste(sample[i], 'is over'))
}


cellchat <- cellchat_list[[3]]

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

p1 + p2

cellchat@data.signaling
pathways.show <- c("APP")
levels(cellchat@idents)   
vertex.receiver = c(1,2,3,4,5,6,7)  
netVisual_aggregate(cellchat, signaling = pathways.show,                      
                    vertex.receiver = vertex.receiver,layout = "hierarchy")

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", 
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 3.5)


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)

netVisual_bubble(cellchat, sources.use = c(1), 
                 targets.use = c(1,2,4,6), remove.isolate = FALSE)

# "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF"
# Take an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(cellchat_list[[2]], pairLR.use = "APP_CD74", point.size = 1.5, 
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                   color.heatmap = "viridis", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))

spatialFeaturePlot(cellchat_list[[5]], pairLR.use = "APP_CD74", point.size = 1.5, 
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                   color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))

spatialFeaturePlot(cellchat_list[[8]], pairLR.use = "APP_CD74", point.size = 1.5, 
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                   color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))
                        
netVisual_aggregate(cellchat_list[[2]], signaling = 'APP', layout = "chord", 
                    edge.width.max = 2, vertex.size.max = 1, 
                    alpha.image = 0.2, vertex.label.cex = 1)


object.list <- list()
object.list[[1]] <- cellchat_list[[3]]
object.list[[2]] <- cellchat_list[[6]]
object.list[[3]] <- cellchat_list[[9]]
names(object.list) <- c('Autism', 'Normal', 'Control')
cellchat.S3 <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

pairLR.use <- extractEnrichedLR(cellchat.S3, signaling = c('NRG', 'PACAP', 'MIF', 'APP', 'OPIOID', 'PSAP', 'SEMA3',
                                                           'THBS', 'SEMA4', 'EPHB', 'NOTCH', 'GAS', 'WNT', 'PDGF', 'CDH', 
                                                           'IGF', 'ENHO', 'GRN', 'FN1'))
p<-netVisual_bubble(cellchat.S3, sources.use = c(1,2,3,4), angle.x = 45,pairLR.use = pairLR.use,
                 targets.use = c(1,2,3,4), remove.isolate = T, comparison = c(1,2))
pdf('./figures/ST/cellchatv2/netVisual_bubble_AN_S3.pdf', width = 7, height = 5)
print(p)
dev.off()



netVisual_bubble(cellchat.S3, sources.use = c(1,2,3,4), angle.x = 45,
                 targets.use = c(1,2,5,6), remove.isolate = T, comparison = c(1,2))

netVisual_bubble(cellchat.ALL, sources.use = c(1,2,3,4), angle.x = 45,
                 targets.use = c(1,2,3,4), remove.isolate = T, comparison = c(7,8))


rankNet(cellchat.ALL, mode = "comparison", comparison = c(7,8), stacked = T, do.stat = TRUE,
        color.use = pal_npg()(3), do.flip = T, measure = 'count')



object.list <- list()
object.list[[1]] <- cellchat_list[[1]]
object.list[[2]] <- cellchat_list[[4]]
object.list[[3]] <- cellchat_list[[7]]
object.list[[4]] <- cellchat_list[[2]]
object.list[[5]] <- cellchat_list[[5]]
object.list[[6]] <- cellchat_list[[8]]
object.list[[7]] <- cellchat_list[[3]]
object.list[[8]] <- cellchat_list[[6]]
object.list[[9]] <- cellchat_list[[9]]
names(object.list) <- c('A1', 'N1', 'C1', 'A2', 'N2', 'C2', 'A3', 'N3', 'C3')
cellchat.ALL <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

p1<-compareInteractions(cellchat.ALL, show.legend = F, group = c(1:9),
                    color.use = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#E64B35FF", "#4DBBD5FF", "#00A087FF",
                                  "#E64B35FF", "#4DBBD5FF", "#00A087FF"))
pdf('./figures/ST/cellchatv2/compareInteractions_number_all.pdf', width = 6, height = 3)
print(p1)
dev.off()


p2<-compareInteractions(cellchat.ALL, show.legend = F, group = c(1:9), measure = "weight",
                    color.use = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#E64B35FF", "#4DBBD5FF", "#00A087FF",
                                  "#E64B35FF", "#4DBBD5FF", "#00A087FF"))
pdf('./figures/ST/cellchatv2/compareInteractions_weight_all.pdf', width = 6, height = 3)
print(p2)
dev.off()


p3<-rankNet(cellchat.ALL, mode = "comparison", comparison = c(1,2,3), stacked = T, do.stat = TRUE,
        color.use = pal_npg()(3), do.flip = T, measure = 'count')
pdf('./figures/ST/cellchatv2/rankNet_Section1.pdf', width = 5, height = 5)
print(p3)
dev.off()

p4<-rankNet(cellchat.ALL, mode = "comparison", comparison = c(4,5,6), stacked = T, do.stat = TRUE,
            color.use = pal_npg()(3), do.flip = T, measure = 'count')
pdf('./figures/ST/cellchatv2/rankNet_Section2.pdf', width = 5, height = 5)
print(p4)
dev.off()

p5<-rankNet(cellchat.ALL, mode = "comparison", comparison = c(7,8,9), stacked = T, do.stat = TRUE,
            color.use = pal_npg()(3), do.flip = T, measure = 'count')
pdf('./figures/ST/cellchatv2/rankNet_Section3.pdf', width = 5, height = 5)
print(p5)
dev.off()

p6 <- netVisual_aggregate(cellchat_list[[2]], signaling = 'APP', layout = "spatial", 
                          edge.width.max = 2, vertex.size.max = 1, 
                          alpha.image = 0.2, vertex.label.cex = 3.5)
p7 <- netVisual_aggregate(cellchat_list[[3]], signaling = 'APP', layout = "spatial", 
                          edge.width.max = 2, vertex.size.max = 1, 
                          alpha.image = 0.2, vertex.label.cex = 3.5)
p8 <- netVisual_aggregate(cellchat_list[[1]], signaling = 'APP', layout = "spatial", 
                          edge.width.max = 2, vertex.size.max = 1, 
                          alpha.image = 0.2, vertex.label.cex = 3.5)

pdf('./figures/ST/cellchatv2/netVisual_aggregate_ST_Autism.pdf', width = 15, height = 5)
print(p6+p7+p8)
dev.off()

p9 <- spatialFeaturePlot(cellchat_list[[2]], pairLR.use = "APP_CD74", point.size = 1.5, 
                     do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                     color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))
p10 <- spatialFeaturePlot(cellchat_list[[5]], pairLR.use = "APP_CD74", point.size = 1.5, 
                         do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                         color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))
p11 <- spatialFeaturePlot(cellchat_list[[8]], pairLR.use = "APP_CD74", point.size = 1.5, 
                          do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                          color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))

pdf('./figures/ST/cellchatv2/spatialFeaturePlot_APP_section2.pdf', width = 15, height = 5)
print(p9+p10+p11)
dev.off()


p12 <- spatialFeaturePlot(cellchat_list[[1]], pairLR.use = "APP_CD74", point.size = 1.5, 
                         do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                         color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))
p13 <- spatialFeaturePlot(cellchat_list[[4]], pairLR.use = "APP_CD74", point.size = 1.5, 
                          do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                          color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))
p14 <- spatialFeaturePlot(cellchat_list[[7]], pairLR.use = "APP_CD74", point.size = 1.5, 
                          do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                          color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))

pdf('./figures/ST/cellchatv2/spatialFeaturePlot_APP_section1.pdf', width = 15, height = 5)
print(p12+p13+p14)
dev.off()

p15 <- spatialFeaturePlot(cellchat_list[[3]], pairLR.use = "APP_CD74", point.size = 1.5, 
                          do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                          color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))
p16 <- spatialFeaturePlot(cellchat_list[[6]], pairLR.use = "APP_CD74", point.size = 1.5, 
                          do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                          color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))
p17 <- spatialFeaturePlot(cellchat_list[[9]], pairLR.use = "APP_CD74", point.size = 1.5, 
                          do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                          color.heatmap = "Reds", direction = 1, color.use = c("#E18727FF", "#20854EFF", "#0072B5FF", 'grey'))

pdf('./figures/ST/cellchatv2/spatialFeaturePlot_APP_section3.pdf', width = 15, height = 5)
print(p15+p16+p17)
dev.off()


netVisual_bubble(cellchat_list[[6]], sources.use = c(1,2,3,4,8,9,10), 
                 targets.use = c(1), remove.isolate = FALSE)


pairLR.use <- extractEnrichedLR(cellchat.S3, signaling = c('NRG', 'PACAP', 'MIF', 'APP', 'OPIOID', 'PSAP', 'SEMA3',
                                                           'THBS', 'SEMA4', 'EPHB', 'NOTCH', 'GAS', 'WNT', 'PDGF', 'CDH', 
                                                           'IGF', 'ENHO', 'GRN', 'FN1'))
p<-netVisual_bubble(cellchat.S3, sources.use = c(1,2,3,4), angle.x = 45,pairLR.use = pairLR.use,
                    targets.use = c(1,2,3,4), remove.isolate = T, comparison = c(1,2))
pdf('./figures/ST/cellchatv2/netVisual_bubble_AN_S3.pdf', width = 7, height = 5)
print(p)
dev.off()


object.list <- list()
object.list[[1]] <- cellchat_list[[1]]
object.list[[2]] <- cellchat_list[[4]]
object.list[[3]] <- cellchat_list[[7]]
names(object.list) <- c('Autism', 'Normal', 'Control')
cellchat.S1 <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

rankNet(cellchat.S1, mode = "comparison", comparison = c(1,2), stacked = T, do.stat = TRUE,
        color.use = pal_npg()(3), do.flip = T, measure = 'count')

pairLR.use <- extractEnrichedLR(cellchat.S1, signaling = c('GRN', 'APP', 'GAS', 'NPY', 'SEMA3', 'JAM', 'ESAM',
                                                           'THBS', 'NOTCH', 'LAMININ', 'EPHB', 'EPHA', 'OPIOID',  
                                                           'MAG'))

p<-netVisual_bubble(cellchat.S1, sources.use = c(1,2,3,4), angle.x = 45,pairLR.use = pairLR.use,
                    targets.use = c(1,2,3,4), remove.isolate = T, comparison = c(1,2))
pdf('./figures/ST/cellchatv2/netVisual_bubble_AN_S1.pdf', width = 7, height = 5)
print(p)
dev.off()

object.list <- list()
object.list[[1]] <- cellchat_list[[2]]
object.list[[2]] <- cellchat_list[[5]]
object.list[[3]] <- cellchat_list[[8]]
names(object.list) <- c('Autism', 'Normal', 'Control')
cellchat.S2 <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

rankNet(cellchat.S2, mode = "comparison", comparison = c(1,2), stacked = T, do.stat = TRUE,
        color.use = pal_npg()(3), do.flip = T, measure = 'count')

pairLR.use <- extractEnrichedLR(cellchat.S2, signaling = c('L1CAM', 'APP', 'SEMA3', 
                                                           'SEMA6', 'ncWNT', 'JAM', 'VEGF', 'FGF', #, 'PASP'
                                                           'SEMA4', 'NGL', 'NT', 'IGF',
                                                            'SEMA7', 'GRN', 'GAS', 'CD39', 'WNT', 'NPR2', 'NPR1', 'PDGF', 'ESAM'))

p<-netVisual_bubble(cellchat.S2, sources.use = c(1,2,3,4), angle.x = 45,pairLR.use = pairLR.use,
                    targets.use = c(1,2,3,4), remove.isolate = T, comparison = c(1,2))
pdf('./figures/ST/cellchatv2/netVisual_bubble_AN_S2.pdf', width = 7, height = 5)
print(p)
dev.off()
#加载包
library(ComplexHeatmap)
library(viridis)
p1 <- netVisual_heatmap(cellchat_list[[1]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p2 <- netVisual_heatmap(cellchat_list[[4]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p3 <- netVisual_heatmap(cellchat_list[[7]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p1 %v% p2 %v% p3
pdf('./figures/ST/cellchatv2/netVisual_heatmap_S1.pdf', width = 5, height = 9)
print(p1 %v% p2 %v% p3)
dev.off()

p4 <- netVisual_heatmap(cellchat_list[[2]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p5 <- netVisual_heatmap(cellchat_list[[5]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p6 <- netVisual_heatmap(cellchat_list[[8]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p4%v%p5%v%p6
pdf('./figures/ST/cellchatv2/netVisual_heatmap_S2.pdf', width = 5, height = 9)
print(p4%v%p5%v%p6)
dev.off()

p7 <- netVisual_heatmap(cellchat_list[[3]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p8 <- netVisual_heatmap(cellchat_list[[6]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p9 <- netVisual_heatmap(cellchat_list[[9]], color.heatmap = viridis(2, begin = 0.3, end = 0.7))
p7+p8
pdf('./figures/ST/cellchatv2/netVisual_heatmap_S31.pdf', width = 5, height = 6.4)
print(p7%v%p8)
dev.off()

pdf('./figures/ST/cellchatv2/netVisual_heatmap_S32.pdf', width = 5, height = 3.6)
print(p9)
dev.off()

####bayes corr####
library(SingleCellExperiment)
library(BayesSpace)
library(viridis)

markers <- list()
m1 <- FindMarkers(ASD.sce.qc, ident.1 = 'Microgila', logfc.threshold = 0.25, only.pos = T)
m2 <- FindMarkers(ASD.sce.qc, ident.1 = 'Endothelial', logfc.threshold = 0.25, only.pos = T)

markers[["SASP"]] <- SASP$SASP
markers[["MG"]] <- row.names(m1 %>% top_n(n = 20, wt = avg_log2FC))
markers[["Endo"]] <- row.names(m2 %>% top_n(n = 20, wt = avg_log2FC))


markers[["Cdkn1a"]] <- 'Cdkn1a'
markers[["Igfbp7"]] <- 'Igfbp7'
markers[["Cd74"]] <- 'Cd74'

ST <- c('A1', 'N1', 'C1',
        'A2', 'N2', 'C2',
        'A3', 'N3', 'C3')

for (i in 1:9) {
  sce <- ASD.Bys[[i]]
  sce.enhanced<-ASD.BysE[[i]]
  sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                  model="xgboost",
                                  feature_names=purrr::reduce(markers, c),
                                  nrounds=0)
  
  sum_counts <- function(sce, features) {
    if (length(features) > 1) {
      colSums(logcounts(sce)[features, ])
    } else {
      logcounts(sce)[features, ]
    }
  }
  spot_expr <- purrr::map(markers, function(xs) sum_counts(sce, xs))
  enhanced_expr <- purrr::map(markers, function(xs) sum_counts(sce.enhanced, xs))
  
  #plot function
  plot_expression <- function(sce, expr, title) {
    featurePlot(sce, expr, color=NA) +ggplot2::scale_fill_gradientn(colours = viridis(50, option = "D"))+
      labs(title=title, fill="Log-normalized\nexpression")
  }
  plot_expression_comparison <- function(cell_type) {
    spot.plot <- plot_expression(sce, 
                                 spot_expr[[cell_type]], 
                                 "Spot")
    enhanced.plot <- plot_expression(sce.enhanced,
                                     enhanced_expr[[cell_type]], 
                                     "Enhanced")
    
    (enhanced.plot) + 
      plot_annotation(title=cell_type,
                      theme=theme(plot.title=element_text(size=18)))
  }
  
  p1<-plot_expression_comparison('SASP')
  pdf(paste("./figures/ST/COR/Figure_",ST[i],"_", "SASP",".pdf",sep = ""),width =4,height = 4)
  print(p1)
  dev.off()

  for (j in names(markers)[-1]) {
    p1<-plot_expression_comparison(j)
    tmp <- data.frame(enhanced_expr$SASP,enhanced_expr[[j]])
    colnames(tmp) <-c('SASP',j)
    p2<- ggscatter( tmp ,x = 'SASP', y = j,
                    add = "reg.line", conf.int = TRUE,color = 'SASP',#col1[1],
                    add.params = list(fill = '#167CB4'),
                    ggtheme = theme_minimal()
    )+gradient_color(viridis(50, option = "G"))+
      stat_cor(method = "pearson",
               color='black',
               p.accuracy = 0.001,
               size=6
      )+theme_bw()+theme(panel.grid.major = element_blank() ,legend.position = "none",
                         axis.text.x = element_text(size = 12,color = "black",angle = 0,hjust = 1,vjust = 1),
                         axis.text.y = element_text(size = 12,color = "black"),
                         axis.title = element_text(size = 15,color = 'black'),
                         legend.text = element_text(size = 15,color = "black"))+labs(x="SASP",y=j,title = ST[i])
    pdf(paste("./figures/ST/COR/Figure_",ST[i],"_",j,".pdf",sep = ""),width =4,height = 4)
    print(p1)
    dev.off()
    pdf(paste("./figures/ST/COR/Figure_",ST[i],"_",j,"_cor_SASP.pdf",sep = ""),width =4,height = 4)
    print(p2)
    dev.off()
  }
  print(paste(ST[i], 'is over', sep = ' '))
}

exp.sce.enhanced <- list()
####
for (i in 1:9) {
  sce <- ASD.Bys[[i]]
  sce.enhanced<-ASD.BysE[[i]]
  sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                  model="xgboost",
                                  feature_names=purrr::reduce(markers, c),
                                  nrounds=0)
  
  sum_counts <- function(sce, features) {
    if (length(features) > 1) {
      colSums(logcounts(sce)[features, ])
    } else {
      logcounts(sce)[features, ]
    }
  }
  spot_expr <- purrr::map(markers, function(xs) sum_counts(sce, xs))
  enhanced_expr <- purrr::map(markers, function(xs) sum_counts(sce.enhanced, xs))
  
  #plot function
  plot_expression <- function(sce, expr, title) {
    featurePlot(sce, expr, color=NA) +ggplot2::scale_fill_gradientn(colours = viridis(50, option = "D"))+
      labs(title=title, fill="Log-normalized\nexpression")
  }
  plot_expression_comparison <- function(cell_type) {
    spot.plot <- plot_expression(sce, 
                                 spot_expr[[cell_type]], 
                                 "Spot")
    enhanced.plot <- plot_expression(sce.enhanced,
                                     enhanced_expr[[cell_type]], 
                                     "Enhanced")
    
    (enhanced.plot) + 
      plot_annotation(title=cell_type,
                      theme=theme(plot.title=element_text(size=18)))
  }
  
  for (j in names(markers)[-3]) {
    # p1<-plot_expression_comparison(j)
    tmp <- data.frame(enhanced_expr$Cd74,enhanced_expr[[j]])
    colnames(tmp) <-c('Cd74',j)
    p2<- ggscatter( tmp ,x = 'Cd74', y = j,
                    add = "reg.line", conf.int = TRUE,color = 'Cd74',#col1[1],
                    add.params = list(fill = '#167CB4'),
                    ggtheme = theme_minimal()
    )+gradient_color(viridis(50, option = "G"))+
      stat_cor(method = "pearson",
               color='black',
               p.accuracy = 0.001,
               size=6
      )+theme_bw()+theme(panel.grid.major = element_blank() ,legend.position = "none",
                         axis.text.x = element_text(size = 12,color = "black",angle = 0,hjust = 1,vjust = 1),
                         axis.text.y = element_text(size = 12,color = "black"),
                         axis.title = element_text(size = 15,color = 'black'),
                         legend.text = element_text(size = 15,color = "black"))+labs(x="Cd74",y=j,title = ST[i])

    pdf(paste("./figures/ST/COR/Addition_",ST[i],"_",j,"_cor_Cd74.pdf",sep = ""),width =4,height = 4)
    print(p2)
    dev.off()
  }
  exp.sce.enhanced[[i]] <- sce.enhanced
  print(paste(ST[i], 'is over', sep = ' '))
}


# MG_A <- c('H2-Aa', 'H2-Ab1', 'H2-Eb1', 'H2-Ea', 'Cd74')
MG_A <- as.list(as.data.frame(c('H2-Aa', 'H2-Ab1', 'H2-Eb1', 'H2-Ea', 'Cd74', 'Cd68')))
names(MG_A) <- 'MG_A'
ASD.brain <- AddModuleScore(ASD.brain, features = MG_A, ctrl = 100, name = 'MG_A')

p<-VlnPlot(ASD.brain, features = "MG_A1", pt.size = 0, group.by = 'group',sort = T,
           cols = pal_npg()(4))+
  stat_compare_means(comparisons = list(c('Autism', 'Normal'), c('Autism', 'Control')),label = 'p.signif')+
  ylim(-0.1, 3)+NoLegend()
pdf('./figures/ST/SASP/Vlnplot_MGA.pdf', width = 4, height = 6)
print(p)
dev.off()




VlnPlot(ASD.brain, features = "Cd68", pt.size = 0, group.by = 'group',sort = T,
        cols = pal_npg()(4))+
  stat_compare_means(comparisons = list(c('Autism', 'Normal')))+
  ylim(-0.1, 3)+NoLegend()


#up
region_heatmap('Gfap')
region_heatmap('Atf3')
region_heatmap('Gap43')
region_heatmap('Sprr1a')
#down
region_heatmap('Rbfox3')
region_heatmap('Snap25')
region_heatmap('Bax')

region_heatmap('App')
region_heatmap('Cd74')

region_heatmap('Cd68')
region_heatmap('MG_A1')


####regionSub####
ASD.brain.CTX <- subset(ASD.brain, Region == 'Isocortex')
cols <- pal_npg()(4)
names(cols) <- names(byscol)[1:4]
p<-SpatialDimPlot(ASD.brain.CTX, cols = cols, pt.size = 1.25, image.alpha = 0, crop = F, group.by = 'RegionAll',
               images = c('slice1.1', 'slice1.4', 'slice1.7'))
pdf('./figures/ST/SASP/SPD_CTX.pdf', width = 12, height = 4)
print(p)
dev.off()


ASD.brain.HIP <- subset(ASD.brain, Region == 'Hippocampal region')
cols <- pal_npg()(3)
names(cols) <- names(byscol)[10:12]
p<-SpatialDimPlot(ASD.brain.HIP, cols = cols, pt.size = 1.25, image.alpha = 0, crop = F, group.by = 'RegionAll',
               images = c('slice1.1', 'slice1.4', 'slice1.7'))
pdf('./figures/ST/SASP/SPD_HIP.pdf', width = 12, height = 4)
print(p)
dev.off()

length(unique(Markers_AN_merge[which(Markers_AN_merge$Region=='CTX'),]$gene))

length(unique(Markers_AC_merge[which(Markers_AN_merge$Region=='CTX'),]$gene))


table(ASD.brain$Sample)


meta1 <- ASD.brain@meta.data[,c(1,2,8)]

result <- meta1[, .(median_value1 = median(nCount_Spatial), median_value2 = median(nFeature_Spatial)), by = Sample]


result <- meta1 %>%
  group_by(Sample) %>%
  summarise(median_value1 = mean(nCount_Spatial), median_value2 = mean(nFeature_Spatial))


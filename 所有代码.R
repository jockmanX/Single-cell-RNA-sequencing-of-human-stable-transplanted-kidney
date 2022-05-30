library(Seurat)
library(SeuratObject)
library(cowplot)
library(patchwork)
library(dplyr)
library(harmony)
library(data.table)
library(ggplot2)
library(monocle)
library(tidyverse)
rm(list = ls())


assays <- dir("./matrix/")
dir <- paste0("./matrix/", assays)

# 按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
samples_name = c('KT1', 'KT2')
### 批量创建seurat对象
scRNAlist1 <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist1[[i]] <- CreateSeuratObject(counts, project = samples_name[i],
                                        min.cells=3, min.features = 200)
  #给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist1[[i]] <- RenameCells(scRNAlist1[[i]], add.cell.id = samples_name[i])
}


scRNA <- merge(scRNAlist1[[1]], scRNAlist1[[2]])



dir.create('./RData')


save(scRNA, file = './RData/1.scRNA.RData')
rm(scRNAlist,scRNAlist1,scRNAlist2,counts,assays,dir,i,samples_name)

#### 2.计算线粒体基因占比 ####
scRNA[['percent.mt']] <- PercentageFeatureSet(scRNA, pattern = 'MT-')
VlnPlot(scRNA, features = c('nFeature_RNA','nCount_RNA','percent.mt'), ncol = 3)

plot1 <- FeatureScatter(scRNA, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 <- FeatureScatter(scRNA, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
CombinePlots(plots = list(plot1, plot2))

#### 3.数据过滤+标准化 ####
#筛选条件：nFeature_RNA > 200 & nFeature_RNA < 2500 &percent.mt < 40
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                  percent.mt < 40)

scRNA <- NormalizeData(scRNA) #这里直接用默认的标准化方法

#筛选可变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = 'vst', nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA),10)
#查看前10个高可变基因
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))


#### 4.整理数据ScaleData ####
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA),   #选择高可变基因作为后续分析，也可以选择all.genes，但是运行会比较慢
                   vars.to.regress = 'percent.mt') #设置忽略线粒体基因

#### 5.PCA线性降维及结果展示 ####
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))

#1.查看PC1-PC5的前5个基因
print(scRNA[['pca']],dims = 1:5, nfeatures = 5)

#2.查看特定PC中的基因和基因对应的贡献度
VizDimLoadings(scRNA, dims = 1:2, reduction = 'pca')

#3.根据PCA结果绘制基因表达量热图
DimHeatmap(scRNA, dims = 1:15, cells = 500, balanced = T)

#4.PCA图
DimPlot(scRNA, reduction = 'pca')

#### 6.存在批次效应，使用harmony整合数据 ####
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony = 20)

DimPlot(scRNA, reduction = 'harmony')
save(scRNA, file = './RData/2.scRNA_after_harmony_mt40.RData')

ElbowPlot(scRNA, ndims = 50)
pc.num=1:20 #根据碎石图确定数字

#### 7.细胞聚类 ####
scRNA <- FindNeighbors(scRNA, dims = pc.num, reduction = 'harmony')
scRNA <- FindClusters(scRNA, resolution = 0.5)
scRNA <- RunUMAP(scRNA, dims = pc.num, resolution = 0.5, reduction = 'harmony')
# RunTSNE 得到tsne图, 个人喜欢UMAP，能把相似的细胞聚在一起

#按cluster，看umap图
DimPlot(scRNA, reduction = 'umap', label = T, label.size = 4)

#按样本，看umap图
DimPlot(scRNA, group.by = "orig.ident", reduction = 'umap')

save(scRNA, file = './RData/3.scRNA_after_umap_mt40.RData')

#### 8.Marker基因筛选，绘制热图 ####
cluster.marker <- FindAllMarkers(scRNA, only.pos = T, min.pct = 0.25,
                                 logfc.threshold = 0.25)
top10 <- cluster.marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(scRNA, features = top10$gene) + NoLegend()


#### 细胞注释 ####

scRNA@meta.data[["celltype"]]<- as.character(scRNA@meta.data$seurat_clusters)
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "0")] <- "PT"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "1")] <- "T cell"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "2")] <- "PT"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "3")] <- "T cell"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "4")] <- "PT"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "5")] <- "EC"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "6")] <- "myeloid cell"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "7")] <- "LOH"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "8")] <- "EC"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "9")] <- "fibrocyte"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "10")] <- "unknown"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "11")] <- "B cell"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters == "12")] <- "mast cell"

DimPlot(scRNA, group.by = "celltype", reduction = 'umap', label = T, label.size = 4)

save(scRNA, file = './RData/3.scRNA_after_celltype_mt40.RData')



#### 细胞数量及比例图 ####
celltype_colors <- c("#FF0000", "#FF44FF", "#4A4AFF", "#00EC00", "#F9F900",'#FF5809','#AD5A5A', '#5CADAD')
Idents(scRNA) <- 'celltype'
scRNA_known <- subset(scRNA, idents = c('PT', 'T cell', 'EC', 'myeloid cell', 'LOH', 'fibrocyte', 'B cell', 'mast cell'))

data_plotC <- table(scRNA_known@meta.data$celltype) %>% melt()
colnames(data_plotC) <- c("CellType","Number")

ggplot(data = data_plotC, aes(x = CellType, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="stack")+
  scale_fill_manual(values=celltype_colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell number")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) 

#出自https://zhuanlan.zhihu.com/p/25234546
#    https://www.cnblogs.com/ljhdo/p/4514106.html

cellnumber <- sum(data_plotC$Number)

ggplot(data = data_plotC, mapping = aes(x = 'content', y = Number, fill = CellType)) +
  geom_bar(stat = "identity",width = 0.5,position = 'stack',size = 5) +
  scale_fill_manual(values = celltype_colors) +  
  coord_polar(theta = "y") +
  labs(x = '', y = '', title = '') + 
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank())+
  geom_text(stat = 'identity', aes(y = Number, label = scales::percent(Number/cellnumber)), size = 4,
            position = position_stack(vjust = 0.5))


#### 细胞周期 ####
#计算周期得分
scRNA <- CellCycleScoring(scRNA, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

plot(scRNA$S.Score,scRNA$G2M.Score,
     col=factor(scRNA$Phase),
     main="CellCycleScoring")
legend("topleft",inset=.05,
       title = "cell cycle",  
       c("G1","S","G2M"), pch = c(1),col=c("black","green","red"))

DimPlot(scRNA, reduction = "umap", 
        group.by = "celltype",
        split.by = "Phase", 
        pt.size = 0.5)


######## T细胞细分 ########
Idents(scRNA) <- 'celltype'
T_cell <- subset(scRNA, idents='T cell')

ElbowPlot(T_cell, ndims = 50)
pc.num=1:15
# resolution = 0.2
T_cell <- FindNeighbors(T_cell, dims = pc.num, reduction = 'harmony')
T_cell <- FindClusters(T_cell, resolution = 0.2)
T_cell <- RunUMAP(T_cell, dims = pc.num, resolution = 0.2, reduction = 'harmony')
DimPlot(T_cell)

allmarkers <- FindAllMarkers(T_cell, only.pos = T, min.pct = 0.25,
                             logfc.threshold = 0.25)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(T_cell, features = top10$gene) + NoLegend()


VlnPlot(T_cell, features = 'CD8A')
VlnPlot(T_cell, features = 'FCGR3A')
VlnPlot(T_cell, features = 'IL7R')

T_cell@meta.data[["celltype"]]<- as.character(T_cell@meta.data$seurat_clusters)
T_cell@meta.data$celltype[which(T_cell@meta.data$seurat_clusters == "0")] <- "CD8 T"
T_cell@meta.data$celltype[which(T_cell@meta.data$seurat_clusters == "1")] <- "NKT"
T_cell@meta.data$celltype[which(T_cell@meta.data$seurat_clusters == "2")] <- "CD4 T"
DimPlot(T_cell, group.by = 'celltype')


######## 内皮细胞细分 ########
Idents(scRNA) <- 'celltype'
EC <- subset(scRNA, idents='EC')

ElbowPlot(EC, ndims = 50)
pc.num=1:15
# resolution = 0.2
EC <- FindNeighbors(EC, dims = pc.num, reduction = 'harmony')
EC <- FindClusters(EC, resolution = 0.2)
EC <- RunUMAP(EC, dims = pc.num, resolution = 0.2, reduction = 'harmony')
DimPlot(EC)

allmarkers <- FindAllMarkers(EC, only.pos = T, min.pct = 0.25,
                             logfc.threshold = 0.25)
top20 <- allmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(EC, features = top20$gene) + NoLegend()


VlnPlot(EC, features = 'EHD3')
VlnPlot(EC, features = 'KCNN3')

EC@meta.data[["celltype"]]<- as.character(EC@meta.data$seurat_clusters)
EC@meta.data$celltype[which(EC@meta.data$seurat_clusters == "0")] <- "EC"
EC@meta.data$celltype[which(EC@meta.data$seurat_clusters == "1")] <- "glomerulus-EC"
EC@meta.data$celltype[which(EC@meta.data$seurat_clusters == "2")] <- "tubule-EC"
DimPlot(EC, group.by = 'celltype')


######## 髓系细胞细分 ########
Idents(scRNA) <- 'celltype'
myeloid <- subset(scRNA, idents='myeloid cell')

ElbowPlot(myeloid, ndims = 50)
pc.num=1:15
# resolution = 0.1
myeloid <- FindNeighbors(myeloid, dims = pc.num, reduction = 'harmony')
myeloid <- FindClusters(myeloid, resolution = 0.1)
myeloid <- RunUMAP(myeloid, dims = pc.num, resolution = 0.1, reduction = 'harmony')
DimPlot(myeloid)

allmarkers <- FindAllMarkers(myeloid, only.pos = T, min.pct = 0.25,
                             logfc.threshold = 0.25)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(myeloid, features = top10$gene) + NoLegend()


VlnPlot(myeloid, features = 'MRC1')
VlnPlot(myeloid, features = 'EREG')
VlnPlot(myeloid, features = 'CD300E')

myeloid@meta.data[["celltype"]]<- as.character(myeloid@meta.data$seurat_clusters)
myeloid@meta.data$celltype[which(myeloid@meta.data$seurat_clusters == "0")] <- "M2"
myeloid@meta.data$celltype[which(myeloid@meta.data$seurat_clusters == "1")] <- "monocyte"
DimPlot(myeloid, group.by = 'celltype')

#### 拟时序分析 ####
data <- as(as.matrix(myeloid@assays$RNA@counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = myeloid@meta.data)
fData <- data.frame(gene_short_name = row.names(data),
                    row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
                              lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

#查看phenodata、featuredata
head(pData(monocle_cds))
head(fData(monocle_cds))

#2.数据过滤+归一化 
HSMM <- monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 3)
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

# 3.降维、聚类 
#选择细胞间异常变异的基因用于聚类
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)

#用DDRTree算法进行降维分析
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')

HSMM <- orderCells(HSMM)

# 6.绘制轨迹 
colnames(pData(HSMM)) #可以查看用哪些参数替换color_by

plot_cell_trajectory(HSMM, color_by = 'State')  

plot_cell_trajectory(HSMM, color_by = 'celltype') 

plot_cell_trajectory(HSMM, color_by = 'Pseudotime') 








#让Pseudotime反方向
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1) {
    T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,'0']
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HSMM_1 <- orderCells(HSMM, root_state = GM_state(HSMM))
plot_cell_trajectory(HSMM_1, color_by = 'Pseudotime')



BEAM_res <- BEAM(HSMM[expressed_genes,], branch_point = 1, cores = 3)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c('gene_short_name', 'pval', 'qval')]

plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,qval<1e-19)),],
                            branch_point = 1, num_clusters = 3, cores = 1,
                            use_gene_short_name = T, show_rownames = T)






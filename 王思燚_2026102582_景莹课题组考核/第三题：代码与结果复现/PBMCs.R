
##  第 0 段：准备工作
# Install package
library(dplyr)
library(Seurat)
library(patchwork)
# presto 为可选项：仅加速 FindAllMarkers 的 Wilcoxon 检验
# Seurat 会自行探测，未安装时回退到标准 wilcox.test
if (requireNamespace("presto", quietly = TRUE)) library(presto)
library(ggplot2)

# 10x CellRanger 输出目录（barcodes.tsv / genes.tsv / matrix.mtx）
# 相对路径，需与本脚本同级；放在别处时改成对应绝对路径
DATA_DIR      <- "filtered_gene_bc_matrices/hg19"
OUT_DIR       <- "output"     

#QC and selecting cells for further analysis
GENE_MIN_CELL <- 3            # 仅保留在 ≥3 个细胞中检出的基因
CELL_MIN_GENE <- 200          # 仅保留检出基因数 ≥200 的细胞
MT_MAX        <- 5            # 线粒体 UMI 占比上限：高占比——膜破损——Low-quality/dying cells
CELL_MAX_GENE <- 2500         # 细胞检出上限：检出基因数 ≥2500 疑为双细胞，即两细胞共用一个液滴

#降维
N_HVG         <- 2000 
N_PC          <- 10           
RESOLUTION    <- 0.5

# 建输出目录
dir.create(OUT_DIR, showWarnings = FALSE)

##  第 1 段：数据读入与 Seurat 对象构建

# 读取barcodes.tsv / genes.tsv / matrix.mtx
# 组装为稀疏表达矩阵：行为基因，列为细胞，元素值为该细胞该基因的 UMI 计数
pbmc.data <- Read10X(data.dir = DATA_DIR)

# 构建 Seurat 对象，同时施加两道初步过滤：
#   min.cells    —— 剔除只在极少数细胞中检出的基因
#   min.features —— 剔除检出基因数过少的细胞
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells    = GENE_MIN_CELL,
                           min.features = CELL_MIN_GENE)

##  第 2 段：Quality Control

# 计算线粒体基因 UMI 占各细胞总 UMI 的比例（percent.mt）
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 绘制三个 QC 指标的分布
qc_plot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(OUT_DIR, "01_质控分布.png"), qc_plot, width = 12, height = 5)

# 绘制指标间相关性散点图：nCount_RNA 与 nFeature_RNA 通常高度正相关，
# 二者共同受测序深度影响
scatter_plot <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(file.path(OUT_DIR, "01b_指标关系.png"), scatter_plot, width = 10, height = 5)

# 施加三重过滤条件：
#   nFeature_RNA 下限 —— 剔除空液滴与低质量细胞碎片
#   nFeature_RNA 上限 —— 剔除疑似双细胞（doublet）
#   percent.mt   上限 —— 剔除濒死或膜破损细胞
pbmc <- subset(pbmc, subset = nFeature_RNA >  CELL_MIN_GENE &
                 nFeature_RNA <  CELL_MAX_GENE &
                 percent.mt   <  MT_MAX)

##  第 3 段：归一化、特征选择与线性降维

# LogNormalize
# 消除细胞间测序深度差异带来的偏差
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 基于VST筛选高变基因
# 即细胞间表达变异最显著的前 N_HVG 个基因，用于后续降维以提升信噪比
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = N_HVG)

cat("第 3 段：挑出", length(VariableFeatures(pbmc)), "个关键基因，前 10 个是\n")
print(head(VariableFeatures(pbmc), 10))

# 可视化高变基因：横轴为平均表达量，纵轴为标准化方差
hvg_top <- head(VariableFeatures(pbmc), 10)
hvg_plot <- LabelPoints(plot = VariableFeaturePlot(pbmc), points = hvg_top, repel = TRUE)
ggsave(file.path(OUT_DIR, "02_高变基因.png"), hvg_plot, width = 10, height = 5)

# 对全部基因做 z-score 标准化
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# 碎石图（Elbow plot）
elbow_plot <- ElbowPlot(pbmc)
ggsave(file.path(OUT_DIR, "03_碎石图.png"), elbow_plot, width = 6, height = 5)

##  第 4 段：聚类、标志基因鉴定与细胞类型注释

# 基于前 N_PC 个主成分构建 KNN 图，并转换为共享最近邻图（SNN graph）
pbmc <- FindNeighbors(pbmc, dims = 1:N_PC)

# 以 Louvain 算法执行社区发现；resolution 控制划分粒度，值越大簇数越多
pbmc <- FindClusters(pbmc, resolution = RESOLUTION)

cat("第 4 段：分成", nlevels(pbmc), "个群，各群细胞数：\n")
print(table(Idents(pbmc)))

# UMAP 非线性降维至二维，用于可视化
pbmc <- RunUMAP(pbmc, dims = 1:N_PC)

# 差异表达分析：以 Wilcoxon 秩和检验识别各簇相对其余细胞显著上调的标志基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# 以已知谱系标志基因（canonical markers）反向检索其在各簇中的富集程度，
# 作为细胞类型注释的直接依据。不宜仅取差异倍数最大的基因——
# 此类基因往往只在极少数细胞中表达，缺乏谱系判别力
known_markers <- c(
  "IL7R", "CCR7",      # Naive CD4 T cell
  "CD14", "LYZ",       # CD14+ monocyte
  "S100A4",            # Memory CD4 T cell
  "MS4A1",             # B cell
  "CD8A",              # CD8 T cell
  "FCGR3A", "MS4A7",   # FCGR3A+ monocyte
  "GNLY", "NKG7",      # NK cell
  "FCER1A", "CST3",    # dendritic cell
  "PPBP"               # platelet
)
check <- pbmc.markers[pbmc.markers$gene %in% known_markers,
                      c("gene", "cluster", "avg_log2FC")]
cat("\n已知标记基因各自在哪一群最高（用来核对命名）：\n")
print(check[order(check$gene, -check$avg_log2FC), ], row.names = FALSE)

# 在 UMAP 投影上可视化标志基因的表达分布，用于交叉验证注释结果
marker_plot <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14",
                                              "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
ggsave(file.path(OUT_DIR, "04_标记基因.png"), marker_plot, width = 14, height = 10)

# 绘制各簇表达量最高的 10 个标志基因热图
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

heatmap_plot <- DoHeatmap(pbmc, features = top_markers$gene) + NoLegend()
ggsave(file.path(OUT_DIR, "05_标志基因热图.png"), heatmap_plot, width = 10, height = 12)

# 依据上述标志基因证据对 9 个簇赋予细胞类型标签
cell_type_names <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(cell_type_names) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, cell_type_names)

# 打印细胞名单
print(table(Idents(pbmc)))


##  第 5 段：结果可视化

# 生成最终 UMAP 图：标注细胞类型、保留图例、调整坐标轴标题与字号
final_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title   = element_text(size = 18),
        legend.text  = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

# 导出图形
ggsave(file.path(OUT_DIR, "pbmc3k_umap.pdf"), final_umap,
       width = 12, height = 7)

# 环境配置
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
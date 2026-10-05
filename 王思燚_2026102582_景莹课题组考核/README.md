# 生物信息学课题组入组考核 · 三题作业

> 考核人：王思燚（学号 2026102582）　申请课题组：景莹课题组

本仓库收录入组考核三道题的完整材料：数据可视化与数据清洗、英文文献批判性阅读、以及 Seurat 单细胞标准流程的代码复现与差异归因。三道题分别考察**数据处理与图形表达**、**文献理解与事实核查**、**可复现的计算流程与调试能力**。

---

## Abstract (English)

This repository contains three assignments for the lab admission assessment: (1) **data visualization** — cleaning a deliberately corrupted bulk expression matrix and answering analytical questions on boxplot / violin plot / log-scale interpretation; (2) **critical reading** — auditing an AI-generated summary of Hanahan (2022) *Hallmarks of Cancer: New Dimensions* for inaccuracies, missing information and overinterpretation; (3) **reproduction** — re-running the Satija Lab PBMC3k guided clustering tutorial under Seurat 5.x, documenting every error/warning encountered and attributing the visual differences from the original tutorial. All analyses are in R; raw data and outputs are included.

---

## 仓库结构

```
.
├── README.md
│
├── 第一题：数据可视化/
│   ├── task1_expression_data.csv    # 原始表达量矩阵（含人为植入的脏数据）
│   ├── visualization.R              # 清洗 + 绘图脚本
│   ├── 问题回答.md                   # 三道分析题的作答
│   ├── boxplot.png                  # 图 1：GENE01 分组箱线图（log10 轴）
│   ├── violin.png                   # 图 2：GENE05 小提琴图
│   └── heatmap.png                  # 图 3：12 基因 × 59 样本 z-score 热图
│
├── 第二题：英文文献阅读.md            # Hanahan (2022) 摘要的 AI 总结审校
│
└── 第三题：代码与结果复现/
    ├── PBMCs.R                      # Seurat 全流程脚本（含中文注释）
    ├── UMAP .pdf                    # 最终 UMAP 细胞类型注释图
    └── 复现记录.md                   # 环境配置、报错排查、与教程图差异归因
```

---

## 题目总览

| 题号 | 主题 | 核心工作 | 主要产出 |
|:--|:--|:--|:--|
| 第一题 | 数据可视化 | 识别并修复 6 类数据质量问题，用三种图形回答分布与尺度问题 | `visualization.R` + 3 张图 + `问题回答.md` |
| 第二题 | 英文文献阅读 | 逐条核查 AI 摘要与原文的事实偏差，辨析 hallmark capability 与 enabling characteristic | `第二题：英文文献阅读.md` |
| 第三题 | 代码与结果复现 | 在 Seurat 5.x 上复现 PBMC3k 教程，排查 5 处报错/警告并解释图形差异 | `PBMCs.R` + UMAP 图 + `复现记录.md` |

---

## 第一题：数据可视化

### 数据与问题

`task1_expression_data.csv` 为 12 个基因 × 61 行的表达量矩阵（字段：`sample_id`, `group`, `batch`, `GENE01`–`GENE12`）。原始文件被植入了 6 类质量问题：

| # | 问题类型 | 具体证据 | 处理方式 |
|:--|:--|:--|:--|
| 1 | 分组名写法混乱 | `control / Control / CTRL`、`Treatment / treatment / TREAT` 共 6 种写法 | `tolower()` + `trimws()` 统一，`ctrl→control`、`treat→treatment` |
| 2 | 重复样本 | `S012` 在文件末尾重复出现 | `!duplicated()` 保留首次 |
| 3 | 分组缺失 | `S005` 的 `group` 为空 | 无法归类，删除该样本 |
| 4 | 字符串混入数值列 | `GENE03` 含 `"n/a"`；`GENE07` 含 `"1,240.5"` | `na.strings` 转 `NA`；`gsub(",", "", x)` 后转 `numeric` |
| 5 | 过高值 | `S044` 的 `GENE07 = 1240.5` 显著离群 | 判为录入错误，置 `NA` |
| 6 | 整批单位错误 | `GENE10` 的 batch B 整体比 batch A 大两个数量级 | batch B 整体除以 100 校正 |

清洗后得到 **59 个样本（control 29 / treatment 30）**。

### 三张图与结论

**图 1 · 箱线图（`boxplot.png`）**：GENE01 按分组中位数从高到低排序，y 轴取 log10。treatment 组中位数 11.93，control 组 8.23，差异清晰可辨。

**图 2 · 小提琴图（`violin.png`）**：GENE05 两组的**中位数几乎相同**（control 7.58，treatment 7.35），仅看箱线图会得出"无差异"的结论。小提琴图暴露出 treatment 组呈**明显双峰分布**——16 个样本落在 5.0–7.6，14 个样本落在 13.0–14.7，中间 7.6–13.0 区间**无任何样本**。这说明 treatment 组内部存在"应答 / 不应答"两个亚群，均值与中位数把这一结构完全掩盖了。

**图 3 · 热图（`heatmap.png`）**：12 基因 × 59 样本，表达值先 `log2` 转换，再**按基因逐行做 z-score 标准化**（用 `apply()` 手写而非 `scale()`，因为 `scale()` 遇到 `NA` 会把整行变成 `NA`）。行按基因、列按 `group` 再按 `batch` 排序，`cluster_rows/cols = FALSE`，避免聚类打乱分组顺序。注释条显示 Group 与 Batch 两个变量。

### 关于对数坐标

- **为什么必须用**：对数轴上，**相同的倍数变化对应相同长度**；线性轴上则是相同差值等距。表达量差异的生物学意义在于倍数关系，故用 log 轴更贴合。
- **不用的后果**：未处理的离群值（如 1240.5）会把线性轴范围拉到极大，其余大部分样本被挤压在图底部，无法分辨组间差异。
- **何时反而误导**：当数据本身都在同一数量级、差异本来就小时，对数轴会把真实差距**视觉上抹平**，让读者低估差异。

---

## 第二题：英文文献阅读

阅读对象：**Hanahan D. Hallmarks of Cancer: New Dimensions. *Cancer Discovery*, 2022.**

任务是审校一段 AI 生成的摘要，找出其中的不准确、信息缺失与过度解读，共 4 条：

1. **不准确——"two new hallmark capabilities" 的表述有误。** 原文将 phenotypic plasticity and disrupted differentiation 表述为一个**离散的（discrete）单一 hallmark capability**，而非两项独立能力。
2. **不准确——混淆了 enabling characteristic 与 hallmark capability。** 原文明确 nonmutational epigenetic reprogramming 与 polymorphic microbiomes **"both constitute distinctive enabling characteristics"**，AI 摘要却称之为"新的 hallmark capabilities"，属于框架中两个层级的混淆。
3. **过度解读——senescent cells 并非"newly confirmed"。** 原文用词是 "may be added to the roster"，属于**提议性**表述，而非已确认的结论。
4. **信息缺失——忽略了框架的启发式与暂定性。** 原文将 hallmarks 概念定位为 "heuristic tool"、产出 "provisional set of underlying principles"，而 AI 摘要把框架描述得过于确定。

第三部分进一步论证了区分 **hallmark capability** 与 **enabling characteristic** 的必要性：前者是肿瘤细胞实际获得的功能能力（如持续增殖、免疫逃逸、表型可塑性），后者是促成这些能力的上游条件或机制（如非突变性表观遗传重编程、多态性微生物组）。这一区分对**因果推理、实验设计与治疗靶点选择**有实际意义；同时该边界并非绝对——因为框架本身是启发式、暂定的。

---

## 第三题：代码与结果复现

复现对象：**Satija Lab · Seurat PBMC3k Guided Clustering Tutorial**（教程基于 Seurat 4.x，本次在 **Seurat 5.5.1** 上运行）。

### 流程概览

`PBMCs.R` 按标准单细胞流程组织，全部关键参数提取为脚本开头的常量，便于调整：

| 阶段 | 操作 | 关键参数 |
|:--|:--|:--|
| 读入 | `Read10X()` + `CreateSeuratObject()` | `min.cells = 3`，`min.features = 200` |
| 质控 | `PercentageFeatureSet("^MT-")` + `subset()` | `200 < nFeature_RNA < 2500`，`percent.mt < 5` |
| 归一化 | `NormalizeData()` LogNormalize | `scale.factor = 10000` |
| 特征选择 | `FindVariableFeatures()` | VST，`nfeatures = 2000` |
| 标准化 / 降维 | `ScaleData()` → `RunPCA()` | 碎石图定主成分数 |
| 聚类 | `FindNeighbors()` → `FindClusters()` | `dims = 1:10`，Louvain，`resolution = 0.5` |
| 可视化 | `RunUMAP()` → `DimPlot()` | 前 10 个主成分 |
| 注释 | `FindAllMarkers()` + 已知标志基因交叉验证 | Wilcoxon 秩和检验，`only.pos = TRUE` |

最终得到 **9 个细胞群**，注释为：Naive CD4 T、CD14+ Mono、Memory CD4 T、B、CD8 T、FCGR3A+ Mono、NK、DC、Platelet。

注释策略上，脚本特意**不直接采用 `avg_log2FC` 最高的基因**——这类基因往往只在极少数细胞中表达，缺乏谱系判别力；而是用一组已知谱系标志基因（`IL7R`/`CCR7`、`CD14`/`LYZ`、`MS4A1`、`CD8A`、`FCGR3A`/`MS4A7`、`GNLY`/`NKG7`、`FCER1A`/`CST3`、`PPBP` 等）反向检索其在各簇中的富集情况，并用 `FeaturePlot()` 做交叉验证。

### 遇到的报错与解决（详见 `复现记录.md`）

| # | 类型 | 现象 | 根因 | 解决 |
|:--|:--|:--|:--|:--|
| 1 | 报错 | `Read10X: Directory provided does not exist` | 教程里的 `/brahms/mollag/practice/...` 是**作者本机路径** | 改成本机数据目录 |
| 2 | 警告 | `Feature names cannot have underscores` | Seurat 5.x 规定基因名不得含 `_` | 正常行为，忽略 |
| 3 | 报错 | `The RStudio Plots window may be too small to show this patchwork` | 提示语**具有误导性**；`print.patchwork` 用 `tryCatch` 把真实错误 `Viewport has zero dimension(s)` 替换掉了 | 绕过面板，改用 `ggsave()` 直接出图 |
| 4 | 报错 | `saveRDS: cannot open the connection` | `saveRDS()` **不会自动创建目录**，`../output/` 不存在 | 改用项目内路径 + `dir.create()` |
| 5 | 疏漏 | 跑完没有任何图片产出 | 教程中 `ggsave()` 处于注释状态 | 取消注释并修正输出路径 |

另外，第一次运行 `FindAllMarkers()` 极慢，Seurat 自身提示可安装 `presto` 加速 Wilcoxon 检验；安装后全流程耗时明显下降（脚本中已写成**可选依赖**：检测到则加载，未安装则自动回退到标准 `wilcox.test`）。

### 与教程原图的差异及归因

差异有两处：**UMAP 整体形状与旋转角度不同**、**本图右侧多出一列细胞类型图例**。

- **右侧图例**：教程中存在两版绘图代码。正文展示的那版末尾带 `+ NoLegend()`，而末尾附录给出的那版（赋值给 `plot`）**保留了图例并调整了字号**，只是没有在网页中展示。本仓库采用后者，因此有图例——属于代码选择差异，而非结果差异。
- **形状差异**：归因于两点。(a) **UMAP 底层实现更换**——Seurat 4 默认调用 Python 的 `umap-learn`，Seurat 5 改为 R 实现的 `uwot`，两者算出的二维坐标本就不同，但**邻接关系保持**，故分群与注释结果不受影响；(b) **UMAP 自带随机性**——`uwot` 使用自身的随机数源，不受 R 的 `set.seed()` 控制，即使参数完全一致也难以复现出逐点相同的图。

---

## 运行环境

```
R 4.6.1 (2026-06-24)
```

| 题号 | 依赖包 |
|:--|:--|
| 第一题 | `readr` 2.2.0 · `dplyr` 1.2.1 · `ggplot2` 4.0.3 · `pheatmap` 1.0.13 |
| 第三题 | `Seurat` 5.5.1 · `dplyr` 1.2.1 · `patchwork` 1.3.2 · `ggplot2` 4.0.3 · `presto`（可选，加速用） |

一次性安装：

```r
install.packages(c("readr", "dplyr", "ggplot2", "pheatmap", "patchwork"))
install.packages("Seurat")
# 可选：显著加速 FindAllMarkers 的 Wilcoxon 检验
install.packages("devtools")
devtools::install_github("immunogenomics/presto")
```

---

## 复现步骤

**第一题**

```bash
cd "第一题：数据可视化"
Rscript visualization.R        # 或 source("visualization.R")
```

数据文件与脚本同目录，脚本内使用相对路径，工作目录设为该文件夹即可。输出为 `boxplot.png`、`violin.png`、`heatmap.png`。

**第三题**

1. 从 [10x Genomics 官网](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) 下载 PBMC3k 数据并解压，得到 `filtered_gene_bc_matrices/hg19/` 目录（含 `barcodes.tsv`、`genes.tsv`、`matrix.mtx`）。
2. 将该目录放在 `PBMCs.R` **同级目录**下（脚本使用相对路径 `DATA_DIR <- "filtered_gene_bc_matrices/hg19"`，放在别处请改成绝对路径）。
3. 运行：

```bash
cd "第三题：代码与结果复现"
Rscript PBMCs.R
```

脚本会自动创建 `output/` 目录，输出质控图、高变基因图、碎石图、标志基因 FeaturePlot、标志基因热图、最终 UMAP（`pbmc3k_umap.pdf`）以及 `sessionInfo.txt`。

---

## 参考

1. Hanahan D. Hallmarks of Cancer: New Dimensions. *Cancer Discovery*, 2022, 12(1): 31–46.
2. Satija Lab. *Seurat Guided Clustering Tutorial (PBMC3k)*. https://satijalab.org/seurat/articles/pbmc3k_tutorial
3. Hao Y, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, 2024.（Seurat 5）
4. Korsunsky I, et al. Presto scales Wilcoxon and auROC analyses to millions of observations. *bioRxiv*, 2019.

---

*本仓库为课题组入组考核提交材料，欢迎批评指正。*

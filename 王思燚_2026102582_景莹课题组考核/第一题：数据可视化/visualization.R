# 装包
library(readr)     # 读 csv
library(dplyr)     # 数据整理
library(ggplot2)   # 画 boxplot / violin
library(pheatmap)  # 画热图

# 数据处理
# 1. 读入数据
dat <- read.csv("task1_expression_data.csv",
                stringsAsFactors = FALSE,
                na.strings = c("n/a", "NA", "", "N/A"))

# 去千分位逗号，转数值（不处理的话 GENE07 是字符型）
dat$GENE07 <- as.numeric(gsub(",", "", dat$GENE07))

# 2. 统一分组名
dat$group <- tolower(trimws(dat$group))          # 全小写、去空格
dat$group[dat$group == "ctrl"]  <- "control"
dat$group[dat$group == "treat"] <- "treatment"

# 3. 删掉没有分组的样本
dat <- dat[!duplicated(dat$sample_id), ]
dat <- dat[!is.na(dat$group) & dat$group != "", ]
# 现在：59 个样本（control 29，treatment 30）

# 4. 修异常值 
dat$GENE07[dat$GENE07 > 100] <- NA                       # 1240.5 -> NA
dat$GENE10[dat$batch == "B"] <- dat$GENE10[dat$batch == "B"] / 100   # 批次单位错误

genes <- paste0("GENE", sprintf("%02d", 1:12))

#图 1-Heatmap

m <- t(as.matrix(dat[, genes]))      # 行=基因，列=样本
colnames(m) <- dat$sample_id
m <- log2(m)
# 每个基因自己减均值除标准差（用 apply 而不是 scale()，因为 scale() 遇到 NA 会让整行变 NA）
m <- (m - apply(m, 1, mean, na.rm = TRUE)) / apply(m, 1, sd, na.rm = TRUE)

o   <- order(dat$group, dat$batch)                 # 按分组（再按批次）排样本
ann <- data.frame(Group = dat$group, Batch = dat$batch)
rownames(ann) <- dat$sample_id

pheatmap(m[, o], annotation_col = ann[o, , drop = FALSE],
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
         main = "12 genes x 59 samples (z-score, batch-corrected)")

#图 2-Boxplot

lev <- names(sort(tapply(dat$GENE01, dat$group, median), decreasing = TRUE))
dat$group <- factor(dat$group, levels = lev)   # treatment(11.93) -> control(8.23)

ggplot(dat, aes(x = group, y = GENE01, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  scale_y_log10() +
  labs(title = "GENE01 (log10 scale)", x = "group (median, high to low)", y = "expression") +
  theme_classic() + theme(legend.position = "none")

#图 3-Violin plot

ggplot(dat, aes(x = group, y = GENE05, fill = group)) +
  geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.7) +
  labs(title = "GENE05", x = "group", y = "expression") +
  theme_classic() + theme(legend.position = "none")








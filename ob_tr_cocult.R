setwd("/Users/akhvorov/lab/cocult/cocult_transcr")

# читаем данные и легенду
library(readxl) 
coldata <- data.frame(read_excel("Cocult_trans_legend.xlsx")) 
counts <- read.delim("featureCounts_NM+PH.txt")

# делаем ID названием строк
rownames(counts) <- counts[,1]
counts <- counts[,-1]
# меняем названия столбиков, чтобы названия образцов совпадали с легендой
colnames(counts) <- coldata[,1]

# делаем легенду и даасет каунтов для остеобластов
ost_coldata <- coldata[coldata$Cell.type == 'Osteo', ]
ost_cts <- counts[, rownames(ost_coldata)]

# убираем из coldata тип клеток и лишний столбец sample
ost_coldata <- ost_coldata[, c(-1, -2)]
# делаем все переменные coldata факторами
ost_coldata$Condition <- as.factor(ost_coldata$Condition)
ost_coldata$Donor <- as.factor(ost_coldata$Donor)
str(ost_coldata)

#keep <- rowSums(ost_cts) >= 10
#cts1 <- ost_cts[keep,]

# DESeq2
library("DESeq2")

dds_ost <- DESeqDataSetFromMatrix(countData = ost_cts,
                              colData = ost_coldata,
                              design= ~  Donor + Condition)
dds_ost

# изменение уровней факторов
dds_ost$Сondition <- factor(dds_ost$Condition, levels = c("Control", "sorting", "separated"))
dds_ost$Donor <- factor(dds_ost$Donor, levels = c("9", "10", "15"))

# минимальная фильтрация ()
keep <- rowSums(counts(dds_ost)) >= 10
dds_ost <- dds_ost[keep,]

# Дифференциальный экспрессионный анализ
dds_ost <- DESeq(dds_ost)
res <- results(dds_ost, alpha=0.05)
resultsNames(dds_ost)

                   
# Независимое взвешивание гипотез    ????????
library("IHW")
resIHW <- results(dds_ost, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.05, na.rm=TRUE)
metadata(resIHW)$ihwResult


#### Condition_sorting_vs_Control ####
res_sort_contr <- results(dds_ost, name="Condition_sorting_vs_Control", alpha=0.05)

# shrink log fold changes association with condition:
library(apeglm)
resLFC_sort_contr <- lfcShrink(dds_ost, coef="Condition_sorting_vs_Control", type="apeglm")
resLFC_sort_contr
resNorm_sort_contr <- lfcShrink(dds_ost, coef="Condition_sorting_vs_Control", type="normal")
resAsh_sort_contr <- lfcShrink(dds_ost, coef="Condition_sorting_vs_Control", type="ashr")

# MA-plot
plotMA(res_sort_contr, ylim=c(-2,2))
plotMA(resLFC_sort_contr, ylim=c(-2,2))

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC_sort_contr, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm_sort_contr, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh_sort_contr, xlim=xlim, ylim=ylim, main="ashr")
dev.off()

# p-значения и скорректированные p-значения
# упорядочить нашу таблицу результатов по наименьшему значению p :
resLFC_ordered_sort_contr <- resLFC_sort_contr[order(resLFC_sort_contr$pvalue),]

### записываем в файл

# полная упорядоченная таблица
#write.csv(as.data.frame(resLFC_ordered_sort_contr), file="Ost_Condition_sorting_vs_Control_trans.csv")
library(xlsx)
#write.xlsx(as.data.frame(resLFC_ordered_sort_contr), file="Ost_Condition_sorting_vs_Control_trans.xlsx")

# упорядоченная отфильтрованная таблица по logFC & padj
resLFC_ordered_subset_sort_contr <- subset(resLFC_ordered_sort_contr, padj < 0.05 & abs(log2FoldChange) >= 1)
# выделяем up & down транскрипты
sort_contr_up <- resLFC_ordered_subset_sort_contr[resLFC_ordered_subset_sort_contr$log2FoldChange > 0,]
sort_contr_down <- resLFC_ordered_subset_sort_contr[resLFC_ordered_subset_sort_contr$log2FoldChange < 0,]

# записываем в файл
#write.xlsx(as.data.frame(sort_contr_up), file="Ost_Condition_sorting_vs_Control_trans_up.xlsx")
#write.xlsx(as.data.frame(sort_contr_down), file="Ost_Condition_sorting_vs_Control_trans_down.xlsx")

# некоторые основные подсчеты
summary(res_sort_contr)
summary(resLFC_sort_contr)

# Сколько скорректированных значений p были меньше 0,05
sum(res_sort_contr$padj < 0.05, na.rm=TRUE)

head(resLFC_ordered_subset_sort_contr)
# counts of reads for a single gene across the groups
plotCounts(dds_ost, gene=which.min(res_sort_contr$padj), intgroup="Condition")
# или любое название гена можно вставить
plotCounts(dds_ost, gene='ENSG00000130702.15', intgroup="Condition")
# построения графика с помощью ggplot
d <- plotCounts(dds_ost, gene=which.min(res_sort_contr$padj), intgroup="Condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(100,200,4000))

mcols(res_sort_contr)$description

#### Condition_separated_vs_Control ####
res_sep_contr <- results(dds_ost, name="Condition_separated_vs_Control", alpha=0.05)

# shrink log fold changes association with condition:
library(apeglm)
resLFC_sep_contr <- lfcShrink(dds_ost, coef="Condition_separated_vs_Control", type="apeglm")
resLFC_sep_contr
resNorm_sep_contr <- lfcShrink(dds_ost, coef="Condition_separated_vs_Control", type="normal")
resAsh_sep_contr <- lfcShrink(dds_ost, coef="Condition_separated_vs_Control", type="ashr")

# MA-plot
plotMA(res_sep_contr, ylim=c(-2,2))
plotMA(resLFC_sep_contr, ylim=c(-2,2))

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC_sep_contr, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm_sep_contr, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh_sep_contr, xlim=xlim, ylim=ylim, main="ashr")
dev.off()

# p-значения и скорректированные p-значения
# упорядочить нашу таблицу результатов по наименьшему значению p :
resLFC_ordered_sep_contr <- resLFC_sep_contr[order(resLFC_sep_contr$pvalue),]

### записываем в файл

# полная упорядоченная таблица
#write.csv(as.data.frame(resLFC_ordered_sep_contr), file="HUV_Condition_seporated_vs_Control_trans.csv")
library(xlsx)
#write.xlsx(as.data.frame(resLFC_ordered_sep_contr), file="Ost_Condition_separated_vs_Control_trans.xlsx")

# упорядоченная отфильтрованная таблица по logFC & padj
resLFC_ordered_subset_sep_contr <- subset(resLFC_ordered_sep_contr, padj < 0.05 & abs(log2FoldChange) >= 1)
# выделяем up & down транскрипты
sep_contr_up <- resLFC_ordered_subset_sep_contr[resLFC_ordered_subset_sep_contr$log2FoldChange > 0,]
sep_contr_down <- resLFC_ordered_subset_sep_contr[resLFC_ordered_subset_sep_contr$log2FoldChange < 0,]

# записываем в файл
#write.xlsx(as.data.frame(sep_contr_up), file="Ost_Condition_separated_vs_Control_trans_up.xlsx")
#write.xlsx(as.data.frame(sep_contr_down), file="Ost_Condition_separated_vs_Control_trans_down.xlsx")

# некоторые основные подсчеты
summary(res_sep_contr)
summary(resLFC_ordered_sep_contr)

# Сколько скорректированных значений p были меньше 0,05
sum(resLFC_ordered_subset_sep_contr$padj < 0.05, na.rm=TRUE)

# counts of reads for a single gene across the groups
head(resLFC_ordered_subset_sep_contr)
# counts of reads for a single gene across the groups
plotCounts(dds_ost, gene=which.min(res_sep_contr$padj), intgroup="Condition")
# или любое название гена можно вставить
plotCounts(dds_ost, gene='ENSG00000105825.14', intgroup="Condition")

# построения графика с помощью ggplot
d <- plotCounts(dds_ost, gene=which.min(res_sep_contr$padj), intgroup="Condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(100,200,4000))

mcols(res_sep_contr)$description



# Извлечение преобразованных значений
vsd_ost <- vst(dds_ost, blind=FALSE)
rld_ost <- rlog(dds_ost, blind=FALSE)
assay(vsd_ost)

# Влияние преобразований на дисперсию
# this gives log2(n + 1)
ntd_ost <- normTransform(dds_ost)
library("vsn")
meanSdPlot(assay(ntd_ost))
meanSdPlot(assay(vsd_ost))
meanSdPlot(assay(rld_ost))

# Тепловая карта матрицы счета
library("pheatmap")
select <- order(rowMeans(counts(dds_ost,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_ost)[,c("Condition","Donor")])

pheatmap(assay(ntd_ost)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd_ost)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld_ost)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Тепловая карта расстояний между образцами
sampleDists <- dist(t(assay(vsd_ost)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_ost$Condition, vsd_ost$Donor, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# График главных компонент образцов
plotPCA(vsd_ost, intgroup=c("Condition", "Donor"))

pcaData <- plotPCA(vsd_ost, intgroup=c("Condition", "Donor"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Donor)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

library(mixOmics)
ost_pca <- pca(t(assay(ntd_ost)), ncomp = 5, center = TRUE)
#tiff('Ost_PCA_group.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')
plotIndiv(ost_pca, comp = c(1, 2), ind.names = F, group = ntd_ost$Condition, 
          legend = TRUE, ellipse = T, title = 'PCA')
dev.off()

#PLS-DA 
ordination.optimum.splsda <- splsda(t(assay(ntd_ost)), ost_coldata$Condition, ncomp = 3, keepX = c(45,45,45))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)


#layout(matrix(c(1, 2), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)

plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")

#tiff('PLSDA_ost.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 2, size.axis = 1.0, size.xlabel = 1.2, size.ylabel = 1.2, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
#dev.off()

#tiff('PLSDA_all_ost_transcr.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()

#tiff('PCA_PLS_ost_transcr.tiff', units="in", width=18, height=8, res=600, compression = 'lzw')
layout(matrix(1:2, ncol = 2))
plotIndiv(ost_pca, comp = c(1, 2), ind.names = F, group = ntd_ost$Condition, 
          style = "graphics", size.title = 1.4, ellipse = T, title = 'PCA', legend=TRUE)
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, abline = TRUE,
          style = "graphics", size.title = 1.4, title = "PLS-DA ordination", legend=TRUE)
#dev.off()


####################

dds_ost <- DESeqDataSetFromMatrix(countData = ost_cts,
                                  colData = ost_coldata,
                                  design= ~  Condition)
dds_ost


# Многофакторные дизайны (пока мимо, с учетом батча) вероятнее всего, я уже учитываю донора выше
colData(dds_ost)
ddsMF <- dds_ost
design(ddsMF) <- formula(~ Donor + Condition)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
head(resMF)
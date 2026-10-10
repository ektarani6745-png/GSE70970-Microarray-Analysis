install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("GEOquery")
library(limma)
getwd()
data_raw <- read.delim("GSE70970_series_matrix.txt.gz",
                       comment.char = "!",
                       stringsAsFactors = FALSE)
rownames(data_raw) <- data_raw[,1]
data <- data_raw[,-1]

data <- as.matrix(data)
mode(data) <- "numeric"
dim(data)
head(data[,1:5])
colnames(data)
group <- factor(c(
  "Control","Control","Control",
  "Treated","Treated","Treated"
))
design <- model.matrix(~group)
design
data_norm <- normalizeBetweenArrays(data)
fit <- lmFit(data_norm, design)
fit <- eBayes(fit)
ncol(data)
length(group)
meta <- readLines("GSE70970_series_matrix.txt.gz")
meta[grep("characteristics_ch1", meta)]
head(meta[grep("characteristics_ch1", meta)], 20)
group_info <- meta[grep("characteristics_ch1", meta)]

group <- ifelse(grepl("control", group_info, ignore.case = TRUE),
                "Control", "Treated")

group <- factor(group)
length(group)
ncol(data)
group_line <- meta[grep("^!Sample_characteristics_ch1", meta)]
group_line
group_line <- group_line[1]
group_info <- unlist(strsplit(group_line, "\t"))
group_info <- group_info[-1]
length(group_info)
ncol(data)
group <- ifelse(grepl("control", group_info, ignore.case = TRUE),
                "Control", "Treated")

group <- factor(group)
design <- model.matrix(~group)

data_norm <- normalizeBetweenArrays(data)

fit <- lmFit(data_norm, design)
fit <- eBayes(fit)

results <- topTable(fit, coef=2, number=Inf)
results <- topTable(fit, coef=2, number=Inf)
head(results)
table(group)
head(group_info, 30)
group <- ifelse(grepl("distant.event: 0", group_info),
                "No_Metastasis", "Metastasis")

group <- factor(group)
table(group)
design <- model.matrix(~group)

data_norm <- normalizeBetweenArrays(data)

fit <- lmFit(data_norm, design)
fit <- eBayes(fit)

results <- topTable(fit, coef=2, number=Inf)
head(results)
sig <- results[results$adj.P.Val < 0.05, ]
nrow(sig)
up <- results[results$logFC > 1 & results$adj.P.Val < 0.05, ]
down <- results[results$logFC < -1 & results$adj.P.Val < 0.05, ]

nrow(up)
nrow(down)
write.csv(results, "All_Genes.csv")
write.csv(up, "Upregulated.csv")
write.csv(down, "Downregulated.csv")
plot(results$logFC, -log10(results$P.Value),
     pch=20, main="Volcano Plot")

with(results, points(logFC[adj.P.Val < 0.05],
                     -log10(P.Value[adj.P.Val < 0.05]),
                     col="red", pch=20))
top_genes <- rownames(results)[1:50]
heatmap(data_norm[top_genes, ])
getwd()
png("volcano_plot.png", width=800, height=600)

plot(results$logFC, -log10(results$P.Value),
     pch=20, main="Volcano Plot")

with(results, points(logFC[adj.P.Val < 0.05],
                     -log10(P.Value[adj.P.Val < 0.05]),
                     col="red", pch=20))

dev.off()
png("heatmap.png", width=800, height=600)

top_genes <- rownames(results)[1:50]
heatmap(data_norm[top_genes, ])

dev.off()
sig <- results[results$adj.P.Val < 0.05, ]
write.csv(sig, "Significant_Genes.csv")
save.image("analysis_workspace.RData")
load("analysis_workspace.RData")
savehistory("analysis_code.Rhistory")

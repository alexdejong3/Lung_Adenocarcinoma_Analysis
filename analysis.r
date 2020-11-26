library(edgeR)
library(stats)
library(cluster)
library(factoextra)
library(devtools)
library(ggfortify)
library(ggplot2)
library(maftools)

# There are 510 diseased samples
data <- read.delim("data/data_RNA_Seq_v2_expression_median.txt")
rownames(data) <- make.names(data$Hugo_Symbol, unique = TRUE)
data <- subset(data, select = -c(Entrez_Gene_Id, Hugo_Symbol))
data <- as.data.frame(t(as.matrix(data)))
data[["Normal"]] = FALSE

# There are 59 normal samples
normal <- read.delim("data/data_RNA_Seq_v2_mRNA_median_normals.txt")
rownames(normal) <- make.names(normal$Hugo_Symbol, unique = TRUE)
normal <- subset(normal, select = -c(Entrez_Gene_Id, Hugo_Symbol))
normal <- as.data.frame(t(as.matrix(normal)))
normal[["Normal"]] = TRUE

# Filter lowly expressed genes
samples <- rbind(data, normal)
normalvector <- samples$Normal
sums <- colSums(samples != 0)
kept <- which(sums > 500)
samples <- samples[,kept]

# Voom
mm <- model.matrix(~0 + normalvector)
voom <- voom(samples_t, mm, plot = T)
samples_t <- voom$E
samples <- as.data.frame(t(as.matrix(voom$E)))

# Kmeans
kmeans <- kmeans(samples, 3)
gene.pca <- prcomp(samples)
autoplot(gene.pca, data = samples, colour=kmeans$cluster)

# Limma
design <- model.matrix(~ normalvector)
fit <- lmFit(samples_t, design)
fit <- eBayes(fit)
topTable(fit)

# Clinical Data
clinicalsample <- read.delim("data/data_clinical_sample.txt")
clinicalpatient <- read.delim("data/data_clinical_patient.txt")
clinical <- cbind(clinicalsample, clinicalpatient)
names(clinical)[names(clinical) == "Sample.Identifier"] <- "Tumor_Sample_Barcode"

# Mutations
mut <- read.maf("data/data_mutations_extended.txt", clinicalData = clinical)
mut

plotmafSummary(maf = mut, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = mut, top = 15)

somaticInteractions(maf = mut, top = 25, pvalue = c(0.05, 0.1))

#dgi = drugInteractions(maf = mut, fontSize = 0.75)

OncogenicPathways(maf = mut)
PlotOncogenicPathways(maf = mut, pathways = "RTK-RAS")





# Load the RLassoCox package and required datasets
BiocManager::install("RLassoCox")
library(RLassoCox)
data(mRNA_matrix)
data(dGMMirGraph)
data(survData) 

# Read and prepare mRNA expression data
library(readr)
TCGA_UVM_htseq_counts <- read_delim("C:/Users/Genmol1/Downloads/TCGA-UVM.htseq_counts.tsv/TCGA-UVM.htseq_counts.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
expression <- as.data.frame(TCGA_UVM_htseq_counts)
rownames(expression) <- expression$Ensembl_ID
expressiondata <- expression[, 2:81]
DE_Z = t(scale(t(expressiondata)))
expresionlista <- t(DE_Z)

# Read and prepare survival data
TCGA_UVM_survival_data <- read_delim("C:/Users/Genmol1/Downloads/TCGA-UVM.survival_data.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     col_types = cols(OS = col_factor(levels = c("0", "1")), OS.time = col_integer()), 
                                     trim_ws = TRUE)
survival <- as.data.frame(TCGA_UVM_survival_data)
survival$status <- ifelse(survival$OS == 0, "FALSE", "TRUE")
survival$time <- survival$OS.time
survivalcompleto <- survival[, c("status", "time")]
rownames(survival) <- survival$sample
survivalcompleto$status <- as.logical(survivalcompleto$status)

# Check intersection between expression matrix and survival table
intersect(colnames(TCGA_UVM_htseq_counts), rownames(survival))

# Run LASSO Cox regression with all genes
mod <- RLassoCox(x=mRNA_matrix, y=survData, globalGraph=dGMMirGraph)

# Read gene list for analysis
library(readxl)
MEC_genes <- read_excel("C:/Users/Genmol1/Downloads/ECM TIMER2 (1).xlsx", sheet = "final")

# Filter expression matrix with gene list
colnames(expresionorganizada) <- gsub("\\..*","",colnames(expresionorganizada))
matrizMEC <- expresionorganizada[, MEC_genes$EnsemblID]                    
all((rownames(survivalcompleto) == rownames(matrizMEC)) == TRUE)
MEC_genes$number <- seq(1:43)
colnames(matrizMEC) <- MEC_genes$number
mod3 <- RLassoCox(x=matrizMEC, y=survivalcompleto, globalGraph=dGMMirGraph)
rownames(mRNA_matrix2) <- rownames(matrizMEC)

# Check results
mod3$PT
head(coef(mod3$glmnetRes, s = 0.2))

# Perform cross-validation with 5 folds
cv.mod <- cvRLassoCox(x=matrizMEC, y=survivalcompleto, globalGraph=dGMMirGraph, nfolds = 5)

# Plot cross-validation results
plot(cv.mod$glmnetRes, xlab = "log(lambda)")
cv.mod$glmnetRes$lambda.min
cv.mod$glmnetRes$lambda.1se

# Identify non-zero coefficients
coef.min <- coef(cv.mod$glmnetRes, s = "lambda.min")
nonZeroIdx <- which(coef.min[,1] != 0)
features <- rownames(coef.min)[nonZeroIdx]
features.coef <- coef.min[nonZeroIdx]
names(features.coef) <- features

# Print features with non-zero coefficients
features.coef

# Merge coefficient data with gene information
coeficientesgenes <- as.data.frame(features.coef)
coeficientesgenes$number <- rownames(coeficientesgenes)
tablafinalLAssocox <- merge(coeficientesgenes, MEC_genes, by="number")
write.table(tablafinalLAssocox, "coeficientesdelasso.txt", row.names = FALSE)

# Save plots
png("grafico1LASSO.png", units="in", width=8, height=8, res=600)
plot(cv.mod$glmnetRes, xlab = "log(lambda)")
dev.off()
png("grafico2LASSO.png", units="in", width=8, height=8, res=600)
plot(mod3$glmnetRes)
dev.off()
png("grafico3LASSO.png", units="in", width=8, height=8, res=600)
plot(mod3$glmnetRes, "lambda")
dev.off()

# Save dataframes
write.table(survivalcompleto, "supervivenciaLASSO.txt")
write.table(MEC_genes, "informaciongenesLASSO.txt", row.names = FALSE)
write.table(matrizMEC, "matriz_numeros_lasso.txt")
write.table(lamatrizdeexpresion, "matriz_ensembl_ID_lasso.txt")

# In the essay
head(coef(mod$glmnetRes, s = 0.2))

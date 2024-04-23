#lasso-cox in R
#Install the package and upload the image for running the algorithm
BiocManager::install("RLassoCox")
library(RLassoCox)
data(mRNA_matrix)
data(dGMMirGraph)
data(survData) 
#read the data, if you have a normalized counts matrix (by Deseq2 or EdgeR)
library(readr)
TCGA_UVM_htseq_counts <- read_delim("C:/Users/Genmol1/Downloads/TCGA-UVM.htseq_counts.tsv/TCGA-UVM.htseq_counts.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
expression <- as.data.frame(TCGA_UVM_htseq_counts)
#add the IDs to rownames in case the structure of the file does not contain that
rownames(expression) <- expression$Ensembl_ID
#remove the ID column
expressiondata <- expression[, 2:81]
#scale the data, necessary for the analysis 
DE_Z=t(scale(t(expressiondata)))
expresionlista <- t(DE_Z)
#expresionorganizada <- expresionlista[rownames(survivalcompleto), 1:100]
# read the survival data
library(readr)
TCGA_UVM_survival_data <- read_delim("C:/Users/Genmol1/Downloads/TCGA-UVM.survival_data.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     col_types = cols(OS = col_factor(levels = c("0", 
                                                                                 "1")), OS.time = col_integer()), 
                                     trim_ws = TRUE)
#View(TCGA_UVM_survival_data)
survival <-  as.data.frame(TCGA_UVM_survival_data)
#convert 0s and 1s to FALSE and TRUE
survival$status <- ifelse(survival$OS == 0, "FALSE", "TRUE")
survival$time <- survival$OS.time
#create a subset of survival data with only outcome and time
survivalcompleto <- survival[, c("status", "time")]
#extract the samples IDs present in the survival dataframe
rownames(survival) <- survival$sample
#convert the outcome to a factor
survivalcompleto$status <- as.logical(survivalcompleto$status)
#estimate the intersection between expression matrix and survival table to check you have all samples in both
intersect(colnames(TCGA_UVM_htseq_counts), rownames(survival))
#you can run the LASSO cox for all the genes in the matrix. 
#mod <- RLassoCox(x=expresionorganizada[,1:100], y=survivalcompleto, globalGraph=dGMMirGraph)
mod <- RLassoCox(x=mRNA_matrix, y=survData, globalGraph=dGMMirGraph)
#read your gene list for the analysis
library(readxl)
MEC_genes <- read_excel("C:/Users/Genmol1/Downloads/ECM TIMER2 (1).xlsx", 
                        sheet = "final")
#View(MEC_genes)
#removing some info from the IDs (optional)
colnames(expresionorganizada) <- gsub("\\..*","",colnames(expresionorganizada))
#filtering the expression matrix to include only our list of genes
matrizMEC <- expresionorganizada[, MEC_genes$EnsemblID]                    
#verify that the surival dataframe and expression matrix have the same individuals
all((rownames(survivalcompleto) == rownames(matrizMEC)) == TRUE)
#changing the IDs of the genes to allow running the algorithm
MEC_genes$number <- seq(1:43)
colnames(matrizMEC) <- MEC_genes$number
mod3 <- RLassoCox(x=matrizMEC, y=survivalcompleto, globalGraph=dGMMirGraph)
rownames(mRNA_matrix2) <- rownames(matrizMEC)
#run the Lasso cox regression
#checking the results
mod3$PT
head(coef(mod3$glmnetRes, s = 0.2))
#performing crossvalidation with 5 folds
cv.mod <- cvRLassoCox(x=matrizMEC, y=survivalcompleto,
                      globalGraph=dGMMirGraph, nfolds = 5)
plot(cv.mod$glmnetRes, xlab = "log(lambda)")
cv.mod$glmnetRes$lambda.min
cv.mod$glmnetRes$lambda.1se
coef.min <- coef(cv.mod$glmnetRes, s = "lambda.min")
coef.min
nonZeroIdx <- which(coef.min[,1] != 0)
features <- rownames(coef.min)[nonZeroIdx]
features
features.coef <- coef.min[nonZeroIdx]
names(features.coef) <- features
coef <- features
features.coef
#verifying the IDs that we changed
lamatrizdeexpresion <- expresionorganizada[, MEC_genes$EnsemblID]
plot(mod3$glmnetRes)
print(mod3$glmnetRes)
coeficientesgenes <- as.data.frame(features.coef)
mod3$PT
#ahora uniendo los datos de los coeficientes y los genes
coeficientesgenes$number <- rownames(coeficientesgenes)
tablafinalLAssocox <- merge(coeficientesgenes, MEC_genes, by="number")
write.table(tablafinalLAssocox, "coeficientesdelasso.txt", row.names = FALSE)
#guardando los graficos en fullHD
png("grafico1LASSO.png", units="in", width=8, height=8, res=600)
plot(cv.mod$glmnetRes, xlab = "log(lambda)")
dev.off()
png("grafico2LASSO.png", units="in", width=8, height=8, res=600)
plot(mod3$glmnetRes)
dev.off()
png("grafico3LASSO.png", units="in", width=8, height=8, res=600)
plot(mod3$glmnetRes, "lambda")
dev.off()
#guardando los dataframes por si algo
write.table(survivalcompleto, "supervivenciaLASSO.txt")
write.table(MEC_genes, "informaciongenesLASSO.txt", row.names = FALSE)
write.table(matrizMEC, "matriz_numeros_lasso.txt")
write.table(lamatrizdeexpresion, "matriz_ensembl_ID_lasso.txt")
#en el ensayo
head(coef(mod$glmnetRes, s = 0.2))

#Here we import the necessary libraries
library(DESeq2)
library(TCGAbiolinks)
library(SummarizedExperiment)

## PART 1: Data Download -----
#Downloading data from the TGCA. We download the data for the Lung Squamous Cell Carcinoma (LUSC) project with the following characteristics:
# - Gene expression data for primary tumor samples
# - Gene expression data for solid tissue normal samples
#The workflow type indicates the method used to generate the data. In this case, the data was generated using the STAR aligner and the counts were quantified.

#Here we create a directory to store the data for the LUSC project
proj <- "TCGA-LUSC"
dir.create(file.path(proj))

#Here we download the gene expression data for primary tumor samples using the GDCquery function
#First, we create a query object specifying the project, data category, data type, workflow type, and sample type
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")
#Next we download the data for the primary tumor samples
GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
#Now, we prepare the data for analysis
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata")
rna.expr.data.C <- assay(rna.data.C)
#It is of interest also to know the genes that are being analyzed
genes.info.C <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))

#Now, we do the same for the solid tissue normal samples
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")
#Next we download the data for the solid tissue normal samples
GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
#Now, we prepare the data for analysis
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata")
rna.expr.data.N <- assay(rna.data.N)
#It is of interest also to know the genes that are being analyzed
genes.info.N <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))



## PART 2: Data Preprocessing -----
#Here we preprocess the data to ensure that it is in a suitable format for analysis
#We make the following steps:

#First, we control for duplicates both in the normal and cancer samples
#ncol(rna.expr.data.N) == length(unique(substr(colnames(rna.expr.data.N), 1,12))) #TRUE
#ncol(rna.expr.data.C) == length(unique(substr(colnames(rna.expr.data.C), 1,12))) #FALSE

#Since we have duplicates in the cancer samples, we remove them
#Here we obtain the patients with unique samples
patients.C <- substr(colnames(rna.expr.data.C), 1,12)
unique.patients.C <- names(which(table(patients.C) == 1))
#We then get the index of these unique patients in the list of patients
idx.unique.pats <- match(unique.patients.C, substr(colnames(rna.expr.data.C), 1,12))

#We only keep the patients with one sample
expr.C <- as.data.frame(rna.expr.data.C[,idx.unique.pats])
expr.N <- as.data.frame(rna.expr.data.N)
#We also rename the patients in a shorter way
colnames(expr.C) <- substr(colnames(expr.C), 1,12)
colnames(expr.N) <- substr(colnames(expr.N), 1,12)

#Afterwards we control samples that do not have a cancerous sample to compare to
#match(setdiff(colnames(expr.N), colnames(expr.C)), colnames(expr.N)) #this is zero
#Since all normal samples have a corresponding cancerous sample, we proceed

#We also control for NaN values but we do not have any
#any(is.na(as.matrix(expr.C))) #FALSE
#any(is.na(as.matrix(expr.N))) #FALSE

#Here we make sure to consider only patients for which we have both normal and cancer samples
expr.C <- expr.C[, colnames(expr.N)]

#Now we combine both datasets
full.data <- cbind(expr.N, expr.C)
full.data <- data.frame(full.data)

## PART 3: Data Normalization ---
#Here we use the DESEQ2 package to normalize the data


metad <- rep("cancer", ncol(full.data))
metad[1:(ncol(full.data)/2)] <- "normal"
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1])

full.data <- cbind(rownames(full.data), full.data)

dds <- DESeqDataSetFromMatrix(countData=full.data, 
                              colData=metad, 
                              design= ~condition,
                              tidy=TRUE)

View(counts(dds))
dim(counts(dds))

# filtering: at least ten counts on 90% of patients? 
( 38*90 )/100
keep <- rowSums(counts(dds) >= 10) >= 34
dds <- dds[keep,]
dim(counts(dds))

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == 38) #no null rows


filtr.expr.n <- as.data.frame(normalized_counts[, 1:19])
filtr.expr.c <- as.data.frame(normalized_counts[, 20:38])
#cancerous sample names were added a ".1" in full.data because  
#they had the same names as the normal samples
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1,12)

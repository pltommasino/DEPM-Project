# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BiocGenerics")
# BiocManager::install("DESeq2")
# BiocManager::install("TCGAbiolinks")
# install.packages("NetworkToolbox")
# install.packages("DT")
# install.packages("psych")
# library(devtools)
# devtools::install_github("briatte/ggnet")

library(BiocGenerics) 
library(DESeq2)
library(psych)
library(NetworkToolbox)
library(ggplot2); library(ggnet)
library(GGally);library(sna);library(network)

# PART 1 

#1: Downloading data from the TGCA -------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)

proj <- "TCGA-LUSC"
dir.create(file.path(proj))

##### Data Primary Tumor
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")

GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata")
rna.expr.data.C <- assay(rna.data.C)

#View(BiocGenerics::as.data.frame(rowRanges(rna.data.C)))
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))

##### Solid Tissue Normal
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata"  )
rna.expr.data.N <- assay(rna.data.N)

genes.info2 <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))
all(na.omit(genes.info2) == na.omit(genes.info))

clinical.query<- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
#write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)

View(clinical.query)
table(clinical.query$ajcc_pathologic_stage)

boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical.query,
        col = "gold", main = "Title", xlab = "", ylab= "age", las=2 )


View(rna.expr.data.N)

dim(rna.expr.data.C)
dim(rna.expr.data.N)
length(unique(clinical.query$submitter_id))



#2: Data cleaning -----

ncol(rna.expr.data.C)
head(colnames(rna.expr.data.C))
head(substr(colnames(rna.expr.data.N), 1,12))

dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) #no duplicates 
dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) #duplicates!

patients.C <- substr(colnames(rna.expr.data.C), 1,12)
sort(table(patients.C)) 

unique.patients.C <- names(which(table(patients.C) == 1))
#let's get their index in the list of patients
idx.unique.pats <- match(unique.patients.C, substr(colnames(rna.expr.data.C), 1,12) )

#I'm only keeping patients with one sample
expr.C <- as.data.frame(rna.expr.data.C[,idx.unique.pats])
expr.N <- as.data.frame(rna.expr.data.N)

#let's rename patients in a shorter way
colnames(expr.C) <- substr(colnames(expr.C), 1,12)
unique(colnames(expr.C))

colnames(expr.N) <- substr(colnames(expr.N), 1,12)
unique(colnames(expr.N))

intersect(colnames(expr.N), colnames(expr.C)) #38 instead of 41
setdiff(colnames(expr.N), colnames(expr.C))
# 3 normal sample do not have a cancerous sample to compare to. Let's remove them
match(setdiff(colnames(expr.N), colnames(expr.C)), colnames(expr.N)) #idx to remove
#???????????????? expr.N <- expr.N[,-c(10,15,20)]

length(intersect(colnames(expr.N), colnames(expr.C)))
#normal samples are a subset of the cancer samples. great. 


#let's check the actual counts
typeof(expr.C[1,1]) #ok
any(is.na(expr.C)) #ok
any(is.nan(as.matrix(expr.C))) #ok

typeof(expr.N[1,1]) #ok
any(is.na(expr.N)) #ok
any(is.nan(as.matrix(expr.N))) #ok

#let's consider only patients for which we have both normal and cancer samples
expr.C <- expr.C[, colnames(expr.N)]
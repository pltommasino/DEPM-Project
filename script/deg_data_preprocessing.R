#Here we import the necessary libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(latex2exp)
library(TCGAbiolinks)
library(SummarizedExperiment)

## PART 1: Data Download -----
#Downloading data from the TGCA. We download the data for the Lung Squamous Cell Carcinoma (LUSC) project with the following characteristics:
# - Gene expression data for primary tumor samples
# - Gene expression data for solid tissue normal samples
#The workflow type indicates the method used to generate the data. In this case, the data was generated using the STAR aligner and the counts were quantified.

#First, we make sure we are in the correct working directory
setwd("~/Documents/Master/DigitalEpidemiology/DEPM-Project")
#Here we create a directory to store the data for the LUSC project
proj <- "TCGA-LUSC"
#dir.create(file.path(proj))

#Here we download the gene expression data for primary tumor samples using the GDCquery function
#First, we create a query object specifying the project, data category, data type, workflow type, and sample type
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")
#Next we download the data for the primary tumor samples
#GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
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
#GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
#Now, we prepare the data for analysis
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata")
rna.expr.data.N <- assay(rna.data.N)
#It is of interest also to know the genes that are being analyzed
genes.info.N <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))

#Finally, we download clinical data in case we need it at a later stage
clinical.query <- GDCquery_clinic(project = proj, type = "clinical", save.csv=FALSE)
#Here we save it
#write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)

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

#Now we combine both datasets. Since patients have the same names in both datasets, we add a suffic to the cancerous samples
colnames(expr.C) <- paste0(colnames(expr.C), ".C")
full.data <- cbind(expr.N, expr.C)
full.data <- data.frame(full.data)

## PART 3: Data Normalization ---
#Here we use the DESEQ2 package to normalize the data. To use this, we create a metadata dataframe with the condition of the samples
sample_number <- ncol(full.data)
metad <- rep("cancer", sample_number)
metad[1:(sample_number/2)] <- "normal"
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1])

#Now we create a DESeqDataSet object
full.data <- cbind(rownames(full.data), full.data)
dds <- DESeqDataSetFromMatrix(countData=full.data, 
                              colData=metad, 
                              design= ~condition,
                              tidy=TRUE)

#Here we filter the data to remove genes with low counts. Specifically, we remove genes with less than 10 counts in 100% of the samples
keep <- rowSums(counts(dds) >= 10) >= sample_number
dds <- dds[keep,]
#Now we normalize the data
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

#We separate the data into normal and cancerous samples
filtr.expr.n <- as.data.frame(normalized_counts[, 1:(sample_number/2)])
filtr.expr.c <- as.data.frame(normalized_counts[, (sample_number/2 +1):sample_number])

#Since cancerous sample names were added a ".C", we remove it
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1,12)

##PART 4: Gene Selection ----
#Now, we just rename the genes so we can interpret our results
genes.c <- intersect(rownames(filtr.expr.c), 
                     genes.info.C[ , "gene_id"])
genes.n <- intersect(rownames(filtr.expr.n),
                     genes.info.N[ , "gene_id"])

#Here we select the genes that are in the list of genes we are interested in
filtr.expr.c <- filtr.expr.c[genes.c,]
filtr.expr.n <- filtr.expr.n[genes.n,]

#Here we obtain their names
rownames(filtr.expr.c) <- genes.info.C[genes.c, "gene_name"]
rownames(filtr.expr.n) <- genes.info.N[genes.n, "gene_name"]

##PART 5: Differentially Expressed Genes (DEGs) ----
#To obtain the DEGs first we obtain the fold change
fc <-  log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n) ) 
names(fc) <- rownames(filtr.expr.c)

#Now we perform a t-test between the normal and cancerous samples to obtain the p-values
pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(as.numeric(filtr.expr.c[i,]), as.numeric(filtr.expr.n[i,]), paired = T ))$p.value)
#We need to adjust the p-values to account for multiple testing.
pval.fc.fdr <- p.adjust(pval.fc, method="fdr")

#Here we create a table with the results of the t-test and fold change of the genes
expr.table <- data.frame(cbind(fc, pval.fc.fdr))
#We select as DEGs those genes with a fold change greater than 1.2 and a p-value less than 0.01
deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.2 & expr.table$pval.fc.fdr <0.01,]) 
#Here we save the list
write.table(expr.table[deg.genes,], file = "data/DEG.csv", sep = ";")
#We can also save the normalized counts for the DEGs
deg.expr.c <- filtr.expr.c[deg.genes, ]
deg.expr.n <- filtr.expr.n[deg.genes, ]
#We save them individually as .RData files
save(deg.expr.c, file = "data/deg_expr_c.RData")
save(deg.expr.n, file = "data/deg_expr_n.RData")

##PART 6: Visualization (Volcano Plot) ----
#Here we create a volcano plot to visualize the DEGs. First we display the genes that are not DEGs
expr.table$diffexpressed <- "No"
expr.table$diffexpressed[expr.table$fc >= 1.2 & expr.table$pval.fc.fdr < 0.01] <- "Up"
expr.table$diffexpressed[expr.table$fc <= -1.2 & expr.table$pval.fc.fdr < 0.01] <- "Down"
expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)

#We plot
ggplot(data=expr.table, aes(x=fc, y=-log10(pval.fc.fdr), col=diffexpressed), label=rownames(expr.table)) + 
  theme_bw() +
  geom_point(alpha=0.8) +
  #Here we add the labels for the DEGs that have a fold change greater than 1.2 and a p-value less than 0.05
  geom_text_repel(data=expr.table[abs(expr.table$fc) >= 4.5& expr.table$pval.fc.fdr <=1e-15,], aes(label=rownames(expr.table[abs(expr.table$fc) >= 4.5 & expr.table$pval.fc.fdr <=1e-15,])), box.padding = 0.5, point.padding = 0.5, segment.size = 0.5, size=3, show.legend = FALSE, max.overlaps = Inf) +
  xlim(c(-7,7)) +
  labs(x=TeX(r"(Fold Change ($\log_2$))"),
       y=TeX(r"(-$\log_{10}(p_{adj})$)"),
       color="Expressed") +
  #Here we add the colors for the differentially expressed genes
  scale_color_manual(values=c("#fc8d59", "gray", "#91bfdb"))+
  geom_hline(yintercept=-log10(0.01), col="red", linetype="dashed")+
  geom_vline(xintercept=1.2, col="red", alpha=0.8, linetype="dashed")+
  geom_vline(xintercept=-1.2, col="red", alpha=0.8, linetype="dashed")

#We save the plot in high resolution as a pdf and tight layout
ggsave("figures/volcano_plot.pdf", dpi=300, device="pdf" , width=7, height=5, useDingbats=FALSE)

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
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)

# PART 1 

#1: Data ----

#Downloading data from the TGCA

setwd("~/Documents/GitHub/DEPM-OLD")
proj <- "TCGA-LUSC"
#dir.create(file.path(proj))

##### Data Primary Tumor
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")

#GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata")
rna.expr.data.C <- assay(rna.data.C)

#View(BiocGenerics::as.data.frame(rowRanges(rna.data.C)))
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))

##### Solid Tissue Normal
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

#GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata")
rna.expr.data.N <- assay(rna.data.N)

genes.info2 <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))
all(na.omit(genes.info2) == na.omit(genes.info))

clinical.query<- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
#write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)

#View(clinical.query)
table(clinical.query$ajcc_pathologic_stage)

boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical.query,
        col = "gold", main = "Title", xlab = "", ylab= "age", las=2 )


#View(rna.expr.data.N)

# dim(rna.expr.data.C)
# dim(rna.expr.data.N)
# length(unique(clinical.query$submitter_id))



#Data cleaning

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



#Normalizing data with Deseq2

#detalied explanation: https://hbctrzaining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#youtube explanation: https://www.youtube.com/watch?v=UFB993xufUU

all(rownames(expr.C) == rownames(expr.N))
full.data <- cbind(expr.N, expr.C)

dim(full.data) #60660    38 -------> 60660   102
full.data <- data.frame(full.data)

metad <- rep("cancer", 102)
metad[1:51] <- "normal"
metad
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1])
full.data <- cbind(rownames(full.data), full.data)

dds <- DESeqDataSetFromMatrix(countData=full.data, 
                              colData=metad, 
                              design= ~condition,
                              tidy=TRUE)

#View(counts(dds))
dim(counts(dds))

# filtering: at least ten counts on 90% of patients? 
( 102*90 )/100
keep <- rowSums(counts(dds) >= 10) >= 91
dds <- dds[keep,]
dim(counts(dds))

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == 102) #no null rows


filtr.expr.n <- as.data.frame(normalized_counts[, 1:51])
filtr.expr.c <- as.data.frame(normalized_counts[, 52:102])
#cancerous sample names were added a ".1" in full.data because  
#they had the same names as the normal samples
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1,12)

#2: Differentially Expressed Genes (DEGs) ----

#Gene selection
genes <- read.csv2("targets.csv", row.names = 1)
genes <- genes[,1]
head(genes) #gene symbols, not ensemble ids

#their enseble ids 
genes.info[ genes.info$gene_name %in% genes ,"gene_id"]
#not all of them might be among our genes!

genes.c <- intersect(rownames(filtr.expr.c), 
                     genes.info[ genes.info$gene_name %in% genes , "gene_id"]   ) 

genes.n <- intersect(rownames(filtr.expr.n),  
                     genes.info2[ genes.info2$gene_name %in% genes , "gene_id"]   )  

setdiff(genes.c, genes.n)

#length(genes)
#length(genes.c)
#length(genes.n)

filtr.expr.n <- filtr.expr.n[genes.n, ]
filtr.expr.c <- filtr.expr.c[genes.c, ]

rownames(filtr.expr.n) <- genes.info2[genes.n, "gene_name"]
rownames(filtr.expr.c) <- genes.info[genes.c, "gene_name"]


#Differentially expressed genes (DEGs) (+ brief enrichment overview)

#what are DEGs?
fc <-  log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n) ) 
names(fc) <- rownames(filtr.expr.c)
head(fc)

#what is the fold change?
pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(as.numeric(filtr.expr.c[i,]), as.numeric(filtr.expr.n[i,]), paired = T ))$p.value)
pval.fc.fdr <- p.adjust(pval.fc, method="fdr")

expr.table <- data.frame(cbind(fc, pval.fc.fdr))
expr.table[,1] <- round(expr.table[,1],2)
#expr.table[,2] <- format(expr.table[,2], digits = 4) nice for printing but it converts to string

deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.2 & expr.table$pval.fc.fdr <=0.05,]) 
deg.genes

head(expr.table[deg.genes,], 10)
# fc  pval.fc.fdr
# GCLC     3.19 5.157093e-06
# ENPP4   -2.12 5.798913e-20
# PRSS22   1.88 5.456300e-04
# ALDH3B1 -2.96 4.012309e-17
# DBF4     1.97 4.667448e-16
# E2F2     2.99 1.872244e-15
# JARID2   1.32 3.755466e-13
# NCAPD2   2.32 3.496074e-14
# BRCA1    2.23 2.484447e-14
# SNAI2    2.11 2.822504e-11

#volcano plot

expr.table$diffexpressed <- "NO";
expr.table$diffexpressed[expr.table$fc >= 1.2 & expr.table$pval.fc.fdr <= 0.05] <- "UP"
expr.table$diffexpressed[expr.table$fc <= -1.2 & expr.table$pval.fc.fdr <= 0.05] <- "DOWN"
head(expr.table)

expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)
#summary(expr.table$diffexpressed)

ggplot(data=expr.table, aes(x=fc, y=-log10(pval.fc.fdr), col=diffexpressed))+  
  geom_point() +
  xlab("fold change (log2)") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_vline(xintercept=1.2, col="red")+
  geom_vline(xintercept=-1.2, col="red")


# print and enrichment 
#cat(deg.genes , sep = "\n")
#enrichR



#3: Co-expression networks -----

#Adjacency matrices of co-expression networks
#### --- Cancer network

#Creating the correlation matrix for the “cancer” group
cor.mat.c <- corr.test(t(filtr.expr.c), use = "pairwise", 
                       method = "spearman",adjust="fdr",ci=FALSE)

#the correlation matrix (\rho).
rho.c <- cor.mat.c$r
diag(rho.c) <- 0

#the matrix of p-values associated with the correlations.
qval.c <- cor.mat.c$p
#qvals are reported on the upper triangle only
qval.c[lower.tri(qval.c)] <- t(qval.c)[lower.tri(qval.c)]

adj.mat.c <- ifelse(abs(rho.c) >= 0.65, 1, 0) * (qval.c <= 1e-4)
#sum(adj.mat.c) #3526


#### --- Normal network 
cor.mat.n <- corr.test(t(filtr.expr.n), use = "pairwise", 
                       method = "spearman",adjust="fdr",ci=FALSE)

rho.n <- cor.mat.n$r
diag(rho.n) <- 0
qval.n <- cor.mat.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]

adj.mat.n <- ifelse(abs(rho.n) >= 0.65, 1, 0) * (qval.n <= 1e-4)


#Co-expression networks

#Cancer network 
net.c <- network(adj.mat.c, matrix.type="adjacency",ignore.eval = FALSE, 
                 names.eval = "weights", directed = F)

network.density(net.c)
#0.001472384
network.size(net.c)
#1548
network.edgecount(net.c) 
#1763
nrow(component.largest(net.c, result = "graph")) 
#551
clustcoeff(adj.mat.c, weighted = FALSE)$CC
#0.1225775

sum(adj.mat.c != 0) / 2
#1763
#how many positive/negative correlations? 
sum(adj.mat.c > 0) / 2
#1763
sum(adj.mat.c < 0) / 2
#0

degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)
head(degree.c,10)
sum(degree.c == 0) 
#814 unconnected nodes 

hist(degree.c)
x <- quantile(degree.c[degree.c>0],0.95) #how big is the degree of the most connected nodes?
x
# 95% 
# 18.35
hist(degree.c)
abline(v=x, col="red")

hubs.c <- degree.c[degree.c>=x]
head(names(hubs.c),10) #

net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c),"hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")
network::set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, "red", "blue"))

#coord.c <- gplot.layout.fruchtermanreingold(net.c, NULL)
#net.c %v% "x" = coord.c[, 1]
#net.c %v% "y" = coord.c[, 2]

ggnet2(net.c, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

#this is extremely dense... what if we lower the pval threshold?
#what if it's even lower?
#adj.mat.c <- rho.c * (qval.c <= 1e-3)
#adj.mat.c <- rho.c * (qval.c <= 1e-4) #too much?


#might be useful to look at negative and postive edges separately:
net.c1 <- network(adj.mat.c* (adj.mat.c > 0),
                  matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

ggnet2(net.c1, color = "deepskyblue3", alpha = 0.7, size = 2, 
       edge.color = "red", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 


net.c2 <- network(adj.mat.c* (adj.mat.c < 0),
                  matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

ggnet2(net.c2, color = "deepskyblue3", alpha = 0.7, size = 2, 
       edge.color = "blue", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 


#Normal network 
net.n <- network(adj.mat.n, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.n)
#0.0178657
network.size(net.n)
#1548
network.edgecount(net.n)
#42784
clustcoeff(adj.mat.n, weighted = FALSE)$CC
#0.3631008
nrow(component.largest(net.n, result = "graph")) 
#1412

sum(adj.mat.n != 0) /2
#21392
#how many positive/negative correlations? 
sum(adj.mat.n > 0) /2
#21392
sum(adj.mat.n < 0) /2
#0

degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = T)
head(degree.n,10)
sum(degree.n == 0) 
#132 unconnected nodes 

hist(degree.n)
y <- quantile(degree.n[degree.n>0],0.95) #how big is the degree of the most connected nodes?
y
# 95% 
# 87 
hist(degree.n)
abline(v=y, col="red")

hubs.n <- degree.n[degree.n>=y]
head(names(hubs.n),10) #

net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n),"hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))

ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

#this is extremely dense... what if we change the pval threshold?
adj.mat.n <- rho.c * (qval.c <= 1e-3) 
adj.mat.n <- rho.c * (qval.c <= 1e-4) #too much?


intersect(names(hubs.c), names(hubs.n))


#Plotting the hub subnetwork
hubs.c
hubs.c.ids <- vector("integer",length(hubs.c))
for (i in 1:length(hubs.c)){hubs.c.ids[i] <- match(names(hubs.c)[i],rownames(adj.mat.c))}
hubs.c.ids

#identifying the neighborhood
hubs.c.neigh <- c()
for (f in hubs.c.ids){
  hubs.c.neigh <- append(hubs.c.neigh, get.neighborhood(net.c, f))
}

hubs.c.neigh <- unique(hubs.c.neigh)
hubs.c.neigh
hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh,])
subnet <- unique(c(names(hubs.c), hubs.c.neigh.names))

#creating the subnetwork
hub.c.adj <- adj.mat.c[subnet, subnet]

head(rownames(hub.c.adj))
head(colnames(hub.c.adj))

net.hub <- network(hub.c.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.hub)
#0.0569854
sum(hub.c.adj > 0 )
#2068
sum(hub.c.adj < 0)
#0

net.hub %v% "type" = ifelse(network.vertex.names(net.hub) %in% names(hubs.c),"hub", "non-hub")
net.hub %v% "color" = ifelse(net.hub %v% "type" == "non-hub", "deepskyblue3", "tomato")
set.edge.attribute(net.hub, "ecolor", ifelse(net.hub %e% "weights" > 0, "red", "blue"))

ggnet2(net.hub,  color = "color",alpha = 0.9, size = 2, 
       edge.color = "ecolor", edge.alpha = 0.9,  edge.size = 0.15, 
       node.label = names(hubs.c), label.color = "black", label.size = 3)+
  guides(size = "none") 


#4: Differential Co-expressed Network ---------------------------------------

#Fisher Z-transformation
##for Cancer
z.cancer <- 0.5 * log((1 + rho.c) / (1 - rho.c))
##for Normal
z.normal <- 0.5 * log((1 + rho.n) / (1 - rho.n))

#Z-scores to evaluate the correlation
##Formula Z
z.diff.formula <- function(zc, zn, fec, fen){
  n.cancer <- ncol(fec)
  n.normal <- ncol(fen)
  out <- (zc - zn) / sqrt((1 / (n.cancer - 3)) + (1 / (n.normal - 3)))
  return(out)
}
##Calculate
z.diff <- z.diff.formula(z.cancer, z.normal, filtr.expr.c, filtr.expr.n)

#Binary adjacency matrix with aij=0 if |Z|<3
adj.diff.mat <- ifelse(abs(z.diff) >= 3, 1, 0)
#diag(adj.diff.mat) <- 0


degree.diff.mat <- rowSums(adj.diff.mat)
####Taken from Damian file
hist(degree.diff.mat, breaks = 30, main = "Degree Distribution (DCN)", xlab = "Degree", ylab = "Frequency")
####

#Differential Network 
net.diff <- network(adj.diff.mat, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.diff)
#0.05530584
network.size(net.diff)
#1548
network.edgecount(net.diff)
#132444
clustcoeff(adj.diff.mat, weighted = FALSE)$CC
#0.1387713
nrow(component.largest(net.diff, result = "graph")) 
#1548

sum(adj.diff.mat != 0) /2
#66222
#how many positive/negative correlations? 
sum(adj.diff.mat > 0) /2
#66222
sum(adj.diff.mat < 0) /2
#0

degree.diff <- rowSums(adj.diff.mat != 0)
names(degree.diff) <- rownames(adj.diff.mat)
degree.diff <- sort(degree.diff, decreasing = T)
head(degree.diff,10)
sum(degree.diff == 0) 
#0 unconnected nodes 

hist(degree.diff)
z <- quantile(degree.diff[degree.diff>0],0.95) #how big is the degree of the most connected nodes?
z
# 95% 
# 174 
hist(degree.diff)
abline(v=z, col="red")

hubs.diff <- degree.diff[degree.diff>=z]
head(names(hubs.diff),10) #

net.diff %v% "type" = ifelse(network.vertex.names(net.diff) %in% names(hubs.diff),"hub", "non-hub")
net.diff %v% "color" = ifelse(net.diff %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.diff, "edgecolor", ifelse(net.diff %e% "weights" > 0, "red", "blue"))

ggnet2(net.diff, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

#this is extremely dense... what if we change the pval threshold?
adj.mat.diff <- rho.c * (qval.c <= 1e-3) 
adj.mat.diff <- rho.c * (qval.c <= 1e-4) #too much?


#intersect(names(hubs.c), names(hubs.n))


#Plotting the hub subnetwork
hubs.diff
hubs.diff.ids <- vector("integer",length(hubs.diff))
for (i in 1:length(hubs.diff)){hubs.diff.ids[i] <- match(names(hubs.diff)[i],rownames(adj.mat.diff))}
hubs.diff.ids

#identifying the neighborhood
hubs.diff.neigh <- c()
for (f in hubs.diff.ids){
  hubs.diff.neigh <- append(hubs.diff.neigh, get.neighborhood(net.diff, f))
}

hubs.diff.neigh <- unique(hubs.diff.neigh)
hubs.diff.neigh
hubs.diff.neigh.names <- rownames(adj.mat.diff[hubs.diff.neigh,])
subnet.diff <- unique(c(names(hubs.diff), hubs.diff.neigh.names))

#creating the subnetwork
hub.diff.adj <- adj.mat.diff[subnet.diff, subnet.diff]

head(rownames(hub.diff.adj))
head(colnames(hub.diff.adj))

net.hub.diff <- network(hub.diff.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
network.density(net.hub)
#0.0569854

sum(hub.diff.adj > 0)
#3306
sum(hub.diff.adj < 0)
#180

net.hub.diff %v% "type" = ifelse(network.vertex.names(net.hub.diff) %in% names(hubs.diff),"hub", "non-hub")
net.hub.diff %v% "color" = ifelse(net.hub.diff %v% "type" == "non-hub", "deepskyblue3", "tomato")
set.edge.attribute(net.hub.diff, "ecolor", ifelse(net.hub.diff %e% "weights" > 0, "red", "blue"))

ggnet2(net.hub.diff,  color = "color",alpha = 0.9, size = 2, 
       edge.color = "ecolor", edge.alpha = 0.9,  edge.size = 0.15, 
       node.label = names(hubs.diff), label.color = "black", label.size = 3)+
  guides(size = "none") 


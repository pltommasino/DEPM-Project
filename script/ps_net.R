library(BiocGenerics) 
library(TCGAbiolinks)
library(latex2exp)
library(ggsurvfit)
library(survminer)
library(survival)
library(multtest)
library(SNFtool)
library(ggplot2)
library(igraph)
library(ggnet)
library(proxy)
library(dplyr)
library(tidyr)
library(stats)


##PART 1: Building the Patient Similarity Networks
#Here we build the patient similarity networks for the LUSC project based on the gene expression data
#First, we load the pre-processed, normalized data we obtained for the DEG analysis
setwd("~/Documents/Master/DigitalEpidemiology/DEPM-Project")
#We load the .RData file containing the cancer data
load("data/deg_expr_c.RData")
#We load the .RData file containing the normal data
load("data/deg_expr_n.RData")

#Here we normalize the columns of the data so the columns sum to 1
normalize <- function(x) x/sum(x)
deg.expr.c <- apply(deg.expr.c, 2, normalize)
deg.expr.n <- apply(deg.expr.n, 2, normalize)

#Now we obtain the similarity matrices for the cancer and normal samples using the cosine similarity
sim.matrix.c <- as.matrix(proxy::simil(deg.expr.c, method = "cosine",  by_rows = FALSE))
sim.matrix.n <- as.matrix(proxy::simil(deg.expr.n, method = "cosine",  by_rows = FALSE))
#We ensure the diagional is 0 to avoid self-loops
diag(sim.matrix.c) <- 0
diag(sim.matrix.n) <- 0

#We construct the adjacency matrix by setting a threshold of 0.8 for the (dis)similarity
adj.matrix.c <- ifelse(sim.matrix.c >= 0.8, 1, 0)
adj.matrix.n <- ifelse(sim.matrix.n >= 0.8, 1, 0)

#We build the patient similarity networks
g.c <- graph_from_adjacency_matrix(adj.matrix.c, mode = "undirected")
g.n <- graph_from_adjacency_matrix(adj.matrix.n, mode = "undirected")

#Now we find the communities in the networks using Louvain's algorithm
communities.c <- cluster_louvain(g.c)
communities.n <- cluster_louvain(g.n)
#Here we add to each vertex the community it belongs to
V(g.c)$community <- communities.c$membership
V(g.n)$community <- communities.n$membership

#We also add the degree of each vertex to the network
V(g.c)$degree <- igraph::degree(g.c)
V(g.n)$degree <- igraph::degree(g.n)

#We add as label of the vertex the last 4 digits of the string
V(g.c)$label <- substr(V(g.c)$name, 9, 12)
V(g.n)$label <- substr(V(g.n)$name, 9, 12)

#We plot the networks
ggraph(g.c, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = community), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the communities
  scale_color_viridis_c(option = "plasma")+
  #We add the labels
  geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = Inf)

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/patient_similarity_network_c.pdf", dpi=300, device="pdf", useDingbats=FALSE)

ggraph(g.n) + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = community), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the communities
  scale_color_viridis_c(option = "plasma")+
  #We add the labels
  geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = Inf)

#We create a dataframe with the community membership of the patients
community.df <- data.frame(patient = gsub("\\.", "-", V(g.c)$name), community = V(g.c)$community)
rownames(community.df) <- community.df$patient

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/patient_similarity_network_n.pdf", dpi=300, device="pdf", useDingbats=FALSE)

##PART 2: Clinical Enrichment of the Cancer Modules
#Now we perform clinical enrichment analysis for the cancer modules
#First we download the clinical data
clinical.data <- GDCquery_clinic(project = "TCGA-LUSC", type = "clinical", save.csv=FALSE)
#We rename the rownames of the clinical data to match the submitter_id
rownames(clinical.data) <- substr(clinical.data$submitter_id, 1, 12)
#We filter the clinical data by the patients in the cancer network
#First we substitute the "." by "-" in the cancer names
cancer.names <- gsub("\\.", "-", V(g.c)$name)
#We filter the clinical data by the patients in the cancer network
clinical.data <- clinical.data[rownames(clinical.data) %in% cancer.names,]
#We add the community membership to the clinical data
clinical.data$community <- community.df[rownames(clinical.data),]$community

#We get the top 3 biggest communities and filter the clinical data by them
top.communities <- table(communities.c$membership) %>% sort(decreasing = TRUE) %>% names() %>% head(3)
#We select the rows of the clinical data that are in the top 3 communities
clinical.data <- clinical.data[clinical.data$community %in% top.communities, ]

#We preprocess some of the columns of the clinical data
#First we focus on important data
filtered.clinical.data <- clinical.data[,c("ajcc_pathologic_stage", "prior_malignancy",
                                           "race", "years_smoked", "gender", "age_at_index", "vital_status", "community")]
#We add a binary column for smoking status: 1 if years_smoked is not NA, 0 otherwise
filtered.clinical.data$smoking_status <- ifelse(is.na(filtered.clinical.data$years_smoked), 0, 1)

#Now we perform clinical enrichment analysis
enrichment_categorical <- function(data, feature, cluster_col = "community") {
  feature_enrichment <- list()
  
  for (cluster in unique(data[[cluster_col]])) {
    # Filter for the current cluster and others
    in_cluster <- data %>% filter(.data[[cluster_col]] == cluster)
    not_in_cluster <- data %>% filter(.data[[cluster_col]] != cluster)
    
    # Contingency table
    contingency <- table(data[[feature]], data[[cluster_col]] == cluster)
    
    # Perform chi-square test
    test <- chisq.test(contingency)
    feature_enrichment[[as.character(cluster)]] <- c(
      "Community" = cluster,
      "Feature" = feature,
      "p_value" = test$p.value
    )
  }
  # Convert to data frame
  feature_enrichment_df <- as.data.frame(do.call(rbind, feature_enrichment))
  return(feature_enrichment_df)
}

# Function for continuous feature enrichment (T-test)
enrichment_continuous <- function(data, feature, cluster_col = "community") {
  feature_enrichment <- list()
  
  for (cluster in unique(data[[cluster_col]])) {
    # Filter for the current cluster and others
    in_cluster <- data %>% filter(.data[[cluster_col]] == cluster)
    not_in_cluster <- data %>% filter(.data[[cluster_col]] != cluster)
    
    # Perform t-test
    test <- t.test(in_cluster[[feature]], not_in_cluster[[feature]])
    feature_enrichment[[as.character(cluster)]] <- c(
      "Community" = cluster,
      "Feature" = feature,
      "p_value" = test$p.value
    )
  }
  # Convert to data frame
  feature_enrichment_df <- as.data.frame(do.call(rbind, feature_enrichment))
  return(feature_enrichment_df)
}

# Analyze categorical features: Stage, Outcome, Gender, SmokingStatus
categorical_features <- c("ajcc_pathologic_stage", "prior_malignancy", "race", "gender", "vital_status", "smoking_status")
cat_results <- lapply(categorical_features, function(feature) {
  enrichment_categorical(filtered.clinical.data, feature)
})

# Analyze continuous feature: Age
age_results <- enrichment_continuous(filtered.clinical.data, "age_at_index")

# Combine results and adjust for multiple testing
all_results <- bind_rows(do.call(rbind, cat_results), as.data.frame(age_results))
all_results$p_adjusted <- p.adjust(as.numeric(all_results$p_value), method = "fdr")

ggplot(all_results, aes(x = -log10(p_adjusted), y = Feature, fill = Community)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = TeX(r"($-\log_{10}(p_{adj})$)"), y = "Feature", fill = "Community") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),   # Adjust the size of y-axis labels
    axis.text.x = element_text(size = 10),   # Adjust the size of x-axis labels
    axis.title = element_text(size = 12),    # Title size for both axes
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.ticks = element_line(size = 0.5)    # Set tick line size
  )+scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"))+
  scale_y_discrete(labels = c(vital_status = "Vital Status",
                              prior_malignancy = "Prior Malignancy",
                              ajcc_pathologic_stage = "AJCC Pathologic Stage",
                              gender ="Gender",
                              smoking_status = "Smoking Status",
                              race = "Race",
                              age_at_index = "Age"))
#We save the plot in high resolution as a pdf
ggsave("figures/clinical_enrichment.pdf", dpi=300, device="pdf", useDingbats=FALSE, width = 5, height = 6)


##PART 3: Building the Fused Patient Similarity Network
#Now we create another patient similarity network using the mutational profiles of the patients
#We download the mutation data from the TCGA-LUSC project
#First we create a query to obtain the mutation data
query <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"  # Workflow type for mutation analysis
)
#We download the data
#GDCdownload(query = query, directory = "GDCdata", method = "api")
#Now, we prepare the data for analysis
mutation.data <- GDCprepare(query = query, directory = "GDCdata")
#Here we create a mutation matrix
mutation_matrix <- mutation.data %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 12)) %>%
  distinct() %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(
    names_from = Hugo_Symbol,
    values_from = present,
    values_fill = 0
)
#We set the rownames to the Tumor_Sample_Barcode
mut.matrix <- as.matrix(mutation_matrix[, -1])
rownames(mut.matrix) <- mutation_matrix$Tumor_Sample_Barcode
#We filter out the patients that are not in the first layer
#First we get the cancer names and substitute the "." by "-"
cancer.names <- gsub("\\.", "-", V(g.c)$name)
#We filter the mutation matrix by the patients in the cancer network
mut.matrix <- mut.matrix[intersect(rownames(mut.matrix), cancer.names),]

#Now we normalize
mut.matrix <- apply(t(mut.matrix), 2, normalize)

#Now we obtain the (distance)similarity matrix using the Cosine similarity
sim.matrix.mut <- as.matrix(proxy::simil(mut.matrix, method = "cosine",  by_rows = FALSE))
#We ensure the diagional is 0 
diag(sim.matrix.mut) <- 0

#Now we alspo have to convert the cosine similarity matrix in a distance matrix
dis.matrix.c <- 1 - sim.matrix.c
dis.matrix.mut <- 1 - sim.matrix.mut
diag(dis.matrix.c) <- 0
diag(dis.matrix.mut) <- 0

#We rename the rows and columns of the cancer matrix to match the mutation matrix
colnames(dis.matrix.c) <-  gsub("\\.", "-", colnames(dis.matrix.c)) 
rownames(dis.matrix.c) <- gsub("\\.", "-", rownames(dis.matrix.c))
#We also filter this matrix by the patients that are in the mutation matrix to ensure they match
dis.matrix.c <- dis.matrix.c[row.names(dis.matrix.mut), colnames(dis.matrix.mut)]

#Now we can construct affinity matrices using these distances
W.c <- affinityMatrix(dis.matrix.c, K = 20, sigma = 0.5)
W.mut <- affinityMatrix(dis.matrix.mut, K = 20, sigma = 0.5)

#We can visualize for example the clusters found in both networks side by side
displayClustersWithHeatmap(W.c, spectralClustering(W.c, 3))
displayClustersWithHeatmap(W.mut, spectralClustering(W.mut, 3))

#Now we merge the two networks using the SNF algorithm
W.fused <- SNF(list(W.c, W.mut), K = 20, t = 20)
labels.fused = spectralClustering(W.fused, 3)
#We visualize the clusters found in the fused network and save it in a pdf
pdf("figures/fused_network_clusters.pdf", width = 5, height = 5.1)
displayClustersWithHeatmap(W.fused, labels.fused, tick.labelsize = 5)
dev.off()

#Here we fuse the rownames of the W.fused matrix to the cluster labels obtained
cluster.labels <- data.frame(cluster = labels.fused, row.names = rownames(W.fused))

##PART 4: Survival Analysis
#Now we perform survival analysis for the patients found in the clusters
#First, we download the clinical data
clinical.query <- GDCquery_clinic(project = "TCGA-LUSC", type = "clinical", save.csv=FALSE)
#Here we save it
#write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)

#First we take the subset of the clinical data that we need and then filter by the patients in the clusters
rownames(clinical.query) <- substr(clinical.query$submitter_id, 1, 12)
#We filter the clinical data by the patients in the clusters
clinical.query <- clinical.query[rownames(clinical.query) %in% rownames(W.fused),]
#We make a binary column for the death status
clinical.query$death <- ifelse(clinical.query$vital_status == "Alive", 0, 1)
#We create an overall survival time column that is:
#days_to_death if the patient is dead
#days_to_last_followup if the patient is alive
clinical.query$overall_survival <- ifelse(clinical.query$vital_status == "Alive",
                                         clinical.query$days_to_last_follow_up,
                                         clinical.query$days_to_death)
#We can also add the cluster labels to the clinical data
clinical.query$cluster <- cluster.labels[rownames(clinical.query),]

#Now we perform survival analysis
#First we create a survival object
surv.obj <- Surv(time = clinical.query$overall_survival, event = clinical.query$death)
#We create a survfit object
surv.fit <- survfit(surv.obj ~ cluster, data = clinical.query)
#We plot the survival curves
ggsurvplot(surv.fit, data = clinical.query,
           size = 1, # Change line size
           palette = c("#66c2a5", "#fc8d62", "#8da0cb"),
           censor.shape = "|", # Change censor shape
           censor.size = 4, # Change censor size
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"), # Change legend labels
           risk.table.height = 0.25, # Height of the risk table
           ggtheme = theme_bw()) # Change ggplot2 theme

#Here we perform a log-rank test
logrank.test <- survdiff(surv.obj ~ cluster, data = clinical.query)

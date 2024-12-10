library(BiocGenerics) 
library(maftools)
library(ggplot2)
library(igraph)
library(ggnet)
library(proxy)
library(dplyr)


##PART 1: Building the Patient Similarity Networks
#Here we build the patient similarity networks for the LUSC project based on the gene expression data
#First, we load the pre-processed, normalized data we obtained for the DEG analysis
setwd("~/Documents/Master/DigitalEpidemiology/DEPM-Project")
#We load the .RData file containing the cancer data
load("data/deg_expr_c.RData")
#We load the .RData file containing the normal data
load("data/deg_expr_n.RData")

#Now we obtain the similarity matrices for the cancer and normal samples using the cosine similarity
sim.matrix.c <- as.matrix(proxy::simil(deg.expr.c, method = "cosine",  by_rows = FALSE))
sim.matrix.n <- as.matrix(proxy::simil(deg.expr.n, method = "cosine",  by_rows = FALSE))
#We ensure the diagional is 0
diag(sim.matrix.c) <- 0
diag(sim.matrix.n) <- 0

#We construct the adjacency matrix by setting a threshold of 0.8 for the similarity
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
V(g.n)$degree <- igraph::degree(gi.n)

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

ggraph(g.n, layout = "fr") + 
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
ggsave("figures/patient_similarity_network_n.pdf", dpi=300, device="pdf", useDingbats=FALSE)

##PART 2: Building the 2-Layer Patient Similarity Network
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
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 9, 12)) %>%
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
mut.matrix <- mut.matrix[intersect(rownames(mut.matrix), V(g.c)$label),]

#Now we obtain the similarity matrix using the Cosine similarity
sim.matrix.mut <- as.matrix(proxy::simil(mut.matrix, method = "jaccard",  by_rows = TRUE))
#We ensure the diagional is 0
diag(sim.matrix.mut) <- 0

#We select as a threshold the top 5% of the similarities
threshold <- quantile(sim.matrix.mut, 0.95)
#We construct the adjacency matrix by setting a threshold of 0.8 for the similarity
adj.matrix.mut <- ifelse(sim.matrix.mut >= threshold, 1, 0)
#We build the patient similarity network
g.mut <- graph_from_adjacency_matrix(adj.matrix.mut, mode = "undirected")
#Now we find the communities in the network using Louvain's algorithm
communities.mut <- cluster_louvain(g.mut)
#Here we add to each vertex the community it belongs to
V(g.mut)$community <- communities.mut$membership
#We also add the degree of each vertex to the network
V(g.mut)$degree <- igraph::degree(g.mut)

#We plot the network
ggraph(g.mut, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = community), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the communities
  scale_color_viridis_c(option = "plasma")+
  #We add the labels
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = Inf)

#We save the plot in high resolution as a pdf and tight layout
ggsave("figures/patient_similarity_network_mut.pdf", dpi=300, device="pdf", useDingbats=FALSE)
















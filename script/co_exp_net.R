#Here we import the necessary libraries
library(psych)
library(igraph)
library(ggraph)
library(ggplot2)
library(latex2exp)
library(gridExtra)

##PART 1: Building the Co-Expression Networks -----
#We are going to build co-expression networks for the LUSC project. We will build two networks: one for the primary tumor samples and another for the solid tissue normal samples.
#We are going to use the Spearman correlation coefficient to measure the strength of the association between the genes.

#First, we load the pre-processed, normalized data we obtained for the DEG analysis
setwd("~/Documents/Master/DigitalEpidemiology/DEPM-Project")
#We load the .RData file containing the cancer data
load("data/deg_expr_c.RData")
#We load the .RData file containing the normal data
load("data/deg_expr_n.RData")

#First, we calculate the Spearman correlation adjusting the p-values for multiple testing
cor.test.c <- corr.test(t(deg.expr.c), use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)
cor.test.n <- corr.test(t(deg.expr.n), use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)

#To build the network, we will use the Spearman correlation to create an Adjacency Matrix
#First, we extract the correlation values
rho.c <- cor.test.c$r
rho.n <- cor.test.n$r
#We set the diagonal to zero since we are not interested in self-correlations
diag(rho.c) <- 0
diag(rho.n) <- 0
#Now we extract the p-values
pval.c <- cor.test.c$p
pval.n <- cor.test.n$p
#Since they are reported only for the upper triangle, we need to fill the lower triangle as well
pval.c[lower.tri(pval.c)] <- t(pval.c)[lower.tri(pval.c)]
pval.n[lower.tri(pval.n)] <- t(pval.n)[lower.tri(pval.n)]

#Now we build the adjacency matrix such that |rho| >= 0.7 and p-value <= 0.01
adj.mat.c <- ifelse((pval.c < 1e-3)&(abs(rho.c) >= 0.7), 1, 0)
adj.mat.n <- ifelse((pval.n < 1e-3) & (abs(rho.n) >= 0.7), 1, 0)

#Here we create a graph object from the adjacency matrix and obtain its components
g.c <- graph.adjacency(adj.mat.c, mode = "undirected")
g.n <- graph.adjacency(adj.mat.n, mode = "undirected")

#Here we obtain the connected components of the networks
components.c <- clusters(g.c, mode = "weak")
components.n <- clusters(g.n, mode = "weak")
#Here we obtain the biggest connected component
biggest.c <- which.max(components.c$csize)
biggest.n <- which.max(components.n$csize)
#We extract the vertices of the biggest connected component
v.idx.c <- V(g.c)[components.c$membership == biggest.c]
v.idx.n <- V(g.n)[components.n$membership == biggest.n]
#We extract the subgraphs
giant.c <- induced.subgraph(g.c, v.idx.c)
giant.n <- induced.subgraph(g.n, v.idx.n)
#We extract their degrees
V(giant.c)$degree <- igraph::degree(giant.c)
V(giant.n)$degree <- igraph::degree(giant.n)
#We add the closeness centrality values
V(giant.c)$closeness <- igraph::closeness(giant.c)
V(giant.n)$closeness <- igraph::closeness(giant.n)

#Now we plot the graphs
ggraph(giant.c, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = closeness), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the closeness centrality
  scale_color_viridis_c(option = "plasma")

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/co_expr_cancer.pdf", dpi=300, device="pdf", useDingbats=FALSE)

ggraph(giant.n, layout = "fr") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = closeness), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the closeness centrality
  scale_color_viridis_c(option = "plasma")

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/co_expr_normal.pdf", dpi=300, device="pdf", useDingbats=FALSE)

##PART 2: Degree Distribution ----
#Here we obtain the degree distribution of the networks, fit a power-law distribution

#First, we obtain the degree distribution for both networks (considering the giant component)
degree.dist.c <- igraph::degree(giant.c)
degree.dist.n <- igraph::degree(giant.n)

#Here we fit a power law using igraphs fit_power_law function
pl.c <- igraph::fit_power_law(degree.dist.c, p.value = TRUE)
pl.n <- igraph::fit_power_law(degree.dist.n, p.value = TRUE)
#Here we print the KS test p-value and test statistic
print(c("Cancer Network:", pl.c$KS.p, pl.c$KS.stat))
print(c("Normal Network:", pl.n$KS.p, pl.n$KS.stat))

#Now we make a log-log plot of the degree distributions side by side using ggplot2
#For the Cancer Network we obtain the frequency of the degrees
degree_counts.c <- table(degree.dist.c)
degree_values.c <- as.numeric(names(degree_counts.c))
degree_frequencies.c <- as.numeric(degree_counts.c)
#For the Normal Network we obtain the frequency of the degrees
degree_counts.n <- table(degree.dist.n)
degree_values.n <- as.numeric(names(degree_counts.n))
degree_frequencies.n <- as.numeric(degree_counts.n)

#Now we plot
degdist_plot.c <- ggplot(data = data.frame(degree = degree_values.c, frequency = degree_frequencies.c), aes(x = degree, y = frequency)) +
  geom_point(shape=1, alpha=0.9, color="#67a9cf") +
  #Here I add log ticks to the x and y axis
  annotation_logticks()+
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=c(1,300),ylim=c(1,200))+
  labs(x = TeX(r"($\log_{10}(k)$)"), y = TeX(r"($\log_{10}(F_k)$)")) +
  theme_bw()+
  #We add a "a)" text to the superior right corner
  annotate("text", x = 200, y = 200, label = "a)", size = 7, color = "black")

degdist_plot.n <- ggplot(data = data.frame(degree = degree_values.n, frequency = degree_frequencies.n), aes(x = degree, y = frequency)) +
  geom_point(shape=5, alpha=0.9, color="#ef8a62") +
  annotation_logticks()+
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=c(1,300),ylim=c(1,200))+
  labs(x = TeX(r"($\log_{10}(k)$)"), y="") +
  theme_bw()+
  #We add a "b)" text to the superior right corner
  annotate("text", x = 200, y = 200, label = "b)", size = 7, color = "black")

comb.plot<- grid.arrange(degdist_plot.c, degdist_plot.n, ncol = 2)
#Here we save the combined plot
ggsave("figures/degree_distribution_coexp.pdf", comb.plot, dpi=300,  width=7, height=5, device="pdf", useDingbats=FALSE)

##PART 3: Network Hubs ----
#Here we obtain the hubs given by the degree and closeness centrality
#Since only the cancer network is scale-free, we will only obtain the hubs for this network
#We obtain the hubs of the network i.e. the top 5% of the nodes with the highest degree
degree.cutoff.c <- quantile(V(giant.c)$degree, 0.95)
hubs.c <- V(giant.c)[V(giant.c)$degree >= degree.cutoff.c]
#We do the same for the closeness centrality
closeness.cutoff.c <- quantile(V(giant.c)$closeness, 0.95)
hubs.c.closeness <- V(giant.c)[V(giant.c)$closeness >= closeness.cutoff.c]

#Here we combine the results ranked by degree and closeness centrality and save them
#We create a dataframe with the hubs names and their degree and closeness centrality
#If one gene is only a hub in one of the measures, we assign a degree of 0 to the closeness centrality and vice versa
hubs.df.c <- data.frame(name = V(giant.c)$name, degree = V(giant.c)$degree, closeness = V(giant.c)$closeness)
#Now we filter out the genes that are not hubs
hubs.df.c <- hubs.df.c[hubs.df.c$name %in% hubs.c$name | hubs.df.c$name %in% hubs.c.closeness$name,]
#Now we add a is_degree column to indicate if the gene is a hub by degree
hubs.df.c$is_degree <- ifelse(hubs.df.c$name %in% hubs.c$name, 1, 0)
#Now we add a is_closeness column to indicate if the gene is a hub by closeness centrality
hubs.df.c$is_closeness <- ifelse(hubs.df.c$name %in% hubs.c.closeness$name, 1, 0)
#We save the hubs
write.table(hubs.df.c, file = "data/hubs_cancer_coexp.csv", sep = ";")

##PART 4: Differential Co-Expression Network ----
#To calculate the differential co-expression network we first have to obtain Fisher's Z transformation to compare the networks
#We calculate the Fisher's Z transformation for the correlation coefficients
z.c <- atanh(rho.c)
z.n <- atanh(rho.n)
#now we obtain the Z-score for the difference between the two networks
z_score <- function(z1, z2, n1, n2){
  z <- (z1 - z2) / sqrt(1/(n1-3) + 1/(n2-3))
  return(z)
}
#Therefore we obtain
z_diff <- z_score(z.n, z.c, ncol(deg.expr.n), ncol(deg.expr.c))

#Now we build the new adjacency matrix for |Z| >= 4.5
adj.mat.diff <- ifelse(abs(z_diff) >= 4.5, 1, 0)
#We also construct the signed networks i.e. Z > 3 or Z < -3
adj.mat.diff.positive <- ifelse(z_diff > 4.5, 1, 0)
adj.mat.diff.negative <- ifelse(z_diff < -4.5, 1, 0)
#We ensure there are no self-connections
diag(adj.mat.diff) <- 0
diag(adj.mat.diff.positive) <- 0
diag(adj.mat.diff.negative) <- 0

#Here we create a graph object from the adjacency matrix and obtain its components
g.diff <- graph.adjacency(adj.mat.diff, mode = "undirected")
g.diff.positive <- graph.adjacency(adj.mat.diff.positive, mode = "undirected")
g.diff.negative <- graph.adjacency(adj.mat.diff.negative, mode = "undirected")

#Here we obtain the connected components of the networks
components.diff <- clusters(g.diff, mode = "weak")
components.diff.positive <- clusters(g.diff.positive, mode = "weak")
components.diff.negative <- clusters(g.diff.negative, mode = "weak")

#Here we obtain the biggest connected component
biggest.diff <- which.max(components.diff$csize)
biggest.diff.positive <- which.max(components.diff.positive$csize)
biggest.diff.negative <- which.max(components.diff.negative$csize)

#We extract the vertices of the biggest connected component
v.idx.diff <- V(g.diff)[components.diff$membership == biggest.diff]
v.idx.diff.positive <- V(g.diff.positive)[components.diff.positive$membership == biggest.diff.positive]
v.idx.diff.negative <- V(g.diff.negative)[components.diff.negative$membership == biggest.diff.negative]


#We extract the subgraphs
giant.diff <- induced.subgraph(g.diff, v.idx.diff)
giant.diff.positive <- induced.subgraph(g.diff.positive, v.idx.diff.positive)
giant.diff.negative <- induced.subgraph(g.diff.negative, v.idx.diff.negative)

#Here we plot the degree distribution of the giant components of the differential networks
degree.dist.diff <- igraph::degree(giant.diff)
degree.dist.diff.positive <- igraph::degree(giant.diff.positive)
degree.dist.diff.negative <- igraph::degree(giant.diff.negative)

#We fit a power-law distribution using igraphs fit_power_law function
pl.diff <- igraph::fit_power_law(degree.dist.diff, p.value = TRUE)
pl.diff.positive <- igraph::fit_power_law(degree.dist.diff.positive, p.value = TRUE)
pl.diff.negative <- igraph::fit_power_law(degree.dist.diff.negative, p.value = TRUE)

#Now we make a log-log plot of the degree distributions side by side using ggplot2
#For the complete Differential Network we obtain the frequency of the degrees
degree_counts.diff <- table(degree.dist.diff)
degree_values.diff <- as.numeric(names(degree_counts.diff))
degree_frequencies.diff <- as.numeric(degree_counts.diff)
#For the negative Differential Network we obtain the frequency of the degrees
degree_counts.diff.negative <- table(degree.dist.diff.negative)
degree_values.diff.negative <- as.numeric(names(degree_counts.diff.negative))
degree_frequencies.diff.negative <- as.numeric(degree_counts.diff.negative)
#For the positive Differential Network we obtain the frequency of the degrees
degree_counts.diff.positive <- table(degree.dist.diff.positive)
degree_values.diff.positive <- as.numeric(names(degree_counts.diff.positive))
degree_frequencies.diff.positive <- as.numeric(degree_counts.diff.positive)

#Now we plot first the individual degree distribution of the complete differential network
ggplot(data = data.frame(degree = degree_values.diff, frequency = degree_frequencies.diff), aes(x = degree, y = frequency)) +
  geom_point(shape=1, alpha=0.9, color="#ef8a62") +
  #Here I add log ticks to the x and y axis
  annotation_logticks()+
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=c(1,300),ylim=c(1,110))+
  labs(x = TeX(r"($\log_{10}(k)$)"), y = TeX(r"($\log_{10}(F_k)$)")) +
  theme_bw()

#Here we save the plot
ggsave("figures/diff_degree_distribution.pdf", dpi=300, device="pdf",width=4, height=5, useDingbats=FALSE)

#Now we plot side by side the degree distributions of the positive and negative differential networks

degdist_plot.diff.positive <- ggplot(data = data.frame(degree = degree_values.diff.positive, frequency = degree_frequencies.diff.positive), aes(x = degree, y = frequency)) +
  geom_point(shape=1, alpha=0.9, color="#67a9cf") +
  #Here I add log ticks to the x and y axis
  annotation_logticks()+
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=c(1,300),ylim=c(1,110))+
  labs(x = TeX(r"($\log_{10}(k)$)"), y="") +
  theme_bw()+
  #We add a "a)" text to the superior right corner
  annotate("text", x = 180, y = 100, label = "a)", size = 7, color = "black")

degdist_plot.diff.negative <- ggplot(data = data.frame(degree = degree_values.diff.negative, frequency = degree_frequencies.diff.negative), aes(x = degree, y = frequency)) +
  geom_point(shape=1, alpha=0.9, color="#ef8a62") +
  #Here I add log ticks to the x and y axis
  annotation_logticks()+
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=c(1,300),ylim=c(1,110))+
  labs(x = TeX(r"($\log_{10}(k)$)"), y="") +
  theme_bw()+
  #We add a "b)" text to the superior right corner
  annotate("text", x = 180, y = 100, label = "b)", size = 7, color = "black")

comb.plot.diff <- grid.arrange(degdist_plot.diff.positive, degdist_plot.diff.negative, ncol = 2)
#Here we save the combined plot
ggsave("figures/diff_degree_distribution_signed.pdf", comb.plot.diff, dpi=300,  width=7, height=5, device="pdf", useDingbats=FALSE)

#Now we plot the graphs
#We extract their degrees
V(giant.diff)$degree <- igraph::degree(giant.diff)
V(giant.diff.positive)$degree <- igraph::degree(giant.diff.positive)
V(giant.diff.negative)$degree <- igraph::degree(giant.diff.negative)

#We add the closeness centrality values
V(giant.diff)$closeness <- igraph::closeness(giant.diff)
V(giant.diff.positive)$closeness <- igraph::closeness(giant.diff.positive)
V(giant.diff.negative)$closeness <- igraph::closeness(giant.diff.negative)

#Now we plot the graphs
ggraph(giant.diff, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = closeness), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the closeness centrality
  scale_color_viridis_c(option = "plasma")

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/diff_co_expr.pdf", dpi=300, device="pdf", useDingbats=FALSE)

ggraph(giant.diff.positive, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = closeness), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the closeness centrality
  scale_color_viridis_c(option = "plasma")

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/diff_co_expr_positive.pdf", dpi=300, device="pdf", useDingbats=FALSE)

ggraph(giant.diff.negative, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = closeness), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the closeness centrality
  scale_color_viridis_c(option = "plasma")

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/diff_co_expr_negative.pdf", dpi=300, device="pdf", useDingbats=FALSE)

#We also plot the subnetwork induced by the hubs in the differential network
#We obtain the hubs of the network i.e. the top 5% of the nodes with the highest degree
#Since only the signed networks are scale-free, we will only obtain the hubs for the signed network
#We do it also for the complete network although we doubt it is scale-free
degree.cutoff.diff <- quantile(V(giant.diff)$degree, 0.95)
degree.cutoff.diff.positive <- quantile(V(giant.diff.positive)$degree, 0.95)
degree.cutoff.diff.negative <- quantile(V(giant.diff.negative)$degree, 0.95)
hubs.diff <- V(giant.diff)[V(giant.diff)$degree >= degree.cutoff.diff]
hubs.diff.positive <- V(giant.diff.positive)[V(giant.diff.positive)$degree >= degree.cutoff.diff.positive]
hubs.diff.negative <- V(giant.diff.negative)[V(giant.diff.negative)$degree >= degree.cutoff.diff.negative]

#We take the indices from the hubs
v.idx.diff.hubs <- V(g.diff)[V(g.diff)$name %in% hubs.diff$name]
v.idx.diff.positive.hubs <- V(g.diff.positive)[V(g.diff.positive)$name %in% hubs.diff.positive$name]
v.idx.diff.negative.hubs <- V(g.diff.negative)[V(g.diff.negative)$name %in% hubs.diff.negative$name]
#We extract the subgraphs
giant.diff.hubs <- induced.subgraph(g.diff, v.idx.diff.hubs)
giant.diff.positive.hubs <- induced.subgraph(g.diff.positive, v.idx.diff.positive.hubs)
giant.diff.negative.hubs <- induced.subgraph(g.diff.negative, v.idx.diff.negative.hubs)

#We extract their degrees and closeness centrality and label them
V(giant.diff.hubs)$degree <- igraph::degree(giant.diff.hubs)
V(giant.diff.hubs)$closeness <- igraph::closeness(giant.diff.hubs)
V(giant.diff.hubs)$label <- V(giant.diff.hubs)$name
V(giant.diff.positive.hubs)$degree <- igraph::degree(giant.diff.positive.hubs)
V(giant.diff.positive.hubs)$closeness <- igraph::closeness(giant.diff.positive.hubs)
V(giant.diff.positive.hubs)$label <- V(giant.diff.positive.hubs)$name
V(giant.diff.negative.hubs)$degree <- igraph::degree(giant.diff.negative.hubs)
V(giant.diff.negative.hubs)$closeness <- igraph::closeness(giant.diff.negative.hubs)
V(giant.diff.negative.hubs)$label <- V(giant.diff.negative.hubs)$name

#We plot the hub subnetworks
ggraph(giant.diff.hubs, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = closeness), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the closeness centrality
  scale_color_viridis_c(option = "plasma")+
  #We add the labels
  geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = Inf)

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/diff_co_expr_hubs.pdf", dpi=300, device="pdf", useDingbats=FALSE)

ggraph(giant.diff.positive.hubs, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = closeness), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the closeness centrality
  scale_color_viridis_c(option = "plasma")+
#We add the labels
geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = Inf)

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/diff_co_expr_positive_hubs.pdf", dpi=300, device="pdf", useDingbats=FALSE)

ggraph(giant.diff.negative.hubs, layout = "fr") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
  geom_edge_link(alpha = 0.2) +
  #The size of the nodes is proportional to the degree.
  geom_node_point(aes(size = degree, color = closeness), alpha=0.8, show.legend = FALSE)+
  #We add a sequential color scale for the closeness centrality
  scale_color_viridis_c(option = "plasma")+
  #We add the labels
  geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = Inf)

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/diff_co_expr_negative_hubs.pdf", dpi=300, device="pdf", useDingbats=FALSE)


##PART 5: Hubs in the Differential Network ----
#We obtain the hubs of the differential network for degree and closeness centrality
#Just for the signed networks since they are scale-free (and the complete network for completeness)

#We obtain the hubs of the network i.e. the top 5% of the nodes with the highest degree
degree.cutoff.diff <- quantile(V(giant.diff)$degree, 0.95)
degree.cutoff.diff.positive <- quantile(V(giant.diff.positive)$degree, 0.95)
degree.cutoff.diff.negative <- quantile(V(giant.diff.negative)$degree, 0.95)
hubs.diff <- V(giant.diff)[V(giant.diff)$degree >= degree.cutoff.diff]
hubs.diff.positive <- V(giant.diff.positive)[V(giant.diff.positive)$degree >= degree.cutoff.diff.positive]
hubs.diff.negative <- V(giant.diff.negative)[V(giant.diff.negative)$degree >= degree.cutoff.diff.negative]
#We do the same for the closeness centrality
closeness.cutoff.diff <- quantile(V(giant.diff)$closeness, 0.95)
closeness.cutoff.diff.positive <- quantile(V(giant.diff.positive)$closeness, 0.95)
closeness.cutoff.diff.negative <- quantile(V(giant.diff.negative)$closeness, 0.95)
hubs.diff.closeness <- V(giant.diff)[V(giant.diff)$closeness >= closeness.cutoff.diff]
hubs.diff.positive.closeness <- V(giant.diff.positive)[V(giant.diff.positive)$closeness >= closeness.cutoff.diff.positive]
hubs.diff.negative.closeness <- V(giant.diff.negative)[V(giant.diff.negative)$closeness >= closeness.cutoff.diff.negative]

#Here we combine the results ranked by degree and closeness centrality and save them
#We create a dataframe with the hubs names and their degree and closeness centrality
#If one gene is only a hub in one of the measures, we assign a degree of 0 to the closeness centrality and vice versa
hubs.df.diff.positive <- data.frame(name = V(giant.diff.positive)$name, degree = V(giant.diff.positive)$degree, closeness = V(giant.diff.positive)$closeness)
#We filter more the genes that are not hubs
hubs.df.diff.positive <- hubs.df.diff.positive[hubs.df.diff.positive$name %in% hubs.diff.positive$name | hubs.df.diff.positive$name %in% hubs.diff.positive.closeness$name,]
#Now we add a is_degree column to indicate if the gene is a hub by degree
hubs.df.diff.positive$is_degree <- ifelse(hubs.df.diff.positive$name %in% hubs.diff.positive$name, 1, 0)
#Now we add a is_closeness column to indicate if the gene is a hub by closeness centrality
hubs.df.diff.positive$is_closeness <- ifelse(hubs.df.diff.positive$name %in% hubs.diff.positive.closeness$name, 1, 0)
#We save the hubs
write.table(hubs.df.diff.positive, file = "data/hubs_diff_positive_coexp.csv", sep = ";")

#We do the same for the negative differential network
hubs.df.diff.negative <- data.frame(name = V(giant.diff.negative)$name, degree = V(giant.diff.negative)$degree, closeness = V(giant.diff.negative)$closeness)
#We filter more the genes that are not hubs
hubs.df.diff.negative <- hubs.df.diff.negative[hubs.df.diff.negative$name %in% hubs.diff.negative$name | hubs.df.diff.negative$name %in% hubs.diff.negative.closeness$name,]
#Now we add a is_degree column to indicate if the gene is a hub by degree
hubs.df.diff.negative$is_degree <- ifelse(hubs.df.diff.negative$name %in% hubs.diff.negative$name, 1, 0)
#Now we add a is_closeness column to indicate if the gene is a hub by closeness centrality
hubs.df.diff.negative$is_closeness <- ifelse(hubs.df.diff.negative$name %in% hubs.diff.negative.closeness$name, 1, 0)

#Finally for the complete network
hubs.df.diff <- data.frame(name = V(giant.diff)$name, degree = V(giant.diff)$degree, closeness = V(giant.diff)$closeness)
#We filter more the genes that are not hubs
hubs.df.diff <- hubs.df.diff[hubs.df.diff$name %in% hubs.diff$name | hubs.df.diff$name %in% hubs.diff.closeness$name,]
#Now we add a is_degree column to indicate if the gene is a hub by degree
hubs.df.diff$is_degree <- ifelse(hubs.df.diff$name %in% hubs.diff$name, 1, 0)
#Now we add a is_closeness column to indicate if the gene is a hub by closeness centrality
hubs.df.diff$is_closeness <- ifelse(hubs.df.diff$name %in% hubs.diff.closeness$name, 1, 0)

#We save the hubs
write.table(hubs.df.diff.negative, file = "data/hubs_diff_negative_coexp.csv", sep = ";")
write.table(hubs.df.diff.negative, file = "data/hubs_diff_negative_coexp.csv", sep = ";")
write.table(hubs.df.diff, file = "data/hubs_diff_coexp.csv", sep = ";")
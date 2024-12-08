#Here we import the necessary libraries
library(psych)
library(igraph)
library(ggraph)
library(ggplot2)
library(poweRlaw)
library(latex2exp)

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
cor.test.c <- corr.test(t(deg.expr.c), use = "pairwise", method = "spearman", adjust = "fdr", ci = FALSE)
cor.test.n <- corr.test(t(deg.expr.n), use = "pairwise", method = "spearman", adjust = "fdr", ci = FALSE)

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

#Now we build the adjacency matrix such that |rho| >= 0.6 and p-value <= 0.01
adj.mat.c <- ifelse(abs(rho.c) >= 0.6 & pval.c <= 0.01, 1, 0)
adj.mat.n <- ifelse(abs(rho.n) >= 0.6 & pval.n <= 0.01, 1, 0)

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

#We obtain the hubs of the network i.e. the top 5% of the nodes with the highest degree
degree.cutoff.c <- quantile(V(giant.c)$degree, 0.95)
degree.cutoff.n <- quantile(V(giant.n)$degree, 0.95)
hubs.c <- V(giant.c)[V(giant.c)$degree >= degree.cutoff]
hubs.n <- V(giant.n)[V(giant.n)$degree >= degree.cutoff]

#We just show the labels of the hubs
V(giant.c)$label <- ifelse(V(giant.c)$name %in% hubs.c$name, V(giant.c)$name, "")
V(giant.n)$label <- ifelse(V(giant.n)$name %in% hubs.n$name, V(giant.n)$name, "")

#Now we plot the graphs
ggraph(giant.c, layout = "fr") + 
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
ggsave("figures/co_expr_cancer.pdf", dpi=300, device="pdf", useDingbats=FALSE)

ggraph(giant.n, layout = "fr") +
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
ggsave("figures/co_expr_normal.pdf", dpi=300, device="pdf", useDingbats=FALSE)

##PART 2: Degree Distribution ----
#Here we obtain the degree distribution of the networks, fit a power-law distribution

#First, we obtain the degree distribution for both networks
degree.dist.c <- igraph::degree(g.c)
degree.dist.n <- igraph::degree(g.n)
#We delete the nodes with degree 0 since we cannot fit a power-law distribution to them
degree.dist.c <- degree.dist.c[degree.dist.c != 0]
degree.dist.n <- degree.dist.n[degree.dist.n != 0]
#We fit a power-law distribution
pl_model.c <- displ$new(degree.dist.c)
pl_model.n <- displ$new(degree.dist.n)
pl_model.c$setXmin(estimate_xmin(pl_model.c))
pl_model.n$setXmin(estimate_xmin(pl_model.n))
#We also fit a log-normal and poisson distribution to see which one fits better
ln_model.c <- dislnorm$new(degree.dist.c)
ln_model.n <- dislnorm$new(degree.dist.n)
ln_model.c$setXmin(estimate_xmin(ln_model.c))
ln_model.n$setXmin(estimate_xmin(ln_model.n))
#Here we fit a poisson distribution
po_model.c <- dispois$new(degree.dist.c)
po_model.n <- dispois$new(degree.dist.n)
po_model.c$setXmin(estimate_xmin(po_model.c))
po_model.n$setXmin(estimate_xmin(po_model.n))

#First we plot the CDF of the cancer degree distribution
pdf("figures/cdf_degree_distribution_coexp.pdf", width = 10, height = 7)
#Here we plot by mimicking the ggplot theme_bw function
par(bg = "white",          # White background
    col.axis = "black",    # Black axis text
    col.lab = "black",     # Black axis labels
    col.main = "black",    # Black title text
    col.sub = "black",     # Black subtitle text
    fg = "black",          # Black plot border
    las = 1,               # Horizontal axis labels
    family = "sans",
    mgp = c(3, 1, 0),
    mfrow = c(1,2)
)
legend_expr.c.pl <-paste("Power Law, $\\alpha$=", round(pl_model.c$pars,2),", $x_{min}$=", round(pl_model.c$xmin,2))
legend_expr.c.ln <-paste("Log-Normal, $\\mu$=", round(ln_model.c$pars[1],2),", $\\sigma$=", round(ln_model.c$pars[2],2),", $x_{min}$=", round(ln_model.c$xmin,2))
legend_expr.c.po <-paste("Poisson, $\\lambda$=", round(po_model.c$pars,2),", $x_{min}$=", round(po_model.c$xmin,2))
plot(pl_model.c, xlab = TeX(r"($k$)"), ylab= TeX(r"(CDF ($k$))"),col = "#0571b0")
lines(pl_model.c, col="#ca0020")
lines(ln_model.c, col="#f4a582")
lines(po_model.c, col="#92c5de")
#Here we add the "a)" text at the bottom right corner
mtext("(a)", side = 3, line =-2, adj = 0.95, cex = 1.5)
legend("bottomleft", legend = c(TeX(legend_expr.c.ln),TeX(legend_expr.c.pl), TeX(legend_expr.c.po)), col = c("#f4a582","#ca0020",  "#92c5de"), lty = 1, cex = 0.65)
grid()

#Now we plot the CDF of the normal degree distribution
legend_expr.n.pl <-paste("Power Law, $\\alpha$=", round(pl_model.n$pars,2),", $x_{min}$=", round(pl_model.n$xmin,2))
legend_expr.n.ln <-paste("Log-Normal, $\\mu$=", round(ln_model.n$pars[1],2),", $\\sigma$=", round(ln_model.n$pars[2],2),", $x_{min}$=", round(ln_model.n$xmin,2))
legend_expr.n.po <-paste("Poisson, $\\lambda$=", round(po_model.n$pars,2),", $x_{min}$=", round(po_model.n$xmin,2))
plot(pl_model.n, xlab = TeX(r"($k$)"), ylab= "", col = "#0571b0")
lines(pl_model.n, col="#ca0020")
lines(ln_model.n, col="#f4a582")
lines(po_model.n, col="#92c5de")
#Here we add the "b)" text at the bottom right corner
mtext("(b)", side = 3, line =-2, adj = 0.95, cex = 1.5)
legend("bottomleft", legend = c(TeX(legend_expr.n.ln),TeX(legend_expr.n.pl), TeX(legend_expr.n.po)), col = c("#f4a582","#ca0020",  "#92c5de"), lty = 1, cex = 0.65)
grid()

#We save the plots in high resolution as a pdf and tight layout
dev.off()

#Now to make sure we are have indeed a power law, we test the power law hypothesis
#We perform a goodness-of-fit test via bootstrapping using poweRlaw package
bs_p.c <- bootstrap_p(pl_model.c, no_of_sims = 1000, threads = 4) #this gives p=0.415, this means that we cannot reject the null hypothesis that the data follows a power-law distribution
bs_p.n <- bootstrap_p(pl_model.n, no_of_sims = 1000, threads = 4) #this gives p=0.011, this means that we can reject the null hypothesis that the data follows a power-law distribution

##PART 3: Network Hubs ----
#Here we obtain the hubs given by the degree and closeness centrality
#We obtain the hubs of the network i.e. the top 5% of the nodes with the highest degree
degree.cutoff.c <- quantile(V(giant.c)$degree, 0.95)
degree.cutoff.n <- quantile(V(giant.n)$degree, 0.95)
hubs.c <- V(giant.c)[V(giant.c)$degree >= degree.cutoff.c]
hubs.n <- V(giant.n)[V(giant.n)$degree >= degree.cutoff.n]
#We do the same for the closeness centrality
closeness.cutoff.c <- quantile(V(giant.c)$closeness, 0.95)
closeness.cutoff.n <- quantile(V(giant.n)$closeness, 0.95)
hubs.c.closeness <- V(giant.c)[V(giant.c)$closeness >= closeness.cutoff.c]
hubs.n.closeness <- V(giant.n)[V(giant.n)$closeness >= closeness.cutoff.n]

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

#We do the same for the normal network although these are not hubs since the degree distribution is not a power-law
hubs.df.n <- data.frame(name = V(giant.n)$name, degree = V(giant.n)$degree, closeness = V(giant.n)$closeness)
#We filter out the genes that are not hubs
hubs.df.n <- hubs.df.n[hubs.df.n$name %in% hubs.n$name | hubs.df.n$name %in% hubs.n.closeness$name,]
#Now we add a is_degree column to indicate if the gene is a hub by degree
hubs.df.n$is_degree <- ifelse(hubs.df.n$name %in% hubs.n$name, 1, 0)
#Now we add a is_closeness column to indicate if the gene is a hub by closeness centrality
hubs.df.n$is_closeness <- ifelse(hubs.df.n$name %in% hubs.n.closeness$name, 1, 0)
#We save the hubs
write.table(hubs.df.n, file = "data/hubs_normal_coexp.csv", sep = ";")

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
#Now we build the new adjacency matrix for |Z| >= 3
adj.mat.diff <- ifelse(abs(z_diff) >= 3, 1, 0)
#We also construct the signed networks i.e. Z > 3 or Z < -3
adj.mat.diff.positive <- ifelse(z_diff > 3, 1, 0)
adj.mat.diff.negative <- ifelse(z_diff < -3, 1, 0)
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

#We extract their degrees
V(giant.diff)$degree <- igraph::degree(giant.diff)
V(giant.diff.positive)$degree <- igraph::degree(giant.diff.positive)
V(giant.diff.negative)$degree <- igraph::degree(giant.diff.negative)

#We add the closeness centrality values
V(giant.diff)$closeness <- igraph::closeness(giant.diff)
V(giant.diff.positive)$closeness <- igraph::closeness(giant.diff.positive)
V(giant.diff.negative)$closeness <- igraph::closeness(giant.diff.negative)

#We obtain the hubs of the network i.e. the top 5% of the nodes with the highest degree
degree.cutoff.diff <- quantile(V(giant.diff)$degree, 0.95)
degree.cutoff.diff.positive <- quantile(V(giant.diff.positive)$degree, 0.95)
degree.cutoff.diff.negative <- quantile(V(giant.diff.negative)$degree, 0.95)
hubs.diff <- V(giant.diff)[V(giant.diff)$degree >= degree.cutoff.diff]
hubs.diff.positive <- V(giant.diff.positive)[V(giant.diff.positive)$degree >= degree.cutoff.diff.positive]
hubs.diff.negative <- V(giant.diff.negative)[V(giant.diff.negative)$degree >= degree.cutoff.diff.negative]

#We just show the labels of the hubs
V(giant.diff)$label <- ifelse(V(giant.diff)$name %in% hubs.diff$name, V(giant.diff)$name, "")
V(giant.diff.positive)$label <- ifelse(V(giant.diff.positive)$name %in% hubs.diff.positive$name, V(giant.diff.positive)$name, "")
V(giant.diff.negative)$label <- ifelse(V(giant.diff.negative)$name %in% hubs.diff.negative$name, V(giant.diff.negative)$name, "")

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
  scale_color_viridis_c(option = "plasma")+
  #We add the labels
  geom_node_text(aes(label = label), repel = TRUE, size = 3, max.overlaps = Inf)

#We save the plots in high resolution as a pdf and tight layout
ggsave("figures/diff_co_expr_negative.pdf", dpi=300, device="pdf", useDingbats=FALSE)

#We also plot the subnetwork induced by the hubs in the differential network
#First we take the indices from the hubs
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

#We plot the subnetworks
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


##PART 5: Degree Distribution (Diff Network) ----
#Here we obtain the degree distribution of the differential network, fit a power-law distribution

#First, we obtain the degree distribution for the differential network
degree.dist.diff <- igraph::degree(g.diff)
degree.dist.diff.positive <- igraph::degree(g.diff.positive)
degree.dist.diff.negative <- igraph::degree(g.diff.negative)

#We delete the nodes with degree 0 since we cannot fit a power-law distribution to them
degree.dist.diff <- degree.dist.diff[degree.dist.diff != 0]
degree.dist.diff.positive <- degree.dist.diff.positive[degree.dist.diff.positive != 0]
degree.dist.diff.negative <- degree.dist.diff.negative[degree.dist.diff.negative != 0]

#We fit a power-law distribution
pl_model.diff <- displ$new(degree.dist.diff)
pl_model.diff.positive <- displ$new(degree.dist.diff.positive)
pl_model.diff.negative <- displ$new(degree.dist.diff.negative)

pl_model.diff$setXmin(estimate_xmin(pl_model.diff))
pl_model.diff.positive$setXmin(estimate_xmin(pl_model.diff.positive))
pl_model.diff.negative$setXmin(estimate_xmin(pl_model.diff.negative))

#We also fit a log-normal and poisson distribution to see which one fits better
ln_model.diff <- dislnorm$new(degree.dist.diff)
ln_model.diff.positive <- dislnorm$new(degree.dist.diff.positive)
ln_model.diff.negative <- dislnorm$new(degree.dist.diff.negative)

ln_model.diff$setXmin(estimate_xmin(ln_model.diff))
ln_model.diff.positive$setXmin(estimate_xmin(ln_model.diff.positive))
ln_model.diff.negative$setXmin(estimate_xmin(ln_model.diff.negative))

#Here we fit a poisson distribution
po_model.diff <- dispois$new(degree.dist.diff)
po_model.diff.positive <- dispois$new(degree.dist.diff.positive)
po_model.diff.negative <- dispois$new(degree.dist.diff.negative)

po_model.diff$setXmin(estimate_xmin(po_model.diff))
po_model.diff.positive$setXmin(estimate_xmin(po_model.diff.positive))
po_model.diff.negative$setXmin(estimate_xmin(po_model.diff.negative))

#First we plot the CDF of the differential degree distribution
pdf("figures/cdf_degree_distribution_diff_coexp.pdf", width = 5, height = 6)
#Here we plot
par(bg = "white",          # White background
    col.axis = "black",    # Black axis text
    col.lab = "black",     # Black axis labels
    col.main = "black",    # Black title text
    col.sub = "black",     # Black subtitle text
    fg = "black",          # Black plot border
    las = 1,               # Horizontal axis labels
    family = "sans",
    mgp = c(3, 1, 0),
    mfrow = c(1,1)
)
legend_expr.diff.pl <-paste("Power Law, $\\alpha$=", round(pl_model.diff$pars,2),", $x_{min}$=", round(pl_model.diff$xmin,2))
legend_expr.diff.ln <-paste("Log-Normal, $\\mu$=", round(ln_model.diff$pars[1],2),", $\\sigma$=", round(ln_model.diff$pars[2],2),", $x_{min}$=", round(ln_model.diff$xmin,2))
legend_expr.diff.po <-paste("Poisson, $\\lambda$=", round(po_model.diff$pars,2),", $x_{min}$=", round(po_model.diff$xmin,2))
plot(pl_model.diff, xlab = TeX(r"($k$)"), ylab= TeX(r"(CDF ($k$))"),col = "#0571b0")
lines(pl_model.diff, col="#ca0020")
lines(ln_model.diff, col="#f4a582")
lines(po_model.diff, col="#92c5de")
legend("bottomleft", legend = c(TeX(legend_expr.diff.ln),TeX(legend_expr.diff.pl), TeX(legend_expr.diff.po)), col = c("#f4a582","#ca0020",  "#92c5de"), lty = 1, cex = 0.65)
grid()

#We save the plots in high resolution as a pdf and tight layout
dev.off()

#Now we plot the CDF of the signed networks side by side
pdf("figures/cdf_degree_distribution_diff_coexp_signed.pdf", width = 10, height = 7)
#Here we plot
par(bg = "white",          # White background
    col.axis = "black",    # Black axis text
    col.lab = "black",     # Black axis labels
    col.main = "black",    # Black title text
    col.sub = "black",     # Black subtitle text
    fg = "black",          # Black plot border
    las = 1,               # Horizontal axis labels
    family = "sans",
    mgp = c(3, 1, 0),
    mfrow = c(1,2)
)
legend_expr.diff.pl.positive <-paste("Power Law, $\\alpha$=", round(pl_model.diff.positive$pars,2),", $x_{min}$=", round(pl_model.diff.positive$xmin,2))
legend_expr.diff.ln.positive <-paste("Log-Normal, $\\mu$=", round(ln_model.diff.positive$pars[1],2),", $\\sigma$=", round(ln_model.diff.positive$pars[2],2),", $x_{min}$=", round(ln_model.diff.positive$xmin,2))
legend_expr.diff.po.positive <-paste("Poisson, $\\lambda$=", round(po_model.diff.positive$pars,2),", $x_{min}$=", round(po_model.diff.negative$xmin,2))
plot(pl_model.diff.positive, xlab = TeX(r"($k$)"), ylab= TeX(r"(CDF ($k$))"),col = "#0571b0")
lines(pl_model.diff.positive, col="#ca0020")
lines(ln_model.diff.positive, col="#f4a582")
lines(po_model.diff.positive, col="#92c5de")
#Here we add the "a)" text at the bottom right corner
mtext("(a)", side = 3, line =-2, adj = 0.95, cex = 1.5)
legend("bottomleft", legend = c(TeX(legend_expr.diff.ln.positive),TeX(legend_expr.diff.pl.positive), TeX(legend_expr.diff.po.positive)), col = c("#f4a582","#ca0020",  "#92c5de"), lty = 1, cex = 0.65)
grid()

legend_expr.diff.pl.negative <-paste("Power Law, $\\alpha$=", round(pl_model.diff.negative$pars,2),", $x_{min}$=", round(pl_model.diff.negative$xmin,2))
legend_expr.diff.ln.negative <-paste("Log-Normal, $\\mu$=", round(ln_model.diff.negative$pars[1],2),", $\\sigma$=", round(ln_model.diff.negative$pars[2],2),", $x_{min}$=", round(ln_model.diff.negative$xmin,2))
legend_expr.diff.po.negative <-paste("Poisson, $\\lambda$=", round(po_model.diff.negative$pars,2),", $x_{min}$=", round(po_model.diff.negative$xmin,2))
plot(pl_model.diff.negative, xlab = TeX(r"($k$)"), ylab= TeX(r"(CDF ($k$))"),col = "#0571b0")
lines(pl_model.diff.negative, col="#ca0020")
lines(ln_model.diff.negative, col="#f4a582")
lines(po_model.diff.negative, col="#92c5de")
#Here we add the "b)" text at the bottom right corner
mtext("(b)", side = 3, line =-2, adj = 0.95, cex = 1.5)
legend("bottomleft", legend = c(TeX(legend_expr.diff.ln.negative),TeX(legend_expr.diff.pl.negative), TeX(legend_expr.diff.po.negative)), col = c("#f4a582","#ca0020",  "#92c5de"), lty = 1, cex = 0.65)
grid()

#We save the plots in high resolution as a pdf and tight layout
dev.off()

#Now we test the power law hypothesis
#We perform a goodness-of-fit test via bootstrapping using poweRlaw package
bs_p.diff <- bootstrap_p(pl_model.diff, no_of_sims = 1000, threads = 4) 
bs_p.diff.positive <- bootstrap_p(pl_model.diff.positive, no_of_sims = 1000, threads = 4) 
bs_p.diff.negative <- bootstrap_p(pl_model.diff.negative, no_of_sims = 1000, threads = 4) 

#All give p>1 meaning we cannot reject the null hypothesis that the data follows a power-law distribution
#We get p=0.256, p=0.2, and p=0,715 respectively

#Now we can calculate the hubs of the differential network for degree and closeness centrality
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
hubs.df.diff <- data.frame(name = V(giant.diff)$name, degree = V(giant.diff)$degree, closeness = V(giant.diff)$closeness)
#Now we filter out the genes that are not hubs
hubs.df.diff <- hubs.df.diff[hubs.df.diff$name %in% hubs.diff$name | hubs.df.diff$name %in% hubs.diff.closeness$name,]
#Now we add a is_degree column to indicate if the gene is a hub by degree
hubs.df.diff$is_degree <- ifelse(hubs.df.diff$name %in% hubs.diff$name, 1, 0)
#Now we add a is_closeness column to indicate if the gene is a hub by closeness centrality
hubs.df.diff$is_closeness <- ifelse(hubs.df.diff$name %in% hubs.diff.closeness$name, 1, 0)
#We save the hubs
write.table(hubs.df.diff, file = "data/hubs_diff_coexp.csv", sep = ";")

#We do the same for the positive differential network
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
#We save the hubs
write.table(hubs.df.diff.negative, file = "data/hubs_diff_negative_coexp.csv", sep = ";")
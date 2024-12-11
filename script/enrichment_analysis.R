library(enrichR)
library(ggplot2)

#First, we load the hubs obtained from the co-expression network analysis
setwd("~/Documents/Master/DigitalEpidemiology/DEPM-Project")
#We extract the cancer hubs from the coexpression network
cancer_data <- read.csv("data/hubs_cancer_coexp.csv", sep = ";", row.names = 1)
#We extract the names of the hubs that are from degree i.e. have is_degree = 1
cancer_genes <-cancer_data[cancer_data$is_degree == 1,]$name

#We do the same for the differential coexpression network
diff_data <- read.csv("data/hubs_diff_coexp.csv", sep = ";", row.names = 1)
diff_genes <-diff_data[diff_data$is_degree == 1,]$name

# We select our databases for enrichment
dbs_go <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
dbs_pw <- c("KEGG_2021_Human")
dbs_dd <- c("PheWeb_2019", "ClinVar_2019")

# #Part 1: GO Enrichment
#For the cancer hubs of the coexpression network
enriched_coexp.go <- enrichr(cancer_genes, databases = dbs_go)
#For the cancer hubs of the differential coexpression network
enriched_diff.go <- enrichr(diff_genes, databases = dbs_go)

# Plotting results for co-expression hubs
plotEnrich(enriched_coexp.go[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value")+
  theme_bw()+
  #Here we change the title of the plot
  labs(title = "Co-Expression Network")
#Here we save the plot
ggsave("figures/enrichment_coexp_GO_2018.pdf", dpi=300, device="pdf",width=6, height=5, useDingbats=FALSE)

# Plotting results for differential co-expression hubs
plotEnrich(enriched_diff.go[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value")+
  theme_bw()+
  #Here we change the title of the plot
  labs(title = "Differential Co-Expression Network")
#Here we save the plot
ggsave("figures/enrichment_diff_GO_2018.pdf", dpi=300, device="pdf",width=6, height=5, useDingbats=FALSE)




# #Part 2: Pathway Enrichment
#For the cancer hubs of the coexpression network
enriched_coexp.pw <- enrichr(cancer_genes, databases = dbs_pw)
#For the cancer hubs of the differential coexpression network
enriched_diff.pw <- enrichr(diff_genes, databases = dbs_pw)

# Plotting results for co-expression hubs
plotEnrich(enriched_coexp.pw[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value") +
  theme_bw()+
  #Here we change the title of the plot
  labs(title = "Co-Expression Network")
#Here we save the plot
ggsave("figures/enrichment_coexp_KEGG_2021_Human.pdf", dpi=300, device="pdf",width=6, height=5, useDingbats=FALSE)

# Plotting results for differential co-expression hubs
plotEnrich(enriched_diff.pw[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value") +
  theme_bw()+
  #Here we change the title of the plot
  labs(title = "Differential Co-Expression Network")
#Here we save the plot
ggsave("figures/enrichment_diff_KEGG_2021_Human.pdf", dpi=300, device="pdf",width=6, height=5, useDingbats=FALSE)

# #Part 3: Disease Enrichment
#For the cancer hubs of the coexpression network
enriched_coexp.dd <- enrichr(cancer_genes, databases = dbs_dd)
#For the cancer hubs of the differential coexpression network
enriched_diff.dd <- enrichr(diff_genes, databases = dbs_dd)

# Plotting results for co-expression hubs
plotEnrich(enriched_coexp.dd[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value") +
  theme_bw()+
  #Here we change the title of the plot
  labs(title = "Co-Expression Network")
#Here we save the plot
ggsave("figures/enrichment_coexp_PheWeb_ClinVar.pdf", dpi=300, device="pdf",width=6, height=5, useDingbats=FALSE)

# Plotting results for differential co-expression hubs
plotEnrich(enriched_diff.dd[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value") +
  theme_bw()+
  #Here we change the title of the plot
  labs(title = "Differential Co-Expression Network")
#Here we save the plot
ggsave("figures/enrichment_diff_PheWeb_ClinVar.pdf", dpi=300, device="pdf",width=6, height=5, useDingbats=FALSE)


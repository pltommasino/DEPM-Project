library(devtools)
install_github("wjawaid/enrichR")
library(enrichR)

#options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
#Load DEG.csv
deg_data <- read.csv("data/DEG.csv", sep = ";", row.names = 1)

# Extract genes names
deg_genes <- rownames(deg_data)

# Find the list of all available databases from Enrichr #####NON FUNZIONA
dbs <- listEnrichrDbs()
print(head(dbs$name)) 

# Esecuzione dell'enrichment analysis
enrich_results <- enrichR::enrichr(deg_genes, databases = c("KEGG_2019_Human", "MSigDB_Hallmark_2020"))

# Visualizza i risultati per KEGG
print(head(enrich_results[["KEGG_2019_Human"]]))

# Visualizza i risultati per MSigDB Hallmark
print(head(enrich_results[["MSigDB_Hallmark_2020"]]))

library(ggplot2)

# Converti i risultati in un data frame
kegg_res <- enrich_results[["KEGG_2019_Human"]]

# Ordina i risultati per p-value
kegg_res <- kegg_res[order(kegg_res$Adjusted.P.value), ]

# Visualizzazione: primi 10 pathway arricchiti
ggplot(kegg_res[1:10, ], aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 KEGG Pathways", x = "Pathway", y = "-log10(Adjusted P-value)") +
  theme_minimal()

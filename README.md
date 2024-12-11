# Digital Epidemiology and Precision Medicine -- Final Project

This is a Github repository created to submit the final project of the **Digital Epidemiology and Precision Medicine** course for the MSc. in Data Science at the Sapienza University of Rome.

--- 
## What's inside this repository?

1. `README.md`: A markdown file that explains the content of the repository.

2. ``script/``: A folder including 4 R scripts used to perform the analysis for the report. The files included are:

    - `deg_data_preprocessing.R`: An R script to preprocess TGCA data in order to obtain the Differentially Expressed Genes (DEGs).

     - `co_exp_net.R`: An R script to obtain gene co-expression networks and differential co-expression networks for cancerous and healthy samples. The objective was to characterize hubs and scale-free behavior.

     - `ps_net.R`: An R script to obtain patient similarity networks with gene expression and mutational profiles and perform PSN and Survival Analysis to characterize clusters of patients.

    - `enrichment_analysis.py`: An R script to perform Gene Set Enrichment Analysis to the hub genes found in previous analysis.


3. `data/`: A folder containing intermediate data used for running the scripts above.

4. `figures/`: A folder containing the figures obtained given the analysis.

5. ``.gitignore``: A predetermined `.gitignore` file that tells Git which files or folders to ignore in an R project.

6. `LICENSE`: A file containing an MIT permissive license.

---

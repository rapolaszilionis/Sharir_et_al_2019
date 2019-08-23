This repository contains Jupyter notebooks with Python code recreating selected single cell RNAseq data analyses described in [1].

### Before you can run the code in this repository
To run the code in the jupyter notebooks, you need to add the following to the repository:
1) The folder "data_from_geo" containing [data from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131204),specifically:
- GSE131204_average_gene_expression_per_population.tsv
- GSE131204_cell_info_8594x25.tsv
- GSE131204_gene_names_alphabetically.txt
- GSE131204_raw_counts_8594x27998.mtx
The following files are found under the [sample GSM3767568 (control)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3767568)
- GSM3767568_control_barcodes.tsv
- GSM3767568_control_genes.tsv
- GSM3767568_control_loom.loom.hdf5
- GSM3767568_control_matrix.mtx
The following files are found under the [sample GSM3767569 (injury)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3767569)
- GSM3767569_injury_barcodes.tsv
- GSM3767569_injury_genes.tsv
- GSM3767569_injury_loom.loom.hdf5
- GSM3767569_injury_matrix.mtx

2) The repository "SPRING_dev-master", which you can download from [here](https://github.com/AllonKleinLab/SPRING).

### Interactive SPRING plots
SPRING[2,3] was used to visualize in 2D and interactively explore single cell RNAseq data. [The first notebook to run](spring_save_counts_for_all_plots.ipynb) saves the count data in a SPRING-compatible format. This is done once for all SPRING subplots, even those used to vizualise only a subset of all cells.
The following tables summarizes the notebooks to locally recreate SPRING plots shown in the manuscript together with links to online versions of the plots.

First shown in the following figure in [1] | Plot online | Jupyter notebook
 --- | --- | ---
Fig. 1d | [control_epithelial](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Sharir_et_al_2019/control_epithelial) | [spring_control_epithelial.ipynb](spring_control_epithelial.ipynb)
Fig. 1i | [control_class1_no_cc](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Sharir_et_al_2019/control_class1_no_cc) | [spring_control_class1_epithelial_no_cell_cycle.ipynb](spring_control_class1_epithelial_no_cell_cycle.ipynb)
Fig. 4b  | [control_epithelial_no_cc](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Sharir_et_al_2019/control_epithelial_no_cc) | [spring_control_epithelial_no_cell_cycle.ipynb](spring_control_epithelial_no_cell_cycle.ipynb)
Fig. 6f | [spring_control_AND_injury_epithelial](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Sharir_et_al_2019/control_and_injury_epithelial) | [spring_control_AND_injury_epithelial.ipynb](spring_control_AND_injury_epithelial.ipynb)
Fig. S5h | [spring_injury_class1_epithelial_no_cell_cycle](https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Sharir_et_al_2019/injury_class1_no_cc) | [spring_injury_class1_epithelial_no_cell_cycle.ipynb](spring_injury_class1_epithelial_no_cell_cycle.ipynb)


### Other analyses

Related to the following figures in [1] | Description | Jupyter notebook
--- | --- | ---
Fig. 1e | Identify population-enriched genes, plot heatmap, save .rnk files for GSEAPreranked tool | [population_enriched_genes_control_and_save_rnk_files.ipynb](population_enriched_genes_control_and_save_rnk_files.ipynb)
Fig. 1f | Score cell cycle, classify cells as "G1", "S", "G2/M" | [score_cell_cycle_in_control_epithelial_cells.ipynb](score_cell_cycle_in_control_epithelial_cells.ipynb)
Fig. 1i | Classify single cells by reference transcriptomes using selected genes | [classify_class1_cells_by_class2_and_3_populations.ipynb](classify_class1_cells_by_class2_and_3_populations.ipynb)
Fig. 2a | Overlay RNA velocity vectors onto chosen 2D SPRING representation | [project_RNA_velocity_on_control_epithelial_spring_plot.ipynb](project_RNA_velocity_on_control_epithelial_spring_plot.ipynb)
Fig. 2c | Order selected cells by pseudotime | [pseudotime_ordering_of_ameloblasts.ipynb](pseudotime_ordering_of_ameloblasts.ipynb)
Fig. 4g-i | Run FateID | [FateID_dataprep_and_plots.ipynb](FateID_dataprep_and_plots.ipynb)
Fig. 6g, S5c,d | Classify single cells by reference transcriptomes using all genes | [classify_injury_cells_by_control_populations.ipynb](classify_injury_cells_by_control_populations.ipynb)
Fig. S5d | Identify population-enriched genes in injury condition and plot heatmap showing similar expression patterns in counterpart populations in control condition | [gene_expression_similarities_in_same_population_ctrl_vs_injury.ipynb](gene_expression_similarities_in_same_population_ctrl_vs_injury.ipynb)


#### References.  
[1] Sharir A, Marangoni P, Zilionis R et al. A large pool of actively cycling progenitors orchestrates self-renewal and injury repair of an ectodermal appendage. Nature Cell Biology. In print.  
[2] Weinreb C, Wolock S, Klein AM. SPRING: a kinetic interface for visualizing high dimensional single-cell expression data. Bioinformatics. 2018 Apr 1;34(7):1246–8.  
[3] https://github.com/AllonKleinLab/SPRING_dev

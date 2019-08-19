This repository contains Jupyter notebooks with Python code recreating selected single cell RNAseq data analyses described in [1].

### Before you can run the code in this repository
To run the code in the jupyter notebooks, you need to add the following to the repository:
1) the folder "data_from_geo" containing data from GEO, instructions for downloading the data will be provided upon public release.
2) the repository "SPRING_dev-master", which you can download from [here](https://github.com/AllonKleinLab/SPRING).

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
Fig. ?? | Identify population-enriched genes, plot heatmap, save .rnk files for GSEAPreranked tool | [population_enriched_genes_control_and_save_rnk_files.ipynb](population_enriched_genes_control_and_save_rnk_files.ipynb)
Fig. ?? | Score cell cycle, classify cells as "G1", "S", "G2/M" | [score_cell_cycle_in_control_epithelial_cells.ipynb](score_cell_cycle_in_control_epithelial_cells.ipynb)
Fig. ?? | Overlay RNA velocity vectors onto chosen 2D SPRING representation | [project_RNA_velocity_on_control_epithelial_spring_plot.ipynb](project_RNA_velocity_on_control_epithelial_spring_plot.ipynb)
Fig. ?? | Order selected cells by pseudotime | [pseudotime_ordering_of_ameloblasts.ipynb](pseudotime_ordering_of_ameloblasts.ipynb)
Fig. ?? | Classify single cells by reference transcriptomes using all genes | [classify_injury_cells_by_control_populations.ipynb](classify_injury_cells_by_control_populations.ipynb)
Fig. ?? | Classify single cells by reference transcriptomes using selected genes | [classify_class1_cells_by_class2_and_3_populations.ipynb](classify_class1_cells_by_class2_and_3_populations.ipynb)
Fig. ?? | Identify population-enriched genes in injury condition and plot heatmap showing similar expression patterns in counterpart populations in control condition | [gene_expression_similarities_in_same_population_ctrl_vs_injury.ipynb](gene_expression_similarities_in_same_population_ctrl_vs_injury.ipynb)
Fig. ?? | Run FateID | [FateID_dataprep_and_plots.ipynb](FateID_dataprep_and_plots.ipynb),\n[FateID_R_code.ipynb](FateID_R_code.ipynb),\n[FateID_get_distance_matrix.ipynb](FateID_get_distance_matrix.ipynb)

#### References.  
[1] Sharir A, Marangoni P, Zilionis R et al. A large pool of actively cycling progenitors orchestrates self-renewal and injury repair of an ectodermal appendage. Nature Cell Biology. Accepted.
[2] Weinreb C, Wolock S, Klein AM. SPRING: a kinetic interface for visualizing high dimensional single-cell expression data. Bioinformatics. 2018 Apr 1;34(7):1246–8.  
[3] https://github.com/AllonKleinLab/SPRING_dev

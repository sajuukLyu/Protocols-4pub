# **Protocols-4pub**

## *----Reusable code collection*

##### 31 March, 2021

Multi-omics analysis protocols concluded by Yulin Lyu (Cheng Li Lab, Peking University, <lvyulin@pku.edu.cn>).

These files including three main parts:

- **Preprocess pipeline** (*etc. bulkRNA, bulkATAC, ChIP, WGBS ...*)

	The R scripts are designed to generate a directory tree containing ordered directories for intermediate data and Linux bash scripts for each step.
	In convenient to submit jobs to different nodes of a computer cluster at same time, commands for multiple samples can be divided to several batches and than submit together.

- **Downstream analysis script** (*etc. DEG, GO, GSEA ...*)

	The framework of different analyses are designed to cope with the most simple scene.

- **Visualization** (*etc. heatmap, PCA plot, Venn plot, track plot ...*)

	Some Plotting functions are designed to generate plots suitable for articles.

### File description

- Preprocess pipeline

| file name     | description                                                |
| ------------- | ---------------------------------------------------------- |
| bulkATACpre.R | Preprocess pipelines for bulk ATAC-seq                     |
| bulkRNApre.R  | Preprocess pipelines for bulk RNA-seq                      |
| ChIPpre.R     | Preprocess pipelines for ChIP-seq                          |
| WGBSpre.R     | Preprocess pipelines for Whole Genome Bisulfite Sequencing |

- Downstream analysis script

| file name                | omics | description                                           |
| ------------------------ | ----- | ----------------------------------------------------- |
| downloadData.R           | any   | Download data from public datasets                    |
| bulkRNAana_1_loadCount.R | RNA   | Load multiple results and convert to count Matrix     |
| bulkRNAana_2a_DESeq.R    | RNA   | Perform DESeq analysis (for sample without replicate) |
| bulkRNAana_2b_DESeq2.R   | RNA   | Perform DESeq2 analysis (for sample with replicates)  |

- Visualization

| file name                | omics | description                                    |
| ------------------------ | ----- | ---------------------------------------------- |
| Visulz_bulkRNA_PCA.R     | RNA   | PCA plot for bulkRNA samples                   |
| Visulz_bulkRNA_MAplot.R  | RNA   | MA plot for DEGs between group of samples      |
| Visulz_bulkRNA_volcano.R | RNA   | Volcano plot for DEGs between group of samples |
| Visulz_bulkRNA_heatmap.R | RNA   | Heatmap for bulkRNA samples                    |

### Pipeline

- bulk RNA-seq

<img src="mermaidPlot\bulkRNApre.svg" align=center>

- bulk ATAC-seq

<img src="mermaidPlot\bulkATACpre.svg" align=center>

- ChIP-seq

<img src="mermaidPlot\ChIPpre.svg" align=center>

- WGBS

<img src="mermaidPlot\WGBSpre.svg" align=center>



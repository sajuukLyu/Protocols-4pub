# Protocols-4pub

## *Reusable code collection*

##### 26 April, 2021

Multi-omics analysis protocols concluded by SajuukLyu <lvyulin@pku.edu.cn>.

---

*声明：该项目仅供学习用途，请共同遵守[开源协议](https://github.com/sajuukLyu/Protocols-4pub/blob/main/LICENSE)，禁止用于商业盈利。未经作者同意请勿转载至其它媒体。*

---

Integrated step-by-step vignettes present on the [wiki](https://github.com/sajuukLyu/Protocols-4pub/wiki).

These files including three main parts:

- **Preprocess pipeline** (*etc. bulkRNA, bulkATAC, ChIP, WGBS ...*)

	These R scripts are designed to generate a directory tree containing ordered directories for intermediate data and Linux bash scripts for each step.
	In convenient to submit jobs to different nodes of a computer cluster at same time, commands for multiple samples can be divided to several batches and than submit together.

- **Downstream analysis script** (*etc. DEG, GO, GSEA ...*)

	These scripts are the framework of different analyses designed to cope with the most simple scene of biology analysis.

- **Visualization script** (*etc. heatmap, PCA plot, Venn plot, track plot ...*)

	These scripts can generate some plots suitable for articles with very few post processing.

### File description

- Preprocess pipeline

| file name                                                    | description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [bulkATACpre.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/a_PreprocessPipeline/bulkATACpre.R) | Preprocess pipelines for bulk ATAC-seq [*(vignettes)*](https://sajuuklyu.github.io/Protocols-4pub/exampleData/ATAC/bulkATACpre.html) |
| [bulkRNApre.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/a_PreprocessPipeline/bulkRNApre.R) | Preprocess pipelines for bulk RNA-seq [*(vignettes)*](https://sajuuklyu.github.io/Protocols-4pub/exampleData/RNA/bulkRNApre.html) |
| [ChIPpre.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/a_PreprocessPipeline/ChIPpre.R) | Preprocess pipelines for ChIP-seq                            |
| [WGBSpre.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/a_PreprocessPipeline/WGBSpre.R) | Preprocess pipelines for Whole Genome Bisulfite Sequencing   |

- Downstream analysis script

| file name                                                    | omics | description                                           |
| ------------------------------------------------------------ | ----- | ----------------------------------------------------- |
| [downloadData.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/downloadData.R) | any   | Download data from public databases                   |
| [bulkRNAana_1_loadCount.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkRNAana_1_loadCount.R) | RNA   | Load multiple results and convert to count matrix     |
| [bulkRNAana_2a_DESeq.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkRNAana_2a_DESeq.R) | RNA   | Perform DESeq analysis (for sample without replicate) |
| [bulkRNAana_2b_DESeq2.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkRNAana_2b_DEseq2.R) | RNA   | Perform DESeq2 analysis (for sample with replicates)  |
| [bulkRNAana_3_GO.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkRNAana_3_GO.R) | RNA   | Perform GO analysis                                   |
| [bulkATACana_1_QC.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkATACana_1_QC.R) | ATAC  | Quality control                                       |
| [bulkATACana_2_loadCount.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkATACana_2_loadCount.R) | ATAC  | Load data and convert to peak matrix                  |
| [bulkATACana_3_annotatePeak.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkATACana_3_annotatePeak.R) | ATAC  | Annotate peaks to nearest genes                       |
| [bulkATACana_4_GO.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkATACana_4_GO.R) | ATAC  | Perform GO analysis                                   |
| [bulkATACana_5_peakTF.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/b_DownstreamAnalysisScript/bulkATACana_5_peakTF.R) | ATAC  | Perform gene regulate network analysis                |

- Visualization script

| file name                                                    | omics | description                                    |
| ------------------------------------------------------------ | ----- | ---------------------------------------------- |
| [Visulz_bulkRNA_PCA.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkRNA_PCA.R) | RNA   | PCA plot for samples                           |
| [Visulz_bulkRNA_MAplot.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkRNA_MAplot.R) | RNA   | MA plot for DEGs between group of samples      |
| [Visulz_bulkRNA_volcano.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkRNA_volcano.R) | RNA   | Volcano plot for DEGs between group of samples |
| [Visulz_bulkRNA_heatmap.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkRNA_heatmap.R) | RNA   | Heatmap of given genes for samples             |
| [Visulz_bulkRNA_GO.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkRNA_GO.R) | RNA   | GO plot for DEGs                               |
| [Visulz_bulkATAC_trackPlot.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkATAC_trackPlot.R) | ATAC  | Track plot for samples                         |
| [Visulz_bulkATAC_PCA.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkATAC_PCA.R) | ATAC  | PCA plot for samples                           |
| [Visulz_bulkATAC_heatmapPeak.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkATAC_heatmapPeak.R) | ATAC  | Heatmap of given peaks for samples             |
| [Visulz_bulkATAC_heatmapTrack.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkATAC_heatmapTrack.R) | ATAC  | Heatmap of given peak tracks for samples       |
| [Visulz_bulkATAC_peakAnnoDisp.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkATAC_peakAnnoDisp.R) | ATAC  | Histogram for peak annotation distribution     |
| [Visulz_bulkATAC_motifEnrich.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkATAC_motifEnrich.R) | ATAC  | Scatter plot for enriched TFs for peak sets    |
| [Visulz_bulkATAC_network.R](https://github.com/sajuukLyu/Protocols-4pub/blob/main/c_VisualizationScript/Visulz_bulkATAC_network.R) | ATAC  | Network plot for peak sets                     |

### Pipeline

- bulk RNA-seq

<img src="mermaidPlot\bulkRNApre.svg" align=center>

- bulk ATAC-seq

<img src="mermaidPlot\bulkATACpre.svg" align=center>

- ChIP-seq

<img src="mermaidPlot\ChIPpre.svg" align=center>

- WGBS

<img src="mermaidPlot\WGBSpre.svg" align=center>



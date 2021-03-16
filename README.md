# **Protocols-4pub**

## *----Reusable code collection*

##### 16 March, 2021

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

1. Preprocess pipeline

| **file name** | **description**                                            |
| ------------- | ---------------------------------------------------------- |
| bulkATACpre.R | Preprocess pipelines for bulk ATAC-seq                     |
| bulkRNApre.R  | Preprocess pipelines for bulk RNA-seq                      |
| ChIPpre.R     | Preprocess pipelines for ChIP-seq                          |
| WGBSpre.R     | Preprocess pipelines for Whole Genome Bisulfite Sequencing |

2. Downstream analysis script

| **file name**            | **Omics** | **description**                                       |
| ------------------------ | --------- | ----------------------------------------------------- |
| bulkRNAana_1_loadCount.R | RNA       | Load multiple results and convert to count Matrix     |
| bulkRNAana_2a_DESeq.R    | RNA       | Perform DESeq analysis (for sample without replicate) |
| bulkRNAana_2b_DESeq2.R   | RNA       | Perform DESeq2 analysis (for sample with replicates)  |

3. Visualization

### Pipeline

1. bulk RNA-seq

[![](https://mermaid.ink/svg/eyJjb2RlIjoiZ3JhcGggVEQ7XG5cdEEoXCJSYXcgZGF0YSAoLmZhc3RxLmd6KVwiKSAtLT58ZmFzdHFjfCBhKFwiUXVhbGl0eSBDb250cm9sXCIpXG5cdEEgLS0-fHRyaW1tb21hdGljfCBCKFwiVHJpbSAoLmZhc3RxLmd6KVwiKVxuXHRCIC0tPnxmYXN0cWN8IGFcblx0QiAtLT58U1RBUnwgQyhcIk1hcCAoLmJhbSlcIilcblx0QyAtLT58UnNlUUN8IEQoXCJJbmZlciBleHBlcmltZW50XCIpXG5cdEMgLS0-fGZlYXR1cmVDb3VudHwgRShcIkNvdW50ICgudHh0KVwiKVxuXHRFIC0tPiBGKFwiRG93biBzdHJlYW0gYW5hbHlzaXNcIikiLCJtZXJtYWlkIjp7InRoZW1lIjoiZGVmYXVsdCJ9LCJ1cGRhdGVFZGl0b3IiOmZhbHNlfQ)



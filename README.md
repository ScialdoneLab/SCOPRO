# SCOPRO
SCOPRO (SCOre PROjection) is an R package that assigns a score projection from 0 to 1 between a given in vivo stage and each single cluster from an in vitro dataset. The score is assigned based on the the fraction of specific markers of the in vivo stage that are conserved in the in vitro clusters.

## Installation

You can install the development version from [GitHub](https://github.com/) with:

```r
devtools::install_github("ScialdoneLab/SCOPRO")
```

## Getting started 
The main function of the package is  **SCOPRO**


### SCOPRO

```r
SCOPRO(norm_vitro, norm_vivo, cluster_vitro, cluster_vivo, name_vivo, marker_stages_filter, threshold = 0.1, number_link = 1, fold_change = 3, threshold_fold_change = 0.1,  marker_stages, selected_stages)
```
requires as input:

1. **norm_vitro**: norm count matrix (n_genes X n_cells) for in vitro dataset
2. **norm_vitro**: norm count matrix (n_genes X n_cells) for in vivo dataset
3. **norm_vitro** Norm count matrix (n_genes X n_cells) for in vitro dataset
4. **norm_vivo** Norm count matrix (n_genes X n_cells) for in vivo dataset
5. **cluster_vitro** cluster for in vitro dataset
6. **cluster_vivo**  cluster for in vivo dataset
7. **name_vivo**  name of the in vivo stage on which SCOPRO is run
8. **marker_stages_filter**  output from the function **filter_in_vitro**
9. **threshold** Numeric value. For a given gene, the jaccard index between the links from the in vivo and in vitro datasets is computed. If the jaccard index is above **threshold**, then the gene is considered to be conserved between the two datasets.
10. **number_link** Numeric value. For a given gene in the in vivo dataset with links above **number_link**, the jaccard index between the links from in vitro and in vivo dataset is computed.
11. **threshold_fold_change** Numeric value. Above **threshold** the fold change between genes is computed. Below **threshold** the difference between genes is computed.
12. **fold_change** Numeric value. For a given gene, the fold change between all the other genes is computed. If fold change is above \emph{fold_change}, then there is a link with weight 1 between the two genes.
13. **marker_stages** Second element of the list given as output by the function **select_top_markers**
14. **selected_stages** In vivo stages for which the markers are computed with the function **select_top_markers**

The mean expression profile of **marker_stages_filter** genes is computed for each cluster in the in vivo and in vitro dataset.
For a given cluster, a connectivity matrix is computed with number of rows and number of columns equal to the length of **marker_stages_filter**. Each entry (i,j)  in the matrix can be 1 if the fold_change between gene i and gene j is above **fold_change**. Otherwise is 0.
Finally the connectivity matrix of a given **name_vivo** stage and all the clusters in the in vitro dataset are compared.
A gene i is considered to be conserved between **name_vivo** and an in vitro cluster if the jaccard index of the links of gene i is above **threshold**.

## Example 
Below an example of input using the development version of **SCOPRO** from GitHub

```r
current_wd <- getwd()
url = "https://hmgubox2.helmholtz-muenchen.de/index.php/s/EHQSnjMJxkR7QYT/download/SCOPRO.zip"
destfile <- paste0(current_wd,"/SCOPRO.zip")
download.file(url, destfile, quiet = FALSE)
unzip(destfile, exdir=current_wd)
```
Load in vitro dataset (single cell RNA seq mouse data from [Iturbe et al., 2021](https://www.nature.com/articles/s41594-021-00590-w))
```r
setwd(paste0(current_wd,"/SCOPRO"))
load(file='mayra_dati_raw_0.Rda')
mayra_seurat_0=cluster_analysis_integrate_rare(mayra_dati_raw_0,"Mayra_data_0",0.1,5,30)
norm_es_vitro=as.matrix(GetAssayData(mayra_seurat_0, slot = "data",assay="RNA"))
cluster_es_vitro=as.vector(mayra_seurat_0$RNA_snn_res.0.1)
```

Load in vivo mouse dataset (bulk RNA seq data from [Deng et al. , 2014](https://www.science.org/doi/10.1126/science.1245316) and [Mohammed et al. , 2017](https://www.sciencedirect.com/science/article/pii/S2211124717309610))

```r
setwd(paste0(current_wd,"/SCOPRO"))
load(file="seurat_genes_published_mouse.Rda")

norm_vivo <- as.matrix(GetAssayData(seurat_genes_published_mouse, slot = "data",assay="RNA"))
```

Compute markers for selected in vivo stages using CIARA function **markers_cluster_seurat** based on package Seurat
```r
DefaultAssay(seurat_genes_published_mouse) <- "RNA"
cluster_mouse_published <- as.vector(seurat_genes_published_mouse$stim)


relevant_stages <- c("Late_2_cell", "epiblast_4.5", "epiblast_5.5", "epiblast_6.5")

DefaultAssay(seurat_genes_published_mouse) <- "RNA"

markers_first_ESC_small <- CIARA::markers_cluster_seurat(seurat_genes_published_mouse[,cluster_mouse_published%in%relevant_stages],cluster_mouse_published[cluster_mouse_published%in%relevant_stages],names(seurat_genes_published_mouse$RNA_snn_res.0.2)[cluster_mouse_published%in%relevant_stages],10)

markers_mouse <- as.vector(markers_first_ESC_small[[3]])
stages_markers <- names(markers_first_ESC_small[[3]])

## Keeping only the genes in common between in vitro and in vivo datasets
stages_markers <- stages_markers[markers_mouse %in% row.names(norm_es_vitro)]

markers_small <- markers_mouse[markers_mouse %in% row.names(norm_es_vitro)]
names(markers_small) <- stages_markers
```
For each in vivo stage, we select only the markers for which the median is above 0.1 and is below 0.1 in all the other stages.

```r
marker_result <- select_top_markers(relevant_stages, cluster_mouse_published, norm_vivo, markers_small, max_number = 100, threshold = 0.1)
marker_all <- marker_result[[1]]
marker_stages <- marker_result[[2]]
```

### Run SCOPRO

We run SCOPRO between the cluster of the mouse ESCs dataset and the in vivo stage "Late 2-cells".

The function SCOPRO first computes the mean expression profile of **marker_stages_filter** genes for each cluster in the in vivo and in vitro dataset.
For a given cluster, a connectivity matrix is computed with number of rows and number of columns equal to the length of **marker_stages_filter**. Each entry (i,j)  in the matrix can be 1 if the fold_change between gene i and gene j is above **fold_change**. Otherwise is 0.
Finally the connectivity matrix of Late 2-cells stage and all the clusters in the in vitro dataset are compared.
A gene i is considered to be conserved between Late 2-cells stage and an in vitro cluster if the jaccard index of the links of gene i is above **threshold**.

There are 25 markers of the Late 2-cells stage that are also expressed in the mouse ESC datasets.
More than 75% of these 25 markers are conserved in the cluster number 2.
This result is expected since cluster 2 is made up by 2CLC, a rare population of cells known to be transcriptionally similar to the late 2 cells-stage in the mouse embryo development (typical markers of 2CLC are the Zscan4 genes, also highly expressed in the late 2 cells-stage).

```r
marker_stages_filter <- filter_in_vitro(norm_es_vitro,cluster_es_vitro ,marker_all, fraction = 0.10, threshold = 0)

analysis_2cell <- SCOPRO(norm_es_vitro,norm_vivo,cluster_es_vitro,cluster_mouse_published,"Late_2_cell",marker_stages_filter, threshold = 0.1, number_link = 1, fold_change = 3, threshold_fold_change = 0.1 ,marker_stages, relevant_stages)

plot_score(analysis_2cell, marker_stages, marker_stages_filter, relevant_stages, "Late_2_cell", "Final score", "Cluster", "Late_2_cell")
```
<img src="https://github.com/ScialdoneLab/SCOPRO/blob/master/figures/SCOPRO_1.png" width="500" height="500">


### Visualization of conserved/ not conserved genes between late 2 cells stage and in vitro clusters

 We can visualize which are the markers of the late 2 cells stage that are conserved/ not conserved in cluter 2.
As expected the Zscan4 family genes are conserved.
```r
common_genes <- select_common_genes(analysis_2cell, marker_stages, relevant_stages, "Late_2_cell", cluster_es_vitro, "2")
no_common_genes <- select_no_common_genes(analysis_2cell, marker_stages, relevant_stages, "Late_2_cell", cluster_es_vitro, "2")




all_genes <- c(no_common_genes[1:4], common_genes[1:10])
all_genes_label <- c(paste0(no_common_genes[1:4], "-no_conserved"), paste0(common_genes[1:10], "-conserved"))


rabbit_plot <- plot_score_genes(all_genes, "Mouse ESC", "Mouse vitro", norm_es_vitro,norm_vivo[ , cluster_mouse_published=="Late_2_cell"],cluster_es_vitro, cluster_mouse_published[cluster_mouse_published == "Late_2_cell"], all_genes_label, 7, 10, "Late_2_cell")
rabbit_plot


```
<img src="https://github.com/ScialdoneLab/SCOPRO/blob/master/figures/SCOPRO_2.png" width="500" height="500">

## Comparison with projection based on Seurat

SCOPRO assigns an absolute score. Commonly used methods for projection generally give as output a relative score normalized to one. This is the case for the projection based on seurat functions FindTransferAnchors and TransferData. Seurat projection is implemented in the function **findCellTypesSeurat**.
The limit of giving a relative score is that in vitro cells will always be assign to an in vivo stage, even if they are not transcriptionally similar.
If the late 2 cells stage is present, then both Seurat and SCOPRO are able to correctly assign cluster 2 to the late 2 cells stage.
Below the result from Seurat projection is shown. The result from SCOPRO is shown in the previous section **Example**
```r
selected_stages = c("Late_2_cell","epiblast_4.5","epiblast_5.5","epiblast_6.5")
cluster_mouse_published_small = cluster_mouse_published[cluster_mouse_published %in% selected_stages]
cluster_mouse_rename_small=revalue(cluster_mouse_published_small, c("Late_2_cell"="Late_2cell","epiblast_4.5"="epi_4.5","epiblast_5.5"="epi_5.5","epiblast_6.5"="epi_6.5"))
seurat_genes_published_mouse_small = seurat_genes_published_mouse[, cluster_mouse_published %in% selected_stages]
```

```r
predictions <- findCellTypesSeurat(queryObj = mayra_seurat_0, referenceObj = seurat_genes_published_mouse_small, k.anchor =5, k.filter = 200, namedLabels <- cluster_mouse_rename_small, k.weight = 20)
assign("mayra_seurat_0", AddMetaData(object = mayra_seurat_0, metadata = predictions))
cluster_vitro_factor <- factor(cluster_mouse_rename_small,levels=c("Late_2cell", "epi_4.5", "epi_5.5", "epi_6.5"))
cellTypesPerClusterBalloonPlot_small(obj = mayra_seurat_0, cluster_vitro_factor,  main = "Projection Mouse ESC on mouse embryos " ,0.7, 0.4, 0.1)
 
```

<img src="https://github.com/ScialdoneLab/SCOPRO/blob/master/figures/SCOPRO_2.png" width="200" height="200">

<img src="https://github.com/ScialdoneLab/SCOPRO/blob/master/figures/SCOPRO_2.png" width="200" height="200">

If only epiblast stages from 4.5 to 6.5 are used, then Seurat will still assign cluster 2 to epiblast 4.5 and epiblast 5.5, although this cluster shares just a few markers with these in vivo stages.
On the other hand SCOPRO assigns a low score (below 0.5) in cluster 2 for both epiblast 4.5 and epiblast 5.5. 


<img src="https://github.com/ScialdoneLab/SCOPRO/blob/master/figures/SCOPRO_2.png" width="400" height="400">      <img src="https://github.com/ScialdoneLab/SCOPRO/blob/master/figures/SCOPRO_2.png" width="400" height="400">      <img src="https://github.com/ScialdoneLab/SCOPRO/blob/master/figures/SCOPRO_2.png" width="400" height="400">

## Vignette

The following vignette is available and completely reproducible. 
In this vignette it is shown the projection performed between single cell RNA seq mouse data from [Iturbe et al., 2021](https://www.nature.com/articles/s41594-021-00590-w) and in vivo mouse datasets from [Deng et al. , 2014](https://www.science.org/doi/10.1126/science.1245316) and [Mohammed et al. , 2017](https://www.sciencedirect.com/science/article/pii/S2211124717309610). 
It can be accessed within R with:
```r
utils::vignette("SCOPRO_vignette")
```




## Contributions and Support
Contributions in the form of feedback, comments, code and bug report are welcome.
* For any contributions, feel free to fork the source code and [submit a pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork).
* Please report any issues or bugs here: https://github.com/ScialdoneLab/SCOPRO/issues.
Any questions and requests for support can also be directed to the package maintainer (gabriele[dot]lubatti[at]helmholtz-muenchen[dot]de).



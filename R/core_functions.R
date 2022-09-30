

#' SCOPRO
#'
#' @param norm_vitro Norm count matrix (n_genes X n_cells) for in vitro dataset
#' @param norm_vivo Norm count matrix (n_genes X n_cells) for in vivo dataset
#' @param cluster_vitro cluster for in vitro dataset
#' @param cluster_vivo  cluster for in vivo dataset
#' @param name_vivo  name of the in vivo stage on which SCOPRO is run
#' @param marker_stages_filter  output from the function \emph{filter_in_vitro}
#' @param threshold Numeric value. For a given gene, the jaccard index between the links from the in vivo and in vitro datasets is computed. If the jaccard index is above \emph{threshold}, then the gene is considered to be conserved between the two datasets.
#' @param number_link Numeric value. For a given gene in the in vivo dataset with links above \emph{number_link}, the jaccard index between the links from in vitro and in vivo dataset is computed.
#' @param threshold_fold_change Numeric value. Above \emph{threshold} the fold change between genes is computed. Below \emph{threshold} the difference between genes is computed.
#' @param fold_change Numeric value. For a given gene, the fold change between all the other genes is computed. If fold change is above \emph{fold_change}, then there is a link with weight 1 between the two genes.
#' @param marker_stages Second element of the list given as output by the function \emph{select_top_markers}
#' @param selected_stages In vivo stages for which the markers where computed with the function \emph{select_top_markers}
#' @description The mean expression profile of \emph{marker_stages_filter} genes is computed for each cluster in the in vivo and in vitro dataset.
#' For a given cluster, a connectivity matrix is computed with number of rows and number of columns equal to the length of \emph{marker_stages_filter}. Each entry (i,j)  in the matrix can be 1 if the fold_change between gene i and gene j is above \emph{fold_change}. Otherwise is 0.
#' Finally the connectivity matrix of a given \emph{name_vivo} stage and all the clusters in the in vitro dataset are compared.
#' A gene i is considered to be conserved between \emph{name_vivo} and an in vitro cluster if the jaccard index of the links of gene i is above \emph{threshold}.
#' @return A list with five elements:
#'
#' \item{common_genes}{Vector with the names of the genes conserved between \emph{name_vivo} and all the clusters in the vitro dataset}
#' \item{no_common_genes}{Vector with the names of the genes not conserved between \emph{name_vivo} and  the clusters in the vitro dataset}
#' \item{genes_kept}{List with the names of the genes conserved between \emph{name_vivo} and each single cluster in the vitro dataset}
#' \item{genes_no_kept}{List with the names of the genes not conserved between \emph{name_vivo} and each single cluster in the vitro dataset}
#' \item{fraction_conserved_markers}{Numeric value, given by the fraction of conserved markers of \emph{name_vivo} and each single cluster in the in vitro dataset}
#' \item{final_score}{fraction of conserved links between the \emph{name_vivo} network and the network
#' from each cluster in the in vitro dataset}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @examples
#' load(system.file("extdata", "norm_es_vitro_small.Rda", package = "SCOPRO"))
#' n_es= norm_es_vitro_small
#' load(system.file("extdata", "norm_vivo_small.Rda", package = "SCOPRO"))
#' n_v = norm_vivo_small
#' load(system.file("extdata", "cluster_es_vitro_small.Rda", package = "SCOPRO"))
#' c_es=cluster_es_vitro_small
#' load(system.file("extdata", "cluster_vivo_small.Rda", package = "SCOPRO"))
#' c_v=cluster_vivo_small
#' load(system.file("extdata", "marker_stages_filter.Rda", package = "SCOPRO"))
#' m_s_f = marker_stages_filter
#' load(system.file("extdata", "marker_stages.Rda", package = "SCOPRO"))
#' m_s = marker_stages
#' stages = c("Late_2_cell","epiblast_4.5","epiblast_5.5","epiblast_6.5")
#' output_SCOPRO = SCOPRO(n_es,n_v,c_es,c_v,"Late_2_cell",m_s_f,0.1,1,3,0.1,m_s,stages)
#' plot_score(output_SCOPRO,m_s,m_s_f,stages,"Late_2_cell","Score","Cluster","2-cells")
#'
#'
#' @export SCOPRO
#'
#'
#' @importFrom magrittr %>%
#' @importFrom grDevices hcl
#' @importFrom stats quantile
#' @importFrom stats median
#' @importFrom stats cor
#' @importFrom stats as.dist
#' @importFrom utils download.file
#' @importFrom utils unzip
#' @importFrom stats sd



SCOPRO <- function (norm_vitro, norm_vivo, cluster_vitro, cluster_vivo,
                      name_vivo, marker_stages_filter, threshold = 0.1, number_link = 1,
                      fold_change = 3, threshold_fold_change = 0.1, marker_stages,
                      selected_stages)
{
  if (sum(selected_stages %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector selected_stages")
  }
  if (sum(cluster_vivo %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector cluster_vivo")
  }
  if (sum(selected_stages %in% cluster_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector cluster_vivo")
  }
  if (length(marker_stages_filter) == 0) {
    stop("Vector marker_stages_filter has 0 length. Please provide a non-zero length vector.")
  }
  mean_2_cell <- mean_stage(norm_vivo, cluster_vivo, name_vivo,
                            marker_stages_filter)
  connectivity_vivo <- find_connectivity(mean_2_cell, fold_change,
                                         threshold_fold_change)
  name_cluster <- levels(factor(cluster_vitro))
  genes_kept <- rep(list(0), length(name_cluster))
  genes_no_kept <- rep(list(0), length(name_cluster))
  fraction_conserved_markers <- rep(list(0), length(name_cluster))
  final_score <- rep(list(0), length(name_cluster))
  for (i in 1:length(name_cluster)) {
    mean_cluster <- mean_stage(norm_vitro, cluster_vitro,
                               name_cluster[i], marker_stages_filter)
    connectivity_vitro <- find_connectivity(mean_cluster,
                                            fold_change, threshold_fold_change)
    connectivity_vitro <- connectivity_vitro[row.names(connectivity_vivo),
                                             colnames(connectivity_vivo)]
    index_specific <- which(selected_stages %in% name_vivo)
    marker_specific <- marker_stages[[index_specific]]
    if (length(marker_specific) == 0) {
      stop(paste0("There are no markers for the stage: ",
                  name_vivo))
    }
    result_comparison <- comparison_vivo_vitro(connectivity_vivo,
                                               connectivity_vitro, threshold, number_link, marker_specific)
    genes_kept[[i]] <- names(result_comparison[[1]])[result_comparison[[2]] ==
                                                       1]
    if (length(genes_kept[[i]]) == 0) {
      warning(paste0("There are not conserved genes between ",
                     name_vivo, " and ", name_cluster[i]))
    }
    genes_no_kept[[i]] <- names(result_comparison[[1]])[result_comparison[[2]] !=
                                                          1]
    names(genes_kept[[i]]) <- rep(name_cluster[i], length(genes_kept[[i]]))
    names(genes_no_kept[[i]]) <- rep(name_cluster[i], length(genes_no_kept[[i]]))
    fraction_conserved_markers[[i]] <- result_comparison[[3]]
    names(fraction_conserved_markers[[i]]) <- name_cluster[i]


    final_score[[i]] <- common_link_detection(connectivity_vivo, connectivity_vitro)
    names(final_score[[i]]) <- name_cluster[i]
  }
  common_genes <- Reduce(intersect, genes_kept)
  if (length(common_genes) == 0) {
    warning(paste0("There are not conserved genes between all the in vitro cluster and ",
                   name_vivo))
  }
  no_common_genes <- Reduce(intersect, genes_no_kept)
  return(list(common_genes, no_common_genes, genes_kept, genes_no_kept,
              fraction_conserved_markers, final_score))
}


#' SCOPRO_random
#' @inheritParams SCOPRO
#' @param iteration Integer. Number of time SCOPRO algorithm should run starting from
#' @description Run SCOPRO function for a number of times equal to \emph{iteration}. Each time the original links in the in vivo and in vitro networks are randomly
#' reassigned.
#' @return A list with length equal to the number of \emph{cluster_vivo}. Each entry is a vetor with length equal to \emph{iteration} with the fraction of conserved links between the \emph{name_vivo} network and the network
#' from each cluster in the in vitro dataset.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export SCOPRO_random


SCOPRO_random <- function (norm_vitro, norm_vivo, cluster_vitro, cluster_vivo,
                           name_vivo, marker_stages_filter, fold_change = 3, threshold_fold_change = 0.1,
                           selected_stages, iteration = 100)
{
  if (sum(selected_stages %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector selected_stages")
  }
  if (sum(cluster_vivo %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector cluster_vivo")
  }
  if (sum(selected_stages %in% cluster_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector cluster_vivo")
  }
  if (length(marker_stages_filter) == 0) {
    stop("Vector marker_stages_filter has 0 length. Please provide a non-zero length vector.")
  }
  mean_2_cell <- mean_stage(norm_vivo, cluster_vivo, name_vivo,
                            marker_stages_filter)
  connectivity_vivo <- find_connectivity(mean_2_cell, fold_change,
                                         threshold_fold_change)
  name_cluster <- levels(factor(cluster_vitro))
  final_score <- rep(list(0), length(name_cluster))
  for (i in 1:length(name_cluster)) {
    mean_cluster <- mean_stage(norm_vitro, cluster_vitro,
                               name_cluster[i], marker_stages_filter)
    connectivity_vitro <- find_connectivity(mean_cluster,
                                            fold_change, threshold_fold_change)
    connectivity_vitro <- connectivity_vitro[row.names(connectivity_vivo),
                                             colnames(connectivity_vivo)]
    names(final_score[[i]]) <- name_cluster[i]
    step <- iteration/2
    tot_len <- (step * iteration) * 2
    seed_vector <- seq(1,tot_len,step)
    g=1
    for ( h in 1:iteration){
      seed_first <- seed_vector[g]
      seed_second <- seed_vector[g+1]
      random_score <- common_link_detection_random(connectivity_vivo, connectivity_vitro, seed_first, seed_second)
      final_score[[i]][h] <- random_score
      g <- g+2
    }
  }
  return((final_score))
}


#' common_link_detection
#' @noRd
common_link_detection <- function (connectivity_vivo, connectivity_vitro)
{


  frac <- rep(0, length(row.names(connectivity_vivo)))
  for (i in 1:length(row.names(connectivity_vivo))) {
    index_vivo <- which(connectivity_vivo[i, ] == 1)
    index_vitro <- which(connectivity_vitro[i, ] == 1)

    frac[i] <- length(intersect(index_vitro, index_vivo))
  }

  names(frac) <- row.names(connectivity_vivo)

  frac_final <- sum(frac)
  link_vivo <- sum(connectivity_vivo == 1)
  link_vitro <- sum(connectivity_vitro == 1)
  frac_total <- frac_final/min(link_vivo,link_vitro)


  return(frac_total)
}

#' common_link_detection_random
#' @noRd
common_link_detection_random <- function (connectivity_vivo, connectivity_vitro, seet_seed_1, seet_seed_2)
{

  connectivity_vivo_random <- make_random(connectivity_vivo, seet_seed_1)
  connectivity_vitro_random <- make_random(connectivity_vitro, seet_seed_2)
  frac <- rep(0, length(row.names(connectivity_vivo)))
  for (i in 1:length(row.names(connectivity_vivo))) {
    index_vivo <- which(connectivity_vivo_random[i, ] == 1)
    index_vitro <- which(connectivity_vitro_random[i, ] == 1)

    frac[i] <- length(intersect(index_vitro, index_vivo))
  }

  names(frac) <- row.names(connectivity_vivo)
  frac_final <- sum(frac)
  link_vivo <- sum(connectivity_vivo==1)
  link_vitro <- sum(connectivity_vitro==1)
  frac_total <- frac_final/min(link_vivo,link_vitro)

  return(frac_total)
}

#' make_random
#' @noRd
make_random <-function(matrix_input, seed = seed){
  set.seed(seed)
  for (i in 1:length(row.names(matrix_input))){
    random_weight <- sample(matrix_input[i,])
    matrix_input[i,] <- as.numeric(random_weight)
  }
  return(matrix_input)
}



#' SCOPRO_single_cell
#'
#' @inheritParams SCOPRO
#' @description The expression profile of \emph{marker_stages_filter} genes is computed for each stage in the in vivo dataset (mean) and each cell in the in vitro dataset.
#' For a given cluster, a connectivity matrix is computed with number of rows and number of columns equal to the length of \emph{marker_stages_filter}. Each entry (i,j)  in the matrix can be 1 if the fold_change between gene i and gene j is above \emph{fold_change}. Otherwise is 0.
#' Finally the connectivity matrix of a given \emph{name_vivo} stage and all the cells in the in vitro dataset are compared.
#' A gene i is considered to be conserved between \emph{name_vivo} and an in vitro cells if the jaccard index of the links of gene i is above \emph{threshold}.
#' @return A list with five elements:
#'
#' \item{common_link}{Vector with the names of the genes conserved between \emph{name_vivo} and all the cells in the vitro dataset}
#' \item{no_common_link}{Vector with the names of the genes not conserved between \emph{name_vivo} and  the cells in the vitro dataset}
#' \item{link_kept}{List with the names of the genes conserved between \emph{name_vivo} and each single cells in the vitro dataset}
#' \item{link_no_kept}{List with the names of the genes not conserved between \emph{name_vivo} and each single cells in the vitro dataset}
#' \item{final_score}{Numeric value, given by the fraction of conserved markers of \emph{name_vivo} and each single cells in the in vitro dataset}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#' @export SCOPRO_single_cell
#'
#'
#'
#'
#'



SCOPRO_single_cell <- function (norm_vitro, norm_vivo, cluster_vivo,
                                name_vivo, marker_stages_filter, threshold = 0.1, number_link = 1,
                                fold_change = 3, threshold_fold_change = 0.1, marker_stages,
                                selected_stages)
{
  if (sum(selected_stages %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector selected_stages")
  }
  if (sum(cluster_vivo %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector cluster_vivo")
  }
  if (sum(selected_stages %in% cluster_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector cluster_vivo")
  }
  if (length(marker_stages_filter) == 0) {
    stop("Vector marker_stages_filter has 0 length. Please provide a non-zero length vector.")
  }
  mean_2_cell <- mean_stage(norm_vivo, cluster_vivo, name_vivo,
                                     marker_stages_filter)
  connectivity_vivo <- find_connectivity(mean_2_cell, fold_change,
                                                  threshold_fold_change)
  name_cells <- colnames(norm_vitro)
  link_kept <- rep(list(0), length(name_cells))
  link_no_kept <- rep(list(0), length(name_cells))
  final_score <- rep(list(0), length(name_cells))
  for (i in 1:length(colnames(norm_vitro))) {
    single_cell_value= norm_vitro[marker_stages_filter,i]
    connectivity_vitro <- find_connectivity(single_cell_value,
                                                     fold_change, threshold_fold_change)
    connectivity_vitro <- connectivity_vitro[row.names(connectivity_vivo),
                                             colnames(connectivity_vivo)]
    index_specific <- which(selected_stages %in% name_vivo)
    marker_specific <- marker_stages[[index_specific]]
    if (length(marker_specific) == 0) {
      stop(paste0("There are no markers for the stage: ", name_vivo))
    }
    result_comparison <- comparison_vivo_vitro(connectivity_vivo,
                                                        connectivity_vitro, threshold, number_link, marker_specific)
    link_kept[[i]] <- names(result_comparison[[1]])[result_comparison[[2]] == 1]
    if (length(link_kept[[i]]) == 0) {
      warning(paste0("There are not conserved genes between ", name_vivo, " and ", colnames(norm_vitro)[i]))}
    link_no_kept[[i]] <- names(result_comparison[[1]])[result_comparison[[2]] != 1]
    names(link_kept[[i]]) <- rep(colnames(norm_vitro)[i], length(link_kept[[i]]))
    names(link_no_kept[[i]]) <- rep(colnames(norm_vitro)[i], length(link_no_kept[[i]]))
    final_score[[i]] <- result_comparison[[3]]
    names(final_score[[i]]) <- colnames(norm_vitro)[i]
  }
  common_link <- Reduce(intersect, link_kept)
  if (length(common_link) == 0) {
    warning(paste0("There are not conserved genes between all the in vitro cells and ", name_vivo))
  }
  no_common_link <- Reduce(intersect, link_no_kept)
  return(list(common_link, no_common_link, link_kept, link_no_kept, final_score))
}






#' mean_stage
#' @noRd
mean_stage <- function(norm_counts, stage, name_stage, marker_stages_filter){
  norm_counts <- norm_counts[,stage == name_stage]
  mean_stage <- apply(norm_counts[marker_stages_filter, ], 1, mean)
  mean_stage <- sort(mean_stage, decreasing = T)
  return(mean_stage)
}

#' find_connectivity
#' @noRd
find_connectivity <- function(mean_stage, fold_change, threshold_fold_change){
  genes_fold_change <- rep(list(0), length(mean_stage))
  for (i in 1:length(mean_stage)){
    if (mean_stage[i] <= 0.1){
      ratio <- mean_stage[i]-mean_stage
      genes_fold_change[[i]][ratio == threshold_fold_change] <- 1
      genes_fold_change[[i]][ratio < threshold_fold_change] <- 0

    }
    else{
      ratio <- mean_stage[i] / mean_stage
      genes_fold_change[[i]][ratio <= fold_change] <- 0
      genes_fold_change[[i]][ratio > fold_change] <- 1
    }
  }

  connectivity <- matrix(unlist(genes_fold_change), byrow=TRUE, nrow=length(genes_fold_change) )
  row.names(connectivity) <- names(mean_stage)
  colnames(connectivity) <- names(mean_stage)
  return(connectivity)}


#' comparison_vivo_vitro
#' @noRd
comparison_vivo_vitro <- function(connectivity_vivo, connectivity_vitro, threshold, number_link, marker_specific){
  frac <- rep(0,length(row.names(connectivity_vivo)))
  for ( i in 1:length(row.names(connectivity_vivo))){
    index_vivo <- which(connectivity_vivo[i,]==1)
    index_vitro <- which(connectivity_vitro[i,]==1)
    if (length(index_vivo) > number_link){
      frac[i] <- length(intersect(index_vitro, index_vivo))/length(union(index_vitro, index_vivo))
    }
    else{frac[i]="Excluded"}
  }

  names(frac) <- row.names(connectivity_vivo)
  frac_1 <- frac[frac != "Excluded"]
  frac_1 <- sort(frac_1,decreasing = TRUE)
  frac_final <- frac_1
  frac_final[frac_1 > threshold] <- 1
  frac_final[frac_1 <= threshold] <- 0
  frac_specific <- sum(frac_1[names(frac_1)%in%marker_specific]>threshold)/length(frac_1[names(frac_1)%in%marker_specific])
  return(list(frac_1,frac_final,frac_specific))
}



#' create_seurat_object
#' @param raw_counts Raw count matrix (n_genes X n_cells).
#' @param project_name Character name of the Seurat project.
#' @param resolution Numeric value specifying the parameter \emph{resolution}
#' used in the Seurat function \emph{FindClusters}.
#' @param neighbors Numeric value specifying the parameter \emph{k.param} in
#' the Seurat function \emph{FindNeighbors}
#' @param max_dimension Numeric value specifying the maximum number of the PCA
#' dimensions used in the parameter \emph{dims} for the Seurat function
#' \emph{FindNeighbors}
#' @param feature_genes vector of features specifying the argument
#' \emph{features} in the Seurat function \emph{RunPCA}.
#' @return Seurat object including raw and normalized counts matrices, UMAP coordinates and cluster result.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso
#' \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindClusters}
#' \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindNeighbors}
#' \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/RunPCA}
#'
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#' @export create_seurat_object
#'



create_seurat_object <- function (raw_counts, project_name, resolution, neighbors, max_dimension, feature_genes = NULL)
{
  if (!(requireNamespace("Seurat", quietly = TRUE))) {
    stop("Package Seurat needed for this function to work. Please install it: install.packages('Seurat')")
  }
  nFeature_RNA <- NULL
  seurat_object <- Seurat::CreateSeuratObject(counts = raw_counts,
                                              project = project_name)
  seurat_object <- subset(seurat_object, subset = nFeature_RNA >
                            0)
  seurat_object <- Seurat::NormalizeData(seurat_object, verbose = FALSE)
  seurat_object <- Seurat::FindVariableFeatures(seurat_object,
                                                selection.method = "vst", nfeatures = 2000)
  seurat_object <- Seurat::ScaleData(seurat_object, verbose = FALSE)
  seurat_object <- Seurat::RunPCA(seurat_object, npcs = max_dimension,
                                  verbose = FALSE, features = feature_genes)
  seurat_object <- Seurat::RunUMAP(seurat_object, reduction = "pca",
                                   dims = 1:20, set.seed = 42)
  seurat_object <- Seurat::FindNeighbors(seurat_object, reduction = "pca",
                                         dims = 1:max_dimension, k.param = neighbors)
  seurat_object <- Seurat::FindClusters(seurat_object, resolution = resolution)
  return(seurat_object)
}

#' significance_value
#' @inheritParams plot_score
#' @param threshold Numeric value.
#' @description Detect the in vitro clusters for which the difference between the final score given by SCOPRO and the median value of the final scores given by the function SCOPRO_random
#' is greater than \emph{threshold}.
#'
#'
#' @export significance_value
significance_value <- function(SCOPRO_output, SCOPRO_output_random, threshold = 0.10){
  median <- rep(list(0),length(SCOPRO_output_random))
  for(i in 1:length(median)){
    median[[i]] <- median(SCOPRO_output_random[[i]])
  }
  median_all <- unlist(median)
  fraction_explained <- (unlist(SCOPRO_output[[6]]) - median_all)/ unlist(SCOPRO_output[[6]])
  significance_clusters <- names(fraction_explained)[fraction_explained > threshold]
  return(significance_clusters)
}




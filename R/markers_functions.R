#' select_top_markers
#'
#' @inheritParams SCOPRO
#' @param selected_stages Character vector with the name of the selected in vivo stages
#' @param markers_small Output given by the function \emph{markers_cluster_seurat} of the package CIARA
#' @param threshold Numeric value.
#' @param max_number Numeric value. Maximum number of top markers to consider for each stage in \emph{selected_stages}

#' @description For each stage in \emph{selected_stages}, starting from the markers given by \emph{markers_cluster_seurat} function of the package CIARA, only the markers
#' with a median above \emph{threshold} in the stage and below \emph{threshold} in all the other stages are kept.

#' @return A list with two elements:
#'
#' \item{marker_all}{Vector with the union of all the \emph{top_number} markers for each stage in \emph{selected_stages}}
#' \item{marker_stages}{List with length equal to number of stages in \emph{selected_stages} . Each element contains the \emph{top_number} markers for a given stage in \emph{selected_stages}  }
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{ https://CRAN.R-project.org/package=CIARA}
#'
#' @export select_top_markers
#' @importFrom stats median




select_top_markers <- function(selected_stages,cluster_vivo, norm_vivo, markers_small, max_number = 100, threshold = 0.1){
  if (length(markers_small) == 0) {
    stop("Length of markers_small is zero. Please provide a non-zero length vector.")
  }
  if (!all(selected_stages %in% cluster_vivo)) {
    stop("One or more stages are not present. Please check that all stages in selected_stages are also present in cluster_vivo")
  }
  marker_stages=rep(list(0),length(selected_stages))
  for (i in 1:length(selected_stages)){
    white_black=white_black_marker(cluster_vivo[cluster_vivo%in%selected_stages],selected_stages[i],norm_vivo[,cluster_vivo%in%selected_stages],markers_small,threshold)
    if (sum(white_black) == 0){
      warning(paste0("For stage ",selected_stages[i]," no markers were selected"))
    }
    marker_stages[[i]] <- names(white_black)[white_black]
    if (sum(white_black) > max_number){
      marker_stages[[i]] <- marker_stages[[i]][1:max_number]
    }
  }
  marker_all <- unlist(marker_stages)
  return(list(marker_all,marker_stages))
}



#' select_common_genes
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @param name_vitro  name of the in vitro stage for which we want to know the conserved markers with the \emph{name_vivo} stage
#' @param SCOPRO_output  output given by function \emph{SCOPRO}


#' @return Character vector with the names of the conserved markers of \emph{name_vivo} stage in the \emph{name_vitro} stage
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export select_common_genes
#'

select_common_genes <- function(SCOPRO_output, marker_stages, selected_stages, name_vivo, cluster_vitro, name_vitro){
  if (sum(selected_stages %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector selected_stages")

  }
  if (sum(cluster_vitro %in% name_vitro) == 0) {
    stop("name_vitro must be one the stages present in the vector cluster_vitro")

  }
  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  index_vitro <- which(levels(factor(cluster_vitro)) %in% name_vitro)
  if (sum(SCOPRO_output[[3]][[index_vitro]] %in% marker_specific) == 0) {
    warning("There are not conserved genes between the in vivo stage and the in vitro cluster")

  }
  common_genes <- SCOPRO_output[[3]][[index_vitro]][SCOPRO_output[[3]][[index_vitro]] %in% marker_specific]
  return(common_genes)
}


#' select_no_common_genes
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams select_common_genes
#' @param name_vitro  name of the in vitro stage for which we want to know the non-conserved markers with the \emph{name_vivo} stage
#' @param SCOPRO_output  output given by function \emph{SCOPRO}


#' @return Character vector with the names of the non-conserved markers of \emph{name_vivo} stage in the \emph{name_vitro} stage
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export select_no_common_genes
#'

select_no_common_genes <- function(SCOPRO_output,marker_stages,selected_stages,name_vivo,cluster_vitro,name_vitro){
  if (sum(selected_stages %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector selected_stages")

  }
  if (sum(cluster_vitro %in% name_vitro) == 0) {
    stop("name_vitro must be one the stages present in the vector cluster_vitro")

  }
  index_specific <- which(selected_stages%in%name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  index_vitro <- which(levels(factor(cluster_vitro)) %in% name_vitro)
  if (sum(SCOPRO_output[[4]][[index_vitro]] %in% marker_specific) == 0) {
    warning("There aren't non-conserved genes between the in vivo stage and the in vitro cluster")

  }
  no_common_genes <- SCOPRO_output[[4]][[index_vitro]][SCOPRO_output[[4]][[index_vitro]] %in% marker_specific]
  return(no_common_genes)
}


#' white_black_marker
#' @noRd

white_black_marker <- function(cluster_vivo, name_vivo, norm_vivo, marker_list, threshold = 0){

  white_black <- apply(norm_vivo[marker_list[names(marker_list) == name_vivo], ], 1, function(x){
    mean_one <- median(x[cluster_vivo == name_vivo])

    mean_different <- rep(0,length(levels(factor(cluster_vivo[cluster_vivo != name_vivo]))))
    for ( i in 1:length(levels(factor(cluster_vivo[cluster_vivo!=name_vivo])))){
      mean_different[i] <- median(x[cluster_vivo == levels(factor(cluster_vivo[cluster_vivo != name_vivo]))[i]])}
    if (mean_one > threshold & (sum(mean_different< threshold) == (length(levels(factor(cluster_vivo)))-1))){
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  })
  return(white_black)


}


#' filter_in_vitro
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams select_common_genes
#' @param marker_all First element of the list given as output by the function \emph{select_top_markers}
#' @param fraction Numeric value.
#' @param threshold Numeric value
#' @description For a given gene in in \emph{marker_all}, if the fraction of cells in one or more clusters with an expression above \emph{threshold} is greater than \emph{fraction}, then the gene
#' is kept

#' @return Character vector with the names of kept genes
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export filter_in_vitro
#'
#'
filter_in_vitro <- function(norm_vitro, cluster_vitro, marker_all, fraction = 0.10, threshold  = 0){
  if (length(marker_all) < 2) {
    stop("Vector marker_all must be at least of length 2")
  }
  livelli <- levels(factor(cluster_vitro))
  result_livelli <- rep(list(0),length(livelli))
  for ( i in 1:length(livelli)){
    filter_gene <- apply(norm_vitro[marker_all, cluster_vitro == livelli[i]],1,function(x){
      y <- x[x > threshold]
      ratio <- length(y) / length(x)
      if (ratio >= fraction){
        return(TRUE)
      }
      else{return(FALSE)}
    })
    result_livelli[[i]] <- filter_gene
  }
  df <- data.frame(matrix(unlist(result_livelli), nrow=length(result_livelli), byrow=TRUE))
  final_logic <- apply(df,2, sum)
  hvg_vivo_final <- marker_all[final_logic > 0]
  return(hvg_vivo_final)
}


#' select_common_genes_sc
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams select_common_genes
#' @param threshold  fraction of cells in the in vitro cluster \emph{name_vitro}.
#' Only the markers conserved in more than \emph{threshold} fraction of cells are given as output.
#' @param SCOPRO_output  output given by function \emph{SCOPRO_single_cell}


#' @return Character vector with the names of the conserved markers of \emph{name_vivo} stage in the cells of the \emph{name_vitro} cluster.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export select_common_genes_sc
#'


select_common_genes_sc <- function (SCOPRO_output, marker_stages, selected_stages, name_vivo, cluster_vitro, name_vitro, threshold = 0.5)
{
  if (sum(name_vivo %in% selected_stages) == 0) {
    stop("name_vivo must be one the stages present in the vector selected_stages")
  }
  if (sum(name_vitro %in% cluster_vitro) == 0) {
    stop("name_vitro must be one the stages present in the vector cluster_vitro")
  }
  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  cells_specific <- SCOPRO_output[[3]][cluster_vitro == name_vitro]
  common_genes_cells <- names(table(unlist(cells_specific)))[table(unlist(cells_specific)) > threshold * sum(cluster_vitro == name_vitro)]
  if (length(common_genes_cells) == 0) {
    warning("There are not conserved genes between the cells in the in vitro cluster")
    return(NULL)
  }
  if (length(common_genes_cells) > 0) {
    common_genes <- common_genes_cells[common_genes_cells %in%
                                         marker_specific]
    if (length(common_genes) == 0) {
      warning("There are not conserved genes between the cells in the in vitro cluster and the in vivo stage")}
    return(common_genes)
  }
  else{
    return(NULL)
  }
}


#' select_no_common_genes_sc
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams select_common_genes
#' @inheritParams select_no_common_genes
#' @param name_vitro  name of the in vitro stage for which we want to know the non-conserved markers with the \emph{name_vivo} stage
#' @param SCOPRO_output output given by function \emph{SCOPRO_single_cell}
#' @param threshold  fraction of cells in the in vitro cluster \emph{name_vitro}.
#' Only the markers not conserved in more than \emph{threshold} fraction of cells are given as output.


#' @return Character vector with the names of the non-conserved markers of \emph{name_vivo} stage in the \emph{name_vitro} cluster
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export select_no_common_genes_sc
#'

select_no_common_genes_sc <- function (SCOPRO_output, marker_stages, selected_stages, name_vivo, cluster_vitro, name_vitro, threshold = 0.5)
{
  if (sum(name_vivo %in% selected_stages) == 0) {
    stop("name_vivo must be one the stages present in the vector selected_stages")
  }
  if (sum(name_vitro %in% cluster_vitro) == 0) {
    stop("name_vitro must be one the stages present in the vector cluster_vitro")
  }
  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  cells_specific <- SCOPRO_output[[4]][cluster_vitro == name_vitro]
  common_genes_cells <- names(table(unlist(cells_specific)))[table(unlist(cells_specific)) > threshold * sum(cluster_vitro == name_vitro)]
  if (length(common_genes_cells) == 0) {
    warning("There are not conserved genes between the cells in the in vitro cluster")
    return(NULL)
  }
  if (length(common_genes_cells) > 0) {
    common_genes <- common_genes_cells[common_genes_cells %in% marker_specific]
    if (length(common_genes) == 0) {
      warning("There are not conserved genes between the cells in the in vitro cluster and the in vivo stage")}
    return(common_genes)
  }
  else{
    return(NULL)
  }
}







#' markers_cluster_seurat
#'
#' The Seurat function \emph{FindMarkers} is used to identify general marker
#' for each cluster (specific cluster vs all other cluster). This list of
#' markers is then filtered keeping only the genes that appear as markers in a
#' unique cluster.
#'
#' @param seurat_object Seurat object as returned by
#' \emph{create_seurat_object}.
#' @param cell_names Vector of length equal to the number of cells, with cell
#' names.
#' @param number_top Integer. Number of top marker genes to keep for each
#' cluster.
#' @param cluster Vector of length equal to the number of cells, with cluster
#' assignment.
#' @return List of three elements. The first is a vector with \emph{number_top}
#' marker genes for each cluster. The second is a vector with \emph{number_top}
#' marker genes and corresponding cluster. The third element is a vector with
#' all marker genes for each cluster.
#'
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/FindMarkers}
#'
#' @export markers_cluster_seurat


markers_cluster_seurat <- function (seurat_object, cluster, cell_names, number_top)
{
  if (sum(!unlist(lapply(list(c("Seurat","Biobase")),requireNamespace, quietly = TRUE)))>0) {
    stop("Package Seurat and Biobase needed for this function to work. Please install them: install.packages('Seurat') and BiocManager::install('Biobase')")
  }
  level <- levels(as.factor(cluster))
  marker_cluster <- as.list(rep(0, length(level)))
  marker_all <- de_seurat_cluster(seurat_object, cluster, cell_names, 0.05)
  for (i in 1:length(level)) {
    marker_cluster[[i]] <- marker_all[[i]]
    message(paste0("Cluster ", level[i]))
  }
  marker_all <- unlist(marker_cluster)
  marker_all <- marker_all[Biobase::isUnique(marker_all)]
  for (i in 1:length(level)) {
    marker_cluster[[i]] <- marker_cluster[[i]][marker_cluster[[i]] %in%
                                                 marker_all]
  }
  marker_to_see <- as.list(rep(0, length(level)))
  marker_complete <- as.list(rep(0, length(level)))
  for (i in 1:length(level)) {
    if (length(marker_cluster[[i]]) >= number_top) {
      marker_to_see[[i]] <- marker_cluster[[i]][1:number_top]
      names(marker_to_see[[i]]) <- rep(level[i], length(marker_to_see[[i]]))
      marker_complete[[i]] <- marker_cluster[[i]]
      names(marker_complete[[i]]) <- rep(level[i], length(marker_complete[[i]]))
    }
    if ((length(marker_cluster[[i]]) > 0) & (length(marker_cluster[[i]]) <
                                             number_top)) {
      marker_to_see[[i]] <- marker_cluster[[i]][1:length(marker_cluster[[i]])]
      names(marker_to_see[[i]]) <- rep(level[i], length(marker_to_see[[i]]))
      marker_complete[[i]] <- marker_cluster[[i]]
      names(marker_complete[[i]]) <- rep(level[i], length(marker_complete[[i]]))
    }

    if (length(marker_cluster[[i]]) == 0) {
      marker_to_see[[i]] <- NULL
      marker_complete[[i]] <- NULL

    }
  }
  marker_top <- unlist(marker_to_see)
  marker_top <- marker_top[marker_top != 0]
  marker_complete <- unlist(marker_complete)
  marker_complete <- marker_complete[marker_complete != 0]
  marker_all_cluster_sing <- as.list(rep(0, length(level)))
  for (i in 1:length(level)) {
    if (sum(names(marker_top) == level[i]) > 0) {
      marker_all_cluster_sing[[i]] <- paste(marker_top[names(marker_top) ==
                                                         level[i]], level[i], sep = "_")
    }

    else {

    }
  }
  marker_all_cluster <- unlist(marker_all_cluster_sing)
  marker_all_cluster <- marker_all_cluster[marker_all_cluster !=
                                             0]
  return(list(marker_top, marker_all_cluster, marker_complete))
}


#' de_seurat_cluster
#' @noRd
de_seurat_cluster <- function(seurat_object, cluster, names_cell, max_p_value) {
  level <- levels(as.factor(cluster))
  final_markers <- vector("list", length(level))
  for (i in 1:length(level)) {
    markers <- Seurat::FindMarkers(seurat_object, ident.1 = names_cell[cluster == level[i]], ident.2 = names_cell[cluster!=level[i]], only.pos = T)
    markers_new <- markers[markers$p_val_adj <= max_p_value, ]
    markers_new <- markers_new[order(markers_new$p_val_adj), ]
    markers_final <- row.names(markers_new)
    final_markers[[i]] <- markers_final
  }
  return(final_markers)
}




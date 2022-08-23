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


#' detect_expressed_genes
#' @noRd
detect_expressed_genes <- function(norm_vitro, marker_stages, selected_stages, name_vivo, threshold = 0){
  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  if (length(marker_specific) == 0) {
    stop(paste0("There are no markers for the stage: ",
                name_vivo))
  }

  marker_specific <- marker_specific[marker_specific %in% row.names(norm_vitro)]
  expressed_genes <- apply(norm_vitro[marker_specific,],2,function(x){
    number_genes <- sum(x > threshold)
    return(number_genes)
  })
  return(expressed_genes)
}

#' detect_names_genes
#' @noRd
detect_names_genes <- function(norm_vitro, marker_stages, selected_stages, name_vivo, threshold = 0){
  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  if (length(marker_specific) == 0) {
    stop(paste0("There are no markers for the stage: ",
                name_vivo))
  }

  marker_specific <- marker_specific[marker_specific %in% row.names(norm_vitro)]

  name_genes <- apply(norm_vitro[marker_specific,],2,function(x){
    index_genes <- which(x > threshold)
    names <- row.names(norm_vitro[marker_specific,])[index_genes]
    return(names)
  })
  return(name_genes)
}


#' detect_names_genes
#' @noRd
names_marker_genes = function (list_intersect, genes_sort, max_number = 5)
{
  text <- rep(0, length(list_intersect))
  for (i in 1:length(list_intersect)) {
    if (length(list_intersect[[i]]) == 0) {
      marker_genes <- "No markers expressed"
      text[i] <- "No markers expressed"
    }
    else {
      marker_genes <- genes_sort[which(genes_sort %in%
                                         list_intersect[[i]])]
      max_genes <- min(max_number, length(marker_genes))
      if (max_genes < max_number) {
        warning("max_number bigger than length of marker_genes. Set max_number = length(marker_genes) ")
      }
      marker_genes <- marker_genes[1:max_genes]
      text[i] <- paste(marker_genes, collapse = "-")
    }
  }
  return(text)
}

#' find_median_cluster
#' @noRd
find_median_cluster <- function(cluster_vivo, selected_stages, norm_vivo, genes){
  number_cluster <- selected_stages
  median_es <- rep(list(0),length(number_cluster))
  for (i in 1:length(number_cluster)){
    median_es[[i]] <- apply(norm_vivo[genes, cluster_vivo == number_cluster[i]],1,'median')
  }
  return(median_es)}

#' find_mean_cluster
#' @noRd
find_mean_cluster <- function(cluster_vivo, selected_stages, norm_vivo, genes){
  number_cluster <- selected_stages
  median_es <- rep(list(0),length(number_cluster))
  for (i in 1:length(number_cluster)){
    median_es[[i]] <- apply(norm_vivo[genes, cluster_vivo == number_cluster[i]],1,'mean')
  }
  return(median_es)}



#' gg_color_hue
#' @noRd
gg_color_hue = function (n)
{
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' create_barplot
#' @noRd

create_barplot <- function(cell_cycle, cluster, title,fraction = TRUE){
  dfForPlot <- data.frame(Phases = cell_cycle, clu = cluster)
  if(fraction == TRUE){
    pnm <- ggplot2::ggplot(dfForPlot, ggplot2::aes(y = cell_cycle, x = cluster)) + ggplot2::geom_bar( ggplot2::aes(fill = cell_cycle, y = (..count..)/sum(..count..)), position = "fill") + ggplot2::ggtitle(title) + ggplot2::ylab("Fraction of cells") + ggplot2::xlab("Cluster") + ggplot2::labs(fill = "In vivo stage")
    pnm}

  if(fraction != TRUE){
    pnm <- ggplot2::ggplot(dfForPlot, ggplot2::aes(y = cell_cycle, x = cluster)) + ggplot2::geom_bar( ggplot2::aes(fill = cell_cycle, y = (..count..), position="fill"))+ ggplot2::ggtitle(title) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Cluster") + ggplot2::labs(fill = "In vivo stage")
    pnm
  }
  return(list(pnm))}



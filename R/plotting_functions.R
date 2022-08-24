
#' plot_score_genes
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @param markers_to_plot Character vector with the names of the genes to plot.
#' @param label_1 Character value. Label for the in vitro dataset
#' @param label_2 Character value. Label for the in vivo dataset
#' @param final_name Character vector with the names of the genes to show in the plot.
#' @param max_size Numeric value, specifying the size of the dot.
#' @param text_size Numeric value, specifying the size of the text in the plot.
#' @param title_name Character value.

#' @return ggplot2::ggplot2 object showing balloon plot.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export plot_score_genes
#'


plot_score_genes <- function(markers_to_plot ,label_1, label_2, norm_vitro,norm_vivo, cluster_vitro,cluster_vivo, final_name, max_size = 9, text_size= 9.5, title_name)
{

  markers_to_plot <- markers_to_plot[markers_to_plot %in% row.names(norm_vitro)]



  all_markers_common_name <- final_name

  label_mouse <- c(rep(label_1,length(colnames(norm_vitro))),rep(label_2,length(colnames(norm_vivo))))
  mouse_norm_common=cbind(norm_vitro[markers_to_plot,],norm_vivo[markers_to_plot,])




  cluster_vitro_name <- paste(label_1, factor(cluster_vitro), sep="-")
  cluster_vivo_name <- paste(label_2, factor(cluster_vivo), sep="-")
  mouse_norm_common <- cbind(norm_vitro[markers_to_plot, ], norm_vivo[markers_to_plot, ])
  label_mouse <- c(cluster_vitro_name, cluster_vivo_name)
  out_2 <- plot_balons(mouse_norm_common, label_mouse,markers_to_plot, length(markers_to_plot), length(colnames(mouse_norm_common)), max_size = max_size,text_size = text_size, all_markers_common_name, keep_order = FALSE)
  out_2 <- out_2 + ggplot2::ggtitle(title_name)

  return(out_2)
}


#' plot_score
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams  select_common_genes
#' @inheritParams plot_score_genes
#' @param y_name Character value
#' @param fill_name Character value.

#' @return ggplot2::ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export plot_score
#'

plot_score <- function(SCOPRO_output, marker_stages, marker_stages_filter, selected_stages, name_vivo, y_name, fill_name, title_name){
  final_score <- c(as.vector(unlist(SCOPRO_output[[5]])))
  label <- c(names(unlist(SCOPRO_output[[5]])))
  df <- data.frame(final_score,label)

  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = label, y = final_score, fill = label)) +
    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_minimal()
  p + ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
             axis.title.x = ggplot2::element_blank()) + ggplot2::ylab(y_name) + ggplot2::labs(fill=fill_name) +  ggplot2::ylim(c(0,1)) + ggplot2::ggtitle(paste0(title_name," (",sum(marker_stages_filter %in% marker_specific)," genes)",sep=" "))
}


#' plot_score_sc
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams  select_common_genes
#' @inheritParams plot_score_genes
#' @inheritParams plot_score
#' @param SCOPRO_output output given by function \emph{SCOPRO_single_cell}

#' @return ggplot2::ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export plot_score_sc
#'

plot_score_sc <- function (SCOPRO_output, cluster_vitro, marker_stages, marker_stages_filter, selected_stages, name_vivo, y_name, fill_name, title_name)
{
  final_score <- c(as.vector(unlist(SCOPRO_output[[5]])))
  label <- cluster_vitro
  df <- data.frame(final_score, label)
  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = label, y = final_score,
                                        fill = label)) + ggplot2::geom_boxplot()+ggplot2::geom_jitter(position=ggplot2::position_jitter(0.2)) +
    ggplot2::theme_minimal()
  p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         vjust = 1, hjust = 1), axis.title.x = ggplot2::element_blank()) +
    ggplot2::ylab(y_name) + ggplot2::labs(fill = fill_name) +
    ggplot2::ylim(c(-0.1, 1)) + ggplot2::ggtitle(paste0(title_name,
                                                        " (", sum(marker_stages_filter %in% marker_specific),
                                                        " genes)", sep = " "))
}




#' plot_balons
#' @noRd
#'
plot_balons <- function(norm_counts, final_cluster, genes_to_plot, max_number, total_cells, max_size=5, text_size=7,label_final, keep_order=TRUE){

  fattore <- factor(final_cluster,levels=unique(final_cluster))
  if (keep_order == FALSE){
    fattore <- factor(final_cluster)}
  livelli <- levels(factor(fattore))



  all_markers <- genes_to_plot

  all_markers <- lapply(all_markers, function(x) {x[x!=0]})


  all_markers_values <- rep(list(0),length(livelli))
  for ( i in 1:length(livelli)){
    gene_values_0 <- apply(norm_counts[genes_to_plot,final_cluster==livelli[i]], 1, mean)
    all_markers_values[[i]] <- gene_values_0
  }


  values <- rep(list(0),length(livelli))
  for ( i in 1:length(livelli)){
    valore=apply(norm_counts[genes_to_plot,final_cluster==livelli[i]],1,function(x){sum(x>0.1)/length(x)})
    values[[i]]=valore
  }

  cluster_label_all=rep(list(0),length(livelli))
  for ( i in 1:length(livelli)){
    cluster_label_all[[i]] <- rep(livelli[i],length(genes_to_plot))
  }


  all_markers_plot <- genes_to_plot
  cluster_label_all <- unlist(cluster_label_all)
  values <- unlist(values)
  all_markers_values <- unlist(all_markers_values)
  balloon_melted <- data.frame(all_markers_plot,cluster_label_all,values,label_final)
  colnames(balloon_melted) <- c("genes_plot","cluster_label_all","values","label_all_final")

  p <- ggplot2::ggplot(balloon_melted, ggplot2::aes(x =factor(cluster_label_all,levels=unique(cluster_label_all)), y = factor(balloon_melted[,4],levels=(unique(label_final)))))
  p + ggplot2::geom_point( ggplot2::aes(size=values,colour=as.numeric(all_markers_values))) + ggplot2::theme(panel.background=ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "blue", fill=NA, size=3))+ggplot2::scale_size_area(max_size=max_size)+ ggplot2::scale_colour_gradient(low = "grey", high = "brown", na.value = NA)+ggplot2::labs(colour="Mean",size="Value")+ggplot2::xlab("Cluster")+ggplot2::ylab("Markers")+ggplot2::theme(text = ggplot2::element_text(size=text_size),axis.text.x = ggplot2::element_text(angle=45,vjust=1,hjust=1),axis.title.x=ggplot2::element_blank())
}


#' cellTypesPerClusterBalloonPlot
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams  select_common_genes
#' @inheritParams plot_score_genes
#' @description All the stages in \emph{cluster_vivo_factor} are showns in the balloon plot
#' @param obj Seurat object given as output from function \emph{findCellTypesSeurat}
#' @param cluster_vivo_factor Factor vector specifyng the cluster partition of the in vivo datataset
#' @param order_label_vivo Character vector specifyng the order of the columns in the balloon plot. The names must be present in \emph{cluster_vivo_factor}
#' @param label_size Numeric value specifyng the label.size parameter in the function \emph{balloonplot} from package \emph{gplots}.
#' @param thresold Numeric value. Cells with a predicted score below or equal to \emph{thresold} are labelled as unassigned. Only used when \emph{method} is Seurat
#' @param method Character name. Must be one of Seurat, scibetR or scmap

#' @return  balloon plot given by function \emph{balloonplot} from package \emph{gplots}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export cellTypesPerClusterBalloonPlot
#' @seealso \url{https://CRAN.R-project.org/package=gplots}

cellTypesPerClusterBalloonPlot <- function (obj, cluster_vivo_factor, order_label_vivo, method = "Seurat", title_name, text_size = 0.4, label_size = 0.4, thresold = 0.5)
{if (method == "Seurat"){
  if (!(requireNamespace("gplots", quietly = TRUE))) {
    stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
  }
  obj@meta.data$predicted.id[obj@meta.data$prediction.score.max <=
                               thresold] <- "Unassigned"
  levels_vitro <- c(levels(cluster_vivo_factor))
  levels_vitro_factor <- factor(levels_vitro, levels = order_label_vivo)
  if (sum(obj@meta.data$predicted.id == "Unassigned") > 0) {
    levels_vitro <- c(levels(cluster_vivo_factor), "Unassigned")
    levels_vitro_factor <- factor(levels_vitro, levels = c(order_label_vivo,
                                                           "Unassigned"))
  }
  cellTypes <- rep(levels(levels_vitro_factor), each = length(unique(as.vector(obj$seurat_clusters))))
  clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(levels(levels_vitro_factor)))
  nElements <- numeric()
  for (cellType in levels(levels_vitro_factor)) {
    nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predicted.id) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
  }
  data <- data.frame(cellTypes, clusters, nElements)
  for (cluster in unique(data$clusters)) {
    data$nElements[which(data$clusters == cluster)] <- 100 *
      data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters == cluster)])
  }
  return(gplots::balloonplot(data$cellTypes, data$clusters, data$nElements, label.size = label_size, text.size = text_size,
                             ylab = "", xlab = "", main = title_name))
}
  if (method == "scibetR"){
    if (!(requireNamespace("gplots", quietly = TRUE))) {
      stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
    }
    levels_vitro <- c(levels(cluster_vivo_factor))
    levels_vitro_factor <- factor(levels_vitro, levels = order_label_vivo)
    cellTypes <- rep(levels(levels_vitro_factor), each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(levels(levels_vitro_factor)))
    nElements <- numeric()
    for (cellType in levels(levels_vitro_factor)) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predictions) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 *
        data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters ==  cluster)])}
    return(gplots::balloonplot(data$cellTypes, data$clusters,
                               data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))
  }
  if (method == "scmap"){
    if (!(requireNamespace("gplots", quietly = TRUE))) {
      stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
    }
    levels_vitro <- c(levels(cluster_vivo_factor), "unassigned")
    levels_vitro_factor <- factor(levels_vitro, levels = c(order_label_vivo,
                                                           "unassigned"))
    cellTypes <- rep(levels(levels_vitro_factor), each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(levels(levels_vitro_factor)))
    nElements <- numeric()
    for (cellType in levels(levels_vitro_factor)) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predictions) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 *
        data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters == cluster)])
    }
    return(gplots::balloonplot(data$cellTypes, data$clusters,
                               data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))

  }
  if (method != "Seurat" & method != "scibetR" & method != "scmap"){
    stop("method must be equal to Seurat, ScibetR or scmap.")
  }

}





#' cellTypesPerClusterBalloonPlot_small
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams  select_common_genes
#' @inheritParams plot_score_genes
#' @inheritParams cellTypesPerClusterBalloonPlot
#' @return  balloon plot given by function \emph{balloonplot} from package \emph{gplots}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#' @description Only the stages in \emph{cluster_vivo_factor} with an assignment greater than zero are showns in the balloon plot
#' @export cellTypesPerClusterBalloonPlot_small
#' @seealso \url{https://CRAN.R-project.org/package=gplots}




cellTypesPerClusterBalloonPlot_small <- function (obj, cluster_vivo_factor, order_label_vivo, title_name, method = "Seurat",text_size = 0.4,
                                                  label_size = 0.7, thresold = 0.5)
{
  if (method == "Seurat"){
    if (!(requireNamespace("gplots", quietly = TRUE))) {
      stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
    }
    obj@meta.data$predicted.id[obj@meta.data$prediction.score.max <= thresold] <- "Unassigned"
    levels_vitro <- c(levels(cluster_vivo_factor))
    levels_vitro_factor <- factor(levels_vitro, levels = c(order_label_vivo))
    if (sum(obj@meta.data$predicted.id == "Unassigned") > 0) {
      levels_vitro <- c(levels(cluster_vivo_factor), "Unassigned")
      levels_vitro_factor <- factor(levels_vitro, levels = c(order_label_vivo,"Unassigned"))
    }
    cellTypes <- rep(levels_vitro[levels_vitro %in% unique(as.vector(obj@meta.data$predicted.id))], each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(unique(as.vector(obj@meta.data$predicted.id))))
    nElements <- numeric()
    for (cellType in levels_vitro[levels_vitro %in%
                                  unique(as.vector(obj@meta.data$predicted.id))]) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predicted.id) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 * data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters == cluster)])
    }
    return(gplots::balloonplot(data$cellTypes, data$clusters, data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))
  }
  if (method == "scibetR"){
    if (!(requireNamespace("gplots", quietly = TRUE))) {
      stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
    }
    cellTypes <- rep(levels(cluster_vivo_factor)[levels(cluster_vivo_factor) %in% unique(as.vector(obj@meta.data$predictions))], each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(unique(as.vector(obj@meta.data$predictions))))
    nElements <- numeric()
    for (cellType in levels(cluster_vivo_factor)[levels(cluster_vivo_factor) %in% unique(as.vector(obj@meta.data$predictions))]) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predictions) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 * data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters == cluster)])
    }
    print(data)
    return(gplots::balloonplot(data$cellTypes, data$clusters,
                               data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))
  }
  if (method == "scmap"){
    levels_vitro <- c(levels(cluster_vivo_factor), "unassigned")
    levels_vitro_factor <- factor(levels_vitro, levels =  c(order_label_vivo,"unassigned"))
    cellTypes <- rep(levels_vitro[levels_vitro %in%
                                    unique(as.vector(obj@meta.data$predictions))], each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(unique(as.vector(obj@meta.data$predictions))))
    nElements <- numeric()
    for (cellType in levels_vitro[levels_vitro %in% unique(as.vector(obj@meta.data$predictions))]) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predictions) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 *
        data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters ==
                                                                                   cluster)])
    }
    return(gplots::balloonplot(data$cellTypes, data$clusters,
                               data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))
  }
  if (method != "Seurat" & method != "scibetR" & method != "scmap"){
    stop("method must be equal to Seurat, scibetR or scmap")
  }
}



#' plot_boxplot
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams select_common_genes
#' @inheritParams plot_score_genes
#' @inheritParams plot_score
#' @inheritParams cellTypesPerClusterBalloonPlot
#' @param norm_vivo Norm count matrix (n_genes X n_cells)
#' @param cluster_vivo cluster partition
#' @param gene_name Character name of the gene of interest
#' @param x_name Character name of x axis.
#' @return  ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#' @description Only the stages in \emph{cluster_vivo_factor} with an assignment greater than zero are showns in the balloon plot
#' @export plot_boxplot
#' @seealso \url{https://CRAN.R-project.org/package=gplots}


plot_boxplot <- function(norm_vivo, gene_name, cluster_vivo, order_label_vivo, title_name, fill_name, x_name, y_name){
  gene_expr <- norm_vivo[gene_name,]
  cluster_vivo <- factor(cluster_vivo,levels = order_label_vivo)
  data_common <- data.frame(cluster_vivo,gene_expr)
  ggplot2::ggplot(data_common, ggplot2::aes(x=cluster_vivo, y=gene_expr,fill= cluster_vivo)) +
    ggplot2::geom_boxplot() + ggplot2::geom_jitter(position=ggplot2::position_jitter(0.2)) + ggplot2::xlab(x_name) + ggplot2::ylab(y_name) + ggplot2::ggtitle(title_name) + ggplot2::labs(fill = fill_name)
}





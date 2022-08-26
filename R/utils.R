globalVariables(c("..count.."))
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

#' find_mean_cluster
#' @noRd
find_mean_cluster <- function (cluster_vivo, selected_stages, norm_vivo, genes, method = "mean")
{
  if (method == "mean"){
    if (length(selected_stages) > 0){
      number_cluster <- selected_stages
      median_es <- rep(list(0), length(number_cluster))
      for (i in 1:length(number_cluster)) {
        median_es[[i]] <- apply(norm_vivo[genes, cluster_vivo == number_cluster[i]], 1, "mean")
      }
      return(median_es)
    }
    if (length(selected_stages) == 0){
      number_cluster <- selected_stages
      median_es <- rep(list(0), 1)
      median_es[[1]] <- mean(norm_vivo[,cluster_vivo == number_cluster])
      return(median_es)
    }
  }
  if (method == "median"){
    if (length(selected_stages) > 0){
      number_cluster <- selected_stages
      median_es <- rep(list(0), length(number_cluster))
      for (i in 1:length(number_cluster)) {
        median_es[[i]] <- apply(norm_vivo[genes, cluster_vivo == number_cluster[i]], 1, "median")
      }
      return(median_es)
    }
    if (length(selected_stages) == 0){
      number_cluster <- selected_stages
      median_es <- rep(list(0), 1)
      median_es[[1]] <- mean(norm_vivo[,cluster_vivo == number_cluster])
      return(median_es)
    }
  }
  if (method != "mean" & method != "median"){
    stop("parameter method from function find_mean_cluster must be equal to mean or median")
  }
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


#' gg_color_hue
#' @noRd
gg_color_hue = function (n)
{
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' preprocess_mouse_published_data
#' @noRd

preprocess_mouse_published_data <- function(raw_count_mohammed, stages_mohammed, raw_count_deng, stages_deng, convert_name_id_mouse){


  raw_count_deng <- raw_count_deng[, (stages_deng != "fibroblast") & (stages_deng != "liver") & (stages_deng != "Zygote") & (stages_deng != "2_cell") & (stages_deng != "Mid_2_cell")]

  stages_deng <- stages_deng[(stages_deng !="fibroblast") & (stages_deng != "liver") & (stages_deng != "Zygote") & (stages_deng != "2_cell") & (stages_deng != "Mid_2_cell")]


  row.names(convert_name_id_mouse) <- convert_name_id_mouse[,1]

  common <- intersect(row.names(raw_count_mohammed),row.names(convert_name_id_mouse))


  name <- convert_name_id_mouse[common, 3]

  common_new <- common[Biobase::isUnique(name)]

  name_new <- convert_name_id_mouse[common_new, 3]

  raw_count_mohammed_new <- raw_count_mohammed[common_new, ]
  row.names(raw_count_mohammed_new) <- name_new


  common_mouse <- Reduce(intersect, list(row.names(raw_count_mohammed_new), row.names(raw_count_deng)))
  common_mouse <- common_mouse
  raw_count_mohammed_new_new <- raw_count_mohammed_new[common_mouse, ]
  raw_count_deng <- raw_count_deng[common_mouse, ]
  raw_counts_mouse_all <- cbind(raw_count_mohammed_new_new, raw_count_deng)
  row.names(raw_counts_mouse_all) <- common_mouse

  cluster_mouse <- c(stages_mohammed, stages_deng)
  return(list(raw_counts_mouse_all, cluster_mouse))

}





#' difference_significance_vitro_single
#' @noRd
#'
difference_significance_vitro_single <- function(norm_vitro, mean_profile_vitro, mean_profile_vivo, p_value, name_vitro, genes_all, cluster_vitro, selected_stages){
  norm_vitro_4 <- norm_vitro[genes_all,cluster_vitro==name_vitro]
  mean_profile_vitro_cluster = mean_profile_vitro[[which(levels(factor(cluster_vitro)) == name_vitro)]]
  norm_vitro_4_4 <- cbind(norm_vitro_4,mean_profile_vitro_cluster)
  dissimilarity_4_4 <- sqrt((1-cor((norm_vitro_4_4), method = "spearman"))/2)
  number_cluster <- selected_stages
  empirical_p_value <- rep(list(0),length(number_cluster))
  for (i in 1:length(number_cluster)){
    norm_vitro_4_2 <- cbind(norm_vitro_4,mean_profile_vivo[[i]])
    norm_vitro_4_2 <- cbind(norm_vitro_4,mean_profile_vivo[[i]])
    dissimilarity_4_2 <- sqrt((1-cor((norm_vitro_4_2), method = "spearman"))/2)
    mean_4_2 <- mean(dissimilarity_4_2[, length(colnames(dissimilarity_4_4))])
    limit <- as.numeric(quantile(dissimilarity_4_4[, length(colnames(dissimilarity_4_4))], probs=c(0, 1-p_value)))
    if (mean_4_2<limit[1] | mean_4_2>limit[2]){
      print("Difference statistically significant")
    }
    else{
      print("Difference not statistically significant")
    }
    empirical_p_value_one <- sum(dissimilarity_4_4[, length(colnames(dissimilarity_4_4))] > mean_4_2) / length(dissimilarity_4_4[, length(colnames(dissimilarity_4_4))])
    empirical_p_value[[i]]=empirical_p_value_one
  }
  return(empirical_p_value)
}



#' download_example_data
#' @noRd
download_example_data <- function(){
  current_wd <- getwd()
  url <- "https://hmgubox2.helmholtz-muenchen.de/index.php/s/EHQSnjMJxkR7QYT/download/SCOPRO.zip"
  destfile <- paste0(current_wd,"/SCOPRO.zip")
  download.file(url, destfile, quiet = FALSE)
  unzip(destfile, exdir = current_wd)
  setwd(paste0(current_wd, "/SCOPRO"))
  load(file = "seurat_genes_published_mouse.Rda")
  setwd(paste0(current_wd,"/SCOPRO"))
  load(file='mayra_dati_raw_0.Rda')
}




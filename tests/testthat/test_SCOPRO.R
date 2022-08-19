# Example with a norm count matrix of 10 genes and 10 samples
# All the genes have the same expression (10) value in all the samples.
norm_vivo = matrix(10, ncol = 10, nrow = 10)
row.names(norm_vivo) = paste0("Gene", seq(1:10))
colnames(norm_vivo) = paste0("Sample", seq(1:10))
selected_stages = c("Stage1")
markers_small = c()
cluster_vivo = c(rep("Stage1",5),rep("Stage2",5))

# select_top_markers
test_that("select_top_markers returns an error when no markers are provided as input", {
  expect_error(select_top_markers(selected_stages,cluster_vivo, norm_vivo, markers_small, max_number = 100, threshold = 0.1), "Length of markers_small is zero. Please provide a non-zero length vector.")
})

norm_vivo = matrix(10, ncol = 10, nrow = 10)
row.names(norm_vivo) = paste0("Gene", seq(1:10))
colnames(norm_vivo) = paste0("Sample", seq(1:10))
selected_stages = c("Stage4","Stage5")
markers_small = row.names(norm_vivo)[1:5]
names(markers_small) = rep("Stage1",length(markers_small))
cluster_vivo = c(rep("Stage1",5),rep("Stage2",5))

test_that("select_top_markers returns an error when one or more stages are not present ", {
  expect_error(select_top_markers(selected_stages,cluster_vivo, norm_vivo, markers_small, max_number = 100, threshold = 0.1), "One or more stages are not present. Please check that all stages in selected_stages are also present in cluster_vivo")
})

# SCOPRO

norm_vivo = matrix(10, ncol = 10, nrow = 10)
row.names(norm_vivo) = paste0("Gene", seq(1:10))
colnames(norm_vivo) = paste0("Sample", seq(1:10))
selected_stages = c("Stage1","Stage2")
name_vivo = c("Stage3")
marker_stages = list(row.names(norm_vivo)[1:5],row.names(norm_vivo)[6:10])
marker_stages_filter = c()

cluster_vivo = c(rep("Stage1",5),rep("Stage2",5))

norm_vitro = matrix(10, ncol = 10, nrow = 10)
row.names(norm_vitro) = paste0("Gene", seq(1:10))
colnames(norm_vitro) = paste0("Sample", seq(1:10))

cluster_vitro = c(rep("Cluster1",5),rep("Cluster2",5))


test_that("SCOPRO returns error if name_vivo is not one the stages present in the vector selected_stages",{expect_error(SCOPRO(norm_vitro, norm_vivo, cluster_vitro, cluster_vivo, name_vivo, marker_stages_filter, threshold = 0.1, number_link = 1, fold_change = 3, threshold_fold_change = 0.1,  marker_stages, selected_stages)
                                  ,"name_vivo must be one the stages present in the vector selected_stages")})





test_that("SCOPRO returns error if name_vivo is not one the stages present in the vector cluster_vivo",{expect_error(SCOPRO(norm_vitro, norm_vivo, cluster_vitro, cluster_vivo, name_vivo, marker_stages_filter, threshold = 0.1, number_link = 1, fold_change = 3, threshold_fold_change = 0.1,  marker_stages, selected_stages)
                                                                                                                        ,"name_vivo must be one the stages present in the vector selected_stages")})

test_that("SCOPRO returns error if name_vivo is not one the stages present in the vector cluster_vivo",{expect_error(SCOPRO(norm_vitro, norm_vivo, cluster_vitro, cluster_vivo, name_vivo, marker_stages_filter, threshold = 0.1, number_link = 1, fold_change = 3, threshold_fold_change = 0.1,  marker_stages, selected_stages)
                                                                                                                     ,"name_vivo must be one the stages present in the vector selected_stages")})

name_vivo = c("Stage2")

test_that("SCOPRO returns error if vector marker_stages_filter has 0 length",{expect_error(SCOPRO(norm_vitro, norm_vivo, cluster_vitro, cluster_vivo, name_vivo, marker_stages_filter, threshold = 0.1, number_link = 1, fold_change = 3, threshold_fold_change = 0.1,  marker_stages, selected_stages),"Vector marker_stages_filter has 0 length. Please provide a non-zero length vector.")
})











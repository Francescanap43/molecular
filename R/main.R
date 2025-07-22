run_full_analysis <- function() {
  task1_load_and_filter_data(gtf_path, matrix_zip_path, matrix_extract_dir)
  task2_plot_genes_over_3_UMIs(sce_object)
  task3_filter_unwanted_genes(counts_matrix, gtf)
  task4_run_pca_and_plot(sce_object)
  task5_run_umap_and_plot(sce_object, n_pcs = 6)
  task6_cluster_cells_and_plot(sce_object, k = 10)
  task7_annotate_cells_with_singleR(sce_object)
  task8_infer_tissue_origin(gtf_data, panglao_db)
}

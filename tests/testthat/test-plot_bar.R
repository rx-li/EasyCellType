test_that("plot_bar", {
  data("gene_pbmc")
  re <- easyct(gene_pbmc, db="cellmarker", species="Human", 
                    tissue=c("Blood", "Peripheral blood", "Blood vessel",
                             "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA",
                    scoretype="pos")
  expected <- plot_bar(test="GSEA", re)
  expect_s3_class(expected[[1]], "ggplot")
})
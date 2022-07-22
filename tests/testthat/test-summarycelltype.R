test_that("summarycelltype", {
  data("gene_pbmc")
  re <- easyct(gene_pbmc, db="cellmarker", species="Human", 
                     tissue=c("Blood", "Peripheral blood", "Blood vessel",
                              "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA",
                     scoretype="pos")
  expected <- summarycelltype(test="GSEA", re, cluster=0)
  expect_type(expected, "list")
})

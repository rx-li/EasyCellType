test_that("easyct", {
  data("gene_pbmc")
  expected <- easyct(gene_pbmc, db="cellmarker", species="Human", 
                    tissue=c("Blood", "Peripheral blood", "Blood vessel",
                    "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA",
                    scoretype="pos")
  expect_type(expected, "list")
  expect_s4_class(expected[[1]], "gseaResult")
})

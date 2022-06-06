#' Fisher exact test used in function 'easyct'
#'
#' @param testgenes A data frame containing query genes and the expression scores. 
#' @param ref The reference data base.
#' @param cols Column names of the input data frame 
#'
#' @return data frame
#' 
#' @import dplyr
#' 
#'
test_fisher <- function(testgenes, ref, cols){
  cell_n <- nrow(ref)
  testgene_n <- nrow(testgenes)
  
  # get the cell types 
  cellname <- unique(ref$celltype)
  
  p_value <- NULL
  score_m <- NULL
  for (i in seq(length(cellname))) {
    target_cell <- ref[ref$celltype == cellname[i], ]
    target_cell_n <- nrow(target_cell)
    
    common_id <- base::intersect(testgenes[, paste(cols[1])], target_cell$entrezid)
    common_id_n <- length(common_id)
    
    # conduct the fisher exact test only if the gene is enriched 
    if (common_id_n / testgene_n <= target_cell_n / cell_n){
      p <- 1
      s <- NA
    } else {
      d <- data_frame(
        "testgene" = c(length(common_id), testgene_n),
        "targetcell" = c(target_cell_n, cell_n)
      )
      res <- fisher.test(d, alternative = "greater")
      p <- as.numeric(res$p.value)
      
      testgenes_f <- testgenes %>%
        dplyr::filter(cols[1] %in% common_id)
      
      score_f <- testgenes_f[, paste(cols[3])]
      s <- mean(score_f, na.rm = TRUE)
    }
    score_m <- c(score_m, s)
    p_value <- c(p_value, p)
  }
  p_adjusted <- p.adjust(p_value, method = "BH")
  out <- data.frame("cellName" = cellname, "p_value" = p_adjusted, "score" = score_m)
  return(out)
}

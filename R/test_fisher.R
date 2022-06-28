#' Fisher exact test used in function 'easyct'
#' 
#' @description This function is used to conduct the modified Fisher's exact test.
#'
#' @param testgenes A data frame containing query genes and the expression scores. 
#' @param ref The reference data base.
#' @param cols Column names of the input data frame 
#' 
#' @importFrom dplyr filter 
#' @importFrom stats p.adjust fisher.test
#' @importFrom magrittr %>%
#' 
#' @return A data frame containg the results of fisher's exact test. 
#'
test_fisher <- function(testgenes, ref, cols){
  if(!is.data.frame(testgenes)){
    stop("Input data should be a data frame")
  }
  if(is.null(ref)){
    stop("Reference database should be specified")
  }
  if(length(cols) != ncol(testgenes)){
    stop("Dimension does not match; Check the input data frame and column names")
  }
    
  cell_n <- nrow(ref)
  testgene_n <- nrow(testgenes)
  
  # get the cell types 
  cellname <- unique(ref$celltype)
       
  calculator <- function(cellname.i){
    target_cell <- ref[ref$celltype == cellname.i, , drop=FALSE]
    target_cell_n <- nrow(target_cell)
    
    common_id <- intersect(testgenes[, cols[1]], target_cell$entrezid)
    common_id_n <- length(common_id)
    
    # conduct the fisher exact test only if the gene is enriched 
    if (common_id_n / testgene_n <= target_cell_n / cell_n){
      p <- 1
      s <- NA
    } else {
      d <- data.frame(
        "testgene" = c(length(common_id), testgene_n),
        "targetcell" = c(target_cell_n, cell_n)
      )
      res <- fisher.test(d, alternative="greater")
      p <- as.numeric(res$p.value)
      idx <- which(testgenes[, cols[1]] %in% common_id)
      testgenes_f <- testgenes[idx, ]
      score_f <- testgenes_f[, paste(cols[3])]
      s <- mean(score_f, na.rm=TRUE)
    }

    return(list(p=p, s=s))
  }
  
  re <- lapply(seq(length(cellname)), function(x) calculator(cellname[x]))
  p_value <- unlist(lapply(re, function(x) x$p))
  score_m <- unlist(lapply(re, function(x) x$s))
  p_adjusted <- p.adjust(p_value, method = "BH")
  
  out <- data.frame("cellName"=cellname, "p_value"=p_value, 
                    "score"=score_m, "p_adjust"=p_adjusted)
  return(out)
}

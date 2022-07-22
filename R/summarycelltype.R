#' Print test results
#'
#' @description This function is used to print summary table of annotation results
#'  for a specific cluster. 
#' 
#' @param test "GSEA" or "fisher".
#' @param results Annotation results.
#' @param cluster Cluster of interest.
#'
#'
#' @examples
#' data(gene_pbmc)
#' result <- easyct(gene_pbmc, db="cellmarker", species="Human", 
#' tissue=c("Blood", "Peripheral blood", "Blood vessel",
#' "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA", scoretype="pos")
#' summarycelltype(test="GSEA", results=result, cluster=0)
#' 
#' @return Test summary table
#' @export summarycelltype
#' 
summarycelltype <- function(test, results, cluster){
  stopifnot("The test should be specified as GSEA or fisher." = 
              test %in% c("GSEA", "fisher"))
  stopifnot("Annotation results should be specified" = length(results) > 0)
  stopifnot("Input cluster is not valid. Please check your annotation results." =
             mean(cluster %in% names(results)) == 1)
  
  idx <- which(names(results) %in% cluster)
  if(test == "GSEA"){  
    enrich.re.l <- lapply(results, function(x) x[seq(dim(x)[1]), drop=FALSE])
    enrich.re.l[idx]
  }else{
    results_sort <- lapply(results, function(x) x[order(-x$p_adjust, abs(x$score), 
                                                        decreasing = TRUE,
                                                        na.last = TRUE), ])
    results_sort[idx]
  }
   
}
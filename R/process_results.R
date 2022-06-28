#' Title Annotate cell types for single cell RNA data
#' 
#' @description This function is used to process the annotation test results. 
#' Processed data will be used to generate plots. 
#'
#' @param test Test used to annotation cell types: "GSEA" or "fisher"
#' @param data Annotation results.
#' @importFrom dplyr bind_rows slice_min group_by mutate 
#' @importFrom magrittr %>%
#'
#' @return A data frame used to generate plots. 
#'
process_results <- function(test, data){
  if(!test %in% c("GSEA", "fisher")){
    stop("The test should be specified as GSEA or fisher.")
  }
  if(is.null(data)){
    stop("Annotation results should be specified.")
  }
  if(test == "GSEA"){
    enrich.re.l <- lapply(data, function(x) x[seq(dim(x)[1]), drop=FALSE])
    enrich.re.d <- do.call(rbind, enrich.re.l)
    enrich.re.d$cluster <- rep(names(enrich.re.l), 
                               lapply(enrich.re.l, nrow))
    enrich.d <- enrich.re.d[, c("ID", "p.adjust", "cluster"), drop=FALSE]
    rownames(enrich.d) <- NULL
    enrich.hard <- enrich.d %>% 
      group_by(cluster) %>%
      slice_min(n = 1, order_by = p.adjust) %>%
      mutate(method = "hard_enrich")
    out <- merge(enrich.d, enrich.hard, by = c("ID", "p.adjust", "cluster"), all.x = TRUE)
    out$method[is.na(out$method)] <- "soft_enrich"
  }else if(test == "fisher"){
    fisher.hard <- lapply(data, 
                          function(x) x[order(-x$p_adjust, abs(x$score), decreasing = TRUE, 
                                              na.last = TRUE), ][1, ]) %>%
      bind_rows(.id = "cluster") %>%
      mutate(method = "hard_fisher") 
    
    fisher.soft <- lapply(data, 
                          function(x) x[order(-x$p_adjust, abs(x$score), decreasing = TRUE,
                                              na.last = TRUE), ][seq.int(5), ]) %>% 
      bind_rows(.id = "cluster")
    
    fisher.d <- merge(fisher.soft, fisher.hard, by = c("cluster", "cellName", "p_value", "score"), all.x = TRUE)
    fisher.d$method[is.na(fisher.d$method)] <- "soft_fisher"
    out <- data.frame(ID = fisher.d$cellName, p.adjust = fisher.d$p_value,
                             cluster = fisher.d$cluster, method = fisher.d$method)
  }
  return(out)
}

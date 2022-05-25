#' Title Annotate cell types for single cell RNA data
#'
#' @param test Test used to annotation cell types: "GSEA" or "fisher"
#' @param data Annotation results.
#'
#' @return data frame
#' @export process_results
#'

process_results <- function(test, data){
  if(test == "GSEA"){
    enrich.re.l <- lapply(data, function(x) x@result)
    enrich.re.d <- do.call(rbind, enrich.re.l)
    enrich.re.d$cluster <- rep(names(enrich.re.l), 
                               lapply(enrich.re.l, nrow))
    enrich.d <- enrich.re.d[, c("ID", "p.adjust", "cluster")]
    rownames(enrich.d) <- NULL
    enrich.hard <- enrich.d %>% 
      dplyr::group_by(cluster) %>%
      slice_min(n = 1, order_by = p.adjust) %>%
      dplyr::mutate(method = "hard_enrich")
    out <- merge(enrich.d, enrich.hard, by = c("ID", "p.adjust", "cluster"), all.x = TRUE)
    out$method[is.na(out$method)] <- "soft_enrich"
  }else if(test == "fisher"){
    fisher.hard <- lapply(data, 
                          function(x) x[order(-x$p_value, x$score, decreasing = TRUE, 
                                              na.last = TRUE), ][1, ]) %>%
      bind_rows(.id = "cluster") %>%
      dplyr::mutate(method = "hard_fisher") 
    
    fisher.soft <- lapply(data, 
                          function(x) x[order(-x$p_value, x$score, decreasing = TRUE,
                                              na.last = TRUE), ][1:5, ]) %>% 
      bind_rows(.id = "cluster")
    
    fisher.d <- merge(fisher.soft, fisher.hard, by = c("cluster", "cellName", "p_value", "score"), all.x = TRUE)
    fisher.d$method[is.na(fisher.d$method)] <- "soft_fisher"
    out <- data.frame(ID = fisher.d$cellName, p.adjust = fisher.d$p_value,
                             cluster = fisher.d$cluster, method = fisher.d$method)
  }
  return(out)
}

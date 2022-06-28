#' Create dot plot for annotation results
#' 
#' @description This function is used to generate a dor plot presenting the top 5
#' candidate cell types for each cluster. 
#'
#' @param test Test used to annotate cell types: "GSEA" or "fisher"
#' @param data Annotation results
#'
#' @importFrom ggplot2 ggplot geom_point theme aes element_blank element_line
#' @importFrom dplyr arrange slice group_by
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data(gene_pbmc)
#' result <- easyct(gene_pbmc, db="cellmarker", species="Human", 
#' tissue=c("Blood", "Peripheral blood", "Blood vessel",
#' "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA", scoretype="pos")
#' plot_dot("GSEA", result)
#' 
#' @return A dot plot showing the top 5 significant cell types for each cluster.
#' 
#' @export plot_dot
#'
plot_dot <- function(test="GSEA", data){
  if(!test %in% c("GSEA", "fisher")){
    stop("The test should be specified as GSEA or fisher.")
  }
  if(is.null(data)){
    stop("Annotation results should be specified.")
  }
  size <- c(7, 5, 3, 1.5, 1)
  if(test == "GSEA"){
    data.f <- process_results("GSEA", data)
    enrich.d <- data.f %>% group_by(cluster) %>% 
      arrange(p.adjust) %>%
      slice(seq.int(5))
    
    enrich.d$size <- unlist(
      lapply(as.numeric(table(enrich.d$cluster)), function(x) size[seq.int(x)])
    )
    ggplot(enrich.d) +
      geom_point(aes(x = cluster, y = ID), color = "#25B3B4", 
                 alpha = 0.8, size = enrich.d$size) +
      theme(axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.title.x = element_blank(),
            # explicitly set the horizontal lines (or they will disappear too)
            panel.grid.major.y = element_line(size=0.05, color="gray"))
  }else if(test == "fisher"){
    data.f <- process_results("fisher", data)
    data.f$size <- unlist(
      lapply(as.numeric(table(data.f$cluster)), function(x) size[seq.int(x)])
    )
    ggplot(data.f) +
      geom_point(aes(x=cluster, y=ID), color="#F38C36", alpha=0.8, 
                 size=data.f$size) +
      theme(axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.title.x = element_blank(),
            # explicitly set the horizontal lines (or they will disappear too)
            panel.grid.major.y = element_line(size=0.05, color="gray"))
  }
}

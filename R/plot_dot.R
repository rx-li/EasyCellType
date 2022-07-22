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
#' @importFrom rlang .data
#' 
#' @examples 
#' data(gene_pbmc)
#' result <- easyct(gene_pbmc, db="cellmarker", species="Human", 
#' tissue=c("Blood", "Peripheral blood", "Blood vessel",
#' "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA", scoretype="pos")
#' plot_dot("GSEA", result)
#' 
#' @return A dot plot showing the top 5 significant cell types for each cluster.
#' @export plot_dot
#'
plot_dot <- function(test="GSEA", data){
  stopifnot("The test should be specified as GSEA or fisher." =
            test %in% c("GSEA", "fisher"))
  stopifnot("Annotation results should be specified." = length(data) > 0)

  size_max <- 7
  
  if(test == "GSEA"){
    data.f <- process_results("GSEA", data)
    enrich.d <- data.f %>% group_by(.data$cluster) %>%
      arrange(.data$pvalue) %>%
      slice(seq.int(5)) %>%
      mutate(group_size_prop = log10(.data$pvalue) / log10(min(.data$pvalue))) 
    
    enrich.d$size <- size_max * enrich.d$group_size_prop
    enrich.d$size[enrich.d$size == "NaN" | enrich.d$size == 0] <- 1
    
    ggplot(enrich.d) +
      geom_point(aes(x = .data$cluster, y = .data$ID), color = "#25B3B4", 
                 alpha = 0.8, size = enrich.d$size) +
      theme(axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.title.x = element_blank(),
            # explicitly set the horizontal lines (or they will disappear too)
            panel.grid.major.y = element_line(size=0.05, color="gray"))
  }else if(test == "fisher"){
    data.f <- process_results("fisher", data) %>% 
      group_by(.data$cluster) %>%
      mutate(group_size_prop = log10(.data$pvalue) / log10(min(.data$pvalue)))
    
    data.f$size <- size_max * data.f$group_size_prop
    data.f$size[data.f$size == "NaN" | data.f$size == 0] <- 1
    
    ggplot(data.f) +
      geom_point(aes(x=.data$cluster, y=.data$ID), color="#F38C36", alpha=0.8, 
                 size=data.f$size) +
      theme(axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.title.x = element_blank(),
            # explicitly set the horizontal lines (or they will disappear too)
            panel.grid.major.y = element_line(size=0.05, color="gray"))
  }
}

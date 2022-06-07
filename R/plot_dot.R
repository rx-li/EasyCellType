#' Create dot plot for annotation results
#'
#' @param test Test used to annotate cell types: "GSEA" or "fisher"
#' @param data Annotation results
#'
#' @import ggplot2
#' 
#' @examples 
#' data(gene_pbmc)
#' result <- easyct(gene_pbmc, db="cellmarker", species="Human", 
#' tissue=c("Blood", "Peripheral blood", "Blood vessel",
#' "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA", scoretype="pos")
#' plot_dot("GSEA", result)
#' 
#' @return Plots
#' 
#' @export plot_dot
#' 
#'
plot_dot <- function(test="GSEA", data){
  size <- c(7, 5, 3, 1.5, 1)
  if(test == "GSEA"){
    data.f <- process_results("GSEA", data)
    enrich.d <- data.f %>% dplyr::group_by(cluster) %>% 
      dplyr::arrange(p.adjust) %>%
      dplyr::slice(1:5)
    
    enrich.d$size <- unlist(
      lapply(as.numeric(table(enrich.d$cluster)), function(x) size[1:x])
    )
    ggplot(enrich.d) +
      geom_point(aes(x = cluster, y = ID), color = "#25B3B4", 
                 alpha = 0.8, size = enrich.d$size) +
      theme(axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.x = element_blank(),
            #axis.text.x=element_text(angle=45, hjust=1),
            axis.title.x = element_blank(),
            # explicitly set the horizontal lines (or they will disappear too)
            panel.grid.major.y = element_line(size=0.05, color="gray"))
  }else if(test == "fisher"){
    data.f <- process_results("fisher", data)
    data.f$size <- unlist(
      lapply(as.numeric(table(data.f$cluster)), function(x) size[1:x])
    )
    ggplot(data.f) +
      geom_point(aes(x=cluster, y=ID), color="#F38C36", alpha=0.8, 
                 size=data.f$size) +
      theme(axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.x = element_blank(),
            #axis.text.x=element_text(angle=45, hjust=1),
            axis.title.x = element_blank(),
            # explicitly set the horizontal lines (or they will disappear too)
            panel.grid.major.y = element_line(size=0.05, color="gray"))
  }
}

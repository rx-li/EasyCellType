#' Create bar plots for each cluster
#'
#' @param test "GSEA" or "fisher"
#' @param data Annotation results 
#'
#' 
#' @import ggplot2
#' @import forcats
#' 
#' @export plot_bar 
#' 
#' @example plot_bar(test="GSEA", data=data)
#'
#' 
plot_bar <- function(test="GSEA", data){
  plot_one_cluster <- function(d, t){

    if (nrow(d) > 10){
      d <- d[1:10, ]
    }
    
    if (test == "GSEA"){
      ID_color <- ifelse(d$method == "hard_enrich", "#009E73", "#38BE90")
    }else if(test == "fisher"){
      ID_color <- ifelse(d$method == "hard_fisher", "#F9A215", "#FEBA4F")
    }
    
    d %>%
      # order for x axis
      mutate(ID = fct_reorder(ID, length)) %>%
      ggplot(aes(x=ID, y=length)) + 
      geom_bar(stat="identity", width=0.5, fill=ID_color) +
      theme_classic() +
      coord_flip() + guides(x="none") +
      labs(title=paste("Cluster", t)) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position="none",
            title=element_text(size=9))

  }
  data.f <- process_results(test, data)
  data.l <- split(data.f, data.f$cluster)
  data.l_sort <- lapply(data.l, function(x) x[order(x$method, x$p.adjust),])
  data.l_sort_bar <- lapply(data.l_sort, function(x) {x$length <- c(nrow(x):1); x})
  lapply(seq(length(data.l_sort_bar)), function(x) 
    plot_one_cluster(data.l_sort_bar[[x]], names(data.l_sort_bar)[x]))
  
}

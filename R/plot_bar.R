#' Create bar plots for each cluster
#' 
#'@description This function is used to generate set of bar plots presenting up 
#'to 10 candidate cell types for each cluster.
#'
#' @param test "GSEA" or "fisher"
#' @param data Annotation results 
#'
#' @importFrom ggplot2 ggplot geom_bar theme_classic coord_flip labs theme aes element_blank element_text
#' @importFrom forcats fct_reorder
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' 
#' @examples 
#' data(gene_pbmc)
#' result <- easyct(gene_pbmc, db="cellmarker", species="Human", 
#' tissue=c("Blood", "Peripheral blood", "Blood vessel",
#' "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA", scoretype="pos")
#' plot_bar("GSEA", result)
#' 
#' @return Bar plots showing show up to 10 candidate cell types for each cluster. 
#' 
#' @export plot_bar 
#' 
plot_bar <- function(test="GSEA", data){
  
  # check the parameters
  if(!test %in% c("GSEA", "fisher")){
    stop("The test should be specified as GSEA or fisher.")
  }
  if(is.null(data)){
    stop("Data should be annotation results")
  }
  
  plot_one_cluster <- function(d, t){

    if (nrow(d) > 10){
      d <- d[seq.int(10), ,drop=FALSE]
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
      coord_flip() +
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
  data.l_sort <- lapply(data.l, function(x) x[order(x$method, x$p.adjust), , drop=FALSE])
  data.l_sort_bar <- lapply(data.l_sort, function(x) {x$length <- c(nrow(x):1); x})
  lapply(seq(length(data.l_sort_bar)), function(x) 
    plot_one_cluster(data.l_sort_bar[[x, drop=FALSE]], names(data.l_sort_bar)[x]))
}

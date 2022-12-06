#' Create bar plots for each cluster
#' 
#'@description This function is used to generate set of bar plots presenting up 
#'to 10 candidate cell types for each cluster.
#'
#' @param test "GSEA" or "fisher"
#' @param data Annotation results 
#' @param cluster Cluster can be specified to print plots.  
#'
#' @importFrom ggplot2 ggplot geom_bar theme_classic coord_flip labs theme aes element_blank element_text scale_y_discrete
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
#' @export plot_bar 
#' 
plot_bar <- function(test="GSEA", data, cluster=NULL){
  
  # check the parameters
  stopifnot("The test should be specified as GSEA or fisher." = 
            test %in% c("GSEA", "fisher"))
  stopifnot("Annotation results should be specified." =
            length(data) > 0)
  
  length_max <- 5
  data.f <- process_results(test, data) %>% 
    group_by(.data$cluster) %>%
    mutate(group_length_prop = log10(.data$pvalue) / log10(min(.data$pvalue)))
  data.f$length <- length_max * data.f$group_length_prop
  data.f$length[data.f$length == "NaN" | data.f$length == 0] <- 1
  
  data.l <- split(data.f, data.f$cluster, drop=FALSE)
  
  # check parameter 'cluster'
  stopifnot("Input cluster is not valid. Please check your annotation results." =
            length(cluster) == 0 | mean(cluster %in% names(data.l)) == 1)

  plot_one_cluster <- function(d, t){

    if (nrow(d) > 10){
      d <- d[seq.int(10), ,drop=FALSE]
    }
    
    if (test == "GSEA"){
      ID_color <- ifelse(d$method == "hard_enrich", "#009E73", "#38BE90")
    }else if(test == "fisher"){
      ID_color <- ifelse(d$method == "hard_fisher", "#F9A215", "#FEBA4F")
    }
    
    lim_p <- round(-log10(rev((d$pvalue))), 2)
    
    d %>%
      # order for x axis
      mutate(ID = fct_reorder(.data$ID, .data$method, .desc = TRUE)) %>%
      ggplot(aes(x=.data$ID, y=.data$length)) + 
      geom_bar(stat="identity", width=0.5, fill=ID_color) +
      theme_classic() +
      coord_flip() +
      labs(title=paste("Cluster", t), y="-Log10(P value)") +
      scale_y_discrete(limits=factor(lim_p)) + 
      theme(#axis.title.x=element_blank(),
            #axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position="none",
            title=element_text(size=9))

  }
  
  data.l_sort <- lapply(data.l, function(x) x[order(x$method, x$pvalue), , drop=FALSE])
  plots <- lapply(seq(length(data.l_sort)), function(x) 
    plot_one_cluster(data.l_sort[[x, drop=FALSE]], names(data.l_sort)[x]))
  names(plots) <- names(data.l_sort)
  if(length(cluster) == 0){
    plots
  } else{
    plots[which(names(plots) %in% cluster)]
  }
}

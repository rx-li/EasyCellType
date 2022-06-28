#' Annotate cell types for scRNA-seq data.
#' 
#' @description This function is used to run the annotation analysis using either 
#' GSEA or a modified Fisher's exact test. We expect users to input a data frame 
#' containing expressed markers, cluster information and the differential score (log fold change). 
#' The gene lists in that data frame should be sorted by their differential score.
#'  
#' @param data Name of a data frame containing the markers' Entrez IDs, cluster, 
#' and expression scores; Marker genes should be sorted in each cluster.
#' Order of the columns should be gene, cluster and score. An example data can be
#' loaded using `data(gene_pbmc)`.
#' @param db Name of the reference database: cellmarker, clustermole or panglaodb;
#' @param species Human or Mouse. Human in default.
#' @param tissue Tissue types can be specified when running the analysis. Length of tissue can be larger than 1. 
#' The possible tissues can be seen using `data(cellmarker_tissue)`, `data(clustermole_tissue)` and `data(panglao_tissue)`.
#' @param p_cut Cutoff of the P value for GSEA.
#' @param test "GSEA" or "fisher"; "GSEA" is used in default. 
#' @param scoretype Argument used for GSEA. Default value is "std". If all scores are positive, then scoretype should be "pos". 
#'
#' @importFrom magrittr %>%
#' @importFrom clusterProfiler GSEA
#' @importFrom dplyr filter select 
#' 
#' @examples 
#' data(gene_pbmc)
#' result <- easyct(gene_pbmc, db="cellmarker", species="Human", 
#' tissue=c("Blood", "Peripheral blood", "Blood vessel",
#' "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA", scoretype="pos")
#' 
#' @return A list containing the test results for each cluster. 
#' 
#' @export easyct
#' 
#'
easyct <- function(data, db="cellmarker",
                  species="Human", tissue=NULL, p_cut=0.5, test="GSEA", scoretype = "std"){
  # check the parameters
  if(!is.data.frame(data)){
    stop("Input should be a data frame")
  }
  if(!db %in% c("cellmarker", "panglao", "clustermole")){
    stop("The database should be chosen from cellmarker, panglao and clustermole.")
  }
  if(!species %in% c("Human", "Mouse")){
    stop("The speices should be either Human or Mouse.")
  }
  if(!test %in% c("GSEA", "fisher")){
    stop("The test should be specified as GSEA or fisher.")
  }
  if(!scoretype %in% c("std", "pos", "neg")){
    stop("The score type should be chosen from std, pos and neg.")
  }
  if(!is.numeric(p_cut)){
    stop("Cutoff of p value should be a number")
  }
  
  # load the databases
  rda <- paste(db, "_db", sep="")
  db.data <- get(rda)
  # extract the species
  db.data.human <- db.data %>% filter(spe == species)
  
  if(length(tissue) > 0 && !tissue %in% unique(db.data.human$organ)){
    stop("Please refer the tissue types using data(cellmarker_tissue), 
         data(clustermole_tissue) or data(panglao_tissue)")
  }

  # data frame to list 
  cols <- names(data)
  data.l <- split(data, f=data[, paste(cols[2])], drop=FALSE)

  # extract the gene id and cell name
  if (length(tissue) > 0) {
    cells <- db.data.human %>%
      filter(organ %in% tissue) %>%
      select(celltype, entrezid)
  } else {
    cells <- db.data.human %>%
      select(celltype, entrezid)
  }
  
  # analysis 
  if (test == "GSEA") {
    input.l <- lapply(data.l, function(x) x[, paste(cols[3])])
    
    input.l.named <- lapply(seq(length(data.l)), function(x) 
      {names(input.l[[x]]) <- data.l[[x]][, paste(cols[1])]; input.l[x][[1]]})
    names(input.l.named) <- names(input.l)

    results <- lapply(
      input.l.named,
      function(x) GSEA(x, TERM2GENE = cells, pvalueCutoff = p_cut, minGSSize = 1, scoreType = scoretype)
      )
    
    results.f <- lapply(seq(length(data.l)), function(i)
      {if (nrow(results[[i]]) == 0){results[[i]]<- NA}; results[[i]]}
    )
    names(results.f) <- names(results)
    
    results.c <- results.f[!is.na(results.f)]
    out <- results.c
    
  } else if(test == "fisher"){
    out <- lapply(data.l, function(x) test_fisher(x, cells, cols))
  }
  
  return(out)
}


#' Annotate cell types for single cell RNA data
#'
#' @param data Name of a data frame containing the markers' Entrez IDs, cluster, and expression scores; Marker genes should be sorted in each cluster.
#'             Order of the columns should be gene, cluster and score.
#' @param db Name of the reference database; CellMarker database is used in default.
#' @param species Human or Mouse. Human in default.
#' @param tissue Tissue types can be specified when running the analysis. Length of tissue can be larger than 1. 
#' @param p_cut Cutoff of the P value for GSEA.
#' @param test "GSEA" or "fisher"; "GSEA" is used in default. 
#' @param scoretype Argument used for GSEA. Default value is "std". If all scores are positive, then scoretype should be "pos". 
#'
#' 
#' @import clusterProfiler
#' @import dplyr
#' 
#' @examples 
#' data(gene_pbmc)
#' result <- annot(gene_pbmc, db="cellmarker", species="Human", 
#' tissue=c("Blood", "Peripheral blood", "Blood vessel",
#' "Umbilical cord blood", "Venous blood"), p_cut=0.3, test="GSEA", scoretype="pos")
#' 
#' @export annot
#' 
#'
annot <- function(data, db="cellmarker",
                  species="Human", tissue=NULL, p_cut=0.5, test="GSEA", scoretype = "std"){
  # loadRData <- function(fileName){
  #   #loads an RData file, and returns it
  #   load(fileName)
  #   get(ls()[ls() != "fileName"])
  # }
  
  # load the database
  rda <- paste(db, "_db", sep="")
  db.data <- get(rda)
  
  # data frame to list 
  cols <- names(data)
  #classes <- as.character(unlist(unique(data[, paste(cols[2])])))
  data.l <- split(data, f=data[, paste(cols[2])])

  # extract the species
  db.data.human <- db.data %>% dplyr::filter(spe == species)
  
  # extract the gene id and cell name
  if (length(tissue) > 0) {
    cells <- db.data.human %>%
      dplyr::filter(organ %in% tissue) %>%
      dplyr::select(celltype, entrezid)
  } else {
    cells <- db.data.human %>%
      dplyr::select(celltype, entrezid)
  }
  
  # analysis 
  if (test == "GSEA") {
    input.l <- lapply(data.l, function(x) x[, paste(cols[3])])
    
    for (i in seq(length(data.l))) {
      names(input.l[[i]]) <- data.l[[i]][, paste(cols[1])]
      }
    
    results <- lapply(
      input.l,
      function(x) GSEA(x, TERM2GENE = cells, pvalueCutoff = p_cut, minGSSize = 1, scoreType = scoretype)
      )

    for (i in seq(length(data.l))) {
      if (nrow(results[[i]]) == 0){
        results[[i]]<- NA
      }
    }
    results.c <- results[!is.na(results)]
    out <- results.c
    
  } else if(test == "fisher"){
    out <- lapply(data.l, function(x) test_fisher(x, cells, cols))
  }
  
  return(out)
}


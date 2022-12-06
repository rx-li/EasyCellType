#' Title Summarize markers contirbuting to the cell type annotation
#'
#' @param test Test used to annotation cell types: "GSEA" or "fisher"
#' @param data Annotation results.
#' @param species "Human" or "Mouse"
#' 
#' @importFrom dplyr group_by distinct select
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db org.Mm.eg.db
#'
#' @return A data frame containing genes contributed to cell annotation
#' @export
#'
#' @examples ## core_markers <- coremarkers("GSEA", data)
#' 
coremarkers <- function(test, data, species){
  # check the parameters
  stopifnot("The test should be specified as GSEA or fisher." = 
              test %in% c("GSEA", "fisher"))
  stopifnot("Annotation results should be specified." =
              length(data) > 0)
  stopifnot("Species should be specified as Human or Mouse." =
              species %in% c("Human", "Mouse"))
  
  data.f <- process_results(test, data)
  data.f <- as.data.frame(data.f)
  genes_entrezid <- lapply(data.f$core_enrichment, function(x) unlist(strsplit(x, "/")))
  
  if(species == "Human"){
    genes_symbol <- lapply(genes_entrezid, function(x)
      mapIds(org.Hs.eg.db,
             keys=x, #Column containing Ensembl gene ids
             column="SYMBOL",
             keytype="ENTREZID",
             multiVals="first")
    )
  } else if(species == "Mouse"){
    genes_symbol <- lapply(genes_entrezid, function(x)
      mapIds(org.Mm.eg.db,
             keys=x, #Column containing Ensembl gene ids
             column="SYMBOL",
             keytype="ENTREZID",
             multiVals="first")
    )
  }

  genes_symbol_str <- lapply(genes_symbol, function(x) paste(x, collapse=","))
  data.f$genes <- unlist(genes_symbol_str)
  out <- data.f %>% 
    group_by(cluster) %>%
    distinct(ID, .keep_all=TRUE) %>%
    select(cluster, ID, genes) 
  out
  return(out)
}



#' Title Convert gene symbol to Entrez ID
#'
#' @description This function is used to convert the gene symbol to Entrez Id. Used in 
#' easyct function.
#' 
#' @param d A data frame where first column contains gene symbols. 
#' @param species "Human" or "Mouse".
#' 
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom AnnotationDbi mapIds
#'
#' @return A data frame containing gene symbols and the corresponding Entrez ID
#'
mapsymbol <- function(d, species){
  stopifnot("Input should be a data frame" = is.data.frame(d))
  stopifnot("The speices should be either Human or Mouse." = 
              species %in% c("Human", "Mouse"))
  
  if(species == "Human"){
    d$entrezid <- mapIds(org.Hs.eg.db,
                         keys=d[, 1], #Column containing Ensembl gene ids
                         column="ENTREZID",
                         keytype="SYMBOL",
                         multiVals="first")
  }else if(species == "Mouse"){
    d$entrezid <- mapIds(org.Mm.eg.db,
                         keys=d[, 1], #Column containing Ensembl gene ids
                         column="ENTREZID",
                         keytype="SYMBOL",
                         multiVals="first")
  }
  d_f <- d[, c(4, 2, 3), drop=FALSE] 
  return(d_f)
}

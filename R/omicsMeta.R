############# omics for metagenome ############
#' split the metaphlan2 dataset
#'
#' @param data  , merge metaphlan2 dataset
#' @param prefix , output name
#'
#' @return list, include all rank data
#' @export
#'
#' @examples
splitMetaphlan2 <- function(data, prefix){

  # this function use to split the metaphlan2 data
  index <- sapply(strsplit(rownames(data), perl = T, split = "\\|"), length)
  maxrank <- max(index)

  datalist <- list()
  # define the rank
  classify <- c("kingdom", "phylum", "class", "order",
                "family", "genus", "species", "strain")
  # output the profile
  for(i in 1:maxrank){
    index2 <- which(index == i)
    data[index2,] -> outdat
    rname <- sapply(strsplit(rownames(outdat), perl = T, split = "\\|"),
                    function(x){y <- x[length(x)]; return(y)})
    rownames(outdat) <- rname
    datalist[[i]] <- outdat
  }
  names(datalist) <- paste0(prefix, "_", classify[1:maxrank])
  return(datalist)

}

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

#' getMetaphlan2rank
#' get the phylum-species rank
#' @param name , the merge metaphlan2 rownames
#'
#' @return
#' @export
#'
#' @examples
getMetahlan2rank <- function(name){

  # get the metaphlan2 rank
  a <- name[grep("s__", rownames(anme))]
  b <- a[-grep("t__", a)]
  c <- unlist(sapply(b, function(x) {strsplit(x, split = "\\|")}))
  d <- as.vector(c)
  e <- matrix(tmp[,1], ncol = 7, byrow = T)

  return(e)

}

############ GEE analysis ########################

#' myGeeAnalysis
#'
#' @param dataset
#' @param metadata
#' @param transformM
#' @param confounder
#' @param timevar
#' @param scale
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
myGeeAnalysis <- function(dataset, metadata, confounder, timevar,
                          scale, IDvar, ...){

  # ready the data
  id <- intersect(rownames(dataset), rownames(metadata))
  dataset <- dataset[id, ]
  metadata <- metadata[id, ]
  print("confirm the sample ID is order by time")
  metadata <- metadata[order(metadata[,timevar]), ]
  dataset <- dataset[rownames(metadata), ]

  confounderindex <- which(colnames(metadata) %in% c(confounder, timevar, IDvar))

  if((length(confounderindex)-2) != length(confounder)){
    stop("please check the confounder variable")
  }

  datacon <- metadata[, confounderindex]
  metadatafilter <- metadata[, -confounderindex]
  # analysis

  result <- matrix(NA, nrow = ncol(metadatafilter), ncol = ncol(dataset)*2)
  result <- as.data.frame(result)
  rownames(result) <- colnames(metadatafilter)

  for(i in 1:c(ncol(metadatafilter))){
    #print(i)
    for(j in 1:ncol(dataset)){
      #print(j)
      dat_com <- data.frame(x = metadatafilter[,i],
                            y = dataset[,j], datacon,
                            PatientID=metadata[,IDvar])

      #confirm the data type

      formula <- formula(paste0("y~x+", paste(confounder, collapse = "+")))
      dat_com <- dat_com[!apply(dat_com, 1, function(x){any(is.na(x))}), ]

      if(scale){
        dat_com$x <- scale(dat_com$x)
        dat_com$y <- scale(dat_com$y)
      }
      # GEE
      geeInd <- geem(formula, id=PatientID, data=dat_com,
                     family=gaussian, corstr="independence")
      tmp <- summary(geeInd)
      # result
      result[i, c((2*j-1):(2*j))] <- c(tmp$beta[2], tmp$p[2])
      colnames(result)[c((2*j-1):(2*j))] <- paste0(colnames(dataset)[j],
                                            c("_estimate", "_p.value"))

    }
  }
  return(result)

}


############ Partial correlation #################

#' myPcor
#'
#' @param dataset
#' @param metadata
#' @param confounder
#'
#' @return
#' @export
#'
#' @examples
myPcor <- function(dataset, metadata, confounder){

  library(ppcor)
  # ready the data
  id <- intersect(rownames(dataset), rownames(metadata))
  dataset <- dataset[id, ]
  metadata <- metadata[id, ]

  confounderindex <- which(colnames(metadata) %in% confounder)
  if(length(confounderindex) != length(confounder)){
    stop("please check the confounder variable")
  }

  datacon <- metadata[, confounderindex]
  metadatafilter <- metadata[, -confounderindex]
  # analysis

  result <- matrix(NA, nrow = ncol(metadatafilter), ncol = ncol(dataset)*2)
  result <- as.data.frame(result)
  rownames(result) <- colnames(metadatafilter)

  for(i in 1:c(ncol(metadatafilter))){

    for(j in 1:ncol(dataset)){

      dat_com <- data.frame(x = metadatafilter[,i],
                            y = dataset[,j], datacon)

      # confirm the data type
      #dat_com <- dat_com[!is.na(dat_com$x) & !is.na(dat_com$y), ]
      dat_com <- dat_com[!apply(dat_com, 1, function(x){any(is.na(x))}), ]
      # pcc
      # check the number is a constant
      op1 <- sum(dat_com$x==dat_com$x[1])==length(dat_com$x)
      op2 <- sum(dat_com$y==dat_com$y[1])==length(dat_com$y)
      if(op1 | op2){
      result[i, c((2*j-1):(2*j))] <- c(0, 1)
      }else{
      tmp <- pcor.test(dat_com$x, dat_com$y, dat_com[,confounder], method="s")
      # result
      result[i, c((2*j-1):(2*j))] <- c(tmp$estimate, tmp$p.value)
      }
      colnames(result)[c((2*j-1):(2*j))] <- paste0(colnames(dataset)[j],
                                                   c("_estimate", "_p.value"))

    }
  }
  return(result)
}







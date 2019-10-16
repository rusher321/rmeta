#' transform a vector to string
#'
#' @param x
#' @param sep
#'
#' @return string
#' @export
#'
#' @examples
pasteP <- function(x, sep){
  # x is vector
  # sep is a  character string to separate the terms
  # join the vector
  len <- length(x)
  y <- c()
  for(i in 1:len){
    y <- paste(y, x[i], sep = sep)
  }
  return(y)
}

#' merge two dataframe ,who have diffrent colname
#'
#' @param dat1
#' @param dat2
#' @param by
#'
#' @return
#' @export
#'
#' @examples
mergeP <- function(dat1, dat2, by){
  coldat1 <- colnames(dat1)
  coldat2 <- colnames(dat2)
  interfeature <- intersect(coldat1, coldat2)
  if(length(interfeature)==1){
    out <- merge(dat1, dat2, by=by, all=T)
  }else{
    diffdat1 <- setdiff(coldat1, interfeature)
    diffdat2 <- setdiff(coldat2, interfeature)
    tmp1 <- dat1[, c(by, diffdat1)]
    tmp2 <- dat2[, c(by, diffdat2)]
    tmp3 <- merge(tmp1, tmp2, by = by, all=T)
    tmpC <- rbind(dat1[, interfeature], dat2[, interfeature])
    out <- cbind(tmpC, tmp3[, -which(colnames(tmp3)==by)])
  }
  return(out)
}

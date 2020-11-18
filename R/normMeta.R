#' cum Norm stat method form ANCOM_BC
#'
#' @param data
#' @param pFlag
#' @param rel
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cumNormStatFast2 <- function (data, pFlag = FALSE, rel = 0.1, ...){

  mat <- data
  smat = lapply(1:ncol(mat), function(i) {
    sort(mat[which(mat[, i] > 0), i], decreasing = TRUE)
  })
  leng = max(sapply(smat, length))
  if (any(sapply(smat, length) == 1))
    stop("Warning sample with one or zero features")
  smat2 = array(NA, dim = c(leng, ncol(mat)))
  for (i in 1:ncol(mat)) {
    smat2[leng:(leng - length(smat[[i]]) + 1), i] = smat[[i]]
  }
  rmat2 = sapply(1:ncol(smat2), function(i) {
    quantile(smat2[, i], p = seq(0, 1, length.out = nrow(smat2)),
             na.rm = TRUE)
  })
  smat2[is.na(smat2)] = 0
  ref1 = rowMeans(smat2)
  ncols = ncol(rmat2)
  diffr = sapply(1:ncols, function(i) {
    ref1 - rmat2[, i]
  })
  diffr1 = matrixStats::rowMedians(abs(diffr))
  if (pFlag == TRUE) {
    plot(abs(diff(diffr1))/diffr1[-1], type = "h",
         ...)
    abline(h = rel)
    axis(1, at = seq(0, length(diffr1), length.out = 5),
         labels = seq(0, 1, length.out = 5))
  }
  x = which(abs(diff(diffr1))/diffr1[-1] > rel)[1]/length(diffr1)
  if (x <= 0.5 | is.na(x)) {
    message("Default value being used.")
    x = 0.5
  }
  return(x)

}


#' Cumulative sum scaling
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
cssMeta <- function(data){

  library(metagenomeSeq)
  library(matrixStats)

  p <- cumNormStatFast2(data)
  x <- data
  xx = x
  xx[x == 0] <- NA
  qs = colQuantiles(xx, probs = p, na.rm = TRUE)
  normFactors <- sapply(1:ncol(xx), function(i) {
    xx = (x[, i] - .Machine$double.eps)
    sum(xx[xx <= qs[i]])
  })
  #names(normFactors) <- colnames(x)
  return(t(t(data)/normFactors))

}


#' upper quartile normalization
#'
#' @param data
#' @param p
#'
#' @return
#' @export
#'
#' @examples
#'
uqMeta <- function (data, p = 0.75){

  UQ <- function(x) {
    quantile(x[x > 0], p)
  }
  uq <- unlist(apply(data, 2, UQ))

  norm_factor <- uq/median(uq)

  return(t(t(data)/norm_factor))

}

#' the effective library size using UQ (upper quartile)
#'
#' @param data
#' @param p
#'
#' @return
#' @export
#'
#' @examples
elibuqMeta <- function (data, p = 0.75){

  UQ <- function(x) {
    quantile(x[x > 0], p)
  }
  uq <- unlist(apply(data, 2, UQ))
  norm_factor <- uq/median(uq)
  ncolsum <- apply(data, 2, sum)
  norm_factor2 <- norm_factor*ncolsum

  return(t(t(data)/norm_factor2))

}


#' median ratio normalization
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#'
medraioMeta <- function(data){

  presudo <- apply(data, 1 , function(x){prod(x)^(1/length(x))})

  ratio <- apply(data, 2, function(x){x/presudo})

  median_ratio <- apply(ratio, 2, median, na.rm=T)

  normlized <- apply(data, 1, function(x){x/median_ratio})

  return(t(normlized))

}




#' Trimmed Mean of M values
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
tmmMeta <- function(data){

 # import edgeR , don't understant the detail
 library(edgeR)

  normfactor <- calcNormFactors(data, method = "TMM")

  normlized <- t(t(data)/normfactor)

  return(normlized)

}

#' Elib-Trimmed Mean of M values
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
elibtmmMeta <- function(data){

  # import edgeR , don't understant the detail
  library(edgeR)

  normfactor <- calcNormFactors(data, method = "TMM")
  ncolsum <- apply(data, 2, sum)
  normfactor2 <- normfactor*ncolsum
  normlized <- t(t(data)/normfactor2)

  return(normlized)

}



#' Total-Sum Scaling
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#'
tssMeta <- function(data){

  apply(data, 2, sum) -> lib.size
  apply(data, 1, function(x){x/lib.size}) -> normlized

  return(t(normlized))

}




#' Based-rank Inverse Normal Transformation
#' The detail inf is (here)[https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html]
#'
#' @param u, vector
#' @param k, 0<= k <=0.5
#'
#' @return
#' out, vector transformed
#' @export
#'
#' @examples
#' x <- c(0.1, 0,3, 0.5 , NA)
#' RankINT(x)
RankINT <- function (u, k = 0.375)
{
  if (!is.vector(u)) {
    stop("A numeric vector is expected for u.")
  }
  if ((k < 0) || (k > 0.5)) {
    stop("Select the offset within the interval (0,0.5).")
  }
  n <- sum(!is.na(u))
  #r <- rank(u)
  out <- qnorm((rank(r, na.last = "keep") - k)/(n - 2 * k + 1))
  return(out)
}









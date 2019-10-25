##### overview for data ##########

#' JSD distance
#' to compute the jsd distance
#' @param inMatrix
#' @param pseudocount
#' @param ...
#'
#' @return distance
#' @export
#'
#' @examples
distJSD <- function(inMatrix, pseudocount=0.00000001, ...) {

  # if have negtive number, transform them
  if(any(inMatrix < 0)){
    mixV <- min(inMatrix)
    inMatrix <- inMatrix+abs(mixV)+pseudocount
  }

  inMatrix <- t(inMatrix)
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  min <- min(inMatrix)
  inMatrix <- inMatrix+min
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) {
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)

}



#' kBest
#' to get the best k in non-supervision method
#' @param data
#' @param dist
#' @param method
#'
#' @return list , include the best cluster number k & original figure & dataframe
#' @export
#'
#' @examples
kBest <- function(data, dist , method = "kmeans"){

  nclusters=NULL
  sil = NULL
  out <- list()
  res <- matrix(NA, 19, ncol(data))

  for (k in 2:20) {
    #print(k)
    switch (method,
            kmeans = { data.cluster_temp <-kmeans(dist, k)$cluster},
            pam = { data.cluster_temp <-pam(dist, k)$clustering},
            fanny = {  data.cluster_temp <- fanny(dist, k)$clustering}
    )
    res[k-1,] <- data.cluster_temp
    nclusters[k-1] <- index.G1(t(data) , data.cluster_temp,  d = dist,
                               centrotypes = "medoids")
    sil[k-1] <- mean(silhouette(data.cluster_temp, dist = dist)[,3])
  }

  best <- which.max(nclusters)+1
  kCluster <- c(2:20)
  CH_index <- nclusters
  Silhouette <- sil
  cluster <- data.frame(kCluster,  CH_index, Silhouette)
  cluster <- melt(cluster, id = "kCluster")
  colnames(cluster) <- c("kCluster", "Index", "value")

  figure <- ggplot(cluster, aes(x=value, y=kCluster))+
    geom_segment(aes(yend=kCluster),xend=0,colour="grey")+
    geom_point(size=3,aes(colour=Index))+
    scale_colour_brewer(palette="Set1",limits=c("CH_index","Silhouette"))+
    theme_bw()+finalTheme+xlab("")+ylab("Number of cluster")+facet_grid(.~Index, scales = "free")

  out <- list(res[best-1,], best, figure)

  return(out)

}


#' kBest2
#'
#' @param data
#' @param dist
#' @param method
#' @param clustNum
#'
#' @return
#' @export
#'
#' @examples
kBest2 <- function(data, dist , method = "kmeans",  clustNum){

  nclusters=NULL
  sil = NULL
  out <- list()
  res <- matrix(NA, clustNum, ncol(data))

  for (k in 1:clustNum) {

    switch (method,
            kmeans = { data.cluster_temp <-kmeans(dist, k)$cluster},
            pam = { data.cluster_temp <- pam(dist, k)$clustering},
            fanny = {  data.cluster_temp <- fanny(dist, k)$clustering}
    )

    res[k,] <- data.cluster_temp
    nclusters[k] <- index.G1(t(data) , data.cluster_temp,  d = dist,
                             centrotypes = "medoids")

  }

  best <- which.max(nclusters)
  print("Best cluster number is ", best)
  CH_index <- nclusters
  df <- data.frame(clusters = as.factor(1:clustNum), y = CH_index)
  df[1,2] <- min(df[,2], na.rm=T)-0.1*min(df[,2], na.rm=T)
  ylab <- "Calinski-Harabasz index"

  p <- ggpubr::ggline(df, x = "clusters", y = "y", group = 1,
                      color = "steelblue", ylab = ylab,
                      xlab = "Number of clusters k",
                      main = "Optimal number of clusters\nCH-index method"
  )

  p <- p + geom_vline(xintercept = which.max(CH_index), linetype=2, color = "steelblue")

  return(p)

}


#' multiAdonis
#' for pair compare for multigroup use the adonis
#' @param dat
#' @param config
#' @param perm
#' @param method
#'
#' @return data.frame
#' @export
#'
#' @examples
multiAdonis <- function(dat, config, perm=999, method="bray"){
  # reorder the sample id of dat and config
  # dat : row is sample
  # perm : permutation number ,default 999
  # method : distance method

  id <- intersect(rownames(dat), rownames(config))
  dat <- dat[id, ]
  config <- config[id,]
  combng <- combn(unique(config), 2)
  res <- matrix(NA, nrow=ncol(combng), ncol = 6)
  for(i in 1:ncol(combng)){
    group <- c(config[which(config==combng[1,i])],
               config[which(config==combng[2,i])])
    pro <- dat[c(which(config==combng[1,i]),
                 which(config==combng[2,i])), ]
    # here use the adonise function
    library(vegan)
    res[i, 1:6] <- as.numeric(adonis(pro~group, permutations = perm, method = method)$aov[1,])
  }
  rownames(res) <- apply(combng, 2, function(x){paste0(x[1], " vs ", x[2])})
  colnames(res) <- c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2", "Pr(>F)")
  return(res)

}


############# outlier ################

#' checkOutlier
#' output the outlier value on a vector by quantile
#' @param x
#' @param coef
#' @param sname
#' @param plot
#' @param vname
#'
#' @return
#' @export
#'
#' @examples
checkOutlier <- function(x, coef = 2, sname ,plot = F, vname = "test"){

  # whether is a conituous data
  if(is.numeric(x) & length(levels(as.factor(x)))>=5){
  # get the outlier index
    x[is.na(x)] <- median(x, na.rm = T)
    quantiles <- quantile(x, probs=c(0.25,0.75), na.rm = T)
    IQR <- quantiles[2]-quantiles[1]
    index <- x < (quantiles[1]-coef*IQR)|x > (quantiles[2]+coef*IQR)
    if(any(index)){
    # print the outlier
    name <- paste0(sname[index], "_", round(x[index], 2))
    res <-  paste(name, collapse = " ")

    #return(res)
   # print the outlier sample on boxplot
    if(plot){
       dat <- data.frame(value = x)
       dat$label <- ifelse(index, sname, "")
       dat$var <- rep(vname, nrow(dat))
       ggplot(dat,aes(var,value))+
         geom_boxplot()+geom_text(aes(label=label),hjust=-0.1)
      }
      return(res)
    }else{
      return(NULL)
    }
  }else{
    stop("the value is not conituous data\n")
  }

}



#' mutivarOutlier
#' on a dataframe to select the outlier
#' @param dat
#' @param scale
#' @param k
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mutivarOutlier <- function(dat, scale=F, k, ...){

  if(scalue){
    dat <- scale(dat)
  }

  sam_num <- nrow(dat)
  row_name <- rownames(dat)

  # library(DDoutlier)
  outlier.scores <- LOF(dat, k = 10 )
  checkOutlier(outlier.scores, ...)


}




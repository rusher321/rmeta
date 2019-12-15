########## Test from mazhanshan <ma@vandals.uidaho.edu> article ###########
    ##### renhuahui@genomics.cn##############

### DAR (diversity–area relationship): Extending classic SAR
### (species–area relationship) for biodiversity and biogeography
### analyses





########### Test from Benoit Lehallier ################
# Undulating changes in human plasma proteome profiles across the lifespan
#' clusterTrajector
#'
#' @param dataset  metabolism
#' @param phemeta  clinical paremeter eg, age
#' @param span  loess paremeter
#'
#' @return
#' @export
#'
#' @examples
clusterTrajector <- function(dataset, phemeta, span){

  # rm the NA
  phemeta <- phemeta[!is.na(phemeta[,1]), ,drop=F]
  id <- intersect(rownames(dataset), rownames(phemeta))
  dataset2 <- dataset[id, ]
  value <- phemeta[id,]
  # to scale the log10 value
  datasetScale <- apply(dataset2, 2, function(x){scale(log10(x+1))})

  # loess estimate
  datasetLoess <- apply(datasetScale, 2,
                        function(x){predict(loess(x~value, span = span))})

  # cluster
  # library(NbClust)
  # here paremete can to modify
  hc <- hclust(dist(t(datasetLoess), method = "euclidean"),
               method = "ward.D2")
  # get the best clust
  cluster <- NbClust(t(datasetLoess), distance = "euclidean", min.nc=2, max.nc=8,
          method = "complete", index = "ch")
  bestC <- as.numeric(cluster$Best.nc[1])
  clustlable <- cutree(hc, bestC)
  hashtmp <- hash(names(clustlable), as.numeric(clustlable))
  # plot
  datasetLoess <- as.data.frame(datasetLoess, check.names=F)
  datasetLoess$var <- sort(value)
  qdat <- melt(datasetLoess, id.vars = "var")
  qdat$group <- sapply(qdat$variable, function(x){values(hashtmp[x])})

  figure <- ggplot(qdat, mapping = aes(x = var, y = value, group=variable))+
      geom_line()+fontTheme+finalTheme+facet_wrap(.~group)

  return(figure)
}

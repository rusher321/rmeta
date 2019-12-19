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
  datasetScale2 <- apply(datasetScale, 2, function(x){tapply(x, value,
                            mean, na.rm=T)})
  # loess estimate
  value <- as.numeric(rownames(datasetScale2))
  datasetLoess <- apply(datasetScale2, 2,
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


#' twopartEta
#' to compute the partial eta or relative risk
#' @param dat
#' @param phe
#' @param number
#'
#' @return
#' @export
#'
#' @examples
twopartEta <- function(dat, phe, number=20){

  # match the sample ID
  id <- intersect(rownames(dat), rownames(phe))
  if(length(id)==0){
    stop("can't match the sample id")
  }
  dat <- dat[id, ]
  phe <- phe[id, ]

  # output
  varnum <- ncol(phe)
  res <- matrix(NA, nrow = ncol(dat), ncol = 6*varnum+2)
  res <- as.data.frame(res)

  colnames(res)[1:(2*varnum)] <- paste0("all_", c(paste0("effect_", colnames(phe)),
                                                  paste0("pvlaue_", colnames(phe))))

  colnames(res)[(2*varnum+1):(4*varnum)] <- paste0("logistic_", c(paste0("effect_", colnames(phe)),
                                                paste0("pvlaue_", colnames(phe))))

  colnames(res)[(4*varnum+1):(6*varnum)] <- paste0("quantify_", c(paste0("effect_", colnames(phe)),
                                                paste0("pvlaue_", colnames(phe))))

  colnames(res)[(6*varnum+1):(6*varnum+2)] <- c("zero number", "no_zero number")
  formula <- as.formula(paste0("x~", paste(colnames(phe), collapse = "+")))

  for(i in 1:ncol(dat)){

    x <- dat[,i]
    # all part
    alldat <- data.frame( x= log(x+1), phe)
    model1 <- lm(formula, data = alldat)
    tmppvalue <- as.numeric(summary(model1)$coefficients[-1,4])

    library(sjstats)

    alleta <- eta_sq(model = model1, partial = T)[,2] # CI

    res[i, 1:(2*varnum)] <- c(alleta, tmppvalue)

    # logistic part
    part1dat <- data.frame( x= x, phe)
    part1dat$x <- with(part1dat, ifelse(x>0, 1, 0))

    if(sum(part1dat$x==0) < number |sum(part1dat$x==1) < number){
      res[i, (2*varnum+1):(4*varnum)] <- c(rep(0, varnum), rep(1, varnum))
    }else{
      model2 <- glm(formula, family = binomial, data = part1dat)
      tmppvalue <- as.numeric(summary(model2)$coefficients[-1,4])
      part1or <- odds_to_rr(model2)[-1,2]
      res[i, (2*varnum+1):(4*varnum)] <- c(part1or, tmppvalue)
    }
    # quantify part
    part2dat <- alldat[alldat$x!=0, ,drop=F]
    part2dat$x <- log(part2dat$x)
    model3  <-  lm(formula, data = part2dat)
    tmppvalue <- as.numeric(summary(model3)$coefficients[-1,4])
    part3eta <- eta_sq(model = model3, partial = T)[,2] # CI

    res[i, (4*varnum+1):(6*varnum)] <- c(part3eta, tmppvalue)

    res[i,(6*varnum+1):(6*varnum+2)]  <- c(sum(part1dat$x==0), sum(part1dat$x==1))

  }
  rownames(res) <- colnames(dat)

  return(res)

}

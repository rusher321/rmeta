
########compare between group#############

#' matchpair ID
#' to match the sample ID for pair wilcox test
#' @param configdat ,dataframe
#' @param ID , sample ID
#' @param Time , time point
#'
#' @return dataframe
#' @export
#'
#' @examples
matchpairID <- function(configdat, ID, Time){

 names(table(configdat[,ID]))[table(configdat[, ID])==2] -> matchname
 configdat[!is.na(match(configdat[, ID], matchname)),] -> outconfig
 # sort by ID
 outconfig[order(outconfig[,ID]), ] -> matchdat
 matchdat[order(matchdat[,Time]), ] -> matchdat

 return(matchdat)

}


#' mywilcox_2t
#' wilcox pair test for 2 time point
#' @param datamatrix ,
#' @param configdata , from function matchpairID
#' @param Time2 , time point
#' @param ratio , use NA or 0 to compute the ratio
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mywilcox_2t <- function(datamatrix, configdata, Time2, ratio="zero",...){

  # to sort the data
  config <- matchID(configdata, ...)
  data <- datamatrix[rownames(config), ,drop=F]
  # to analysis

  out <- matrix(NA, nrow = ncol(data), ncol = 9)
  config <- as.factor(as.character(config[, Time2]))
  nlevels = levels(droplevels(config))

  for(i in 1:ncol(data)){

    tmp <- as.numeric(data[,i])
    g1 <- tmp[config == nlevels[1]]
    g2 <- tmp[config == nlevels[2]]

    wilcox_sign <- pvalue(wilcoxsign_test(g1~g2))
    effect <- wilcoxonPairedR(x <- tmp, g <- config)

    if(ratio=="zero"){
      or <- tapply(tmp, config, function(x){sum(x!=0, na.rm=T)/length(x)})
    }else{
      or <- tapply(tmp, config, function(x){sum(!is.na(x))/length(x)})
    }
    mean_abun <- tapply(tmp, config, mean, na.rm=T)
    median_abun <-  tapply(tmp, config,  median, na.rm=T)
    z_score <- statistic(wilcoxsign_test(g1~g2))
    out[i, 1:9] <- c(wilcox_sign, or, mean_abun, median_abun, z_score, effect)
  }

  # out
  rownames(out) <- colnames(data)
  colnames(out) <- c("sign_p.value",paste0(rep(nlevels,3),
                    rep(c("_ratio", "_mean", "_median"), each=2)),
                    "z_score", "effect_size")
  out <- as.data.frame(out)

  out$p.adjust <- p.adjust(out$sign_p.value, method = "BH")
  out$enrich <- ifelse(out$p.adjust<0.05, ifelse(out$z_score>0,
                      nlevels[1], nlevels[2]), "none")
  return(out)
}


#' comparePair
#' wilcox pair test for multi time point
#' @param data
#' @param config
#' @param group
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
comparePair <- function(data, config, group, ...){

  # wilcox test on multigroup
  nlevle <- levels(config[, group])
  if(length(nlevle) >2 ){
    combng <- combn(nlevle, 2)
    out <- list()
    for(i in 1:ncol(combng)){
     # print(i)
      subconfig <- config[config[,group] %in% combng[,i], ]
      out[[i]] <- mywilcox_2t(data, subconfig, ...)

    }
    return(out)
  }else{

    res <- mywilcox_2t(data, config, ...)
    return(out)
  }


}

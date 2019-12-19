
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
  config <- matchpairID(configdata, ...)
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

    out <- mywilcox_2t(data, config, ...)
    return(out)
  }


}



#' kwmeta
#'
#' @param pr , dataframe
#' @param config , dataframe
#'
#' @return
#' @export
#'
#' @examples
kwmeta <- function(pr, config){

  # set the alpha
  alpha <- 0.05  # 0.05
  # match names of files
  inter <- intersect(rownames(config), rownames(pr))
  pr <- t(pr[inter, ])
  config <-config[inter,,drop=F]

  sum_name <- sum(ifelse(rownames(config)==colnames(pr),1,0))
  if(sum_name==nrow(config)){print ("All sample matched")}
  print(paste0("the sample number is ",length(inter)))

  num <- nrow(pr)
  fr <- as.factor(config[, 1])
  group <- levels(fr)
  print(group)
  #output
  len <- length(group)
  num2 <- len*len/2+3.5*len + 2
  out <- matrix(NA, num, num2)

  ktp <- apply(pr, 1, function(x){kruskal.test(x ~ fr)$p.value})
  #post hoc dunn test
  library(PMCMR)
  for (i in 1:num) {
    pr1 <- as.numeric(pr[i,])
    index <- is.na(pr1)
    pr1 <- pr1[!index]
    fr1 <- fr[!index]
    rk  <- rank(pr1)
    res <- c(ktp[i], tapply(pr1, fr1, median),tapply(pr1,fr1,mean), tapply(pr1>0, fr1,mean))
    dtp <- posthoc.kruskal.dunn.test(pr1, fr1, p.adjust.method = "BH")$p.value
    dtp <- cbind(dtp, rep(NA, len - 1))
    dtp <- rbind(rep(NA, len), dtp)
    dtp[upper.tri(dtp)] <- t(dtp)[upper.tri(dtp)]
    rownames(dtp)[1] <- colnames(dtp)[1]
    colnames(dtp)[len] <- rownames(dtp)[len]

    mean_rank <- tapply(rank(pr1),fr1,mean)
    res <- c(res,dtp[lower.tri(dtp)], mean_rank)
    #
    conclude <- rep(0,2*len-1)
    or <- order(mean_rank)
    conclude[2*(1:len)-1] <- group[or]
    op <- rep(1,len-1)
    for(j in 1:(len-1)){op[j] <- dtp[or[j],or[j+1]]}
    symbol <- rep("=",len-1)
    symbol[!is.na(op) & op <= alpha] <- "<"
    symbol[is.na(op) | op == 1] <- "<=>"
    for(x in 1:(len-1)){
      if(symbol[x]=="<"){
        p_tmp <- c()
        for(y in 1:x){
          for(z in (x+1):len){
            p_tmp <- c(p_tmp,dtp[or[y],or[z]])
          }
        }
        if(any(p_tmp>0.05)){symbol[x] <- "="}
      }
    }

    conclude[(1:(len - 1)) * 2] <- symbol
    res <- c(res, paste(conclude, collapse = " "))
    if(length(res)==ncol(out)){out[i, ] <- res}else{print (res)}
  }
  rownames(out) <- rownames(pr)
  cn <- c("kw.p", paste("median", group[1:len], sep = "_"))
  cn <- c(cn,paste("mean", group[1:len], sep = "_"))
  cn <- c(cn, paste("or", group[1:len], sep = "_"))
  cn <- c(cn, paste0("p", "_", group[row(dtp)[lower.tri(dtp)]], "_", group[col(dtp)[lower.tri(dtp)]]))
  cn <- c(cn, paste("mean_rank",group[1:len], sep = "_"))
  cn <- c(cn, "nearby")
  colnames(out) <- cn

  return(as.data.frame(out, check.names=F))
}




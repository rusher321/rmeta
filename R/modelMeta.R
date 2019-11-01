########### model for metagenome ############

#' smote.exs2
#' this function from package DMwR smote.exs
#' @param data
#' @param tgt
#' @param N
#' @param k
#'
#' @return
#' @export
#'
#' @examples
smote.exs2 <- function(data,tgt,N,k)
  # INPUTS:
  # data are the rare cases (the minority "class" cases)
  # tgt is the name of the target variable
  # N is the percentage of over-sampling to carry out;
  # and k is the number of nearest neighours to use for the generation
  # OUTPUTS:
  # The result of the function is a (N/100)*T set of generated
  # examples with rare values on the target
{
  nomatr <- c()
  T <- matrix(nrow=dim(data)[1],ncol=dim(data)[2]-1)
  for(col in seq.int(dim(T)[2]))
    if (class(data[,col]) %in% c('factor','character')) {
      T[,col] <- as.integer(data[,col])
      nomatr <- c(nomatr,col)
    } else T[,col] <- data[,col]

  if (N < 100) { # only a percentage of the T cases will be SMOTEd
    nT <- NROW(T)
    idx <- sample(1:nT,as.integer((N/100)*nT))
    T <- T[idx,]
    N <- 100
  }

  p <- dim(T)[2]
  nT <- dim(T)[1]

  ranges <- apply(T,2,max)-apply(T,2,min)

  nexs <-  as.integer(N/100) # this is the number of artificial exs generated
  # for each member of T
  new <- matrix(nrow=nexs*nT,ncol=p)    # the new cases

  for(i in 1:nT) {

    # the k NNs of case T[i,]
    xd <- scale(T,T[i,],ranges)
    for(a in nomatr) xd[,a] <- xd[,a]==0
    dd <- drop(xd^2 %*% rep(1, ncol(xd)))
    kNNs <- order(dd)[2:(k+1)]

    for(n in 1:nexs) {
      # select randomly one of the k NNs
      neig <- sample(1:k,1)

      ex <- vector(length=ncol(T))

      # the attribute values of the generated case
      difs <- T[kNNs[neig],]-T[i,]
      new[(i-1)*nexs+n,] <- T[i,]+runif(1)*difs
      for(a in nomatr)
        new[(i-1)*nexs+n,a] <- c(T[kNNs[neig],a],T[i,a])[1+round(runif(1),0)]

    }
  }
  newCases <- data.frame(new)
  for(a in nomatr)
    newCases[,a] <- factor(newCases[,a],levels=1:nlevels(data[,a]),labels=levels(data[,a]))

  newCases[,tgt] <- factor(rep(data[1,tgt],nrow(newCases)),levels=levels(data[,tgt]))
  colnames(newCases) <- colnames(data)
  newCases
}

#' rmetaSMOTE
#'
#' @param form
#' @param data
#' @param perc.over
#' @param k
#' @param learner
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
rmetaSMOTE <- function(form,data,
                  perc.over=200,k=5,
                  #perc.under=200,
                  learner=NULL,...
){
  # INPUTS:
  # form a model formula
  # data the original training set (with the unbalanced distribution)
  # minCl  the minority class label
  # per.over/100 is the number of new cases (smoted cases) generated
  #              for each rare case. If perc.over < 100 a single case
  #              is generated uniquely for a randomly selected perc.over
  #              of the rare cases
  # k is the number of neighbours to consider as the pool from where
  #   the new examples are generated
  # perc.under/100 is the number of "normal" cases that are randomly
  #                selected for each smoted case
  # learner the learning system to use.
  # ...  any learning parameters to pass to learner
  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  minCl <- levels(data[,tgt])[which.min(table(data[,tgt]))]

  # get the cases of the minority class
  minExs <- which(data[,tgt] == minCl)

  # generate synthetic cases from these minExs
  if (tgt < ncol(data)) {
    cols <- 1:ncol(data)
    cols[c(tgt,ncol(data))] <- cols[c(ncol(data),tgt)]
    data <-  data[,cols]
  }
  newExs <- smote.exs2(data[minExs,],ncol(data),perc.over,k)
  if (tgt < ncol(data)) {
    newExs <- newExs[,cols]
    data <- data[,cols]
  }

  # get the undersample of the "majority class" examples
  # selMaj <- sample((1:NROW(data))[-minExs],
  #                 as.integer((perc.under/100)*nrow(newExs)),
  #                 replace=T)

  selMaj <- (1:NROW(data))[-minExs]

  # the final data set (the undersample+the rare cases+the smoted exs)
  newdataset <- rbind(data[selMaj,],data[minExs,],newExs)

  # learn a model if required
  if (is.null(learner)) return(newdataset)
  else do.call(learner,list(form,newdataset,...))
}

#' smoteMatch
#'
#' @param data
#' @param group
#'
#' @return
#' @export
#'
#' @examples
smoteMatch <- function(data, group){
  # match the sample id
  #library(DMwR)
  id <- intersect(rownames(data), rownames(group))
  dataX <- data[id, ]
  dataY <- group[id, , drop=F]
  # to get the SMOTE function parameter perc.over , perc.under
  less <- min(table(dataY[,1]))
  over <- max(table(dataY[,1]))


  perc.over <- ((over/less)-1)*100
  #perc.under <- over*100*100/(less*perc.over)

  # generate the new data
  dataX$group <- dataY[,1]

  newData <- rmetaSMOTE(group ~ ., dataX, perc.over = perc.over)

  return(newData)

}


#### from the fujingyuan #####
#### two part explain model #######

#' twopartModel
#' From The Gut Microbiome Contributes to a Substantial Proportion of the Variation in Blood Lipids
#' @param dat , microbiome data row is sample id, col is variable
#' @param phe , metadata row is sample id ,col is variable
#' @param response , the variable to explain
#' @param cutoff , ensure to detect or undetect
#' @param number , when sample number is limited , default 10
#'
#' @return dataframe
#' @export
#'
#' @examples
twopartModel <- function(dat , phe, response, cutoff, number=10){
  # match the sample ID
  id <- intersect(rownames(dat), rownames(phe))
  if(length(id)==0){
    stop("can't match the sample id")
  }
  dat <- dat[id, ]
  y <- phe[id, response]
  # part1 transform the matrix to 0-1
  dat2 <- dat
  dat2[dat2 >= cutoff] <- 1
  dat2[dat2 < cutoff ] <- 0

  # filter all zero variable

  out <- matrix(NA, nrow = ncol(dat), ncol = 3+4+4+2+2)
  for(i in 1:ncol(dat)){
    #print(i)
    x <- dat2[,i]
    out[i,1:3] <- c(sum(x==0), sum(x==1), median(dat[, i]))
    if(sum(x==0) < number |sum(x==1) < number){   # here set the cutoff 20, maybe need the sample size to adjust
      out[i,4:7] <- c(0, 0 , 0, 1)
    }else{
      res <- glm(y~x)
      out[i,4:7] <- as.numeric(summary(res)$coefficients[2,])
    }

  }

  # part2 transform the matrix log
  dat3 <- dat
  dat3[dat3 < cutoff] <- 0
  dat3 <- apply(dat3, 2, log10)
  for(i in 1:ncol(dat)){
   # print(i)
    x <- dat3[,i]
    id2 <- !is.infinite(x)
    x2 <- x[id2]
    y2 <- y[id2]
    if(length(x2) < number){   # here set the cutoff 20, maybe need the sample size to adjust
      out[i,8:11] <- c(0, 0 , 0, 1)
    }else{
      res <- glm(y2~x2)
      out[i,8:11] <- as.numeric(summary(res)$coefficients[2,])
    }

    # meta-analysis
    tmp_meta <- data.frame(beta = out[i, c(4,8)], se = out[i, c(5,9)])
    tmp_meta$se <- ifelse(tmp_meta$se == 0, 0.001, tmp_meta$se) # need to discuss for if se==0, can't get meta result
    res_meta <- rma(yi = beta, data = tmp_meta, sei = se, method = "DL")
     out[i, 12:13] <- c(res_meta$zval, res_meta$pval)

    # association value
    pvalue <- c(out[i, c(7,11,13)])
    estimate <- c(out[i, c(4,8,12)])

    out[i ,15] <- min(pvalue)
    zscore <- qnorm(1-(min(pvalue)/2))
    out[i, 14] <- ifelse(estimate[which.min(pvalue)] >0, zscore, -zscore)

  }

  # meta-analysis


  rownames(out) <- colnames(dat)
  colnames(out) <- c("No.absent", "No.Present", "medianAbundance",
                paste0("binary", c("_estimate", "_se", "_tvalue", "_p")),
                paste0("quantitative", c("_estimate", "_se", "_tvalue", "_p")),
                paste0("Meta", c("_zvalue", "_p")),
                paste0("Asso", c("_zvalue", "_p")))

  return(out)

  }



#' r2Twopartmodel
#'
#' From The Gut Microbiome Contributes to a Substantial Proportion of the Variation in Blood Lipids
#' @param dat , microbiome data row is sample id, col is variable
#' @param phe , metadata row is sample id ,col is variable
#' @param response , the variable to expain
#' @param cutoff , ensure  to detect or undetect
#' @param number , when sample number is limited , default 10
#' @param cutoffp , significant feature
#' @param repeatN , repeat number
#'
#' @return vector
#' @export
#'
#' @examples
r2Towpartmodel <- function(dat , phe, response, cutoff, number=10, cutoffp=0.01, repeatN=100, confunder=F){
    # match the sample ID
    id <- intersect(rownames(dat), rownames(phe))
    if(length(id)==0){
      stop("can't match the sample id")
    }
    dat <- dat[id, ]
    phe <- phe[id, ]

    R2 <- c()

    for(m in 1:repeatN){
     print(m)
    # split the data 80/20 percent , 80% is discovery data ,20% is validation data
    sampNum <- nrow(dat)
    discSampindex <- sample(1:sampNum, size = round(sampNum*0.8), replace = NULL)
    valiSampindex <- c(1:sampNum)[-discSampindex]

    discDatax <- dat[discSampindex, ]
    discPhe <- phe[discSampindex, ]
    valiDatax <- dat[valiSampindex, ]
    valiPhe <- phe[valiSampindex, response]

    # twopart Model to get the effect size
    twopartRes <- twopartModel(dat = discDatax, phe = discPhe, response = response,
                               cutoff = cutoff, number = number)

    featurelist <- rownames(twopartRes)[twopartRes[,15] <= cutoffp]
    if(length(featurelist)==0){
      R2[m] <- 0
      stop(paste0("pvalue cutoff :", cutoffp, " can't get any feature"))
    }else{
      print(paste0("pvalue cutoff :", cutoffp, " get the significant feature ",
                   length(featurelist)))
    }
    # get the additive modele
    twopartRes2 <- twopartRes[featurelist, ]
    valiDatax2 <- valiDatax[, featurelist]
    risk <- c()
    for(i in 1:nrow(valiDatax2)){
      riskvalue <- c()
      for(j in 1:ncol(valiDatax2)){
        beta1 <- twopartRes2[, 4]
        beta2 <- twopartRes2[, 8]
        b <- ifelse(valiDatax2[i,j] < cutoff, 0, 1)
        q <- ifelse(valiDatax2[i, j] < cutoff, 0, log10(valiDatax2[i,j]))
        riskvalue[j] <- beta1+b+beta2*q  # can't not ecsure why add b
      }
      risk[i] <- sum(riskvalue)
    }
    # get the R square
    if(!is.null(confunder)){
      # risk
      tmp <- phe[valiSampindex, c(response, confunder)]
      tmp$risk <- risk
      formula <- as.formula(paste0(response, "~."))
      lmmode <- summary(lm(formula, data = tmp))

      # no risk
      tmp2 <- tmp[,-ncol(tmp)]
      formula2 <- as.formula(paste0(response, "~."))
      lmmode2 <- summary(lm(formula2, data = tmp2))

      R2[m] <- lmmode$adj.r.squared-lmmode2$adj.r.squared

    }else{
      lmmode <- summary(lm(valiPhe~risk))
      R2[m] <- lmmode$adj.r.squared
    }
    }
    return(R2)

  }







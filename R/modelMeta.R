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
r2Twopartmodel <- function(dat , phe, response, cutoff, number=10, cutoffp=0.01, repeatN=100, confounder=NULL){
    # match the sample ID
    id <- intersect(rownames(dat), rownames(phe))
    if(length(id)==0){
      stop("can't match the sample id")
    }
    dat <- dat[id, ]
    phe <- phe[id, ]

    R2 <- c()

    for(m in 1:repeatN){
    # print(m)
    # split the data 80/20 percent , 80% is discovery data ,20% is validation data
    sampNum <- nrow(dat)
    discSampindex <- sample(1:sampNum, size = round(sampNum*0.8), replace = F)
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
    twopartRes2 <- twopartRes[featurelist, ,drop=F]
    valiDatax2 <- valiDatax[, featurelist,drop=F]
    risk <- c()
    for(i in 1:nrow(valiDatax2)){
      riskvalue <- c()
      for(j in 1:ncol(valiDatax2)){
        beta1 <- twopartRes2[, 4]
        beta2 <- twopartRes2[, 8]
        b <- ifelse(valiDatax2[i,j] <= cutoff, 0, 1)
        q <- ifelse(valiDatax2[i, j] <= cutoff, 0, log10(valiDatax2[i,j]))
        riskvalue[j] <- beta1+b+beta2*q  # can't not ecsure why add b
      }
      risk[i] <- sum(riskvalue)
    }
    # get the R square
    if(!is.null(confounder)){
      # risk
      tmp <- phe[valiSampindex, c(response, confounder)]
      tmp$risk <- risk
      formula <- as.formula(paste0(response, "~."))
      lmmode <- summary(lm(formula, data = tmp))
      #
      # no risk
      tmp2 <- tmp[,-ncol(tmp)]
      formula2 <- as.formula(paste0(response, "~."))
      lmmode2 <- summary(lm(formula2, data = tmp2))

      #　from the rsq.partial in the rsq package
      R2[m] <- 1-((1-lmmode$r.squared)/(1-lmmode2$r.squared))*(lmmode2$df[2]/lmmode$df[2])

      #R2[m] <- lmmode$adj.r.squared-lmmode2$adj.r.squared

    }else{
      lmmode <- summary(lm(valiPhe~risk))
      R2[m] <- lmmode$adj.r.squared
    }
    }
    return(R2)

  }


#' r2Twopartmodelcv
#' add the cross validation process
#' From The Gut Microbiome Contributes to a Substantial Proportion of the Variation in Blood Lipids
#' @param dat , microbiome data row is sample id, col is variable
#' @param phe , metadata row is sample id ,col is variable
#' @param response , the variable to expain
#' @param cutoff , ensure  to detect or undetect
#' @param number , when sample number is limited , default 10
#' @param cutoffp , significant feature
#' @param repeatN , repeat number
#' @param fold , cross validation fold
#' @return vector
#' @export
#'
#' @examples
r2Twopartmodelcv <- function(dat , phe, response, cutoff, number=10, cutoffp=0.01, repeatN=100, fold,confounder=NULL){
  # match the sample ID
  id <- intersect(rownames(dat), rownames(phe))
  if(length(id)==0){
    stop("can't match the sample id")
  }
  dat <- dat[id, ]
  phe <- phe[id, ]

  R2 <- c()

  for(m in 1:repeatN){
    # print(m)
    # cross vlidation
    sampNum <- nrow(dat)
    foldlist <- createFolds(1:sampNum, k = fold)
    risk <- rep(NA, sampNum)

    for(n in 1:fold){

      discSampindex <- c(1:sampNum)[-foldlist[[n]]]
      valiSampindex <- foldlist[[n]]

      discDatax <- dat[discSampindex, ]
      discPhe <- phe[discSampindex, ]
      valiDatax <- dat[valiSampindex, ]
      valiPhe <- phe[valiSampindex, response]

    # twopart Model to get the effect size
      twopartRes <- twopartModel(dat = discDatax, phe = discPhe, response = response,
                               cutoff = cutoff, number = number)

      featurelist <- rownames(twopartRes)[twopartRes[,15] <= cutoffp]
    if(length(featurelist)==0){
      risk[foldlist[[n]]] <- 0  # if no feature, the rish set zero
      next
      print(paste0("pvalue cutoff :", cutoffp, " can't get any feature"))
    }else{
      print(paste0("pvalue cutoff :", cutoffp, " get the significant feature ",
                   length(featurelist)))
    }
    # get the additive modele
    twopartRes2 <- twopartRes[featurelist, ,drop=F]
    valiDatax2  <- valiDatax[, featurelist, drop=F]

    for(i in 1:nrow(valiDatax2)){
      riskvalue <- c()
      for(j in 1:ncol(valiDatax2)){
        beta1 <- twopartRes2[, 4]
        beta2 <- twopartRes2[, 8]
        b <- ifelse(valiDatax2[i, j] <= cutoff, 0, 1)
        q <- ifelse(valiDatax2[i, j] <= cutoff, 0, log10(valiDatax2[i,j]))
        riskvalue[j] <- beta1+b+beta2*q  # can't not ecsure why add b
      }
      risk[foldlist[[n]][i]] <- sum(riskvalue)
      }
    }

    # get the R square
    foldindex <- as.vector(unlist(foldlist))
    if(!is.null(confounder)){
      # risk

      tmp <- phe[foldindex, c(response, confounder)]
      tmp$risk <- risk
      formula <- as.formula(paste0(response, "~."))
      lmmode <- summary(lm(formula, data = tmp))

      # no risk
      tmp2 <- tmp[,-ncol(tmp)]
      formula2 <- as.formula(paste0(response, "~."))
      lmmode2 <- summary(lm(formula2, data = tmp2))

      #　from the rsq.partial in the rsq package
      R2[m] <- 1-((1-lmmode$r.squared)/(1-lmmode2$r.squared))*(lmmode2$df[2]/lmmode$df[2])

      #R2[m] <- lmmode$adj.r.squared-lmmode2$adj.r.squared

    }else{
      y <- phe[foldindex, response]
      print(y)
      print(risk)
      lmmode <- summary(lm(y~risk))
      R2[m] <- lmmode$r.squared
    }
   }
  return(R2)
}


#' CvLasso
#' selsect feature on diffrent group
#' @param metadata , response
#' @param dataset , metabolism or microbe
#' @param feturetop
#'
#' @return
#' @export
#'
#' @examples
CvLassomulti <- function(metadata , dataset, response, feturetop=NULL){
  # combind data
  id <- intersect(rownames(metadata), rownames(dataset))
  print(paste0("the sample size is ", length(id)))
  tmp <- dataset[id, ]


  #　to tranform the data
  #  inverse-quantile normalized
  invt <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
  tmp2 <- as.data.frame(apply(tmp, 2, invt))
  mdat <- data.frame(metadata[id, response] , tmp2)

  # select feature
  library(glmnet)
  set.seed(123)
  nlev <- length(levels(droplevels(as.factor(mdat[,1]))))
  if(length(nlev)>2){
    lasso <- cv.glmnet(x=as.matrix(mdat[,-1]),
                       y=as.factor(mdat[,1]),
                       family='multinomial',
                       nfolds = 10,
                       alpha = 1,
                       nlambda = 100)
  }else{
    lasso <- cv.glmnet(x=as.matrix(mdat[,-1]),
                       y=as.factor(mdat[,1]),
                       family='binomial',
                       nfolds = 10,
                       alpha = 1,
                       nlambda = 100)

  }
  library(dplyr)
  library(tibble)
  lasso.mk <- coef(lasso, lasso$lambda.min)

  if(length(nlev)>2){
    featurescore <- list()
    for(i in 1:length(lasso.mk)){

      featurescore[[i]]<- data.frame(as.matrix(lasso.mk[[i]])) %>% setNames("Score") %>%
        rownames_to_column("Type") %>%
        slice(-c(1:2)) %>%
        filter(Score!=0) %>% mutate(dir = ifelse(Score >0 ,"pos", "neg")) %>%
        arrange(abs(Score))

    }

    # get the max & min
    minmax <- c()
    for(i in 1:length(featurescore)){
      tmpminmax <- c(max(featurescore[[i]][,2]), min(featurescore[[i]][,2]))
      minmax <- c(minmax, tmpminmax)
    }
    xmin <- min(minmax)
    xmax <- max(minmax)

    # plot feature impotance

    qplot <- function(qdat, top= T){

      if(is.null(top)){
        order_Type <- as.character(qdat$Type)
        qdat$Type <- factor(qdat$Type, levels = order_Type)
      }else{
        qdat <- tail(qdat, 10)
        qdat <- arrange(qdat, Score)
        order_Type <- as.character(qdat$Type)
        qdat$Type <- factor(qdat$Type, levels = order_Type)
      }

      qdat$dir <- factor(qdat$dir, levels = c("pos", "neg"))
      p <- ggplot(qdat, aes(y=Type,x=Score,color=dir)) +
        scale_x_continuous(limits =c(xmin,xmax)) +
        geom_segment(xend=0,aes(yend=Type),size=5)  +
        geom_vline(xintercept = 0) +
        scale_color_manual(values = c("#ef3b2c", "#2171b5"))+
        xlab("Coefficient estimate")+ylab("")+fontTheme+finalTheme


      return(p)
    }
    qlist <- list()
    for(i in 1:length(featurescore)){
      qlist[[i]] <- qplot(featurescore[[i]])
    }

    out <- list(featurescore, qlist)
  }else{
    featurescore <- data.frame(as.matrix(lasso.mk)) %>% setNames("Score") %>%
      rownames_to_column("Type") %>%
      slice(-c(1:2)) %>%
      filter(Score!=0) %>% mutate(dir = ifelse(Score >0 ,"pos", "neg")) %>%
      arrange(abs(Score))
    out <- featurescore


  }
  return(out)

}


#' CvLassoResponse
#' select feature using lasso method , to continous data
#' @param tag , reponse name
#' @param dataset1 , phenotype or response dataset
#' @param dataset2 ,
#'
#' @return
#' @export
#'
#' @examples
CvLassoResponse <- function(tag, dataset1, dataset2, transformM = "IQN"){
  # combind data
  id <- intersect(rownames(dataset1), rownames(dataset2))
  print(paste0("the sample size is ", length(id)))
  mdat <- data.frame(cbind(dataset1[id, tag, drop =F] , dataset2[id, ]))
  # rm the na row
  mdat <- mdat[!is.na(mdat[,1]), ]
  if(nrow(dat) < 10){
    return(NULL)
  }else{

  #　to tranform the data
  #  inverse-quantile normalized
  mdat2 <- as.data.frame(model.matrix(~., mdat))
  if(transformM == "IQN"){
    invt <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
    mdat3 <- as.data.frame(apply(mdat2, 2, invt))
  }else if(transforM == "AST"){
    AST <- function(x) {return(sign(x) * asin(sqrt(abs(x))))}
    mdat3 <- as.data.frame(apply(mdat2, 2, AST))
  }else{
    mdat3 <- mdat2
  }

  # select feature
  library(glmnet)
  set.seed(123)
  lasso <- cv.glmnet(x=as.matrix(mdat3[,-c(1:2)]),
                     y=mdat3[,2],
                     family='gaussian',
                     nfolds = 10,
                     alpha = 1,
                     nlambda = 100)
  library(dplyr)
  library(tibble)
  lasso.mk <- data.frame(as.matrix(coef(lasso, lasso$lambda.min))) %>%
    setNames("Score") %>%
    rownames_to_column("Type") %>%
    slice(-c(1:2)) %>%
    filter(Score!=0)
  if(nrow(lasso.mk) == 0){
    return(NULL)
  }else{
    data_all <- mdat3 %>% dplyr::select(tag, lasso.mk$Type)
    # cross validation of feature
    library(caret)
    set.seed(123)
    method <- "glm";
    num <- nrow(data_all)
    folds <- createFolds(y=data_all[,1], k=num )
    colnames(data_all)[1] <- "y"
    res <- data.frame()
    for (i in 1:num) {
      train_cv <- data_all[-folds[[i]], ]
      test_cv <- data_all[folds[[i]], ]
      if(method=="glm"){
        fit <- glm(y~., data=train_cv, family = "gaussian")
      } else if(method=="rf"){
        library(randomForest)
        fit <- randomForest(y~., data = train_cv, mtry=3, importance=T)
      } else if(method=="gbm"){
        library(gbm)
        fit <- gbm(y~., data = train_cv,
                   distribution = "gaussian",
                   n.trees = 1000,
                   shrinkage = 0.01,
                   n.minobsinnode = 1,
                   bag.fraction = 1,
                   interaction.depth = 8,
                   cv.folds = 5)
      }
      pred <- predict(fit, test_cv) %>% data.frame()
      pred_res <- cbind(test_cv$y, pred) %>%
        setNames(c("True", "Predict"))
      res <- rbind(res, pred_res)
    }
    return(list(res, lasso.mk))
  }
  }
}


# summary the result

summary_pred <- function(response, datalist, dataname, responsename, plotname){
  # need the libray
  library(reshape)
  library(pheatmap)
  library(glmnet)

  # generate the result
  out <- data.frame(True = 1, Predict = 1, type1 = "", type2 = "")
  out2 <- data.frame(Type = "", Score = 1, type1 = "", type2 = "")

  Rresult <- matrix(NA, nrow=length(responsename), ncol=length(datalist))
  colnames(Rresult) <-  dataname
  rownames(Rresult) <- responsename

  for(i in 1:length(responsename)){
    for(j in 1:length(dataname)){
      print(j)
      res <- CvLassoCluster(responsename[i], dataset1 = response, dataset2 = datalist[[j]])
      if(is.null(res)){
        next
        Rresult[i, j] <- NA
      }else{
        res[[1]]$type1 <- rep(responsename[i], nrow(res[[1]]))
        res[[1]]$type2 <- rep(dataname[j], nrow(res[[1]]))
        res[[2]]$type1 <- rep(responsename[i], nrow(res[[2]]))
        res[[2]]$type2 <- rep(dataname[j], nrow(res[[2]]))
        Rresult[i, j] <- cor.test(res[[1]][,1], res[[1]][,2], method="s")$estimate
        out <- rbind(out, res[[1]])
        out2 <- rbind(out2, res[[2]])
      }
    }
  }


  # plot the result
  feature_plot <- out2[-1, ]
  recast(feature_plot, Type+type2~type1) -> qdat
  annotation_row <- as.data.frame(qdat[, "type2", drop = F])
  rownames(annotation_row) <- paste0("test", 1:nrow(annotation_row))
  qdat2 <- qdat[, 3:10]
  rownames(qdat2) <- rownames(annotation_row)
  qdat2[is.na(qdat2)] <- 0
  qdat2[qdat2<0] <- -1
  qdat2[qdat2>0] <- 1

  # remove the only once
  index <- apply(abs(qdat2), 1, function(x){sum(x)>=2})
  annotation_row <- annotation_row[index,, drop=F]
  qdat3 <- as.data.frame(qdat2[index, ])

  pdf(paste0(plotname, ".feature.lasso.pdf"), width = 20, height = 12)
  pheatmap(t(qdat3), labels_col  = qdat$Type[index], cellheight = 8 , cellwidth = 8,legend = F,
           annotation_col = annotation_row)
  dev.off()

  return(Rresult)
}












## figure for metagenome
mytheme <- function(tmp){

  finalTheme <- theme_set(theme_bw()) +
        theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
        )

  fontTheme <- theme(
    axis.title=element_text(size=14, face = "bold", colour = "black"),
    text=element_text(size=12),
    legend.text = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 14, face = "bold", colour = "black"),
    legend.title = element_text(size = 14, face = "bold", colour = "black"),
    strip.text = element_text(face = "bold", size = 14)
  )

}

#' plot volcano
#' this function to plot volcano based ggplot2
#' @param qdat include column name Fdr & fc
#' @param fdcut  the foldchange cutoff
#'
#' @return figure
#' @export
#'
#' @examples
myvolcano <- function(qdat, fdcut){
  # to generate the enrich
  qdat$fc <- log2(qdat$fc)
  qdat$sig <- NA
  qdat$sig[(qdat$Fdr > 0.05|qdat$Fdr=="NA")|(qdat$fc < fdcut | qdat$fc > - fdcut)] <- "no"
  qdat$sig[qdat$Fdr <= 0.05 & qdat$fc >= fdcut] <- "up"
  qdat$sig[qdat$Fdr <= 0.05 & qdat$fc <= -fdcut] <- "down"
  qdat$lable <- ifelse(qdat$sig!="no", qdat$species, "")

  x_lim <- max(qdat$fc, -qdat$fc)
  library(ggplot2)
  library(RColorBrewer)

  p <- ggplot(qdat, aes(fc, -1*log10(Fdr),
                        color = sig))+geom_point(size = 1)+xlim(-x_lim,x_lim)+
    labs(x="log2(FoldChange)",y="-log10(FDR)")

  p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
    geom_hline(yintercept=-log10(0.05),linetype=4)+
    geom_vline(xintercept=c(-fdcut, fdcut),linetype=4)
  p <- p +theme(panel.grid =element_blank())+
    theme(axis.line = element_line(size=0))+ylim(0,20)
  #p <- p + geom_label_repel(aes(fc, -1*log10(Fdr), label=lable))
  p <- p + guides(colour = FALSE)
  p <- p + theme(axis.text=element_text(size=10),axis.title=element_text(size=10))
  p <- p +facet_grid(.~clust)

  return(p)

}


volcanoEffect <- function(qdat, effvar, effcut, showlabel=F, labelsize=8, ymax=10){

  # to generate the enrich
  qdat$fc <- qdat[, effvar]
  qdat$sig <- NA
  qdat$sig[(qdat$Fdr > 0.05|qdat$Fdr=="NA")|(qdat$fc < effcut | qdat$fc > -effcut)] <- "no"
  qdat$sig[qdat$Fdr <= 0.05 & qdat$fc >= effcut] <- "up"
  qdat$sig[qdat$Fdr <= 0.05 & qdat$fc <= -effcut] <- "down"
  qdat$label <- ifelse(qdat$sig!="no", qdat$label, "")

  x_lim <- max(qdat$fc, -qdat$fc)
  library(ggplot2)
  library(RColorBrewer)

  p <- ggplot(qdat, aes(fc, -1*log10(Fdr),
                        color = sig))+geom_point(size = 1)+xlim(-x_lim,x_lim)+
    labs(x="effect size", y="-log10(FDR)")

  p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
    geom_hline(yintercept=-log10(0.05),linetype=4)+
    geom_vline(xintercept=c(-effcut, effcut),linetype=4)
  p <- p +theme(panel.grid =element_blank())+
    theme(axis.line = element_line(size=0))+ylim(0,ymax)
  if(showlabel){
    p <- p + geom_text_repel(aes(fc, -1*log10(Fdr), label=label),size=labelsize)
  }
  p <- p + guides(colour = FALSE)
  p <- p + theme_classic()
  p <- p +facet_grid(.~clust)

  return(p)

}



#' Multiplot
#' this function to plot volcano based ggplot2
#' @param ...
#' @param plotlist
#' @param file
#' @param cols
#' @param layout
#'
#' @return
#' @export
#'
#' @examples
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' pcoaFig
#'
#' @param data.dist
#' @param cluster
#'
#' @return figure
#' @export
#'
#' @examples
pcoaFig <- function(data.dist, cluster){

  #library(RColorBrewer)

  # plot the pcoa
  obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
  var1 <- round((obs.pcoa$eig[1]/sum(obs.pcoa$eig))*100,2)
  var2 <- round((obs.pcoa$eig[2]/sum(obs.pcoa$eig))*100,2)
  minX <- min(obs.pcoa$li[,1])
  maxX <- max(obs.pcoa$li[,1])
  minY <- min(obs.pcoa$li[,2])
  maxY <- max(obs.pcoa$li[,2])

  plot(0,0, main = "Pcoa", type = "n",
       xlab=paste("Pco1 (",var1,"%)"),ylab=paste("Pco2 (",var2,"%)"),
       xlim=c(minX-10^floor(log10(abs(minX))),maxX+10^floor(log10(abs(maxX)))),
       ylim=c(minY-10^floor(log10(abs(minY))),maxY+10^floor(log10(abs(maxY)))),
       frame=TRUE, cex=0.5, add=T)
  s.class(obs.pcoa$li, fac=as.factor(cluster), cell =2 ,
          csta = 0 , col = brewer.pal(8,'Set2'), grid=F, add.plot = T)


}


#' nmdsFig
#'
#' @param data
#' @param method
#' @param config
#'
#' @return figure
#' @export
#'
#' @examples

nmdsFig <- function(data, method = "bray", config){

  mds <- metaMDS(data, distance = method, k = 2, trymax = 50)

  dat2 <- data.frame(mds$point[,1:2])

  minX <- min(dat2[,1])
  maxX <- max(dat2[,1])
  minY <- min(dat2[,2])
  maxY <- max(dat2[,2])

  plot(0,0, main = "Pcoa", type = "n",
       xlab=paste("Pco1 (",var1,"%)"),ylab=paste("Pco2 (",var2,"%)"),
       xlim=c(minX-5^floor(log10(abs(minX))),maxX+5^floor(log10(abs(maxX)))),
       ylim=c(minY-5^floor(log10(abs(minY))),maxY+5^floor(log10(abs(maxY)))),
       frame=TRUE, cex=0.5, add=T)
  s.class(dat2, fac=as.factor(config), cell =2 ,
          csta = 0 , col = brewer.pal(8,'Set2'), grid=F, add.plot = T)


}


#' nmdsFigEx
#' plot the nmds figure
#' @param data
#' @param method
#' @param config
#' @param color
#'
#' @return
#' @export
#'
#' @examples
nmdsFigEx <- function(data, method = "bray", config, color = c("#66C2A5","#E78AC3")){

  mds <- metaMDS(data, distance = method, k = 2, trymax = 50)

  dat2 <- data.frame(mds$point[,1:2])
  colnames(dat2) <- c("NMDS1", "NMDS2")
  # to confirm the config plot
  dat2$group <- as.factor(config[,1])
  dat2$group2 <- as.factor(config[,2])
  plot <- ggplot(dat2,aes(NMDS1, NMDS2,color=group))+geom_point()+
    stat_ellipse(level = 0.7, show.legend = F)+
    ggtitle(paste0("NMDS","_",method))+
    theme(plot.title = element_text(hjust = 0.5,size = 14),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 13),
          panel.grid = element_blank(),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 13))+
    scale_color_manual(values = color)+facet_wrap(.~group2)

  return(plot)

}


#' topTax
#' plot the high abundance tax
#' @param metadata , tax profile row is tax col is sample
#' @param K , top number
#' @param rmUnclass , if to remove the unclass tax
#'
#' @return figure
#' @export
#'
#' @examples
topTax <- function(metadata, K = 20, rmUnclass = F){

  # generate the top20 tax profile
  if(rmUnclass){
    # rm the unclass & renorm
    unclassindex <- which(rownames(metadata) %in% "unclassed")
    metadata <- metadata[-unclassindex, ]
    metadata <- apply(metadata, 2, function(x){y <- x/sum(x); return(y)})
  }
  # order the tax
  top <- names(head(sort(apply(metadata , 1, mean), decreasing = T), K))
  topdata <- metadata[top, ]
  lessdata <- metadata[-which(rownames(metadata) %in% top), ]
  otherdata <- t(data.frame(apply(lessdata, 2, sum)))
  rownames(otherdata) <- "Others"
  qdat <- rbind(topdata, otherdata)
  # plot the result
  naOrder <- rownames(qdat)
  idOrder <- colnames(qdat)[order(qdat[1,], decreasing = T)]
  qdat$sample <- rownames(qdat)
  qdat2 <- melt(qdat)
  colnames(qdat2) <- c("Tax", "Sample", "value")
  qdat2$Tax <- factor(qdat2$Tax, levels = rev(naOrder))
  qdat2$Sample <- factor(qdat2$Sample, levels = idOrder)
  # ggplot
  ggplot(qdat2, aes(Sample, value, fill=Tax))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 1)+
    theme_classic()+
    #mytheme+
    theme(axis.text.x = element_blank(),
          #legend.position = c(0.01,0.99),
          #legend.justification = c(0,1),
          axis.ticks = element_blank())+
    xlab("")+ylab("Relative Abundance")+
    scale_y_continuous(expand = c(0,0),breaks = c(0,.2,.4,.6,.8,1))

}



########## corplot ###############



#' corPlot
#' pheatmap to figure the correlation
#' @param corres
#' @param cutoff
#' @param adjust
#' @param tr
#'
#' @return
#' @export
#'
#' @examples
corPlot <- function(corres , cutoff, adjust, tr){


  #　sub function
  trans <- function(x){
    if(x <= 0.05 & x>0.01){
      out <-"*"
    }else if(x <= 0.01 & x>0.001){
      out <- "**"
    }else if(x <= 0.001){
      out<- "***"
    }else{
      out <- " "
    }
    return(out)
  }
  #  ready the data
  xname <- rownames(corres)
  corres <- apply(corres, 2, as.numeric)
  sp.corr.t <- corres
  rownames(sp.corr.t) <- xname
  # split the dataset
  index <- 2*c(1:(ncol(sp.corr.t)/2))
  dat.pvalue <- sp.corr.t[, index]
  dat.cor <- sp.corr.t[, -index]
  colnames(dat.cor) <- gsub("_p.value", "", colnames(dat.pvalue))
  # adjust
  if(adjust){
    dat.pvalue <- matrix(p.adjust(as.vector(as.matrix(dat.pvalue)),method = "BH"),
                         nrow = nrow(dat.pvalue))
  }
  pvalue.index <- apply(dat.pvalue, 2, function(x) any(x < cutoff))
  pvalue.index2 <- apply(dat.pvalue, 1, function(x) any(x < cutoff))
  dat.cor.cle <- dat.cor[pvalue.index2, pvalue.index]
  dat.pva.cle <- dat.pvalue[pvalue.index2, pvalue.index]
  num<- matrix(NA,nrow = nrow(dat.pva.cle), ncol = ncol(dat.pva.cle))

  for(i in 1:ncol(dat.pva.cle)){
    num[,i] <- mapply(trans, dat.pva.cle[,i])
  }
  colt<-c("#4C38CB","#9191C8","#DADAEC","#F0C1C1","#E28383","#D44545","#CD2626")

  if(tr){
    dat.cor.cle <- t(dat.cor.cle)
    num <- t(num)
  }
  pheatmap(dat.cor.cle,
           treeheight_row=43,
           treeheight_col=23,
           cellwidth=20,
           cellheight=8,
           cluster_cols=T,
           cluster_rows=T,
           fontsize_row=8,
           fontsize_col=13,
           show_colnames=T,
           display_numbers=num,
           color =colt,
           breaks=seq(-0.6,0.6,0.2),
           legend_breaks = c(-0.6,-0.3, 0, 0.3,0.6),
           legend_labels = c("-0.6","-0.3","0", "0.3", "0.6"),
           number_color = "black")

}









## myfigure for metagenome

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




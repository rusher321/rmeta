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
       xlab= "NMDS1",
       ylab= "NMDS2",
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
#' @param sample_order, vector of sample id
#' @param tax_order, vector of taxnomy
#' @param tax_colour, vector of colour
#' @return figure
#' @export
#'
#' @examples
topTax <- function(metadata, K = 20, rmUnclass = F, sample_order =NULL ,
                   tax_order = NULL, tax_colour=NULL){

  # generate the top tax profile
  if(rmUnclass){
    # rm the unclass & renorm
    unclassindex <- which(rownames(metadata) %in% "unclassed")
    metadata <- metadata[-unclassindex, ]
    metadata <- apply(metadata, 2, function(x){y <- x/sum(x); return(y)})
  }

  # order the tax
  if(!is.null(tax_order)){
    data_sub <- metadata[tax_order, ]
    lessdata <- metadata[-which(rownames(metadata) %in% tax_order), ]
    otherdata <- t(data.frame(apply(lessdata, 2, sum)))
    rownames(otherdata) <- "Others"
    qdat <- rbind(lessdata, otherdata)
  }else{
    top <- names(head(sort(apply(metadata , 1, mean), decreasing = T), K))
    topdata <- metadata[top, ]
    lessdata <- metadata[-which(rownames(metadata) %in% top), ]
    otherdata <- t(data.frame(apply(lessdata, 2, sum)))
    rownames(otherdata) <- "Others"
    qdat <- rbind(topdata, otherdata)
  }
  naOrder <- rownames(qdat)

  # order the sample
  if(!is.null(sample_order)){
    idOrder = sample_order
  }else{
    idOrder <- colnames(qdat)[order(qdat[1,], decreasing = T)]
  }
  qdat <- as.data.frame(qdat, check.names=F)
  qdat$sample <- rownames(qdat)
  qdat2 <- melt(qdat)
  colnames(qdat2) <- c("Tax", "Sample", "value")
  qdat2$Tax <- factor(qdat2$Tax, levels = rev(naOrder))
  qdat2$Sample <- factor(qdat2$Sample, levels = idOrder)
  # ggplot
  p = ggplot(qdat2, aes(Sample, value, fill=Tax))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 1)+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank())+
    xlab("")+ylab("Relative Abundance")+
    scale_y_continuous(expand = c(0,0),breaks = c(0,.2,.4,.6,.8,1))
  if(!is.null(tax_colour)){
    p+scale_fill_manual(values = tax_colour)
  }

  return(p)
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



# centor compare

#' centroComp
#' compare the center distance between groups
#' @param pro metagenome profile ,row is sample ID
#' @param method distance method
#' @param config group info
#' @param color set group color in plot
#'
#' @return
#' @export
#'
#' @examples
centorComp <- function(pro, method , config, color){

  id <- intersect(rownames(pro), rownames(config))
  pro <- pro[id,]
  config <- config[id, ]
  # compute the distance

  prodis <- vegan::vegdist(pro, method = method)
  mod <- betadisper(prodis, config)

  qdata <- data.frame(dis = mod$distance, label = config)

  # plot
  my_comparisons = list()
  num <- combn(length(unique(config)),2)
  for(i in 1:ncol(num)){my_comparisons[[i]] <- num[,i]}

  p <- ggboxplot(qdata, x="label", y = "dis", color = "label" ,add = "jitter",alpha=0.6,size = 0.5)+
    stat_compare_means(comparisons = my_comparisons)+
    scale_color_manual(values=color)+
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),legend.position = "none")+xlab("")+
    ylab("Distance to centroid")

  return(p)

}


#' correlation between pairwise
#'
#' @param data
#' @param method
#' @param cor_matrix
#' @param nbreaks
#' @param digits
#' @param name
#' @param low
#' @param mid
#' @param high
#' @param midpoint
#' @param palette
#' @param geom
#' @param min_size
#' @param max_size
#' @param label
#' @param label_alpha
#' @param label_color
#' @param label_round
#' @param label_size
#' @param limits
#' @param drop
#' @param layout.exp
#' @param legend.position
#' @param legend.size
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggcorr <- function(
  data,
  method = c("pairwise", "pearson"),
  cor_matrix = NULL,
  nbreaks = NULL,
  digits = 2,
  name = "",
  low = "#3B9AB2",
  mid = "#EEEEEE",
  high = "#F21A00",
  midpoint = 0,
  palette = NULL,
  geom = "tile",
  min_size = 2,
  max_size = 6,
  label = FALSE,
  label_alpha = FALSE,
  label_color = "black",
  label_round = 1,
  label_size = 4,
  limits = TRUE,
  drop = !limits,
  layout.exp = 0,
  legend.position = "right",
  legend.size = 9,
  ...) {

  # -- required packages -------------------------------------------------------

  require(ggplot2, quietly = TRUE)
  require(reshape2, quietly = TRUE)

  # -- check geom argument -----------------------------------------------------

  if (length(geom) > 1 || !geom %in% c("blank", "circle", "text", "tile")) {
    stop("incorrect geom value")
  }

  # -- correlation method ------------------------------------------------------

  if (length(method) == 1) {
    method = c(method, "pearson") # for backwards compatibility
  }

  # -- check data columns ------------------------------------------------------

  if (!is.null(data)) {

    if (!is.data.frame(data)) {
      data = as.data.frame(data)
    }

    x = which(!sapply(data, is.numeric))

    if (length(x) > 0) {

      warning(paste("data in column(s)",
                    paste0(paste0("'", names(data)[x], "'"), collapse = ", "),
                    "are not numeric and were ignored"))

      data = data[, -x ]

    }

  }

  # -- correlation matrix ------------------------------------------------------

  if (is.null(cor_matrix)) {
    cor_matrix = cor(data, use = method[1], method = method[2])
  }

  m = cor_matrix
  colnames(m) = rownames(m) = gsub(" ", "_", colnames(m)) # protect spaces

  # -- correlation data.frame --------------------------------------------------

  m = data.frame(m * lower.tri(m))
  m$.ggally_ggcorr_row_names = rownames(m)
  m = reshape2::melt(m, id.vars = ".ggally_ggcorr_row_names")
  names(m) = c("x", "y", "coefficient")
  m$coefficient[ m$coefficient == 0 ] = NA

  # -- correlation quantiles ---------------------------------------------------

  if (!is.null(nbreaks)) {

    x = seq(-1, 1, length.out = nbreaks + 1)

    if (!nbreaks %% 2) {
      x = sort(c(x, 0))
    }

    m$breaks = cut(m$coefficient, breaks = unique(x), include.lowest = TRUE,
                   dig.lab = digits)

  }

  # -- gradient midpoint -------------------------------------------------------

  if (is.null(midpoint)) {

    midpoint = median(m$coefficient, na.rm = TRUE)
    message(paste("Color gradient midpoint set at median correlation to",
                  round(midpoint, 2)))

  }

  # -- plot structure ----------------------------------------------------------

  m$label = round(m$coefficient, label_round)
  p = ggplot(na.omit(m), aes(x, y))

  if (geom == "tile") {

    if (is.null(nbreaks)) {

      # -- tiles, continuous ---------------------------------------------------

      p = p +
        geom_tile(aes(fill = coefficient), color = "white")

    } else {

      # -- tiles, ordinal ------------------------------------------------------

      p = p +
        geom_tile(aes(fill = breaks), color = "white")

    }

    # -- tiles, color scale ----------------------------------------------------

    if (is.null(nbreaks) && limits) {

      p = p +
        scale_fill_gradient2(name, low = low, mid = mid, high = high,
                             midpoint = midpoint, limits = c(-1, 1))

    } else if (is.null(nbreaks)) {

      p = p +
        scale_fill_gradient2(name, low = low, mid = mid, high = high,
                             midpoint = midpoint)

    } else if (is.null(palette)) {

      x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))

      p = p +
        scale_fill_manual(name, values = x, drop = drop)

    } else {

      p = p +
        scale_fill_brewer(name, palette = palette, drop = drop)

    }

  } else if (geom == "circle") {

    p = p +
      geom_point(aes(size = abs(coefficient) * 1.25), color = "grey50") # border

    if (is.null(nbreaks)) {

      # -- circles, continuous -------------------------------------------------

      p = p +
        geom_point(aes(size = abs(coefficient), color = coefficient))

    } else {

      # -- circles, ordinal ----------------------------------------------------

      p = p +
        geom_point(aes(size = abs(coefficient), color = breaks))

    }

    p = p +
      scale_size_continuous(range = c(min_size, max_size)) +
      guides(size = FALSE)

    r = list(size = (min_size + max_size) / 2)

    # -- circles, color scale --------------------------------------------------

    if (is.null(nbreaks) && limits) {

      p = p +
        scale_color_gradient2(name, low = low, mid = mid, high = high,
                              midpoint = midpoint, limits = c(-1, 1))

    } else if (is.null(nbreaks)) {

      p = p +
        scale_color_gradient2(name, low = low, mid = mid, high = high,
                              midpoint = midpoint)

    } else if (is.null(palette)) {

      x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))

      p = p +
        scale_color_manual(name, values = x, drop = drop) +
        guides(color = guide_legend(override.aes = r))

    } else {

      p = p +
        scale_color_brewer(name, palette = palette, drop = drop) +
        guides(color = guide_legend(override.aes = r))

    }

  } else if (geom == "text") {

    if (is.null(nbreaks)) {

      # -- text, continuous ----------------------------------------------------

      p = p +
        geom_text(aes(label = label, color = coefficient), size = label_size)

    } else {

      # -- text, ordinal -------------------------------------------------------

      p = p +
        geom_text(aes(label = label, color = breaks), size = label_size)

    }

    # -- text, color scale ----------------------------------------------------

    if (is.null(nbreaks) && limits) {

      p = p +
        scale_color_gradient2(name, low = low, mid = mid, high = high,
                              midpoint = midpoint, limits = c(-1, 1))

    } else if (is.null(nbreaks)) {

      p = p +
        scale_color_gradient2(name, low = low, mid = mid, high = high,
                              midpoint = midpoint)

    } else if (is.null(palette)) {

      x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))

      p = p +
        scale_color_manual(name, values = x, drop = drop)

    } else {

      p = p +
        scale_color_brewer(name, palette = palette, drop = drop)

    }

  }

  # -- coefficient labels ------------------------------------------------------

  if (label) {

    if (isTRUE(label_alpha)) {

      p = p +
        geom_text(aes(x, y, label = label, alpha = abs(coefficient)),
                  color = label_color, size = label_size,
                  show_guide = FALSE)

    } else if (label_alpha > 0) {

      p = p +
        geom_text(aes(x, y, label = label, show_guide = FALSE),
                  alpha = label_alpha, color = label_color, size = label_size)

    } else {

      p = p +
        geom_text(aes(x, y, label = label),
                  color = label_color, size = label_size)

    }

  }

  # -- horizontal scale expansion ----------------------------------------------

  l = levels(m$y)

  if (!is.numeric(layout.exp) || layout.exp < 0) {
    stop("incorrect layout.exp value")
  } else if (layout.exp > 0) {
    l = c(rep(NA, as.integer(layout.exp)), l)
  }

  p = p  +
    geom_text(data = m[ m$x == m$y & is.na(m$coefficient), ],
              aes(label = x), ...) +
    scale_x_discrete(breaks = NULL, limits = l) +
    scale_y_discrete(breaks = NULL, limits = levels(m$y)) +
    labs(x = NULL, y = NULL) +
    coord_equal() +
    theme(
      panel.background = element_blank(),
      legend.key = element_blank(),
      legend.position = legend.position,
      legend.title = element_text(size = legend.size),
      legend.text = element_text(size = legend.size)
    )

  return(p)

}







## myfigure for metagenome

#' this function to plot volcano based ggplot2
#'
#' @param qdat include column name Fdr & fc
#' @param fdcut  the foldchange cutoff
#'
#' @return a figure
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

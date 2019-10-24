########### network for species ###########


#' renorm
#' normalization the compsition data
#' @param dat , dataframe row is variable , col is sample ID
#'
#' @return dataframe
#' @export
#'
#' @examples
#'
renorm <- function(dat){
  # normlization function
  trans <- function(x){
    y <- x/(sum(x))
    return(y)
  }
  # tran the dataframe
  dat2 <- t(dat)
  dat3 <- t(apply(dat2, 2, trans))

  return(dat3)
}

#' sparccNet
#'
#' @param cutoff ,numberic 0.3
#' @param main , str topic of figure
#' @param dat , dataframe  rawdata
#' @param layout , the lay style of node layout
#'
#' @return figure
#' @export
#'
#' @examples

sparccNet <- function(dat, cutoff, main, layout = "layout.circle"){

  # library(igraph)
  # library(Matrix)
  # library(SpiecEasi)
  datRenorm <- round(renorm(dat)*10^5)
  basesparcc <- sparcc(datRenorm)
  graph <- sparcc$Cor
  num.row <- nrow(graph)

  # tran the corelation to adj matrox

  for(i in 1:num.row){
    for(j in 1:num.row){
      a <- graph[i, j]
      graph[i,j] <- ifelse(abs(a) >= cutoff, a, 0)
    }
  }

  diag(graph) <- 0
  igraph <- adj2igraph(Matrix(graph, sparse=TRUE))
  vsize.g1 <- rowMeans(clr(dat, 1))+6

  # set edge color，postive correlation
  # 设定为red, negative correlation设定为blue
  igraph.weight = E(igraph)$weight
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#fc9272", ifelse(E.color<0, "#31a354","grey"))
  E(igraph)$color = as.character(E.color)

  # add the edge width
  E(igraph)$width = abs(igraph.weight)*4
  # add the node color
  V(igraph)$color <- "#636363"
  plot(igraph,layout=layout.circle, vertex.size=vsize.g1,
       vertex.label=NA, main = main )

}




#' pairNet
#'
#' @param Cor1 , corelation matrix 1
#' @param Cor2 , corelation matrix 2
#' @param cutoff, the eage cutoff
#' @param species, retain the species
#'
#' @return
#' @export
#'
#' @examples
pairNet <- function(Cor1, Cor2, cutoff, species){

  # generate adj matrix on species

  bacindex <- which(rownames(Cor1) == species)
  retainbac <- names(Cor1[,bacindex][abs(Cor1[,bacindex])>=cutoff])
  retainnet <- Cor1[retainbac, retainbac]
  which(retainbac== species) -> speindex
  retainnet[-speindex, -speindex] <- 0
  retainnet[speindex, speindex] <- 0

  bacindex2 <- which(rownames(Cor2) == species)
  retainbac2 <- names(Cor2[,bacindex][abs(Cor2[,bacindex])>=cutoff])
  retainnet2 <- Cor2[retainbac2, retainbac2]
  which(retainbac2 == species) -> speindex2
  retainnet2[-speindex2, -speindex2] <- 0
  retainnet2[speindex2, speindex2] <- 0
  rownames(retainnet2)[speindex] <- colnames(retainnet2)[speindex] <- paste0(species,"_C2")
  retainbac2[speindex2] <- paste0(species,"_C2")
  # combine the two network

  comBac <- unique(c(retainbac, retainbac2))
  out <- matrix(0, nrow=length(comBac), ncol=length(comBac))
  rownames(out) <- colnames(out) <- comBac

  # input the corresponde value
  which(comBac %in% retainbac) -> netindex1
  netindex1[speindex] -> netspindex1
  out[netindex1, netspindex1] <- retainnet[, speindex]
  out[netspindex1, netindex1] <- retainnet[speindex, ]

  which(comBac %in% retainbac2) -> netindex2
  netindex2[speindex2] -> netspindex2
  out[netindex2, netspindex2] <- retainnet2[, speindex2]
  out[netspindex2, netindex2] <- retainnet2[speindex2, ]

  return(out)

}

#' simpleNet
#'
#' @param adjmatrix
#' @param data
#' @param main
#' @param layout
#'
#' @return
#' @export
#'
#' @examples
#'
simplenet <- function(adjmatrix, data, main, layout = "layout.circle"){
  #library(igraph)
  #library(Matrix)

  graph <- adjmatrix
  num.row <- nrow(graph)
  igraph <- adj2igraph(Matrix(graph, sparse=TRUE))
  vsize.g1 <- rowMeans(clr(data, 2))+6
  # set edge color，postive correlation 设定为red, negative correlation设定为blue
  igraph.weight = E(igraph)$weight
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#fc9272",ifelse(E.color<0, "#31a354","grey"))
  E(igraph)$color = as.character(E.color)

  # add the edge width
  E(igraph)$width = abs(igraph.weight)*4
  # add the node color
  #V(igraph)$color <- "#636363"
  plot(igraph,layout=layout, vertex.size=vsize.g1,
       vertex.label=colnames(data), vertex.label.cex = .6,
       vertex.label.color = "black",
       main = main )

}

########### model for metagenome ############

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

  perc.under <- ((over/less)-1)*100
  perc.over <- over*100*100/(less*perc.under)

  # generate the new data
  data$group <- group[,1]
  newData <- SMOTE(group ~ ., data, perc.over = perc.over,
                   perc.under=perc.under)

  return(newData)

}



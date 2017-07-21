#' Predict if a transcript is subject to nonsense mediated decay
#'
#' @param orfs orf information data.frame produced by GeneStructureTools::get_orfs()
#' @return vector with predicted class
#' @export
#' @import caret
#' @import gbm
#' @examples
#' @author Beth Signal
predictNMD <- function(orfs){
    m <- match(notNMDnames, colnames(orfs))
    if(all(!is.na(m))){
        orfs <- orfs[,m]
        orfs_p <- predict(preProcValues, orfs)

        # p <- predict(model_gbm, orfs_p, type="prob")[,1]
        n <- predict(model_gbm, orfs_p)
        return(n)

    }else{
        message("please provide the correct input")
    }
}

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
        keep <- which(apply(orfs, 1, function(x) all(!is.na(x))))
        orfs_p <- orfs[keep,]

        orfs_p <- predict(preProcValues, orfs_p)

        # p <- predict(model_gbm, orfs_p, type="prob")[,1]
        n <- predict(model_gbm, orfs_p)

        all_n <- rep(NA, nrow(orfs))
        all_n[keep] <- as.character(n)

        return(all_n)

    }else{
        message("please provide the correct input")
    }
}

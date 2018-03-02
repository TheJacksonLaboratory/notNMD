#' Predict if a transcript is subject to nonsense mediated decay
#'
#' @param orfs orf information data.frame produced by GeneStructureTools::getOrfs()
#' @param output class or prob - output a class prediction or probability of nonsense mediated decay?
#' @param model use base or lncRNA model?
#' @return vector with predicted class or probaility for nonsense mediated decay
#' @export
#' @import caret
#' @import gbm
#' @examples
#' @author Beth Signal
predictNMD <- function(orfs, output="class", model="base"){
    m <- match(notNMDnames, colnames(orfs))
    if(all(!is.na(m))){
        orfs <- orfs[,m]
        keep <- which(apply(orfs, 1, function(x) all(!is.na(x))))
        orfs_p <- orfs[keep,]

        orfs_p <- predict(preProcValues, orfs_p)
        all_n <- rep(NA, nrow(orfs))

        if(output=="prob"){
            n <- predict(model_gbm.long, orfs_p, type="prob")[,1]
            all_n[keep] <- as.numeric(n)
            if(model == "lncRNA"){
                n <- predict(model_gbm.long.lnc, orfs_p, type="prob")[,1]
                all_n[keep] <- as.numeric(n)
            }

        }else{
            n <- predict(model_gbm.long, orfs_p)
            all_n[keep] <- as.character(n)
            if(model == "lncRNA"){
                n <- predict(model_gbm.long.lnc, orfs_p)
                all_n[keep] <- as.character(n)
            }
        }
        return(all_n)
    }else{
        message("please provide the correct input")
    }
}

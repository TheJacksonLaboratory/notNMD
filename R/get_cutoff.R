#' Get min/max cutoff value for a performance metric of the notNMD model
#'
#' @param value value for cutoff
#' @param metric what performance metric the value relates to
#' @param model use base or lncRNA model?
#' @return min or max (or both) cutoff value where the metric is above the value
#' @export
#' @author Beth Signal
#' @examples
#' get_cutoff(0.9, "Pos Pred Value")
#' get_cutoff(0.9, "Neg Pred Value")
#' get_cutoff(0.8, "Accuracy")
get_cutoff <- function(value, metric, model="base"){

    if(model == "base"){
        conf_v <- conf_v.base
    }else if(model == "lncRNA"){
        conf_v <- conf_v.lncRNA
    }else{
        stop("please specify model as 'base' or 'lncRNA'")
    }

    if(metric %in% colnames(conf_v)){

        cutoffs <- c(min(conf_v$val[conf_v[,which(colnames(conf_v) == metric)] > value], na.rm=T),
                     max(conf_v$val[conf_v[,which(colnames(conf_v) == metric)] > value], na.rm=T))

        minmax <- c(min(conf_v$val[!is.na(conf_v[,which(colnames(conf_v) == metric)])]),
                    max(conf_v$val[!is.na(conf_v[,which(colnames(conf_v) == metric)])]))

        if(any(cutoffs %in% minmax)){
            cutoffs <- cutoffs[!(cutoffs%in% minmax)]
        }

        return(cutoffs)
    }else{
        message("please specify metric as one of the following:")
        message(paste(colnames(conf_v)[-19], collapse=", "))
    }

}

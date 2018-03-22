library(caret)

model_perf <- read.csv("source_data/base_model_performance.csv")
vals <- seq(0,1,by=0.01)
c_list <- list()
for(v in seq_along(vals)){

    model_perf$pred_class = ifelse(model_perf$prob > vals[v], "nonsense_mediated_decay","not_nmd")
    c_list[[v]] <- caret::confusionMatrix(model_perf$pred_class, model_perf$Class)
}

conf_v = cbind(as.data.frame(matrix(unlist(lapply(c_list, function(x) x$byClass)), ncol=11, byrow = T)),
               as.data.frame(matrix(unlist(lapply(c_list, function(x) x$overall)), ncol=7, byrow = T)))
colnames(conf_v) = c(names(c_list[[1]]$byClass),names(c_list[[1]]$overall))
conf_v$val = vals
conf_v.base <- conf_v

model_perf <- read.csv("source_data/lnc_model_performance.csv")
vals <- seq(0,1,by=0.01)
c_list <- list()
for(v in seq_along(vals)){

    model_perf$pred_class = ifelse(model_perf$prob > vals[v], "nonsense_mediated_decay","not_nmd")
    c_list[[v]] <- caret::confusionMatrix(model_perf$pred_class, model_perf$Class)
}

conf_v = cbind(as.data.frame(matrix(unlist(lapply(c_list, function(x) x$byClass)), ncol=11, byrow = T)),
               as.data.frame(matrix(unlist(lapply(c_list, function(x) x$overall)), ncol=7, byrow = T)))
colnames(conf_v) = c(names(c_list[[1]]$byClass),names(c_list[[1]]$overall))
conf_v$val = vals
conf_v.lncRNA <- conf_v

load("R/sysdata.rda")
save(model_gbm.long, model_gbm.long.lnc, preProcValues, notNMDnames, conf_v.base, conf_v.lncRNA, file="R/sysdata.rda")

setwd("/Volumes/Terra/Kunoichi_Files/Projects/notNMD")

library(plyr)
library(GenomicRanges)
library(caret)
devtools::install("~/Documents/Projects/GeneStructureTools/")
library(GeneStructureTools)
library(gbm)

source("~/Documents/Projects/mouse_splicing/scripts/support.R")

# Gencode annotations
gtf <- rtracklayer::import("~/Documents/Projects/resources/genomes/gencode.v21.annotation.gtf")
gtf.exons <- gtf[gtf$type=="exon"]
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

load("source_data/core_model_data.Rdata")
load("source_data/model_gbm_long.RData")

################################################################################

orfs_train <- orfs_all[orfs_all$transcript_type == "protein_coding" |
                           orfs_all$transcript_type == "nonsense_mediated_decay" |
                           orfs_all$transcript_type2=="lncRNA",]

#longest ORFs for protein coding genes
orfs_train <- arrange(orfs_train, plyr::desc(orf_length))
orfs_train <- orfs_train[!duplicated(orfs_train$id),]

rm <- which(apply(orfs_train, 1, function(x) any(is.na(x))))
if(length(rm) > 0){
    orfs_train <- orfs_train[-rm,]
}

set.seed(2)
index.lnc <- sample(which(orfs_train$transcript_type2 == "lncRNA"), ((4500*5)*0.15) + 1000*0.15)
train_index.lnc <- index.lnc[1:((4500*5)*0.15)]
test_index.lnc <- index.lnc[-(1:((4500*5)*0.15))]

train_ids.lnc <- orfs_train$id[train_index.lnc]
test_ids.lnc <- orfs_train$id[test_index.lnc]

set.seed(3)
pc_ids <- orfs_train$id[orfs_train$transcript_type2=="protein_coding" & orfs_train$id %in% train_ids]
pc_ids <- sample(pc_ids, (4500*5)-length(train_ids.lnc))
nmd_ids <- orfs_train$id[orfs_train$transcript_type2=="nonsense_mediated_decay" & orfs_train$id %in% train_ids]
train_ids.lnc <- c(train_ids.lnc, pc_ids, nmd_ids)

set.seed(4)
pc_ids <- orfs_train$id[orfs_train$transcript_type2=="protein_coding" & orfs_train$id %in% test_ids]
pc_ids <- sample(pc_ids, (1000)-length(test_ids.lnc))
nmd_ids <- orfs_train$id[orfs_train$transcript_type2=="nonsense_mediated_decay" & orfs_train$id %in% test_ids]
test_ids.lnc <- c(test_ids.lnc, pc_ids, nmd_ids)

train <- orfs_train[orfs_train$id %in% train_ids.lnc,]
train$id <- NULL
train$orf_sequence <- NULL
train$start_site <- NULL
train$stop_site <- NULL
train$frame <- NULL
train$gene_id <- NULL
train$transcript_type2 <- NULL
train$min_dist_to_junction_a <- NULL

test <- orfs_train[orfs_train$id %in% test_ids.lnc,]
test$orf_sequence <- NULL
test$Class <- "nonsense_mediated_decay"
test$Class[test$transcript_type != "nonsense_mediated_decay"] <- "not_nmd"

preProcValues.long.lnc <- preProcess(train, method = c("center", "scale"))
train_p <- predict(preProcValues.long.lnc, train)
train_p$transcript_type[train_p$transcript_type != "nonsense_mediated_decay"] <- "not_nmd"

test_p <- predict(preProcValues.long.lnc, test)

fitControl <- trainControl(method = "repeatedcv",number = 3,repeats = 10,classProbs=TRUE)

model_weights <- ifelse(train_p$transcript_type == "nonsense_mediated_decay",
                        (1/table(train_p$transcript_type)[1])*0.5,
                        (1/table(train_p$transcript_type)[2])*0.5)

if(!file.exists("source_data/model_gbm_long.lnc.RData") | remake_models == TRUE){
    model_gbm.long.lnc <- train(transcript_type ~ .,
                            data=train_p,
                            method = "gbm",
                            trControl=fitControl,
                            weights=model_weights,
                            verbose=F)
    save(model_gbm.long.lnc, preProcValues.long.lnc, file="source_data/model_gbm_long.lnc.RData")
}else{

}
p <- predict(model_gbm.long.lnc, test_p, type="prob")[,1]
test$prob=p
n <- predict(model_gbm.long.lnc, test_p)
c_long.lnc <- confusionMatrix(test$Class, n)

test.lnc <- test[,c(1:2,21, 24,25)]
write.csv(test.lnc, "~/Documents/Projects/notNMD/source_data/lnc_model_performance.csv", row.names = F,quote=F)

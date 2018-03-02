library(plyr)
library(GenomicRanges)
library(caret)
devtools::install("~/Documents/Projects/GeneStructureTools/")
library(GeneStructureTools)

load("source_data/orfs_normal.Rdata")
load("source_data/orfs_lnc.Rdata")
load("source_data/orfs_other.Rdata")
load("source_data/orfs_pseudo.Rdata")

orfs_all <- rbind(orfs_normal,orfs_lncRNA, orfs_other,orfs_pseudo)

load("source_data/upstream_orfs.Rdata")

m = match(orfs_all$id, uorfs_bytrans$id)
orfs_all = cbind(orfs_all, uorfs_bytrans[m,-c(1)])

for(i in 17:ncol(orfs_all)){
    orfs_all[which(is.na(orfs_all[,i])),i] <- 0
}

orfs_all$orf_length <- as.numeric(orfs_all$orf_length)
orfs_all$gene_id <- gtf.exons$gene_id[match(orfs_all$id, gtf.exons$transcript_id)]
orfs_all$transcript_type <- gtf.exons$transcript_type[match(orfs_all$id, gtf.exons$transcript_id)]
orfs_all$exons <- exons_per_transcripts$Freq[match(orfs_all$id, exons_per_transcripts$Var1)]
orfs_all$exon_a_from_start = orfs_all$exon_a_from_start +1
orfs_all$transcript_type2 <- orfs_all$transcript_type
orfs_all$transcript_type2[orfs_all$transcript_type2 %in% c("3prime_overlapping_ncrna", "antisense", "lincRNA","sense_intronic","sense_overlapping")] <- "lncRNA"
orfs_all$transcript_type2[orfs_all$transcript_type2 %in% c("processed_pseudogene", "transcribed_processed_pseudogene",
                                                           "transcribed_unprocessed_pseudogene","translated_processed_pseudogene","unitary_pseudogene","unprocessed_pseudogene")] <- "pseudogene"
orfs_all$transcript_type2[orfs_all$transcript_type2 %in% c("IG_V_pseudogene", "transcribed_processed_pseudogene",
                                                           "polymorphic_pseudogene","TR_V_pseudogene")] <- "remove"


orfs_train <- orfs_all[orfs_all$transcript_type == "protein_coding" | orfs_all$transcript_type == "nonsense_mediated_decay",]

#longest ORFs for protein coding genes
orfs_train <- arrange(orfs_train, plyr::desc(orf_length))
orfs_train <- orfs_train[!duplicated(orfs_train$id),]

rm <- which(apply(orfs_train, 1, function(x) any(is.na(x))))
if(length(rm) > 0){
    orfs_train <- orfs_train[-rm,]
}

index_1 <- which(orfs_train$transcript_type == "nonsense_mediated_decay")
index_2 <- which(orfs_train$transcript_type != "nonsense_mediated_decay")

set.seed(1)

train_index <- c(sample(index_1, 4500),
                 sample(index_2, 4500*5))

test_index <- c(sample(index_1[!(index_1 %in% train_index)], 1000),
                sample(index_2[!(index_2 %in% train_index)],1000))

train_ids <- orfs_train$id[train_index]
test_ids <- orfs_train$id[test_index]

train <- orfs_train[train_index,]
train$id <- NULL
train$orf_sequence <- NULL
train$start_site <- NULL
train$stop_site <- NULL
train$frame <- NULL
train$gene_id <- NULL
train$transcript_type2 <- NULL
train$seq_length <- NULL
train$exons=NULL
train$orf_length=NULL
train$seq_length_nt=NULL
train$total_uorfs=NULL

test <- orfs_train[test_index,]
test$orf_sequence <- NULL
test$Class <- "nonsense_mediated_decay"
test$Class[test$transcript_type != "nonsense_mediated_decay"] <- "not_nmd"

preProcValues <- preProcess(train, method = c("center", "scale"))
train_p <- predict(preProcValues, train)
train_p$transcript_type[train_p$transcript_type != "nonsense_mediated_decay"] <- "not_nmd"
test_p <- predict(preProcValues, test)

fitControl <- trainControl(method = "repeatedcv",number = 3,repeats = 10,classProbs=TRUE)

model_weights <- ifelse(train_p$transcript_type == "nonsense_mediated_decay",
                        (1/table(train_p$transcript_type)[1])*0.5,
                        (1/table(train_p$transcript_type)[2])*0.5)
model_gbm.long <- train(transcript_type ~ .,
                        data=train_p,
                        method = "gbm",
                        trControl=fitControl,
                        weights=model_weights,
                        verbose=F)

preProcValues.long <- preProcValues
save(model_gbm.long, train_ids, test_ids, preProcValues.long, file="source_data/model_gbm_long.RData")

############### +lncrna model ##################

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
train$seq_length <- NULL
train$exons=NULL
train$orf_length=NULL
train$seq_length_nt=NULL
train$total_uorfs=NULL

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

model_gbm.long.lnc <- train(transcript_type ~ .,
                            data=train_p,
                            method = "gbm",
                            trControl=fitControl,
                            weights=model_weights,
                            verbose=F)

p <- predict(model_gbm.long.lnc, test_p, type="prob")[,1]
test$prob=p
n <- predict(model_gbm.long.lnc, test_p)
c_long.lnc <- confusionMatrix(test$Class, n)

save(model_gbm.long.lnc, preProcValues.long.lnc, c_long.lnc, file="source_data/model_gbm_long.lnc.RData")
test.base <- test[,c(1,2,21,24,25)]
write.csv(test.base, "source_data/base_model_performance.csv", row.names = F,quote=F)

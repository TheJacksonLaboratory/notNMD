library(plyr)
library(GenomicRanges)
library(caret)
devtools::install("~/Documents/Projects/GeneStructureTools/")
library(GeneStructureTools)

# Gencode annotations
gtf <- rtracklayer::import("~/Documents/Projects/resources/genomes/gencode.v21.annotation.gtf")
gtf.exons <- gtf[gtf$type=="exon"]
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# protein coding / nmd with TSL > 3
normal_transcripts <- gtf.exons[which(gtf.exons$transcript_type %in%
                                   c("protein_coding", "nonsense_mediated_decay") &
                                   gtf.exons$transcript_support_level < 3)]

orfs_normal <- getOrfs(normal_transcripts, g, returnLongestOnly = FALSE,all_frames = TRUE)

save(orfs_normal, normal_transcripts, file = "orfs_normal.Rdata")
load("orfs_normal.Rdata")

lncRNA_transcripts <- gtf.exons[which(gtf.exons$transcript_type %in%
                                          c("lincRNA", "antisense", "sense_intronic","sense_overlapping") &
                                          gtf.exons$transcript_support_level < 3)]

orfs_lncRNA <- getOrfs(lncRNA_transcripts, g, returnLongestOnly = FALSE,all_frames = TRUE)
save(orfs_lncRNA, lncRNA_transcripts, file = "orfs_lnc.Rdata")

orfs_normal <- rbind(orfs_lncRNA, orfs_normal)

orfs_normal$orf_length <- as.numeric(orfs_normal$orf_length)
orfs_normal$gene_id <- gtf.exons$gene_id[match(orfs_normal$id, gtf.exons$transcript_id)]
orfs_normal$transcript_type <- gtf.exons$transcript_type[match(orfs_normal$id, gtf.exons$transcript_id)]

orfs_pc <- orfs_normal[orfs_normal$transcript_type != "nonsense_mediated_decay",]
orfs_nmd_long <- orfs_normal[orfs_normal$transcript_type == "nonsense_mediated_decay" & orfs_normal$orf_length >= 100,]
orfs_nmd_short <- orfs_normal[orfs_normal$transcript_type == "nonsense_mediated_decay" & orfs_normal$orf_length < 100,]
orfs_nmd <- rbind(orfs_nmd_long, orfs_nmd_short)

#longest ORFs for protein coding genes
orfs_pc <- arrange(orfs_pc, plyr::desc(orf_length))
orfs_pc <- orfs_pc[!duplicated(orfs_pc$id),]

orfs_train <- rbind(orfs_pc, orfs_nmd)
rm <- which(apply(orfs_train, 1, function(x) any(is.na(x))))
if(length(rm) > 0){
    orfs_train <- orfs_train[-rm,]
}
index_1 <- which(orfs_train$transcript_type == "nonsense_mediated_decay")
index_2 <- which(orfs_train$transcript_type != "nonsense_mediated_decay")

set.seed(1)

train_index <- c(sample(index_1, 5000),
                 sample(index_2, 5000))

test_index <- c(sample(index_1[!(index_1 %in% train_index)], 5000),
                sample(index_2[!(index_2 %in% train_index)], 5000))

train <- orfs_train[train_index,]
train$id <- NULL
train$orf_sequence <- NULL
train$start_site <- NULL
train$stop_site <- NULL
train$frame <- NULL
train$gene_id <- NULL

test <- orfs_train[test_index,]
test$orf_sequence <- NULL
test$Class <- "nonsense_mediated_decay"
test$Class[test$transcript_type != "nonsense_mediated_decay"] <- "not_nmd"

preProcValues <- preProcess(train, method = c("center", "scale"))

train_p <- predict(preProcValues, train)
train_p$transcript_type[train_p$transcript_type != "nonsense_mediated_decay"] <- "not_nmd"

test_p <- predict(preProcValues, test)

fitControl <- trainControl(method = "repeatedcv",number = 3,repeats = 10,classProbs=TRUE)
model_gbm <- train(transcript_type ~ ., data=train_p, method = "gbm", trControl=fitControl, verbose=F)

p <- predict(model_gbm, test_p, type="prob")[,1]
plot(density(p))

n <- predict(model_gbm, test_p)

k <- which(test$orf_length > 100)
confusionMatrix(n[k], test$Class[k])
confusionMatrix(n, test$Class)

notNMDnames <- colnames(train_p)
notNMDnames <- notNMDnames[-which(notNMDnames == "transcript_type")]

devtools::use_data(preProcValues, model_gbm, notNMDnames, overwrite = TRUE)
save(preProcValues, model_gbm, notNMDnames, file="R/sysdata.rda")

#
cutoffs <- seq(0.5, 1, 0.01)
c_list <- list()
for(i in seq_along(cutoffs)){
    bottom_c <- 1-cutoffs[i]
    keep <- which(p <= bottom_c | p > cutoffs[i])

    c_list[[i]] <- confusionMatrix(n[keep], test$Class[keep])
    message(i)
}

plot(unlist(lapply(c_list, function(x) as.numeric(x$overall[1]))),unlist(lapply(c_list, function(x) sum(x$table))))


#reset not nmdscore to 0 for

keep <- which(gtf.exons$transcript_type[match(orfs_normal$id, gtf.exons$transcript_id)] %in% c("nonsense_mediated_decay","protein_coding"))

orfs_normal$orf_length <- as.numeric(orfs_normal$orf_length)
orfs_normal$gene_id <- gtf.exons$gene_id[match(orfs_normal$id, gtf.exons$transcript_id)]
orfs_normal$transcript_type <- gtf.exons$transcript_type[match(orfs_normal$id, gtf.exons$transcript_id)]

orfs_pc <- orfs_normal[orfs_normal$transcript_type == "protein_coding",]
orfs_nmd_long <- orfs_normal[orfs_normal$transcript_type == "nonsense_mediated_decay" & orfs_normal$orf_length >= 100,]
orfs_nmd_short <- orfs_normal[orfs_normal$transcript_type == "nonsense_mediated_decay" & orfs_normal$orf_length < 100,]
orfs_nmd <- rbind(orfs_nmd_long, orfs_nmd_short)
orfs_nmd <- arrange(orfs_nmd, plyr::desc(orf_length))
orfs_nmd <- orfs_nmd[!duplicated(orfs_nmd$id),]

#longest ORFs for protein coding genes
orfs_pc <- arrange(orfs_pc, plyr::desc(orf_length))
orfs_pc <- orfs_pc[!duplicated(orfs_pc$id),]

orfs_train <- rbind(orfs_pc, orfs_nmd)
rm <- which(apply(orfs_train, 1, function(x) any(is.na(x))))
orfs_train <- orfs_train[-rm,]

index_1 <- which(orfs_train$transcript_type == "nonsense_mediated_decay")
index_2 <- which(orfs_train$transcript_type != "nonsense_mediated_decay")

set.seed(1)

train_index <- c(sample(index_1, 2000),
                 sample(index_2, 2000))

test_index <- c(sample(index_1[!(index_1 %in% train_index)], 2000),
                sample(index_2[!(index_2 %in% train_index)], 2000))

train <- orfs_train[train_index,]
train$id <- NULL
train$orf_sequence <- NULL
train$start_site <- NULL
train$stop_site <- NULL
train$frame <- NULL
train$gene_id <- NULL

test <- orfs_train[test_index,]
test$orf_sequence <- NULL
test$Class <- "nonsense_mediated_decay"
test$Class[test$transcript_type != "nonsense_mediated_decay"] <- "not_nmd"

preProcValues <- preProcess(train, method = c("center", "scale"))

train_p <- predict(preProcValues, train)
train_p$transcript_type[train_p$transcript_type != "nonsense_mediated_decay"] <- "not_nmd"

test_p <- predict(preProcValues, test)

fitControl <- trainControl(method = "repeatedcv",number = 3,repeats = 10,classProbs=TRUE)
model_gbm <- train(transcript_type ~ ., data=train_p, method = "gbm", trControl=fitControl, verbose=F)

p <- predict(model_gbm, test_p, type="prob")[,1]
plot(density(p))

n <- predict(model_gbm, test_p)

k <- which(test$orf_length > 100)
confusionMatrix(n[k], test$Class[k])
confusionMatrix(n, test$Class)

varImp(model_gbm)

n2 <- n
n2[which(test$min_dist_to_junction_b < 50 & test$exon_b_from_final != 0)] <- "not_nmd"
n2[which(test$min_dist_to_junction_b > 50 & test$exon_b_from_final != 0)] <- "nonsense_mediated_decay"

confusionMatrix(n2[test$exon_b_from_final !=0], test$Class[test$exon_b_from_final !=0])
confusionMatrix(n[test$exon_b_from_final ==0], test$Class[test$exon_b_from_final ==0])



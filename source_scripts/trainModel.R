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

orfs_normal <- get_orfs(normal_transcripts, g, returnLongestOnly = FALSE,all_frames = TRUE)

save(orfs_normal, normal_transcripts, file = "orfs_normal.Rdata")
load("orfs_normal.Rdata")


orfs_normal$orf_length <- as.numeric(orfs_normal$orf_length)
orfs_normal$gene_id <- gtf.exons$gene_id[match(orfs_normal$id, gtf.exons$transcript_id)]
orfs_normal$transcript_type <- gtf.exons$transcript_type[match(orfs_normal$id, gtf.exons$transcript_id)]

orfs_pc <- orfs_normal[orfs_normal$transcript_type == "protein_coding",]
orfs_nmd_long <- orfs_normal[orfs_normal$transcript_type == "nonsense_mediated_decay" & orfs_normal$orf_length >= 100,]
orfs_nmd_short <- orfs_normal[orfs_normal$transcript_type == "nonsense_mediated_decay" & orfs_normal$orf_length < 100,]
orfs_nmd <- rbind(orfs_nmd_long, orfs_nmd_short)

#longest ORFs for protein coding genes
orfs_pc <- arrange(orfs_pc, plyr::desc(orf_length))
orfs_pc <- orfs_pc[!duplicated(orfs_pc$id),]

orfs_train <- rbind(orfs_pc, orfs_nmd)
rm <- which(apply(orfs_train, 1, function(x) any(is.na(x))))
orfs_train <- orfs_train[-rm,]

index_1 <- which(orfs_train$transcript_type == "nonsense_mediated_decay")
index_2 <- which(orfs_train$transcript_type == "protein_coding")

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

preProcValues <- preProcess(train, method = c("center", "scale"))

train_p <- predict(preProcValues, train)
test_p <- predict(preProcValues, test)

fitControl <- trainControl(method = "repeatedcv",number = 3,repeats = 10,classProbs=TRUE)
model_gbm <- train(transcript_type ~ ., data=train_p, method = "gbm", trControl=fitControl, verbose=F)

p <- predict(model_gbm, test_p, type="prob")[,1]
plot(density(p))

n <- predict(model_gbm, test_p)

k <- which(test$orf_length > 100)
confusionMatrix(n[k], test$transcript_type[k])
confusionMatrix(n, test$transcript_type)

notNMDnames <- colnames(train_p)
notNMDnames <- notNMDnames[-which(notNMDnames == "transcript_type")]

devtools::use_data(preProcValues, model_gbm, notNMDnames, overwrite = TRUE)

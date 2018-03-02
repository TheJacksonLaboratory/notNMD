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


################################################################################

orfs_testall <- orfs_all

#longest ORFs for protein coding genes
orfs_testall <- arrange(orfs_testall, plyr::desc(orf_length))
orfs_testall <- orfs_testall[!duplicated(orfs_testall$id),]

testall_p <- predict(preProcValues.long.lnc, orfs_testall)
p <- predict(model_gbm.long.lnc, testall_p, type="prob")[,1]
orfs_testall$prob.lnc=p
n <- predict(model_gbm.long.lnc, testall_p)
orfs_testall$pred_class.lnc=n

testall_p <- predict(preProcValues.long, orfs_testall)
p <- predict(model_gbm.long, testall_p, type="prob")[,1]
orfs_testall$prob=p
n <- predict(model_gbm.long, testall_p)
orfs_testall$pred_class=n

orfs_testall$inTestTrain = "no"
orfs_testall$inTestTrain[orfs_testall$id %in% c(train_ids, train_ids.lnc, test_ids, test_ids.lnc)] <- "yes"
orfs_testall$inTestTrain[orfs_testall$id %in% c(train_ids, train_ids.lnc)] <- "train"
orfs_testall$inTestTrain[orfs_testall$id %in% c(train_ids, train_ids.lnc)] <- "test"

table(orfs_testall$transcript_type2, orfs_testall$pred_class.lnc)
table(orfs_testall$transcript_type2, orfs_testall$pred_class)


library(RColorBrewer)

palette_type2 <- brewer.pal(n=12, "Set3")[c(1,6,8,5,10,12,4)]

table(orfs_testall$transcript_type2)
orfs_testall$transcript_type2_factor <- gsub("_"," ",orfs_testall$transcript_type2)
orfs_testall$transcript_type2_factor <- factor(orfs_testall$transcript_type2_factor, levels=c("retained intron","processed transcript","pseudogene","remove",
                                                                                       "protein coding","lncRNA","nonsense mediated decay"))

orfs_testall$used_in_model <- gsub("no" , "not used in model training/testing" , gsub("yes", "used in model training/testing", orfs_testall$inTestTrain))

ggplot(orfs_testall[orfs_testall$transcript_type2 != "remove",],
       aes(x=prob, y=prob.lnc, col=transcript_type2_factor)) + geom_point(alpha=0.5, stroke=0) +
    facet_wrap(~inTestTrain) + theme_figure + scale_color_manual(values=palette_type2[c(7,3,5,4,1,2)]) +
    scale_y_continuous(name="notNMD (with lncRNA) probability") +
    scale_x_continuous(name="notNMD probability")

ggplot(orfs_testall[orfs_testall$transcript_type2 != "remove",],
       aes(x=prob, y=prob.lnc, col=transcript_type2_factor)) + geom_point(alpha=0.5, stroke=0) +
    facet_wrap(~transcript_type2_factor) + theme_figure + scale_color_manual(values=palette_type2[c(7,3,5,4,1,2)]) +
    scale_y_continuous(name="notNMD (with lncRNA) probability") +
    scale_x_continuous(name="notNMD probability")

p_classed <- table(orfs_testall$transcript_type[orfs_testall$inTestTrain=="no"], orfs_testall$pred_class[orfs_testall$inTestTrain=="no"])
percent_classed <- as.data.frame(p_classed/rowSums(p_classed))
percent_classed$n <- as.data.frame(p_classed)[,3]

percent_classed <- arrange(percent_classed, Var2, Freq)
percent_classed$Var1 <- gsub("_", " ",percent_classed$Var1)
percent_classed$Var1_fac <- factor(percent_classed$Var1, levels=percent_classed$Var1[percent_classed$Var2=="nonsense_mediated_decay"])

keep <- aggregate(n ~ Var1, percent_classed, sum)
keep <- keep$Var1[keep$n >50]


plot.transcript_types <- ggplot(percent_classed[percent_classed$Var1 %in% keep,], aes(x=Var1_fac, y=Freq, fill=Var2)) +
    geom_bar(stat="identity", aes(alpha=log10(n)))+
    geom_bar(stat="identity", aes(alpha=0, col=Var2)) +
    geom_text(aes(label=n, y=ifelse(Var2=="not_nmd", 0.02, 0.98), hjust=ifelse(Var2=="not_nmd", 0, 1)), size=2) +
    theme_figure + coord_flip() +
    scale_y_continuous(name="proportion of total transcripts") +
    scale_x_discrete(name="transcript biotype") +
    scale_fill_manual(values=palette_type2[c(2,4)], name="predicted\nclass", labels=c("NMD","not NMD")) +
    scale_color_manual(values=palette_type2[c(2,4)], name="predicted\nclass", labels=c("NMD","not NMD"))

################################################################################

p_classed.lnc <- table(orfs_testall$transcript_type[orfs_testall$inTestTrain=="no"], orfs_testall$pred_class.lnc[orfs_testall$inTestTrain=="no"])
percent_classed.lnc <- as.data.frame(p_classed.lnc/rowSums(p_classed.lnc))
percent_classed.lnc$n <- as.data.frame(p_classed.lnc)[,3]

percent_classed.lnc <- arrange(percent_classed.lnc, Var2, Freq)
percent_classed.lnc$Var1 <- gsub("_", " ",percent_classed.lnc$Var1)
percent_classed.lnc$Var1_fac <- factor(percent_classed.lnc$Var1, levels=percent_classed$Var1[percent_classed$Var2=="nonsense_mediated_decay"])

keep <- aggregate(n ~ Var1, percent_classed.lnc, sum)
keep <- keep$Var1[keep$n >50]


plot.transcript_types.lnc <- ggplot(percent_classed.lnc[percent_classed.lnc$Var1 %in% keep,], aes(x=Var1_fac, y=Freq, fill=Var2)) +
    geom_bar(stat="identity", aes(alpha=log10(n)))+
    geom_bar(stat="identity", aes(alpha=0, col=Var2)) +
    geom_text(aes(label=n, y=ifelse(Var2=="not_nmd", 0.02, 0.98), hjust=ifelse(Var2=="not_nmd", 0, 1)), size=2) +
    theme_figure + coord_flip() +
    scale_y_continuous(name="proportion of total transcripts") +
    scale_x_discrete(name="transcript biotype") +
    scale_fill_manual(values=palette_type2[c(2,4)], name="predicted\nclass", labels=c("NMD","not NMD")) +
    scale_color_manual(values=palette_type2[c(2,4)], name="predicted\nclass", labels=c("NMD","not NMD"))


library(cowplot)

ggdraw()+
    draw_plot(plot.transcript_types + theme(legend.position = "none"), x=0,y=0,width=0.55, height=1) +
    draw_plot(plot.transcript_types.lnc + theme(axis.title.y = element_blank(), axis.text.y = element_blank()), x=0.55,y=0,width=0.45, height=1)


################################################################################

variable_importance <- varImp(model_gbm.long, scale = F)$importance
variable_importance <- as.data.frame(variable_importance)
variable_importance$variable <- rownames(variable_importance)

variable_importance.lnc <- varImp(model_gbm.long.lnc, scale=F)$importance
variable_importance$lnc_imp <- variable_importance.lnc$Overall[match(variable_importance$variable, rownames(variable_importance.lnc))]


ggplot(variable_importance, aes(x=Overall, y=lnc_imp, label=variable)) +
    geom_point() + geom_label() + theme_figure + geom_abline(slope=1)

################################################################################
orfs_testall$prob_set <- 0.1
orfs_testall$prob_set[orfs_testall$prob > 0.1] <- 0.2
orfs_testall$prob_set[orfs_testall$prob > 0.2] <- 0.3
orfs_testall$prob_set[orfs_testall$prob > 0.3] <- 0.4
orfs_testall$prob_set[orfs_testall$prob > 0.4] <- 0.5
orfs_testall$prob_set[orfs_testall$prob > 0.5] <- 0.6
orfs_testall$prob_set[orfs_testall$prob > 0.6] <- 0.7
orfs_testall$prob_set[orfs_testall$prob > 0.7] <- 0.8
orfs_testall$prob_set[orfs_testall$prob > 0.8] <- 0.9
orfs_testall$prob_set[orfs_testall$prob > 0.9] <- 1.0


orfs_testall$prob_set.lnc <- 0.1
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.1] <- 0.2
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.2] <- 0.3
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.3] <- 0.4
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.4] <- 0.5
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.5] <- 0.6
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.6] <- 0.7
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.7] <- 0.8
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.8] <- 0.9
orfs_testall$prob_set.lnc[orfs_testall$prob.lnc > 0.9] <- 1.0




orfs_testall$exon_b_from_final_set <- "stop in final exon"
orfs_testall$exon_b_from_final_set[orfs_testall$exon_b_from_final > 0] <- "stop not in final exon"


ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding"),],
       aes(x=factor(prob_set), y=orf_length, fill=transcript_type2)) + geom_boxplot() + scale_y_log10() +theme_figure

ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding"),],
       aes(x=factor(prob_set), y=start_site_nt, fill=transcript_type2)) + geom_boxplot() + scale_y_log10() +theme_figure



library(virdis)


ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding")& orfs_testall$inTestTrain!="yes",],
       aes(x=factor(prob_set),fill=factor(exon_b_from_final))) + geom_bar()  +
    scale_fill_manual(values=(viridis(75,option="magma"))[c(31, 41,51:75)], name="stop codon exon\nfrom final exon") +
    theme_figure +
    theme(axis.text.x = element_text(hjust=1, angle=90)) +
    scale_x_discrete(name="notNMD probability score", labels=c(paste0("0.", c(0:8), "-0.", c(1:9)), "0.9-1.0")) +
    facet_wrap(~transcript_type2_factor, scales="free_y")

ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding") & orfs_testall$inTestTrain!="yes",],
       aes(x=factor(prob_set), y=start_site_nt)) +
    theme_figure +
    scale_color_manual(values=palette_type2[c(4,2)], name="transcript biotype") +
    theme(axis.text.x = element_text(hjust=1, angle=90)) +
    scale_y_log10(name="5'UTR length") +
    geom_jitter(alpha=0.2,aes(col=transcript_type2_factor)) + geom_boxplot(alpha=0) +
    scale_x_discrete(name="notNMD probability score", labels=c(paste0("0.", c(0:8), "-0.", c(1:9)), "0.9-1.0")) +
    facet_wrap(~exon_b_from_final_set)

ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding") & orfs_testall$inTestTrain!="yes",],
       aes(x=factor(prob_set), y=min_dist_to_junction_b)) +
    theme_figure +
    scale_color_manual(values=palette_type2[c(4,2)], name="transcript biotype") +
    theme(axis.text.x = element_text(hjust=1, angle=90)) +
    scale_y_log10(name="junction B distance") +
    geom_jitter(alpha=0.2,aes(col=transcript_type2_factor)) + geom_boxplot(alpha=0) +
    geom_hline(yintercept = 50, linetype=2, col="grey30") +
    scale_x_discrete(name="notNMD probability score", labels=c(paste0("0.", c(0:8), "-0.", c(1:9)), "0.9-1.0")) +
    facet_wrap(~exon_b_from_final_set)

################################################################################

ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding", "lncRNA")& orfs_testall$inTestTrain!="yes",],
       aes(x=factor(prob_set.lnc),fill=factor(exon_b_from_final))) + geom_bar()  +
    scale_fill_manual(values=(viridis(75,option="magma"))[c(31, 41,51:75)], name="stop codon exon\nfrom final exon") +
    theme_figure +
    theme(axis.text.x = element_text(hjust=1, angle=90)) +
    scale_x_discrete(name="notNMD probability score", labels=c(paste0("0.", c(0:8), "-0.", c(1:9)), "0.9-1.0")) +
    facet_wrap(~transcript_type2_factor, scales="free_y")

ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding", "lncRNA") & orfs_testall$inTestTrain!="yes",],
       aes(x=factor(prob_set.lnc), y=start_site_nt)) +
    theme_figure +
    scale_color_manual(values=palette_type2[c(4,1,2)], name="transcript biotype") +
    theme(axis.text.x = element_text(hjust=1, angle=90)) +
    scale_y_log10(name="5'UTR length") +
    geom_jitter(alpha=0.2,aes(col=transcript_type2_factor)) + geom_boxplot(alpha=0) +
    scale_x_discrete(name="notNMD probability score", labels=c(paste0("0.", c(0:8), "-0.", c(1:9)), "0.9-1.0")) +
    facet_wrap(~exon_b_from_final_set)

ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding","lncRNA") & orfs_testall$inTestTrain!="yes",],
       aes(x=factor(prob_set.lnc), y=min_dist_to_junction_b)) +
    theme_figure +
    scale_color_manual(values=palette_type2[c(4,1,2)], name="transcript biotype") +
    theme(axis.text.x = element_text(hjust=1, angle=90)) +
    scale_y_log10(name="junction B distance") +
    geom_jitter(alpha=0.2,aes(col=transcript_type2_factor)) + geom_boxplot(alpha=0) +
    geom_hline(yintercept = 50, linetype=2, col="grey30") +
    scale_x_discrete(name="notNMD probability score", labels=c(paste0("0.", c(0:8), "-0.", c(1:9)), "0.9-1.0")) +
    facet_wrap(~exon_b_from_final_set)

################################################################################














ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding") & orfs_testall$inTestTrain!="yes",],
       aes(x=factor(prob_set), y=utr3_length)) +
    theme_figure +
    scale_color_manual(values=palette_type2[c(4,2)], name="transcript biotype") +
    theme(axis.text.x = element_text(hjust=1, angle=90)) +
    scale_y_log10(name="3'UTR length") +
    geom_jitter(alpha=0.2,aes(col=transcript_type2_factor)) + geom_boxplot(alpha=0) +
    scale_x_discrete(name="notNMD probability score", labels=c(paste0("0.", c(0:8), "-0.", c(1:9)), "0.9-1.0")) +
    facet_wrap(~exon_b_from_final_set)

ggplot(orfs_testall[orfs_testall$transcript_type2 %in% c("nonsense_mediated_decay", "protein_coding"),],
       aes(fill=prob_set, x=exon_b_from_final)) + geom_bar(aes(fill=prob_set))

ggplot(orfs_testall, aes(x=prob_set, fill=factor(exon_b_from_final))) + geom_bar()


################################################################################
# cutoffs --> actual preds

load("source_data/core_model_data.Rdata")

cutoff_c <- list()
vals <- seq(0,1, by=0.01)
for(v in seq_along(vals)){
    test$pred_class <- ifelse(test$prob > vals[v], "nonsense_mediated_decay", "not_nmd")
    cutoff_c[[v]] <- confusionMatrix(test$pred_class, test$Class)
}

cutoff_df <- as.data.frame(matrix(unlist(lapply(cutoff_c, function(x) c(x$overall, x$byClass))), byrow = T, ncol=18))
colnames(cutoff_df) <- gsub(" ", "_", names(c(cutoff_c[[1]]$overall, cutoff_c[[1]]$byClass)))
cutoff_df$val <- vals

ggplot(cutoff_df, aes(x=val, y=Accuracy)) + geom_point()
ggplot(cutoff_df, aes(x=val, y=Pos_Pred_Value)) + geom_point()
ggplot(cutoff_df, aes(x=val, y=Neg_Pred_Value)) + geom_point()

### min max cutoffs for changes
min(cutoff_df$val[cutoff_df$Pos_Pred_Value > 0.9], na.rm=TRUE)
min(cutoff_df$val[cutoff_df$Pos_Pred_Value > 0.95], na.rm=TRUE)

max(cutoff_df$val[cutoff_df$Neg_Pred_Value > 0.9], na.rm=TRUE)
max(cutoff_df$val[cutoff_df$Neg_Pred_Value > 0.95], na.rm=TRUE)




cutoff_c[[50]]
cutoff_c[[100]]

plot(unlist(lapply(cutoff_c, function(x) x$overall["Accuracy"])))

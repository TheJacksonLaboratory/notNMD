library(plyr)
library(GenomicRanges)
library(caret)
devtools::install_github("betsig/GeneStructureTools/")
library(GeneStructureTools)
library(gbm)

# Gencode annotations
gtf <- rtracklayer::import("gencode.v21.annotation.gtf")
gtf.exons <- gtf[gtf$type=="exon"]
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# protein coding / nmd with TSL > 3
normal_transcripts <- gtf.exons[which(gtf.exons$transcript_type %in%
                                        c("protein_coding", "nonsense_mediated_decay") &
                                        gtf.exons$transcript_support_level < 3)]

orfs_normal <- getOrfs(normal_transcripts, g, returnLongestOnly = FALSE,allFrames = TRUE,uORFs = TRUE)

save(orfs_normal, normal_transcripts, file = "source_data/orfs_normal.Rdata")
load("source_data/orfs_normal.Rdata")

lncRNA_transcripts <- gtf.exons[which(gtf.exons$transcript_type %in%
                                        c("lincRNA", "antisense", "sense_intronic","sense_overlapping") &
                                        gtf.exons$transcript_support_level < 3)]

orfs_lncRNA <- getOrfs(lncRNA_transcripts, g, returnLongestOnly = FALSE,allFrames = TRUE,uORFs = TRUE)
save(orfs_lncRNA, lncRNA_transcripts, file = "source_data/orfs_lnc.Rdata")
load("source_data/orfs_lnc.Rdata")

# !lnc/pc/nmd with TSL > 3
other_transcripts <- gtf.exons[which(gtf.exons$transcript_type %in%
                                        c("3prime_overlapping_ncrna","polymorphic_pseudogene","processed_transcript",
                                          "retained_intron","unprocessed_pseudogene") &
                                        gtf.exons$transcript_support_level < 3)]

orfs_other <- getOrfs(other_transcripts, g, returnLongestOnly = FALSE,allFrames = T,uORFs = TRUE)
save(orfs_other, other_transcripts, file = "source_data/orfs_other.Rdata")
load("source_data/orfs_other.Rdata")

psuedo_names <- gtf.exons$transcript_id[which(grepl("pseudo", gtf.exons$transcript_type) &
                                                gtf.exons$level==1 & gtf.exons$exon_number==2)]

pseudo_transcripts <- gtf.exons[which(gtf.exons$transcript_id %in% psuedo_names)]

orfs_pseudo <- getOrfs(psuedo_transcripts, g, returnLongestOnly = FALSE,allFrames = TRUE,uORFs = TRUE)
save(orfs_pseudo, pseudo_transcripts, file = "source_data/orfs_pseudo.Rdata")
load("source_data/orfs_pseudo.Rdata")


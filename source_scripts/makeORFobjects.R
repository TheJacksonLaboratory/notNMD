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

orfs_normal <- getOrfs(normal_transcripts, g, returnLongestOnly = FALSE,allFrames = TRUE)

save(orfs_normal, normal_transcripts, file = "orfs_normal.Rdata")
load("orfs_normal.Rdata")

lncRNA_transcripts <- gtf.exons[which(gtf.exons$transcript_type %in%
                                        c("lincRNA", "antisense", "sense_intronic","sense_overlapping") &
                                        gtf.exons$transcript_support_level < 3)]

orfs_lncRNA <- getOrfs(lncRNA_transcripts, g, returnLongestOnly = FALSE,allFrames = TRUE)
save(orfs_lncRNA, lncRNA_transcripts, file = "orfs_lnc.Rdata")
load("orfs_lnc.Rdata")

# !lnc/pc/nmd with TSL > 3
other_transcripts <- gtf.exons[which(gtf.exons$transcript_type %in%
                                        c("3prime_overlapping_ncrna","polymorphic_pseudogene","processed_transcript",
                                          "retained_intron","unprocessed_pseudogene") &
                                        gtf.exons$transcript_support_level < 3)]

orfs_other <- getOrfs(other_transcripts, g, returnLongestOnly = FALSE,allFrames = T)
save(orfs_other, other_transcripts, file = "orfs_other.Rdata")
load("source_data/orfs_other.Rdata")

psuedo_names <- gtf.exons$transcript_id[which(grepl("pseudo", gtf.exons$transcript_type) &
                                                gtf.exons$level==1 & gtf.exons$exon_number==2)]

pseudo_transcripts <- gtf.exons[which(gtf.exons$transcript_id %in% psuedo_names)]

orfs_pseudo <- getOrfs(psuedo_transcripts, g, returnLongestOnly = FALSE,allFrames = TRUE)
save(orfs_pseudo, pseudo_transcripts, file = "orfs_pseudo.Rdata")
load("orfs_pseudo.Rdata")

orfs_all <- rbind(orfs_normal,orfs_lncRNA, orfs_other,orfs_pseudo)

# annotate upstream ORFs

upstream_orfs = getUOrfs(c(normal_transcripts, lncRNA_transcripts,other_transcripts,pseudo_transcripts),
                         BSgenome = g,
                         orfs = orfs_all,
                         findExonB = TRUE)

uorfs_bytrans = aggregate(overlaps_main_ORF ~ id, upstream_orfs, function(x) length(x))
colnames(uorfs_bytrans)[2] = "total_uorfs"

uorfs_bytrans2 = aggregate(overlaps_main_ORF ~ id, upstream_orfs, function(x) length(x[which(x=="upstream")]))
uorfs_bytrans$upstream_count = uorfs_bytrans2[match(uorfs_bytrans$id, uorfs_bytrans2$id),2]

uorfs_bytrans2 = aggregate(overlaps_main_ORF ~ id, upstream_orfs, function(x) length(x[which(x=="downstream")]))
uorfs_bytrans$downstream_count = uorfs_bytrans2[match(uorfs_bytrans$id, uorfs_bytrans2$id),2]

uorfs_bytrans2 = aggregate(uorf_length ~ id, upstream_orfs, function(x) max(x))
uorfs_bytrans$max_uorf = uorfs_bytrans2[match(uorfs_bytrans$id, uorfs_bytrans2$id),2]

uorfs_bytrans2 = aggregate(min_dist_to_junction_b ~ id, upstream_orfs[upstream_orfs$exon_b_from_final !=0,], function(x) max(x))
uorfs_bytrans$uorf_maxb = uorfs_bytrans2[match(uorfs_bytrans$id, uorfs_bytrans2$id),2]

save(upstream_orfs, uorfs_bytrans, file="source_data/upstream_orfs.Rdata")
load("source_data/upstream_orfs.Rdata")


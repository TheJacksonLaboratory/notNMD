#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(GeneStructureTools)
library(notNMD)

load(args[1])

if(exists("orfsX")){
orfsX$nmd_prob <- predictNMD(orfsX, "prob")
orfsX$nmd_class <- predictNMD(orfsX)
}

if(exists("orfsY")){
orfsY$nmd_prob <- predictNMD(orfsY, "prob")
orfsY$nmd_class <- predictNMD(orfsY)
}

if(exists("orfAllGenes")){
    orfAllGenes$nmd_prob <- predictNMD(orfAllGenes, "prob")
    orfAllGenes$nmd_class <- predictNMD(orfAllGenes)
}
save.image(args[1])


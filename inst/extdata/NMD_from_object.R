#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(GeneStructureTools)
library(notNMD)

load(args[1])

if(is.null(useModel)){useModel="base"}

if(exists("orfsX")){
orfsX$nmd_prob <- predictNMD(orfsX, "prob", model = useModel)
orfsX$nmd_class <- predictNMD(orfsX, model = useModel)
}

if(exists("orfsY")){
orfsY$nmd_prob <- predictNMD(orfsY, "prob", model = useModel)
orfsY$nmd_class <- predictNMD(orfsY, model = useModel)
}

if(exists("orfAllGenes")){
    orfAllGenes$nmd_prob <- predictNMD(orfAllGenes, "prob", model = useModel)
    orfAllGenes$nmd_class <- predictNMD(orfAllGenes, model = useModel)
}
save.image(args[1])


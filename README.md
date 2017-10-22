# notNMD

**R package for prediction of nonsense mediated decay potential in transcripts**

## Installation
notNMD can be installed:
```
library(devtools)
install_github("betsig/notNMD")
```
After installation, the package can be loaded into R.
```
library(notNMD)
```

## Usage

Please check the [vignette](http://htmlpreview.github.io/?https://github.com/betsig/notNMD/blob/master/vignette.html) for usage details.

## Training and Performance

notNMD can be used for prediction of NMD potential in transcripts supplied as GRanges.
Scripts used for training are available in source_scripts/
Briefly, we used GeneStructureTools to annotate open reading frame details in transcripts classed as protein coding or lncRNAs (not NMD targets), and transcripts classed as nonsense mediated decay from Human Gencode v21 annotations.
We then trained a GBM model on a training set of these transcripts, and tested performance on an independant testing set:

![](https://github.com/betsig/notNMD/blob/master/performance.png)






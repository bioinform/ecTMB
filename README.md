# ecTMB: estimation and classification of TMB

ecTMB is a powerful and flexible statistical framework for TMB estimation and classification. It uses an explicit background mutation mdoel for more robust and consistent TMB prediction. The backgournd mutation model takes account of unknown as well as known mutational heterogeneous factors, including tri-nucleotide context, sample mutational burden, gene expression level and replication timing by utilization of a Bayesian framework. The discovery of three TMB-based subtypes, including one novel subtype TMB-extreme, enable ecTMB to classify samples to biological and clinically relavent TMB subtypes.

## Table of Contents
**[Dependency](#dependency)**<br>
**[Installation](#installation)**<br>
**[Example Usage](#example-usage)**<br>
**[License](#license)**<br>

## Dependency
ecTMB has been sucessfully tested on Intel(R) Xeon(R) CPU E5-2680 v4 Machine with 28 cores.

ecTMB import following R packages: ggplot2, limma, reshape2, dplyr, R6, MASS, GenomicRanges, data.table, parallel, mixtools

ecTMB also depends on bedtools 2.27.1 and R = 3.5.1

You can install these packages using [anaconda](https://www.anaconda.com/download)/[miniconda](https://conda.io/miniconda.html) :
```
conda install bedtools=2.27.1 r=3.5.1
```
Then you can export the conda paths as:
```
export PATH="/PATH/TO/CONDA/bin:$PATH"
export LD_LIBRARY_PATH="/PATH/TO/CONDA/lib:$LD_LIBRARY_PATH"
```

## Installation
```
install.packages("devtools")
library(devtools);
devtools::install_github("bioinform/ecTMB");
```
## Download Example and Reference Data
```
#Example file download from URL: https://www.dropbox.com/s/knpgl73samhdtvg/ecTMB_data.tar.gz?dl=1
URL <- "https://www.dropbox.com/s/knpgl73samhdtvg/ecTMB_data.tar.gz?dl=1"
download.file(URL,destfile = "ecTMB.example.tar.gz")
untar("./ecTMB.example.tar.gz")
```


## Example Usage
* **Load ecTMB package and genome annotation reference files**
```
library(ecTMB)
load("./ecTMB_data/example/UCEC.rda")
extdataDir             = "./ecTMB_data/references"
exomef                 = file.path(extdataDir, "exome_hg38_vep.Rdata" )  #### hg38 exome file
covarf                 = file.path(extdataDir,"gene.covar.txt")   ### gene properties
mutContextf            = file.path(extdataDir,"mutation_context_96.txt" )  ### 96 mutation contexts
TST170_panel           = file.path(extdataDir,"TST170_DNA_targets_hg38.bed" )  ### 96 mutation contexts
ref                    = file.path(extdataDir,"GRCh38.d1.vd1.fa" )

```
* **Set random 70% as training and rest as test set**
```
set.seed(1002200)
SampleID_all   = UCEC_cli$sample
SampleID_train = sample(SampleID_all, size = round(2 * length(SampleID_all)/3), replace = F)
SampleID_test  = SampleID_all[!SampleID_all %in% SampleID_train]
```
* **Generate train and test data object**
```
## mutations which are inconsistent with reference annotation files will be removed.
## train data
trainData      = UCEC_mafs[UCEC_mafs$Tumor_Sample_Barcode %in% as.character(SampleID_train),]
trainset       = readData(trainData, exomef, covarf, mutContextf, ref)

## test data for panel TST 170
sample         = data.frame(SampleID = SampleID_test, BED = TST170_panel, stringsAsFactors = FALSE)
testData       = UCEC_mafs[UCEC_mafs$Tumor_Sample_Barcode %in% as.character(SampleID_test),]
testset_panel  = readData(testData, exomef, covarf, mutContextf, ref, samplef = sample)
testset_WES    = readData(testData, exomef, covarf, mutContextf, ref)  ## to calculate WES-TMB for test samples
```
* **Background mutation model training** 
---
**NOTE**

This step takes up to ~12 mins when 24 parallel processes are used. You can skip 
and use the pre-loaded parameters defined from training data set.

---
```
MRtriProb_train= getBgMRtri(trainset)
trainedModel   = fit_model(trainset, MRtriProb_train, cores = 24)
```

* **Predict TMB for TST170 panel**
```
TMBs          = pred_TMB(testset_panel, WES = testset_WES,
                        params = trainedModel, mut.nonsil = T, gid_nonsil_p = trainset$get_nonsil_passengers(0.95))
                        
## plot the prediction.    
library(dplyr)
library(ggplot2)

TMBs %>% melt(id.vars = c("sample","WES_TMB")) %>% 
  ggplot( aes(x = WES_TMB, y = value,  color = factor(variable, levels = c("ecTMB_panel_TMB",  "count_panel_TMB")), 
              group = factor(variable))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='log2') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) +
  labs(x = "TMB defined by WES", y = sprintf("Predicted TMB from panel: TST170"))
```

* **Classify sample to 3 subtypes**
```
Sutypes      = assignClass(TMBs$ecTMB_panel_TMB, prior = prior_bs)
```

## License
ecTMB is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.



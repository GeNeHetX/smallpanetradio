

## code to prepare `DATASET` dataset goes here

# source("/Users/remy.nicolle/Workspace/PDAC_aux/TNERadio/smalltneradio/data-raw/DATASET.R")

#usethis::use_data(DATASET, overwrite = TRUE)
setwd("/Users/remy.nicolle/Workspace/PDAC_aux/TNERadio/smalltneradio")
library(devtools)

bothmodel <- readRDS(paste0("../radiopredmodels/bothmodel_pcaglmnet_0.55_20_1_TRUE.rds"))
artmodel <- readRDS(paste0("../radiopredmodels/artmodel_pcaglmnet_0.6_20_1_TRUE.rds"))
portmodel <- readRDS(paste0("../radiopredmodels/portmodel_pcaglmnet_0.6_20_1_TRUE.rds"))

refclin <-  readRDS(paste0("../radiopredmodels/cleanclindf.rds"))


# clindf <- readRDS("radiopredmodels/cleanclindf.rds")
refartdata <- readRDS("../radiopredmodels/cleanartdf.rds")
refportdata <- readRDS("../radiopredmodels/cleanportdf.rds")



lvar=list(both=rownames(bothmodel$PCA$var$coord),art=rownames(artmodel$PCA$var$coord),port=rownames(portmodel$PCA$var$coord))




MODELS=list(bothmodel=bothmodel,artmodel=artmodel,portmodel=portmodel,
bothvars=rownames(bothmodel$PCA$var$coord),
artvars=rownames(artmodel$PCA$var$coord),
portvars=rownames(portmodel$PCA$var$coord),
allvars=unique(unlist(lvar))
# unamevar=unique(sub("^port|^art","",allvar))
)

use_data(MODELS,refartdata,refportdata,refclin,overwrite = TRUE)

#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 12012017_Apodization_wrapper.R
## Manon Martin
## manon.martin@uclouvain.be


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Use of runExampleL ?
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

runExampleL <- FALSE
# 
# if(runExampleL) {
#   ##------------------------------
#   ## Example of arguments
#   ##------------------------------
#   argLs <- list(StudyDir = "Tlse_BPASourisCerveau",
#                 upper = "10.0",
#                 lower = "0.50",
#                 bucket.width = "0.01",
#                 exclusion = "TRUE",
#                 exclusion.zone = list(c(6.5,4.5)),
#                 graph="Overlay")
#   
#   argLs <- c(argLs,
#              list(dataMatrixOut = paste(directory,"_NmrBucketing_dataMatrix.tsv",sep=""),
#                   sampleMetadataOut = paste(directory,"_NmrBucketing_sampleMetadata.tsv",sep=""),
#                   variableMetadataOut = paste(directory,"_NmrBucketing_variableMetadata.tsv",sep=""),
#                   graphOut = paste(directory,"_NmrBucketing_graph.pdf",sep=""),
#                   logOut = paste(directory,"_NmrBucketing_log.txt",sep="")))
# }



##------------------------------
## Options
##------------------------------
strAsFacL <- options()$stringsAsFactors
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# stringsAsFactors = FALSE utilis?? o?? ??
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options(stringsAsFactors = FALSE)




##------------------------------
## Libraries loading
##------------------------------
# For parseCommandArgs function
library(batch) 

# R script call
source_local <- function(fname)
{
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}
#Import the different functions
source_local("Apodization_Manon.R")

##------------------------------
## Errors ?????????????????????
##------------------------------


##------------------------------
## Constants
##------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# interviennent ou ?
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

topEnvC <- environment()
flagC <- "\n"


##------------------------------
## Script
##------------------------------
if(!runExampleL)
  argLs <- parseCommandArgs(evaluate=FALSE) # If FALSE valeurs par d??faut d??finies dans le software


## Parameters Loading
##-------------------
## Inputs
# data

# --- Galaxy_preformatted data
# data <- read.table(argLs[["dataMatrix"]],check.names=FALSE,header=TRUE,sep="\t")
# rownames(data) <- data[,1]
# data <- data[,-1]
# --- 

# --- data from ReadFids
Fid_data <- read.table(argLs[["dataMatrixOut"]],check.names=FALSE,header=TRUE,sep="\t")
Fid_info <- read.table(argLs[["sampleOut"]],check.names=FALSE,header=TRUE,sep="\t")
# --- 

# other inputs (cf. XML wrapper)
DT =  argLs[["DT"]]
type.apod = argLs[["ApodizationMethod"]]

# set default values for optional arguments
phase=0
rectRatio=1/2
gaussLB=1
expLB=1

# change the default values
if (type.apod %in% c("cos2", "hanning", "hamming")) {
  phase = argLs[["phase"]]
} 

if (type.apod == "blockexp") {
  rectRatio = argLs[["rectRatio"]] 
  expLB = argLs[["expLB"]]
}

if (type.apod == "blockcos2") {
  rectRatio = argLs[["rectRatio"]]
}
                     
if (type.apod == "gauss") {
gaussLB = argLs[["gaussLB"]]
}

if (type.apod == "exp") {
  expLB = argLs[["expLB"]]
}

plotWindow = FALSE
returnFactor = FALSE



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Utility of Outputs ??
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Outputs
dataMatrixOut <- argLs[["dataMatrixOut"]] # Names from Saving



## Checking arguments
##-------------------
error.stock <- "\n"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# error.stock utilis?? ou ?
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
if(length(error.stock) > 1)
  stop(error.stock)


## Computation
##------------
outputs <- Apodization(Fid_data = Fid_data, Fid_info = Fid_info, DT = DT, 
                        type.apod = type.apod, phase = phase, rectRatio = rectRatio, 
                        gaussLB = gaussLB, expLB = expLB, plotWindow = plotWindow, returnFactor = returnFactor) 
  
data_matrix <- outputs # Data matrix



## Saving
##-------
# Data matrix
data_matrix <- cbind(rownames(data_matrix),data_matrix)
colnames(data_matrix) <- c("Sample",colnames(data_matrix)[-1])
write.table(data_matrix,file=argLs$dataMatrixOut,quote=FALSE,row.names=FALSE,sep="\t")



## Ending
##---------------------

cat("\nEnd of 'Apodization' Galaxy module call: ", as.character(Sys.time()), sep = "")

## sink(NULL)

options(stringsAsFactors = strAsFacL)

rm(list = ls())
#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 170116_NmrPreprocessing.R
## Manon Martin and Marie Tremblay-Franco

##======================================================
##======================================================
# Preamble
##======================================================
##======================================================

runExampleL <- FALSE

##------------------------------
## Options
##------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

##------------------------------
## Libraries laoding
##------------------------------
# For parseCommandArgs function
ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
# print(ip, row.names=FALSE)

library(batch)
library(ptw)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)


# R script call
source_local <- function(fname)
{
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	source(paste(base_dir, fname, sep="/"))
}
#Import the different functions
source_local("NmrPreprocessing_script.R")
source_local("DrawFunctions.R")

##------------------------------
## Errors ?????????????????????
##------------------------------


##------------------------------
## Constants
##------------------------------
topEnvC <- environment()
flagC <- "\n"


##------------------------------
## Script
##------------------------------
if(!runExampleL)
    argLs <- parseCommandArgs(evaluate=FALSE)

## Starting
cat("\nStart of 'Preprocessing' Galaxy module call: ", as.character(Sys.time()), sep = "")


##======================================================
##======================================================
## Parameters Loading
##======================================================
##======================================================

# 1rst order phase correction ------------------------
  # Inputs
	## Data matrix
Fid_data0 <- read.table(argLs[["dataMatrixFid"]],header=TRUE, check.names=FALSE, sep='\t')
Fid_data0 <- Fid_data0[,-1]
Fid_data0 <- as.matrix(Fid_data0)


	## Samplemetadata
samplemetadataFid <- read.table(argLs[["sampleMetadataFid"]],check.names=FALSE,header=TRUE,sep="\t")
samplemetadataFid <- as.matrix(samplemetadataFid)


# water and solvent(s) correction ------------------------
  # Inputs
lambda <- argLs[["lambda"]]
ptwSS1 <- argLs[["ptwSS"]]
ptwSS <- FALSE
if (ptwSS1=="YES")
 	ptwSS <- TRUE
	

# apodization -----------------------------------------
  # Inputs
phase=0
rectRatio=1/2
gaussLB=1
expLB=1
apodization <- argLs[["apodizationMethod"]]
if (apodization=='exp'){
  expLB <- argLs[["expLB"]]
  } else if (apodization=='cos2'){
  phase <- argLs[["phase"]]
  } else if (apodization=='hanning'){
  phase <- argLs[["phase"]]
  } else if (apodization=='hamming'){
  phase <- argLs[["phase"]]
  } else if (apodization=='blockexp'){
  rectRatio <- argLs[["rectRatio"]]
  expLB <- argLs[["expLB"]]
  } else if (apodization=='blockcos2'){
  rectRatio <- argLs[["rectRatio"]]
  } else if (apodization=='gauss'){
  rectRatio <- argLs[["rectRatio"]]
  gaussLB <- argLs[["gaussLB"]]
  }		


# Fourier transform ----------------------------------
  # Inputs
FTGraph <- argLs[["FTGraph"]]

# Shift referencering ----------------------------------
  # Inputs
ShiftReferencing <- argLs[["shift_referencing.choice"]]
print(ShiftReferencing)
if (ShiftReferencing=="YES")
{
	shiftTreshold = 2 # c
	ppmxaxis = TRUE
	shiftReferencingRangeList = NULL  # fromto.TMSP
	pctNear0 = 0.02 # pc 
	rowindex_graph = NULL
	shiftReferencingMethod <- argLs[["shiftReferencingMethod"]]
	if (shiftReferencingMethod == "thres")
	{
		shiftTreshold <- argLs[["shiftTreshold"]]
	}
	shiftReferencingRange <- argLs[["shiftReferencingRange"]]
	print(shiftReferencingRange)
	if (shiftReferencingRange == "near0"){
		pctNear0 <- argLs[["pctNear0"]]
	}
	 if (shiftReferencingRange == "window"){
			shiftReferencingRangeList <- NULL
			shiftReferencingRangeLeft <- argLs[["shiftReferencingRangeLeft"]]
			shiftReferencingRangeRight <- argLs[["shiftReferencingRangeRight"]]
			shiftReferencingRangeList <- list(shiftReferencingRangeList,c(shiftReferencingRangeLeft,shiftReferencingRangeRight))
	 }
	shiftHandling <- argLs[["shiftHandling"]]
	shiftReferencing <- argLs[["shiftReferencing"]]
	if (shiftReferencingRange == "near0"){
		  pctNear0 <- argLs[["pctNear0"]]
	}
		
	if (shiftReferencingRange == "window"){
		shiftReferencingRangeList <- NULL
		shiftReferencingRangeLeft <- argLs[["shiftReferencingRangeLeft"]]
		shiftReferencingRangeRight <- argLs[["shiftReferencingRangeRight"]]
		shiftReferencingRangeList <- list(shiftReferencingRangeList,c(shiftReferencingRangeLeft,shiftReferencingRangeRight))
	}
	shiftHandling <- argLs[["shiftHandling"]]
}

# Zero Order Phase Correction -------------------------------
  # Inputs
plot_rms = NULL
returnAngle = FALSE
createWindow = TRUE
Angle = NULL
plot_spectra = FALSE
quant = 0.95
freq = TRUE
fromto.0OPC = NULL
zeroOrderPhaseMethod <- argLs[["zeroOrderPhaseMethod"]]
if (zeroOrderPhaseMethod=='rms')
{
  quant <- argLs[["quant"]]
}
if (zeroOrderPhaseMethod=='max')
{
  angle <- argLs[["angle"]]
}

searchZoneZeroPhase <- argLs[["searchZoneZeroPhase.choice"]]
if (searchZoneZeroPhase == "YES") {
  searchZoneZeroPhaseList <- NULL
  searchZoneZeroPhaseLeft <- argLs[["shiftReferencingRangeLeft"]]
  searchZoneZeroPhaseRight <- argLs[["shiftReferencingRangeRight"]]
  searchZoneZeroPhaseList <- list(searchZoneZeroPhaseList,c(searchZoneZeroPhaseLeft,searchZoneZeroPhaseRight))
  }


# Baseline Correction -------------------------------
  # Inputs
ptwBc <- as.logical(argLs[["ptwBc"]])
maxIter <- argLs[["maxIter"]] 
lambdaBc <- argLs[["lambdaBc"]] 
pBc <- argLs[["pBc"]] 
epsilon <- argLs[["epsilon"]] 


# transformation of negative values -------------------------------
  # Inputs
NegativetoZero <- argLs[["NegativetoZero"]]


  # Outputs
nomGraphe <- argLs[["graphOut"]]
dataMatrixOut <- argLs[["dataMatrixOut"]]
log <- argLs[["logOut"]]



## Checking arguments
##-------------------
error.stock <- "\n"

if(length(error.stock) > 1)
  stop(error.stock)
  

##======================================================
##======================================================  
## Computation
##======================================================
##======================================================
pdf(nomGraphe, onefile = TRUE, width = 13, height = 13)

# FirstOrderPhaseCorrection ---------------------------------
Fid_data <- FirstOrderPhaseCorrection(Fid_data0, Fid_info = samplemetadataFid, group_delay = NULL)


# SolventSuppression ---------------------------------
Fid_data <- SolventSuppression(Fid_data, lambda.ss = lambda, ptw.ss = ptwSS, plotSolvent = F, returnSolvent = F)
	

# Apodization ---------------------------------	
Fid_data <- Apodization(Fid_data, Fid_info = samplemetadataFid, DT = NULL, 
                        type.apod = apodization, phase = phase, rectRatio = rectRatio, gaussLB = gaussLB, expLB = expLB, plotWindow = F, returnFactor = F)

# FourierTransform ---------------------------------
Spectrum_data <- FourierTransform(Fid_data, Fid_info = samplemetadataFid, reverse.axis = TRUE)

if (FTGraph == "YES") {
	title = "Fourier transformed spectra"
	DrawSignal(Spectrum_data, subtype = "stacked",
               ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
               xlab = "Frequency", num.stacked = 4, 
               main.title = title, createWindow=FALSE)
}


# InternalReferencing ---------------------------------
if (ShiftReferencing=="YES")
{
	Spectrum_data <- InternalReferencing(Spectrum_data, samplemetadataFid, method = shiftReferencingMethod, range = shiftReferencingRange, shiftHandling = shiftHandling,
										c = shiftTreshold, fromto.TMSP = shiftReferencingRangeList, pc = pctNear0)
	
	title = "Spectra after Shift referencing"
	DrawSignal(Spectrum_data, subtype = "stacked",
               ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
               xlab = "Frequency", num.stacked = 4, 
               main.title = title, createWindow=FALSE)
}


# ZeroOrderPhaseCorrection ---------------------------------
Spectrum_data  <- ZeroOrderPhaseCorrection(Spectrum_data, method = zeroOrderPhaseMethod, plot_rms = plot_rms, returnAngle = returnAngle, createWindow = createWindow, 
											Angle = Angle, plot_spectra = plot_spectra, quant = quant, freq = freq, fromto.0OPC = fromto.0OPC)

title = "Spectra after Zero Order Phase Correction"
DrawSignal(Spectrum_data, subtype = "stacked",
           ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
           xlab = "Frequency", num.stacked = 4, 
           main.title = title, createWindow=FALSE)


# BaselineCorrection ---------------------------------			
Spectrum_data <- BaselineCorrection(Spectrum_data, ptw.bc = ptwBc, maxIter = maxIter, lambda.bc = lambdaBc, p.bc = pBc, eps = epsilon, returnBaseline = F) 

title = "Spectra after Baseline Correction"
DrawSignal(Spectrum_data, subtype = "stacked",
			ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
			xlab = "Frequency", num.stacked = 4, 
			main.title = title, createWindow=FALSE)
 	

# NegativeValuesZeroing ---------------------------------
if (NegativetoZero=="YES") {
    Spectrum_data <- NegativeValuesZeroing(Spectrum_data)
    title = "Spectra after Negative Values Zeroing"
    DrawSignal(Spectrum_data, subtype = "stacked",
               ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
               xlab = "Frequency", num.stacked = 4, 
               main.title = title, createWindow=FALSE)
}

invisible(dev.off())


##======================================================
##======================================================
## Saving
##======================================================
##======================================================
# Data Matrix
write.table(Spectrum_data,file=argLs$dataMatrixOut, quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)


## Ending
cat("\nEnd of 'Preprocessing' Galaxy module call: ", as.character(Sys.time()), sep = "")
options(stringsAsFactors = strAsFacL)
rm(list = ls())

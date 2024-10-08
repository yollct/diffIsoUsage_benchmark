library(tidyverse)
library(dplyr)
library(VGAM)
library(Boom)
library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
library(scDesign3)
library(useful)
library(gridExtra)
library(ggplot2)
library(scran)
library(splatter)
library(Matrix)

readfilesdir <- "/nfs/home/students/chit/is_benchmark/singlecells/rawdata/sm3_salmon_count.h5ad" # directory to rsem output (as input)
######read single cell data and define clusters
sce <- readH5AD(readfilesdir)
geneid <- read.csv("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/sm3_salmon_count.h5ad.genetx.txt",sep="\t")
cellid <- read.csv("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/sm3_cellid.out", header=F)
#celltypes <- 
pbmc <- CreateSeuratObject(counts=assay(sce,"X"), project="pbmcsm3")
row.names(pbmc) <- geneid$isoform_id
colnames(pbmc) <- cellid$V1
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
Idents(pbmc) %>% head
colnames(pbmc) %>% head
#######

countmat <- readMM("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/sm3_salmon_count.mtx")
countspl <- array(countmat, c(5924, 251121))
sce <- SingleCellExperiment(list(counts=t(countspl)), colData=DataFrame(cellid=cellid), rowData=DataFrame(isoform_id=geneid$isoform_id))
counts(sce) <- assay(sce)
####### splatter 
params <- splatEstimate(sce)
######
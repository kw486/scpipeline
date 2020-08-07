install.packages("devtools")
library(devtools)

install.packages("dplyr")
install.packages("igraph") 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
BiocManager::install(c('scater', 'scran', 'uwot'))

install_github("markowetzlab/SCperturb")
library(SCperturb)
SCperturb::listprojects()

devtools::install_github("bitmask/labnetmet") 
library(labnetmet)
devtools::install_github("bitmask/NEMpipeline")
library(NEMpipeline)


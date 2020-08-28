############ BOOTSTRAPPING ###############
#Test robustness. Resample 80% Egenes and rerun network inference on new Egene subset.

############ FUNCTIONS & LIBRARIES ###############
library(dplyr)
install.packages("dgof") 
library("dgof")
devtools::install_github("bitmask/labnetmet") ; 
library(labnetmet)

run_nems <- function(nem_method, expr_data, selected.genes) {
  contr <- c(0.15,0.05)
  type<-"mLL"
  hyper <- nem::set.default.parameters(selected.genes, 
                                       para=contr, 
                                       type = type)
  b <- nem::nem.bootstrap(as.matrix(expr_data),
                          inference=nem_method,
                          control=hyper, 
                          verbose=F,
                          nboot=1 # 
  )
  return(b)
}

################## LOAD ##########################

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## might need to change this working directory if/when pipeline is rerun to make new outputs 
setwd("~/Google Drive File Stream/My Drive/NEMs Project MPhil/Data/seurat_prepared")
e.genes.YPD <-loadRData("e.genes_YPD.Rdata") 
e.genes.RAPA <-loadRData("e.genes_RAPA.Rdata") 
e.genes.DIAUXY <-loadRData("e.genes_DIAUXY.Rdata") 
e.genes.CSTARVE <-loadRData("e.genes_CSTARVE.Rdata") 
e.genes.NLIM <-loadRData("e.genes_NLIM.Rdata") 

###################################################

sce.name<-"CSTARVE"
e.genes<-e.genes.CSTARVE

##Bootstraps
bootstraps<-list()

for (i in 1:200){
    
  sampled.e.genes<-sample_frac(e.genes, size = 0.8, replace = FALSE) # sample 80% of the e.genes
  bootstrap.nem<-run_nems("triples",sampled.e.genes,colnames(sampled.e.genes)) #create a new nem
  bootstraps[[i]]<-as(bootstrap.nem$graph, "matrix")
}
#check distances between these bootstrap nems. Mean will be used to compare to null
bootstrap.distances<-labnetmet::generate_distances(bootstraps, trans_dist) 
sp<-labnetmet:: plot_dist(bootstraps, trans_dist)
bootstrap.distances.long<-labnetmet:: generate_distances_long(bootstraps, trans_dist)

sp + geom_hline(yintercept=2000)

## Null Model
nullmodels<-list()

for (i in 1:200){
  
  null.e.genes<-e.genes
  colnames(null.e.genes)<- sample(colnames(e.genes)) # shuffle just the column names (perturbations) but leave the rest of the data the same
  null.nem<-run_nems("triples",null.e.genes,colnames(null.e.genes))  #create a nem
  nullmodels[[i]]<-as(null.nem$graph, "matrix")
  
}

#Check the distances of the null nems
null.distances<-labnetmet::generate_distances(nullmodels, trans_dist)
labnetmet::plot_dist(nullmodels, trans_dist)
null.distances.long<-labnetmet::generate_distances_long(nullmodels, trans_dist)

#Kolmogorov–Smirnov test

ks<-dgof::ks.test(bootstrap.distances.long$distance, null.distances.long$distance)



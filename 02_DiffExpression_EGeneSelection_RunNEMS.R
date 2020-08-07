############ CALCULATE DIFFERENTIAL EXPRESSION USING SEURAT ###############
#' cycling through the 11 different perterbed genes, extracting their counts, 
#' extract counts for wildtype, do pairwise diff expression

library(dplyr)
library(Seurat)


diffexp_seurat <- function(countsdata, rowmetadata, colmetadata, sce.name, print=TRUE) {
  
  seuratoutput<-list()
  
  perturbed.genes <- unique(colmetadata$gene)
  seuratoutput[["perturbed.genes"]]<-perturbed.genes
  
  #find perturbation with smallest number of cells - use to set k.filter
  dataset.sizes<-data.frame("Num Cells" = rep(0,length(perturbed.genes))) ; rownames(dataset.sizes)<-as.character(perturbed.genes)
  for (gene in as.character(perturbed.genes)) {
    dataset.sizes[gene,]<-sum(colmetadata$gene == gene)
  }
  
  seuratoutput[["dataset.sizes"]]<-dataset.sizes
  
  kfilter<-min(200,min(dataset.sizes))
  
  diffexplist <- list()
  
  start<- Sys.time()
  for (gene in as.character(perturbed.genes)) {
    if (startsWith(gene, "WT(ho)")) {
      next
    }
    creb1 <- which(startsWith(as.character(colmetadata$gene), gene))
    ctrl <- which(startsWith(as.character(colmetadata$gene), "WT(ho)"))
    
    ctrl_counts <- countsdata[,ctrl]
    ctrl_colmetadata <- colmetadata[ctrl,]
    creb1_counts <- countsdata[,creb1]
    creb1_colmetadata <- colmetadata[creb1,]
    
    colnames(ctrl_counts) <- ctrl_colmetadata$Cell_Barcode
    colnames(creb1_counts) <- creb1_colmetadata$Cell_Barcode
    rownames(ctrl_counts) <-rowmetadata$ensembl_gene_id 
    rownames(creb1_counts) <- rowmetadata$ensembl_gene_id
    
    # Set up control object
    ctrl <- Seurat::CreateSeuratObject(counts = ctrl_counts, project = "WT(ho)", min.cells=3)
    ctrl$stim <- "WT(ho)"
    ctrl <- subset(x = ctrl, subset = nFeature_RNA > 50)
    ctrl <- Seurat::NormalizeData(object = ctrl, verbose = FALSE)
    ctrl <- Seurat::FindVariableFeatures(object = ctrl, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
    
    # Set up stimulated object
    stim <- Seurat::CreateSeuratObject(counts = creb1_counts, project = "STIM", min.cells=3)
    stim$stim <- "STIM"
    stim <- subset(x = stim, subset = nFeature_RNA > 50)
    stim <- Seurat::NormalizeData(object = stim, verbose = FALSE)
    stim <- Seurat::FindVariableFeatures(object = stim, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
    
    # Integration of control (WT) vs perturbed (gene) 
    # Identify anchors between the two datasets. Will return an object which holds an integrated 
    # (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
    Seurat :: FindTransferAnchors(reference=stim, query=ctrl, dims=1:20)
    anchors <- Seurat::FindIntegrationAnchors(object.list = list(ctrl, stim),  dims = 1:20, k.filter = kfilter) # max.features = 26
    combined <- Seurat::IntegrateData(anchorset = anchors, dims = 1:20)
    Seurat::DefaultAssay(object = combined) <- "integrated"
    
    # Run the standard workflow for visualization and clustering
    combined <- Seurat::ScaleData(object = combined, verbose = FALSE)
    combined <- Seurat::RunPCA(object = combined, npcs = 20, verbose = FALSE)
    # t-SNE and Clustering
    combined <- Seurat::RunTSNE(object = combined, reduction = "pca", dims = 1:20)
    combined <- Seurat::FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
    combined <- Seurat::FindClusters(combined, resolution = 0.5)
    
    #TSNE_plot<- SeuratmPlot(combined, reduction="tsne", split.by="stim")
    #setwd("~/Google Drive File Stream/My Drive/NEMs Project MPhil/Data/TSNE")
    #png(paste(gene,"_TSNE", ".png", sep=""),width=660, height=480)
    #Seurat::DimPlot(combined, reduction="tsne", split.by="stim")
    #dev.off()
    
    combined$celltype.stim <- paste(Seurat::Idents(combined), combined$stim, sep = "_")
    combined$celltype <- Seurat::Idents(combined)
    nclusters <- length(levels(combined$celltype))
    Seurat::Idents(combined) <- "celltype.stim"
    
    allmarkers <- Seurat::FindAllMarkers(combined)
    allmarkers$sgene <- gene
    allmarkers$egene <- allmarkers$gene
    rownames(allmarkers) <- NULL
    diffexplist[[gene]] <- allmarkers
    
    print(gene)
    
  }
  seuratoutput[["time"]]<-Sys.time() - start
  
  setwd("~/Google Drive File Stream/My Drive/NEMs Project MPhil/Data/seurat")
  output_file_name <- paste0("diffexplist_",sce.name,".Rdata")
  save(diffexplist, file=output_file_name)
  
  seuratoutput[["diffexplist"]]<-diffexplist
  
  return(seuratoutput)
  
}

############ PREPARE DIFF EXP DATA FOR NEMS ###############
#' Uses output of DiffExp_Seurat.R

########## FILTER & BINARISE ###########

setwd("~/Google Drive File Stream/My Drive/NEMs Project MPhil/Data/seurat_prepared")

### FILTERING BY LOG FOLD CHANGE ###
diffexplist<-seurat.output$diffexplist
perturbed.genes<-seurat.output$perturbed.genes

#create new data structure to hold log fold change filtered results
diffexplist_lff<-vector(mode="list", length=12) ; names(diffexplist_lff)<-c(as.character(perturbed.genes))
#Convert between natural log (Seurat output) and log base 2 (to compare with Jackson) & filter
for (gene in as.character(perturbed.genes)){
  diffexplist_lff[[gene]]<-diffexplist[[gene]][diffexplist[[gene]]$avg_logFC * log2(exp(1)) >= 1.5,] 
}

### BINARISE AND FILTER BY ADJ P VALUE ####

binary <- function(shrink.list, adjusted_pvalue_cutoff) {
  # generate discretized matrix
  # (will be needed for downstream modeling of the NEMs)
  sig.names <- c()
  for(i in 1:length(shrink.list)){
    sig.names <- c(sig.names,shrink.list[[i]]$gene[which(shrink.list[[i]]$p_val_adj < adjusted_pvalue_cutoff)])
  }
  sig.names <- unique(sig.names)
  # discrete
  Seurat.disc.mat <- matrix(0L, 
                            nrow = length(sig.names), 
                            ncol= length(names(shrink.list))
  )
  colnames(Seurat.disc.mat) <- names(shrink.list)
  rownames(Seurat.disc.mat) <- sig.names
  # significantly DEGs with padj < 0.05 are set to 1
  for(i in 1:length(colnames(Seurat.disc.mat))){
    Seurat.disc.mat[shrink.list[[i]]$gene[which(shrink.list[[i]]$p_val_adj < adjusted_pvalue_cutoff)],i] <- 1
  }
   return(Seurat.disc.mat)
}

prepared<-binary(diffexplist_lff, 0.05)

########### FIND E GENES ###########
# Cut off for any genes diff expressed (up or down) twice or more, check histograms
df.prepared<-as.data.frame(prepared) ; print(paste0("Num Differentially Expressed Genes:",nrow(df.prepared)))
hist(rowSums(df.prepared), breaks=15, xlim=c(0,12))
e.genes<-  df.prepared[rowSums(df.prepared) >=2, ]  ; print(paste("eGenes Diff Expressed >2x:",nrow(e.genes), sep=" "))
hist(rowSums(e.genes), breaks=15, xlim=c(0,12))

setwd("~/Google Drive File Stream/My Drive/NEMs Project MPhil/Data/seurat_prepared")
output_file_name <- paste0("e.genes_",sce.name,".Rdata")
save(e.genes, file=output_file_name)

output_file_name <- paste0("binary_",sce.name,".Rdata")
save(prepared, file=output_file_name)


############ RUN NEMS ###############

devtools::install_github("bitmask/labnetmet") ; library(labnetmet)
devtools::install_github("bitmask/NEMpipeline") ; library(NEMpipeline)
install.packages("igraph") 

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
                          nboot=1 
  )
  return(b)
}

first.nem<-run_nems("nem.triplet",e.genes,colnames(e.genes))
new.first.nem<-nem::transitive.reduction(first.nem$graph)

### Plotting ###
plot(as(new.first.nem, "graphNEL"))

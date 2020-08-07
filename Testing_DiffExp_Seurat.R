library(dplyr)
library(Seurat)


#' Calculate differential expression using Seurat
#'
#' @param countsdata Sparse matrix containing expression data
#' @param rowmetadata Dataframe describing the rows of countsdata (reporter genes)
#' @param colmetadata Dataframe describing the columns of countsdata (cell barcode and which genes perturbed)
#' @return Dataframe containing the differential expression
#' @export
#' 
#' https://satijalab.org/seurat/v3.0/integration.html << notes on what's happening here!
#' 
#' cycling through the 11 different perterbed genes, extracting their counts, extract counts for wildtype
#' do pairwise diff expression
#' 
#' 


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
    rownames(ctrl_counts) <-rowmetadata$ensembl_gene_id #all perturbed genes are still here at this stage!
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


sce.name<-"NLIM"
seurat.output<-diffexp_seurat(counts(sce2),rowData(sce2),colData(sce2), sce.name)





sce2<-sce_NLIM


MMETOH
NLIMGLN
NLIMNH4
NLIMPRO
NLIMUREA
YPETOH

sce2<-sce_YPETOH
sce2<-sce_CSTARVE
sce2<-sce_NLIM
sce2<-sce_YPD

countsdata<- counts(sce2)
rowmetadata<-rowData(sce2)
colmetadata<-colData(sce2)


countsdata<- counts(sce_RAPA)
rowmetadata<-rowData(sce_RAPA)
colmetadata<-colData(sce_RAPA)

diffexpl<-diffexp_seurat(countsdata, rowmetadata, colmetadata)




#' To restore R object
# load(seurat_output.RData)



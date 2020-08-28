##################### Quality Control - Test & Check ########################
# Originally designed as a QC pipeline. For this specific data, this pipeline wasn't appropriate 
# (data had been pre-QCed by original authors/mitochondrial QC was too stringent for low cell 
# counts).
##########################################################################################

library(scater)
library(SCperturb)

#Read in data & convert to single cell experiment format
data(counts.jackson2019_RAPA,
     colmetadata.jackson2019_RAPA,
     rowmetadata.jackson2019)

sce2<-SingleCellExperiment(assays = list(counts = counts.jackson2019_RAPA), 
                           colData = colmetadata.jackson2019_RAPA, 
                           rowData = rowmetadata.jackson2019)

#Find mitochondrial genes based on the fact they always start with "Q"
is.mito <- grepl("^Q", rownames(sce2)); genes<-rownames(sce2); genes[is.mito]

#Use single step scater functions to perform QC analysis -
qcstats <- perCellQCMetrics(sce2, subsets=list(Mito=is.mito)); qcstats
colSums(as.matrix(filtered)) # Summarize the number of cells proposed to be removed for each reason

# QuickCellQC is just a wrapper for `isOutlier`, so check the thresholds it was pulling out:
qc.lib2 <- isOutlier(qcstats$sum, log=TRUE, nmads=3, type="lower")
attr(qc.lib2, "thresholds")
qc.nexprs2 <- isOutlier(qcstats$detected, nmads=3, log=TRUE, type="lower")
attr(qc.nexprs2, "thresholds")
qc.mito2 <- isOutlier(qcstats$subsets_Mito_percent, nmads=3, type="higher")
attr(qc.mito2, "thresholds") 

#Examine the plots
hist(qcstats$total/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(qcstats$detected, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(qcstats$subsets_Mito_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

batch <- paste0(sce2$grna)
batch.reasons <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent", batch=batch)
colSums(as.matrix(batch.reasons))
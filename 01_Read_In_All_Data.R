##################### Read in data & convert to single cell format ########################
#  Does this but for any condition i.e. "RAPA" substitute: 
#  data(counts.jackson2019_RAPA,
#      colmetadata.jackson2019_RAPA,
#       rowmetadata.jackson2019)
#  sce2<-SingleCellExperiment(assays = list(counts = counts.jackson2019_RAPA), 
#                             colData = colmetadata.jackson2019_RAPA, 
#                             rowData = rowmetadata.jackson2019)
# Loads from original data repository
##########################################################################################

data(counts.jackson2019_RAPA,
       colmetadata.jackson2019_RAPA,
            rowmetadata.jackson2019)

sce_RAPA<-SingleCellExperiment(assays = list(counts = counts.jackson2019_RAPA), 
                                                      colData = colmetadata.jackson2019_RAPA, 
                                                      rowData = rowmetadata.jackson2019)


library(SCperturb)
library(SingleCellExperiment)
data_dir<- "~/Google Drive File Stream/My Drive/NEMs Project MPhil/Data/sce"

listofnames<- c("CSTARVE", "NLIMUREA", "NLIMNH4", "NLIMPRO", "NLIMGLN", "MMETOH", "MMD", "YPETOH", "RAPA", "DIAUXY", "YPD")

for(i in 1:length(listofnames)){
  first = paste0("counts.jackson2019_",listofnames[i])
  second = paste0("colmetadata.jackson2019_",listofnames[i])
  third = paste0("rowmetadata.jackson2019")
  
  eval(
    str2expression(
      paste0("data(",first,",",second,",",third,")") 
    )
  )
  
  eval(
    str2expression(
      paste0("data(",paste0("counts.jackson2019_",listofnames[i]),",",second,",",third,")") 
    )
  )
  
  
  # assign allows you to create and assign to a variable name from paste
  assign(paste0("sce_",listofnames[i]), SingleCellExperiment(assays = list(counts = eval(str2expression(first))), 
                                                             colData = eval(str2expression(second)), 
                                                             rowData = eval(str2expression(third))) 
  )
  
  # SAVE #
  #read.data<-eval(str2expression(paste0("sce_",listofnames[i])))
  #output_file_name <- paste0("sce_",listofnames[i],".Rdata")
  #save(read.data, file=output_file_name)
  
  rm(first,second, third)
  rm(list=paste0("counts.jackson2019_",listofnames[i]))
  rm(list=paste0("colmetadata.jackson2019_",listofnames[i]))
  rm(list=paste0("rowmetadata.jackson2019"))
  rm(i); rm(read.data) ; rm(output_file_name)
}

################## COMBINE THE NLIM SCEs ##############
sce_NLIM<-cbind(sce_NLIMGLN,sce_NLIMNH4,sce_NLIMPRO,sce_NLIMUREA, deparse.level=1)

################## LOAD ##########################
# Use if you want to manually load single data files that you 
# previously stored locally

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("~/Google Drive File Stream/My Drive/NEMs Project MPhil/Data/sce")
#...for softcoded single celled object...
sce2 <-loadRData("sce_NLIM.Rdata") 

#..or for hard coded single cell objects
sce_YPD<-loadRData("sce_YPD.Rdata") 
sce_CSTARVE<-loadRData("sce_CSTARVE.Rdata")
sce_DIAUXY<-loadRData("sce_DIAUXY.Rdata")
sce_RAPA<-loadRData("sce_RAPA.Rdata")
sce_MMD<-loadRData("sce_MMD.Rdata") #
sce_YPETOH<-loadRData("sce_YPETOH.Rdata") ;sce.name<-"YPETOH"
sce_MMETOH<-loadRData("sce_MMETOH.Rdata") ;sce.name<-"MMETOH"



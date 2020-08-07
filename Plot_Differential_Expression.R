### Create Tile Plots of Number of Differentially Expressed Genes in Each Condition ###

## Prefilter / prepared
mm<-data.frame(condition=rep(c("NLIM","CSTARVE","RAPA","DIAUXY","YPD"), each=11), gene=rep(0,55), number=rep(0,55))

for (i in 1:11){
  mm[i+44,3]<-sum(prepared[,i])
  mm[i+44,2]<-colnames(prepared)[i]
}

gg<- ggplot(mm, aes(condition, gene))
gg+ 
  geom_tile(aes(fill = number), colour = "white") + 
  geom_text(aes(label=number)) +
  scale_fill_gradient(low = "white", high = "dodgerblue", trans="log", na.value = "white") +
  labs(x = "",y = "") +
  theme(legend.position="none")

## E genes

ee<-data.frame(condition=rep(c("NLIM","CSTARVE","RAPA","DIAUXY","YPD"), each=11), gene=rep(0,55), number=rep(0,55))
for (i in 1:11){
  ee[i+44,3]<-sum(e.genes[,i])
  ee[i+44,2]<-colnames(e.genes)[i]
}

ge<- ggplot(ee, aes(condition, gene))
ge+ 
  geom_tile(aes(fill = number), colour = "white") + 
  geom_text(aes(label=number)) +
  scale_fill_gradient(low = "white", high = "dodgerblue", trans="log", na.value = "white") +
  labs(x = "",y = "") +
  theme(legend.position="none")



        
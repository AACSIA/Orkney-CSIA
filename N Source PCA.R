## Set working directory and load data sheet
getwd()
setwd("C:/Users/annie/OneDrive/R")
library(giscoR)
library(ggplot2)
library("ggspatial")
library(corrplot)
library(FactoMineR)
library(factoextra)
humanreldat <- read.csv("PCADataSource.csv")

aminoacids <- humanreldat[,3:6]
plot(aminoacids,pch=20)
dev.off()
x <- cor(aminoacids)
corrplot(x,diag=F)
dev.off()
corrplot(x, order = "FPC", diag=F)
dev.off()
res.pca <- PCA(aminoacids,scale.unit=TRUE,graph=FALSE)

png(file = "SOURCEPCAsplit1.png", height = 1500, width = 2000, units = "px", res=300)
par(mar=c(4,4,1,0))
plot(res.pca,choix="var")+
  labs(title ="", x = "PC1", y = "PC2")
dev.off()

png(file = "SOURCEPCAsplit2.png", height = 1500, width = 2000, units = "px", res=300)
par(mar=c(4,4,1,0))
plot(res.pca)+
  labs(title ="", x = "PC1", y = "PC2")
dev.off()

png(file = "SOURCEPCA.png", height = 1500, width = 2000, units = "px", res=300)
par(mar=c(4,4,1,0))
fviz_pca_biplot(res.pca, label="var",geom="point", pointsize = 2.5,col.var="black") +
  labs(title ="", x = "PC1 (54.81%)", y = "PC2 (25.66%)")
dev.off()

res.pca$eig

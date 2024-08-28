getwd()
setwd("C:/Users/annie/OneDrive/R")
library(giscoR)
library(ggplot2)
library("ggspatial")
library(corrplot)
library(FactoMineR)
library(factoextra)
humanreldat <- read.csv("PCADataLimpets.csv")

aminoacids <- humanreldat[,2:14]
plot(aminoacids,pch=20)
dev.off()
x <- cor(aminoacids)
corrplot(x,diag=F)
dev.off()
corrplot(x, order = "FPC", diag=F)
dev.off()
res.pca <- PCA(aminoacids,scale.unit=TRUE,graph=FALSE)

png(file = "PCAsplit1LIMPETS.png", height = 1500, width = 2000, units = "px", res=300)
par(mar=c(4,4,1,0))
plot(res.pca,choix="var")+
  labs(title ="", x = "PC1", y = "PC2")
dev.off()

png(file = "PCAsplit2LIMPETS.png", height = 1500, width = 2000, units = "px", res=300)
par(mar=c(4,4,1,0))
plot(res.pca)+
  labs(title ="", x = "PC1", y = "PC2")
dev.off()

png(file = "PCAL&B.png", height = 1500, width = 2000, units = "px", res=300)
par(mar=c(4,4,1,0))
fviz_pca_biplot(res.pca, label="var",geom="point", pointsize = 2.5,col.var="black") +
  labs(title ="", x = "PC1 (79.65%)", y = "PC2 (9.91)")
dev.off()

res.pca$eig

png(file = "PCAL&B.png", height = 1500, width = 2000, units = "px", res = 300)
par(mar = c(4, 4, 1, 0))

# Create the PCA biplot without gridlines
fviz_pca_biplot(res.pca, label = "var", geom = "point", pointsize = 2.5, col.var = "black") +
  labs(title = "", x = "PC1 (79.65%)", y = "PC2 (9.91)") +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines

dev.off()

png(file = "NitrogenPCA.png", height = 1500, width = 2000, units = "px", res = 300)
par(mar = c(4, 4, 1, 0))

# Create the PCA biplot without gridlines and central axis lines, but with an outer border
fviz_pca_biplot(res.pca, label = "var", geom = "point", pointsize = 2.5, col.var = "black") +
  labs(title = "", x = "PC1 (79.65%)", y = "PC2 (9.91)") +
  theme(panel.grid.major = element_blank(),    # Remove major gridlines
        panel.grid.minor = element_blank(),    # Remove minor gridlines
        axis.line = element_blank(),           # Remove central axis lines
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))  # Add a black border around the plot

dev.off()




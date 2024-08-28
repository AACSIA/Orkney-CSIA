## Set working directory and load data sheet
getwd()
setwd("C:/Users/annie/OneDrive/R")
library(dunn.test)
library(sf)
library(giscoR)
library(ggplot2)
library("ggspatial")
library(corrplot)
library(FactoMineR)
library(factoextra)
humanreldat <- read.csv("AAdata.csv")

humanreldat$colour <- NA
HumanIndex <- which(humanreldat$species %in% "Human")
SaitheIndex <- which(humanreldat$species %in% "Saithe")
CodIndex <- which(humanreldat$species %in% "Cod")
CattleIndex <- which(humanreldat$species %in% "Cow")
humanreldat$colour[HumanIndex] <- "red"
humanreldat$colour[SaitheIndex] <- "navy"
humanreldat$colour[CodIndex] <- "pink"
humanreldat$colour[CattleIndex] <- "orange"

# N glu phe trophic position plot
png(file = "Ngluphe.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$phe_d15N,y=humanreldat$glu_d15N, xlab=expression("δ"^15*"N"[Phe]* "(‰)"),ylab=expression("δ"^15*"N"[Glu]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4, xlim = c(0, 15),
     ylim = c(0, 30),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(0)),
                 round(max(15)), by = 5))
axis(2, at = seq(round(min(0)),
                 round(max(30)), by = 5))
abline(coef = c(-8.4,1),col="red",lty = "dashed")
abline(coef = c(-0.8,1),col="red",lty = "dashed")
abline(coef = c(6.8,1),col="red",lty = "dashed")
abline(coef = c(11,1),col="darkblue",lty = "dashed")
abline(coef = c(18.6,1),col="darkblue",lty = "dashed")
abline(coef = c(26.2,1),col="darkblue",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3)
dev.off()

# N glu phe trophic position plot with Nielsen changes
png(file = "NglupheNielsen.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$phe_d15N,y=humanreldat$glu_d15N, xlab=expression("δ"^15*"N"[Phe]* "(‰)"),ylab=expression("δ"^15*"N"[Glu]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4, xlim = c(0, 15),
     ylim = c(0, 30),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=3)
axis(1, at = seq(round(min(0)),
                 round(max(15)), by = 5))
axis(2, at = seq(round(min(0)),
                 round(max(30)), by = 5))
abline(coef = c(-8.4,1),col="red",lty = "dashed")
abline(coef = c(-1.8,1),col="red",lty = "dashed")
abline(coef = c(4.8,1),col="red",lty = "dashed")
abline(coef = c(9.5,1),col="darkblue",lty = "dashed")
abline(coef = c(16.1,1),col="darkblue",lty = "dashed")
abline(coef = c(22.7,1),col="darkblue",lty = "dashed")
abline(coef = c(29.3,1),col="darkblue",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3)
dev.off()

# N glu phe trophic position plot with Nielsen changes AND BARLEY
png(file = "NglupheNielsen.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$phe_d15N,y=humanreldat$glu_d15N, xlab=expression("δ"^15*"N"[Phe]* "(‰)"),ylab=expression("δ"^15*"N"[Glu]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4, xlim = c(0, 15),
     ylim = c(0, 30),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=3)
axis(1, at = seq(round(min(0)),
                 round(max(15)), by = 5))
axis(2, at = seq(round(min(0)),
                 round(max(30)), by = 5))
abline(coef = c(-8.4,1),col="red",lty = "dashed")
abline(coef = c(-1.8,1),col="red",lty = "dashed")
abline(coef = c(4.8,1),col="red",lty = "dashed")
abline(coef = c(9.5,1),col="darkblue",lty = "dashed")
abline(coef = c(16.1,1),col="darkblue",lty = "dashed")
abline(coef = c(22.7,1),col="darkblue",lty = "dashed")
abline(coef = c(29.3,1),col="darkblue",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3)
dev.off()


# N glu-phe by bulk N
png(file = "NgluphebyNcoll.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$d15N,y=humanreldat$gluphe_d15N, xlab=expression("δ"^15*"N"[Coll]* "(‰)"),ylab=expression("δ"^15*"N"[Glu-Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4, xlim = c(0, 20),
     ylim = c(-5, 30),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=4)
axis(1, at = seq(round(min(0)),
                 round(max(20)), by = 5))
axis(2, at = seq(round(min(-5)),
                 round(max(30)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3,bty = "n")
dev.off()
cor.test(x=humanreldat$d15N,y=humanreldat$gluphe_d15N)

# Set up the plot with the desired dimensions and margins
png(file = "NgluphebyNcoll.png", height = 720, width = 720) # Adjust width to match height for a squarer plot
par(mar = c(5, 6, 4, 2))

# Create the plot with outlined points
plot(x = humanreldat$d15N, y = humanreldat$gluphe_d15N, 
     xlab = expression("δ"^15*"N"[Coll]* "(‰)"),
     ylab = expression("δ"^15*"N"[Glu-Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(0, 18), ylim = c(-5, 30), xaxt = "n", yaxt = "n", 
     cex.lab = 2, cex.axis = 4)

# Add the axes
axis(1, at = seq(0, 18, by = 5), cex.axis = 1.5)
axis(2, at = seq(-5, 30, by = 5), cex.axis = 1.5)

# Add a dashed line
abline(coef = c(0, 1), col = "black", lty = "dashed")

# Add the legend
legend(x = "bottomright", legend = c("Human", "Saithe", "Cod", "Cattle"), 
       col = "black", pt.bg = c("red", "navy", "pink", "orange"), pch = 21, 
       cex = 2, pt.cex = 3, bty = "n")

# Finish the plot
dev.off()


# N glu-phe by bulk C compared with only bulk
png(file = "gluphebyCcoll.png", height = 720, width = 960)
par(mar=c(5, 6, 4, 15), xpd=TRUE)
plot(x=humanreldat$d13C,y=humanreldat$gluphe_d15N, xlab="δ13CColl (‰)",ylab="δ15NGlu-Phe (‰)", col=humanreldat$colour,  pch=20, cex=4, xlim = c(-25, -10),
     ylim = c(-5, 30),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-25)),
                 round(max(-10)), by = 5))
axis(2, at = seq(round(min(-5)),
                 round(max(30)), by = 5))
legend("topright", inset=c(-0.23, 0.36),legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,bty = "n",pt.cex=3)
dev.off()

#bulk N by bulk C
png(file = "CCollNcoll.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$d13C,y=humanreldat$d15N, xlab=expression("δ"^13*"C"[Coll]* "(‰)"),ylab=expression("δ"^15*"N"[Coll]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4, xlim = c(-25, -10),
     ylim = c(-5, 30),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-25)),
                 round(max(-10)), by = 5))
axis(2, at = seq(round(min(-5)),
                 round(max(30)), by = 5))
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3,bty = "n")
dev.off()

# Assuming humanreldat is a data frame with columns 'd13C', 'd15N', and 'colour'
# If it's not, please replace 'humanreldat' with the actual data frame name.

# Set up the plot with the desired dimensions and margins
png(file = "CCollNcoll.png", height = 720, width = 960)
par(mar = c(5, 6, 4, 2))

# Create the plot with outlined points
plot(x = humanreldat$d13C, y = humanreldat$d15N, 
     xlab = expression("δ"^13*"C"[Coll]* "(‰)"),
     ylab = expression("δ"^15*"N"[Coll]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(-25, -10), ylim = c(-5, 30), xaxt = "n", yaxt = "n", 
     cex.lab = 2, cex.axis = 2)

# Add the axes
axis(1, at = seq(-25, -10, by = 5), cex.axis = 1.5)
axis(2, at = seq(-5, 30, by = 5), cex.axis = 1.5)

# Add the legend
legend(x = "bottomright", legend = c("Human", "Saithe", "Cod", "Cattle"), 
       col = "black", pt.bg = c("red", "navy", "pink", "orange"), pch = 21, 
       cex = 2, pt.cex = 3, bty = "n")

# Finish the plot
dev.off()


#1x2 bulk n against gluphe n
png(file = "BULKNvsGLUPHEN.png", height = 420, width = 960)
par(mfrow=c(1,2))
par(mar=c(5,6,4,2))
plot(x=humanreldat$d13C,y=humanreldat$d15N, xlab=expression("δ"^13*"C"[Coll]* "(‰)"),ylab=expression("δ"^15*"N"[Coll]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4, xlim = c(-25, -10),
     ylim = c(-5, 30),xaxt = "n",yaxt = "n",cex.lab=1.6, cex.axis=2)
axis(1, at = seq(round(min(-25)),
                 round(max(-10)), by = 5))
axis(2, at = seq(round(min(-5)),
                 round(max(30)), by = 5))
plot(x=humanreldat$d13C,y=humanreldat$gluphe_d15N, xlab=expression("δ"^13*"C"[Coll]* "(‰)"),ylab=expression("δ"^15*"N"[Glu-Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4, xlim = c(-25, -10),
     ylim = c(-5, 30),xaxt = "n",yaxt = "n",cex.lab=1.6, cex.axis=2)
axis(1, at = seq(round(min(-25)),
                 round(max(-10)), by = 5))
axis(2, at = seq(round(min(-5)),
                 round(max(30)), by = 5))
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=1.5,pt.cex=2.5,bty = "n")
dev.off()


# Set up the plot with the desired dimensions and margins
png(file = "BULKNvsGLUPHEN2.png", height = 420, width = 960)
par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 2))

# First plot
plot(x = humanreldat$d13C, y = humanreldat$d15N, 
     xlab = expression("δ"^13*"C"[Coll]* "(‰)"),
     ylab = expression("δ"^15*"N"[Coll]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3, 
     xlim = c(-25, -10), ylim = c(-5, 30), xaxt = "n", yaxt = "n", 
     cex.lab = 1.6, cex.axis = 2)
axis(1, at = seq(-25, -10, by = 5), cex.axis = 1.5)
axis(2, at = seq(-5, 30, by = 5), cex.axis = 1.5)

# Second plot
plot(x = humanreldat$d13C, y = humanreldat$gluphe_d15N, 
     xlab = expression("δ"^13*"C"[Coll]* "(‰)"),
     ylab = expression("δ"^15*"N"[Glu-Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3, 
     xlim = c(-25, -10), ylim = c(-5, 30), xaxt = "n", yaxt = "n", 
     cex.lab = 1.6, cex.axis = 2)
axis(1, at = seq(-25, -10, by = 5), cex.axis = 1.5)
axis(2, at = seq(-5, 30, by = 5), cex.axis = 1.5)

# Add the legend
legend(x = "bottomright", legend = c("Human", "Saithe", "Cod", "Cattle"), 
       col = "black", pt.bg = c("red", "navy", "pink", "orange"), pch = 21, 
       cex = 1.5, pt.cex = 2.5, bty = "n")

# Finish the plot
dev.off()



# C Val by C phe
png(file = "CVALbyCPHE.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$val_d13C,y=humanreldat$phe_d13C, xlab=expression("δ"^13*"C"[Val]* "(‰)"),ylab=expression("δ"^13*"C"[Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(-30, -15),
     ylim = c(-30, -20),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-30)),
                 round(max(-15)), by = 5))
axis(2, at = seq(round(min(-30)),
                 round(max(-20)), by = 5))
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,bty="n",pt.cex=3)
dev.off()

# Assuming humanreldat is a data frame with columns 'val_d13C', 'phe_d13C', and 'colour'
# If it's not, please replace 'humanreldat' with the actual data frame name.

# Set up the plot with the desired dimensions and margins
png(file = "CVALbyCPHE2.png", height = 720, width = 960)
par(mar = c(6, 7, 4, 2), mgp = c(4, 1, 0))  # Moved x-axis label slightly down

# Create the plot with outlined points
plot(x = humanreldat$val_d13C, y = humanreldat$phe_d13C, 
     xlab = expression("δ"^13*"C"[Val]* "(‰)"),
     ylab = expression("δ"^13*"C"[Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(-30, -15), ylim = c(-30, -20), xaxt = "n", yaxt = "n", 
     cex.lab = 2.5, cex.axis = 2)

# Add the axes
axis(1, at = seq(-30, -15, by = 5), cex.axis = 1.5)
axis(2, at = seq(-30, -20, by = 2), cex.axis = 1.5)

# Add the legend
legend(x = "bottomright", legend = c("Human", "Saithe", "Cod", "Cattle"), 
       col = "black", pt.bg = c("red", "navy", "pink", "orange"), pch = 21, 
       cex = 2, pt.cex = 3, bty = "n")

# Finish the plot
dev.off()


# N PCA
aminoacid<-data.frame(humanreldat$ala_d15N,humanreldat$asp_d15N,humanreldat$glu_d15N,humanreldat$gly_d15N,humanreldat$hyp_d15N,humanreldat$ile_d15N,humanreldat$leu_d15N,humanreldat$lys_d15N,humanreldat$phe_d15N,humanreldat$pro_d15N,humanreldat$ser_d15N,humanreldat$thr_d15N,humanreldat$val_d15N,humanreldat$nle_d15N)
plot(aminoacid,pch=20,col=humanreldat$colour)
dev.off()
x <- cor(aminoacid)
corrplot(x,diag=F)
dev.off()
corrplot(x, order = "FPC", diag=F)
dev.off()
res.pca <- PCA(aminoacid,scale.unit=TRUE,graph=FALSE)
par(mar=c(4,4,1,0))
plot(res.pca,choix="var")
dev.off()


#hyp-pro C
png(file = "CHYPPRO.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$pro_d13C,y=humanreldat$hyp_d13C, xlab=expression("δ"^13*"C"[Pro]* "(‰)"),ylab=expression("δ"^13*"C"[Hyp]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(-25, -10),
     ylim = c(-25, -10),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-25)),
                 round(max(-10)), by = 5))
axis(2, at = seq(round(min(-25)),
                 round(max(-10)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=2.5)
dev.off()
cor.test(x=humanreldat$pro_d13C,y=humanreldat$hyp_d13C)

#hyp-pro N
png(file = "NHYPPRO.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$pro_d15N,y=humanreldat$hyp_d15N, xlab=expression("δ"^15*"N"[Pro]* "(‰)"),ylab=expression("δ"^15*"N"[Hyp]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(5, 20),
     ylim = c(5, 20),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(5)),
                 round(max(20)), by = 5))
axis(2, at = seq(round(min(5)),
                 round(max(20)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=2.5)
dev.off()
cor.test(x=humanreldat$pro_d15N,y=humanreldat$hyp_d15N)

#est by obs C
png(file = "Cestbyobs.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$d13C,y=humanreldat$d13C_est, xlab=expression("δ"^13*"C"[Obs]* "(‰)"),ylab=expression("δ"^13*"C"[Est]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(-25, -5),
     ylim = c(-25, -5),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-25)),
                 round(max(-5)), by = 5))
axis(2, at = seq(round(min(-25)),
                 round(max(-5)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=2.5)
dev.off()
cor.test(x=humanreldat$d13C,y=humanreldat$d13C_est)

#est by obs N
png(file = "Nestbyobs.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$d15N,y=humanreldat$d15N_est, xlab=expression("δ"^15*"N"[Obs]* "(‰)"),ylab=expression("δ"^15*"N"[Est]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(0, 20),
     ylim = c(0, 20),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(0)),
                 round(max(20)), by = 5))
axis(2, at = seq(round(min(0)),
                 round(max(20)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=2.5)
dev.off()
cor.test(x=humanreldat$d15N,y=humanreldat$d15N_est)

#4x4 QC plots
png(file = "QCPLOTS.png", height = 720, width = 960)
par(mfrow=c(2,2))
par(mar=c(5,6,4,2))
plot(x=humanreldat$d13C,y=humanreldat$d13C_est, xlab=expression("δ"^13*"C"[Obs]* "(‰)"),ylab=expression("δ"^13*"C"[Est]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(-25, -5),
     ylim = c(-25, -5),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-25)),
                 round(max(-5)), by = 5))
axis(2, at = seq(round(min(-25)),
                 round(max(-5)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
plot(x=humanreldat$d15N,y=humanreldat$d15N_est, xlab=expression("δ"^15*"N"[Obs]* "(‰)"),ylab=expression("δ"^15*"N"[Est]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(0, 20),
     ylim = c(0, 20),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(0)),
                 round(max(20)), by = 5))
axis(2, at = seq(round(min(0)),
                 round(max(20)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
plot(x=humanreldat$pro_d13C,y=humanreldat$hyp_d13C, xlab=expression("δ"^13*"C"[Pro]* "(‰)"),ylab=expression("δ"^13*"C"[Hyp]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(-25, -10),
     ylim = c(-25, -10),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-25)),
                 round(max(-10)), by = 5))
axis(2, at = seq(round(min(-25)),
                 round(max(-10)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
plot(x=humanreldat$pro_d15N,y=humanreldat$hyp_d15N, xlab=expression("δ"^15*"N"[Pro]* "(‰)"),ylab=expression("δ"^15*"N"[Hyp]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(5, 20),
     ylim = c(5, 20),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(5)),
                 round(max(20)), by = 5))
axis(2, at = seq(round(min(5)),
                 round(max(20)), by = 5))
abline(coef = c(0,1),col="black",lty = "dashed")
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3,bty = "n")
dev.off()

# Set up the plot with the desired dimensions and margins
png(file = "QCPLOTS.png", height = 720, width = 960)
par(mfrow = c(2, 2))
par(mar = c(5, 6, 4, 2))

# Plot 1
plot(x = humanreldat$d13C, y = humanreldat$d13C_est, 
     xlab = expression("δ"^13*"C"[Obs]* "(‰)"),
     ylab = expression("δ"^13*"C"[Est]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(-25, -5), ylim = c(-25, -5), xaxt = "n", yaxt = "n", 
     cex.lab = 2, cex.axis = 2)
axis(1, at = seq(-25, -5, by = 5), cex.axis = 2)
axis(2, at = seq(-25, -5, by = 5), cex.axis = 2)
abline(coef = c(0, 1), col = "black", lty = "dashed")

# Plot 2
plot(x = humanreldat$d15N, y = humanreldat$d15N_est, 
     xlab = expression("δ"^15*"N"[Obs]* "(‰)"),
     ylab = expression("δ"^15*"N"[Est]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(0, 20), ylim = c(0, 20), xaxt = "n", yaxt = "n", 
     cex.lab = 2, cex.axis = 2)
axis(1, at = seq(0, 20, by = 5), cex.axis = 2)
axis(2, at = seq(0, 20, by = 5), cex.axis = 2)
abline(coef = c(0, 1), col = "black", lty = "dashed")

# Plot 3
plot(x = humanreldat$pro_d13C, y = humanreldat$hyp_d13C, 
     xlab = expression("δ"^13*"C"[Pro]* "(‰)"),
     ylab = expression("δ"^13*"C"[Hyp]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(-25, -10), ylim = c(-25, -10), xaxt = "n", yaxt = "n", 
     cex.lab = 2, cex.axis = 2)
axis(1, at = seq(-25, -10, by = 5), cex.axis = 2)
axis(2, at = seq(-25, -10, by = 5), cex.axis = 2)
abline(coef = c(0, 1), col = "black", lty = "dashed")

# Plot 4
plot(x = humanreldat$pro_d15N, y = humanreldat$hyp_d15N, 
     xlab = expression("δ"^15*"N"[Pro]* "(‰)"),
     ylab = expression("δ"^15*"N"[Hyp]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(5, 20), ylim = c(5, 20), xaxt = "n", yaxt = "n", 
     cex.lab = 2, cex.axis = 2)
axis(1, at = seq(5, 20, by = 5), cex.axis = 2)
axis(2, at = seq(5, 20, by = 5), cex.axis = 2)
abline(coef = c(0, 1), col = "black", lty = "dashed")

# Add the legend
legend(x = "bottomright", legend = c("Human", "Saithe", "Cod", "Cattle"), 
       col = "black", pt.bg = c("red", "navy", "pink", "orange"), pch = 21, 
       cex = 2, pt.cex = 3, bty = "n")

# Finish the plot
dev.off()


#val-phe by lys-phe
png(file = "Cvalphebylysphe.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$lysphe_d13C,y=humanreldat$valphe_d13C, xlab=expression("δ"^13*"C"[Lys-Phe]* "(‰)"),ylab=expression("δ"^13*"C"[Val-Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(6, 12),
     ylim = c(-3, 5),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(6)),
                 round(max(12)), by = 2))
axis(2, at = seq(round(min(-3)),
                 round(max(5)), by = 2))
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3,bty="n")
dev.off()

# Assuming humanreldat is a data frame with columns 'lysphe_d13C', 'valphe_d13C', and 'colour'
# If it's not, please replace 'humanreldat' with the actual data frame name.

# Set up the plot with the desired dimensions and margins
png(file = "Cvalphebylysphe2.png", height = 720, width = 960)
par(mar = c(6, 7, 4, 2), mgp = c(4, 1, 0))  # Moved x-axis label slightly down

# Create the plot with outlined points
plot(x = humanreldat$lysphe_d13C, y = humanreldat$valphe_d13C, 
     xlab = expression("δ"^13*"C"[Lys-Phe]* "(‰)"),
     ylab = expression("δ"^13*"C"[Val-Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(6, 12), ylim = c(-3, 5), xaxt = "n", yaxt = "n", 
     cex.lab = 2.5, cex.axis = 2)

# Add the axes
axis(1, at = seq(6, 12, by = 2), cex.axis = 1.5)
axis(2, at = seq(-3, 5, by = 2), cex.axis = 1.5)

# Add the legend
legend(x = "bottomright", legend = c("Human", "Saithe", "Cod", "Cattle"), 
       col = "black", pt.bg = c("red", "navy", "pink", "orange"), pch = 21, 
       cex = 2, pt.cex = 3, bty = "n")

# Finish the plot
dev.off()


#val-phe C by glu-phe N
png(file = "CvalphebyglupheN.png", height = 720, width = 960)
par(mar=c(5,6,4,2))
plot(x=humanreldat$valphe_d13C,y=humanreldat$gluphe_d15N, xlab=expression("δ"^13*"C"[Val-Phe]* "(‰)"),ylab=expression("δ"^15*"N"[Glu-Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(-3, 5),
     ylim = c(-5, 25),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-3)),
                 round(max(5)), by = 1))
axis(2, at = seq(round(min(-5)),
                 round(max(25)), by = 5))
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3,bty="n")
dev.off()

png(file = "valphebygluphe.png", height = 720, width = 960)
par(mar = c(6, 7, 4, 2), mgp = c(4, 1, 0))

# Create the plot with outlined points
plot(x = humanreldat$valphe_d13C, y = humanreldat$gluphe_d15N, 
     xlab = expression("δ"^13*"C"[Val-Phe]* "(‰)"),
     ylab = expression("δ"^15*"N"[Glu-Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 3.5, 
     xlim = c(-3, 5), ylim = c(-5, 25), xaxt = "n", yaxt = "n", 
     cex.lab = 2.5, cex.axis = 2)

# Add the axes
axis(1, at = seq(-3, 5, by = 1), cex.axis = 1.5)
axis(2, at = seq(-5, 25, by = 5), cex.axis = 1.5)

# Add the legend
legend(x = "bottomright", legend = c("Human", "Saithe", "Cod", "Cattle"), 
       col = "black", pt.bg = c("red", "navy", "pink", "orange"), pch = 21, 
       cex = 2, pt.cex = 3, bty = "n")

# Finish the plot
dev.off()

# c thr by c gly-phe
png(file = "Cthrbyglyphe.png", height = 720, width = 960)
par(mar=c(8,6,4,2))
plot(x=humanreldat$glyphe_d13C,y=humanreldat$thr_d13C, xlab="",ylab=expression("δ"^13*"C"[Thr]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(5, 21),
     ylim = c(-15, -1),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(5)),
                 round(max(21)), by = 2))
axis(2, at = seq(round(min(-15)),
                 round(max(-1)), by = 2))
title(xlab=expression("δ"^13*"C"[Gly-Phe]* "(‰)"), cex.lab=2, line=5.2)
legend(x="bottomright",legend=c("Human", "Saithe", "Cod","Cattle"), col=c("red","navy","pink","orange"),pch=c(20,20),cex=2,pt.cex=3,bty="n")
dev.off()


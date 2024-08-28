## Set working directory and load data sheet
getwd()
setwd("C:/Users/annie/OneDrive/R")
humanreldat <- read.csv("ControlData.csv")

humanreldat$colour <- NA
FreshwaterHumanIndex <- which(humanreldat$Species %in% "FHuman")
C3HumanIndex <- which(humanreldat$Species %in% "Human")
C4HumanIndex <- which(humanreldat$Species %in% "C4Human")
MarineHumanIndex <- which(humanreldat$Species %in% "MHuman")
SaitheIndex <- which(humanreldat$Species %in% "Saithe")
CodIndex <- which(humanreldat$Species %in% "Cod")
CattleIndex <- which(humanreldat$Species %in% "Cattle")
HumanOrkIndex <- which(humanreldat$Species %in% "HumanOrk")
DeerIndex <- which(humanreldat$Species %in% "Deer")
HareIndex <- which(humanreldat$Species %in% "Hare")
MooseIndex <- which(humanreldat$Species %in% "Moose")
SealIndex <- which(humanreldat$Species %in% "Seal")
SeaLionIndex <- which(humanreldat$Species %in% "Sea lion")
DolphinIndex <- which(humanreldat$Species %in% "Dolphin")
PorpoiseIndex <- which(humanreldat$Species %in% "Porpoise")

humanreldat$colour[C3HumanIndex] <- "green"
humanreldat$colour[C4HumanIndex] <- "darkgreen"
humanreldat$colour[MarineHumanIndex] <- "purple"
humanreldat$colour[FreshwaterHumanIndex] <- "salmon"
humanreldat$colour[SaitheIndex] <- "navy"
humanreldat$colour[CodIndex] <- "navy"
humanreldat$colour[CattleIndex] <- "orange"
humanreldat$colour[HumanOrkIndex] <- "red"
humanreldat$colour[DeerIndex] <- "yellow"
humanreldat$colour[HareIndex] <- "yellow"
humanreldat$colour[MooseIndex] <- "yellow"
humanreldat$colour[SealIndex] <- "turquoise"
humanreldat$colour[SeaLionIndex] <- "turquoise"
humanreldat$colour[DolphinIndex] <- "turquoise"
humanreldat$colour[PorpoiseIndex] <- "turquoise"

#val by phe
png(file = "CVALbyCPHEControl.png", height = 720, width = 960)
par(mar=c(7, 6, 4, 15), xpd=TRUE)
plot(x=humanreldat$val_d13C,y=humanreldat$phe_d13C, xlab="",ylab=expression("δ"^13*"C"[Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(-32, -14),
     ylim = c(-32, -14),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(-32)),
                 round(max(-14)), by = 2))
axis(2, at = seq(round(min(-32)),
                 round(max(-14)), by = 2))
title(xlab=expression("δ"^13*"C"[Val]* "(‰)"), cex.lab=2, line=4.6)
legend("topright", inset=c(-0.31, 0.32),legend=c("Cairns Humans", "C3 Humans","C4 Humans", "Marine Humans","Freshwater Humans", "Terrestrial Fauna","Orkney Cattle", "Orkney Fish","Marine Mammals" ), col=c("red","green","darkgreen", "purple","salmon","yellow", "orange","navy","turquoise"),pch=c(20,20),cex=1.5,bty = "n", pt.cex=2.5)
dev.off()




#val-phe by lys-phe
png(file = "CvalphebylyspheControl.png", height = 720, width = 960)
par(mar=c(7, 6, 4, 15), xpd=TRUE)
plot(x=humanreldat$lysphe_d13C,y=humanreldat$valphe_d13C, xlab="",ylab=expression("δ"^13*"C"[Val-Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(2, 12),
     ylim = c(-4.5, 8),xaxt = "n",yaxt = "n",cex.lab=2, cex.axis=2)
axis(1, at = seq(round(min(2)),
                 round(max(12)), by = 2))
axis(2, at = seq(round(min(-4.5)),
                 round(max(8)), by = 2))
title(xlab=expression("δ"^13*"C"[Lys-Phe]* "(‰)"), cex.lab=2, line=4.6)
legend("topright", inset=c(-0.31, 0.32),legend=c("Cairns Humans", "C3 Humans","C4 Humans", "Marine Humans","Freshwater Humans", "Terrestrial Fauna","Orkney Cattle", "Orkney Fish","Marine Mammals" ), col=c("red","green","darkgreen", "purple","salmon","yellow", "orange","navy","turquoise"),pch=c(20,20),cex=1.5,bty = "n", pt.cex=2.5)
dev.off()

#2x1
png(file = "2x1Control.png", height = 900, width = 1782)
par(mfrow=c(1,2))
par(mar=c(8,8,4,2))
plot(x=humanreldat$val_d13C,y=humanreldat$phe_d13C, xlab="",ylab=expression("δ"^13*"C"[Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(-32, -14),
     ylim = c(-32, -14),xaxt = "n",yaxt = "n",cex.lab=3, cex.axis=2)
axis(1, at = seq(round(min(-32)),
                 round(max(-14)), by = 2))
axis(2, at = seq(round(min(-32)),
                 round(max(-14)), by = 2))
title(xlab=expression("δ"^13*"C"[Val]* "(‰)"), cex.lab=3, line=5.2)
plot(x=humanreldat$lysphe_d13C,y=humanreldat$valphe_d13C, xlab="",ylab=expression("δ"^13*"C"[Val-Phe]* "(‰)"), col=humanreldat$colour,  pch=20, cex=4,xlim = c(2, 12),
     ylim = c(-4.5, 8),xaxt = "n",yaxt = "n",cex.lab=3, cex.axis=2)
axis(1, at = seq(round(min(2)),
                 round(max(12)), by = 2))
axis(2, at = seq(round(min(-4.5)),
                 round(max(8)), by = 2))
title(xlab=expression("δ"^13*"C"[Lys-Phe]* "(‰)"), cex.lab=3, line=5.2)
legend("bottomright",legend=c("Cairns Humans", "C3 Humans","C4 Humans", "Marine Humans","Freshwater Humans", "Terrestrial Fauna","Orkney Cattle", "Orkney Fish","Marine Mammals" ), col=c("red","green","darkgreen", "purple","salmon","yellow", "orange","navy","turquoise"),pch=c(20,20),cex=1.5, pt.cex=2.5)
dev.off()

png(file = "2x1Control.png", height = 900, width = 1782)
par(mfrow = c(1, 2))
par(mar = c(8, 8, 4, 2))

# Plot 1
plot(x = humanreldat$val_d13C, y = humanreldat$phe_d13C, 
     xlab = "", ylab = expression("δ"^13*"C"[Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 4, 
     xlim = c(-32, -14), ylim = c(-32, -14), xaxt = "n", yaxt = "n", 
     cex.lab = 3, cex.axis = 2)
axis(1, at = seq(round(-32), round(-14), by = 2), cex.axis = 2)
axis(2, at = seq(round(-32), round(-14), by = 2), cex.axis = 2)
title(xlab = expression("δ"^13*"C"[Val]* "(‰)"), cex.lab = 3, line = 5.2)

# Plot 2
plot(x = humanreldat$lysphe_d13C, y = humanreldat$valphe_d13C, 
     xlab = "", ylab = expression("δ"^13*"C"[Val-Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 4, 
     xlim = c(2, 12), ylim = c(-4.5, 8), xaxt = "n", yaxt = "n", 
     cex.lab = 3, cex.axis = 2)
axis(1, at = seq(round(2), round(12), by = 2), cex.axis = 2)
axis(2, at = seq(round(-4.5), round(8), by = 2), cex.axis = 2)
title(xlab = expression("δ"^13*"C"[Lys-Phe]* "(‰)"), cex.lab = 3, line = 5.2)

# Add legend
legend("bottomright", legend = c("Cairns Humans", "C3 Humans", "C4 Humans", 
                                 "Marine Humans", "Freshwater Humans", 
                                 "Terrestrial Fauna", "Orkney Cattle", 
                                 "Orkney Fish", "Marine Mammals"), 
       col = "black", pt.bg = c("red", "green", "darkgreen", "purple", 
                                 "salmon", "yellow", "orange", "navy", 
                                 "turquoise"), pch = 21, cex = 1.5, 
       pt.cex = 2.5, bty = "n")

dev.off()

png(file = "2x1Control.png", height = 1800, width = 1000)  # Increased height for better vertical fit
par(mfrow = c(2, 1))  # Plots one on top of the other
par(mar = c(8, 8, 4, 2))  # Margins for the plots

# Plot 1
plot(x = humanreldat$val_d13C, y = humanreldat$phe_d13C, 
     xlab = "", ylab = expression("δ"^13*"C"[Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 4, 
     xlim = c(-32, -14), ylim = c(-32, -14), xaxt = "n", yaxt = "n", 
     cex.lab = 3, cex.axis = 2)
axis(1, at = seq(round(-32), round(-14), by = 2), cex.axis = 2)
axis(2, at = seq(round(-32), round(-14), by = 2), cex.axis = 2)
title(xlab = expression("δ"^13*"C"[Val]* "(‰)"), cex.lab = 3, line = 5.2)

# Plot 2
plot(x = humanreldat$lysphe_d13C, y = humanreldat$valphe_d13C, 
     xlab = "", ylab = expression("δ"^13*"C"[Val-Phe]* "(‰)"), 
     col = "black", bg = humanreldat$colour, pch = 21, cex = 4, 
     xlim = c(2, 12), ylim = c(-4.5, 8), xaxt = "n", yaxt = "n", 
     cex.lab = 3, cex.axis = 2)
axis(1, at = seq(round(2), round(12), by = 2), cex.axis = 2)
axis(2, at = seq(round(-4.5), round(8), by = 2), cex.axis = 2)
title(xlab = expression("δ"^13*"C"[Lys-Phe]* "(‰)"), cex.lab = 3, line = 5.2)

# Add legend to the second plot
legend("bottomright", legend = c("Cairns Humans", "C3 Humans", "C4 Humans", 
                                 "Marine Humans", "Freshwater Humans", 
                                 "Terrestrial Fauna", "Orkney Cattle", 
                                 "Orkney Fish", "Marine Mammals"), 
       col = "black", pt.bg = c("red", "green", "darkgreen", "purple", 
                                "salmon", "yellow", "orange", "navy", 
                                "turquoise"), pch = 21, cex = 2, 
       pt.cex = 3, bty = "n")

dev.off()


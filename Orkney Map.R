## Set working directory and load data sheet
getwd()
setwd("C:/Users/annie/OneDrive/R/Summative")
library(dunn.test)
library(sf)
library(giscoR)
library(ggplot2)
library("ggspatial")
humanreldat <- read.csv("HumanData.csv")

### Create a map of sites where the human data is from
#Save sites with x and y coordinates as a data frame
Sites <- data.frame(ID=c("The Cairns","BNKS","Knap of Howar","Ness of Brodgar"), x=c(-2.945228,-2.9813,-2.91085,-3.22333244), y=c(58.76611,59.2591,59.34935,59.0))
#Download map of Orkney
Orkney <- gisco_get_nuts(epsg = 4326, resolution = "01", nuts_id = "UKM65")
#Save map as a png 
png(file = "Orkneymapdiss3.png", height = 2000, width = 2000, units = "px", res=300)
#Make a ggplot map with site points, labels and a title
ggplot(Orkney) + geom_sf(fill="navajowhite2") + annotation_scale(location = "bl", width_hint = 0.5) + annotation_north_arrow(location = "bl", which_north = "true", pad_y = unit(1, "cm"),style = north_arrow_fancy_orienteering) + coord_sf(xlim = c(-3.6, -2.3), ylim = c(59.45, 58.6)) + geom_point(data=Sites, aes(x=x, y=y), cex=2) + geom_text(label="The Cairns", x=-2.75,y=58.76611, col="black",cex=5) + geom_text(label="BNKS", x=-3.1,y=59.263, col="black", cex=5) + geom_text(label="Knap of Howar", x=-3.15,y=59.34935, col="black",cex=5) + geom_text(label="Ness of Brodgar", x=-2.97,y=59.0, col="black",cex=5) + ggtitle(label="") + labs(x=NULL,y=NULL) + theme(panel.background=element_rect('lightblue')) + theme(panel.grid.major=element_line(colour='lightblue'))
dev.off()
library(googlesheets4)
library(ggplot2)
library (devtools)
library(tidyverse)
library(dplyr)

getwd()
setwd("C:/Users/annie/OneDrive/R")
Larsen <- read.csv("Larsen.csv")

# Load colour palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

Larsen_pp <- filter(Larsen, PP == "Y")
#Select groups and amino acids
Larsen_pp <- Larsen_pp %>%  dplyr::select(Group, Ala, Asp, Glx, Gly, Ile, Thr, Leu, Lys, Val, Phe)
Larsen_pp <- na.omit(Larsen_pp) #remove N/As

all_AA <- Larsen_pp %>%  dplyr::select(Ala, Asp, Glx, Gly, Ile, Thr, Leu, Lys, Val, Phe)
PCA <- prcomp(all_AA, scale. = TRUE)
# Extract PC axes for plotting
PCAvalues <- data.frame(group = Larsen_pp$Group, PCA$x)
# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
# Plot PC1 vs PC2
p <- ggplot (data=PCAvalues, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size=3) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10),
                                       yend = (PC2*10)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  annotate ("text", x = (PCAloadings$PC1*11), y = (PCAloadings$PC2*11),
            label = PCAloadings$Variables, )+
  xlim(5,-5)+
  ylim(5,-5)+
  coord_fixed (ratio = 1)
p
dev.off()
#Calculate the AAn values which equals each values variance from the mean and repeat the analysis.

#First organise the data and normalise to row mean.
Larsen_pp <- filter(Larsen, PP == "Y")
Larsen_pp <- Larsen_pp %>%  dplyr::select(Group, Ala, Asp, Glx, Gly, Ile, Thr, Leu, Lys, Val, Phe)
Larsen_pp <- na.omit(Larsen_pp) #
#Create a new coloumn containing the row mean of the key amino acid for each row. .
Larsen_pp <- Larsen_pp %>% rowwise() %>% mutate(row_mean = mean(c(Ala, Asp, Glx, Gly, Ile, Thr, Leu, Lys, Val, Phe)))#
#Subtract the row mean from each of the row values to normalize
Larsen_pp_n <-Larsen_pp %>% rowwise() %>%
  mutate(across(c(Ala, Asp, Glx, Gly, Ile, Thr, Leu, Lys, Val, Phe), ~ .x - row_mean))


#10. Then perform the PCA and visualise the results. This creates Figure 1 in Larsen et al.


#extract the amino acid values
all_AA <- Larsen_pp_n %>%  dplyr::select(Ala, Asp, Glx, Gly, Ile, Thr, Leu, Lys, Val, Phe))
PCA <- prcomp(all_AA, scale. = TRUE)
# Extract PC axes for plotting
PCAvalues <- data.frame(group = Larsen_pp_n$Group, PCA$x)
# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
# Visuqlise the PCA
p <- ggplot (data=PCAvalues, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size=3) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10),
                                       yend = (PC2*10)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  annotate ("text", x = (PCAloadings$PC1*11), y = (PCAloadings$PC2*11),
            label = PCAloadings$Variables, )+
  xlim(5,-5)+
  ylim(5,-5)+
  coord_fixed (ratio = 1)
p
dev.off()

p2 <- ggplot(data = PCAvalues, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 3) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1 * 10), yend = (PC2 * 10)), 
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  annotate("text", x = (PCAloadings$PC1 * 11), y = (PCAloadings$PC2 * 11), 
           label = PCAloadings$Variables) +
  scale_color_manual(values = c(
    "Cattle" = "orange", 
    "Human" = "red", 
    "Saithe" = "darkblue", 
    "Cod" = "pink", 
    "Terr_Plants" = "darkgreen", 
    "Seagrasses" = "green", 
    "Bacteria" = "brown", 
    "Algae" = "blue"
  )) +
  xlim(-5, 5) +
  ylim(-5, 5) +
  coord_fixed(ratio = 1)
p2
dev.off()

#11. Repeat PCA just for source aminos

Larsen_ppS <- filter(Larsen, PP == "Y")
Larsen_ppS <- Larsen_pp %>%  dplyr::select(Group, Ile, Leu, Lys, Val, Phe)
Larsen_ppS <- na.omit(Larsen_pp) #
#Create a new coloumn containing the row mean of the key amino acid for each row. .
Larsen_ppS <- Larsen_ppS %>% rowwise() %>% mutate(row_mean = mean(c(Ile, Leu, Lys, Val, Phe)))#
#Subtract the row mean from each of the row values to normalize
Larsen_pp_Sn <-Larsen_pp %>% rowwise() %>%
  mutate(across(c(Ile, Leu, Lys, Val, Phe), ~ .x - row_mean))
#Conduct the PCA
all_AA <- Larsen_pp_Sn %>%  dplyr::select(Ile, Leu, Lys, Val, Phe)
PCA1 <- prcomp(all_AA, scale. = TRUE)
# Extract PC axes for plotting
PCAvalues1 <- data.frame(group = Larsen_pp_Sn$Group, PCA1$x)
# Extract loadings of the variables
PCAloadings1 <- data.frame(Variables = rownames(PCA1$rotation), PCA1$rotation)
# Visuqlise the PCA
p3 <- ggplot (data=PCAvalues1, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size=3) +
  geom_segment(data = PCAloadings1, aes(x = 0, y = 0, xend = (PC1*10),
                                        yend = (PC2*10)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  annotate ("text", x = (PCAloadings1$PC1*11), y = (PCAloadings1$PC2*11),
            label = PCAloadings1$Variables, )+
  xlim(10,-10)+
  ylim(10,-10)+
  coord_fixed (ratio = 1)
p3
dev.off()
png(file = "LarsenPCAEA.png", height = 2000, width = 2500, units = "px", res=300)
par(mar=c(4,4,1,0))
ggplot(data = PCAvalues1, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 3) +
  geom_segment(data = PCAloadings1, aes(x = 0, y = 0, xend = (PC1 * 10), yend = (PC2 * 10)), 
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  annotate("text", x = (PCAloadings1$PC1 * 11), y = (PCAloadings1$PC2 * 11), 
           label = PCAloadings1$Variables) +
  scale_color_manual(values = c(
    "Cattle" = "orange", 
    "Human" = "red", 
    "Saithe" = "black", 
    "Cod" = "hotpink", 
    "Terr_Plants" = "darkgreen", 
    "Seagrasses" = "green", 
    "Bacteria" = "brown", 
    "Algae" = "blue"
  )) +
  xlim(-10, 10) +
  ylim(-10, 10) +
  coord_fixed(ratio = 1)
dev.off()

# Calculate the variance explained by each principal component
eigenvalues <- PCA1$sdev^2

# Calculate the percentage variance explained by each PC
percent_variance_explained <- (eigenvalues / sum(eigenvalues)) * 100

# Print the percentage variance explained by each PC
print(percent_variance_explained)

# Conduct the PCA
all_AA <- Larsen_pp_Sn %>%  dplyr::select(Ile, Leu, Lys, Val, Phe)
PCA1 <- prcomp(all_AA, scale. = TRUE)

# Calculate the variance explained by each principal component
eigenvalues <- PCA1$sdev^2

# Calculate the percentage variance explained by each PC
percent_variance_explained <- (eigenvalues / sum(eigenvalues)) * 100

# Print the percentage variance explained by each PC
print(percent_variance_explained)

# Extract PC axes for plotting
PCAvalues1 <- data.frame(group = Larsen_pp_Sn$Group, PCA1$x)

# Extract loadings of the variables
PCAloadings1 <- data.frame(Variables = rownames(PCA1$rotation), PCA1$rotation)

# Visualize the PCA
p3 <- ggplot(data = PCAvalues1, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 3) +
  geom_segment(data = PCAloadings1, aes(x = 0, y = 0, xend = (PC1 * 10), yend = (PC2 * 10)), 
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  annotate("text", x = (PCAloadings1$PC1 * 11), y = (PCAloadings1$PC2 * 11), 
           label = PCAloadings1$Variables) +
  scale_color_manual(values = c(
    "Cattle" = "orange", 
    "Human" = "red", 
    "Saithe" = "black", 
    "Cod" = "hotpink", 
    "Terr_Plants" = "darkgreen", 
    "Seagrasses" = "green", 
    "Bacteria" = "brown", 
    "Algae" = "blue"
  )) +
  xlim(-10, 10) +
  ylim(-10, 10) +
  coord_fixed(ratio = 1)

p3
dev.off()

# Calculate the variance explained by each principal component
eigenvalues <- PCA1$sdev^2
percent_variance_explained <- (eigenvalues / sum(eigenvalues)) * 100

# Create axis labels with variance explained
x_label <- paste0("PC1 (", round(percent_variance_explained[1], 2), "%)")
y_label <- paste0("PC2 (", round(percent_variance_explained[2], 2), "%)")

# Visualize the PCA with variance in labels
png(file = "LarsenPCAEA.png", height = 2000, width = 2500, units = "px", res=300)
par(mar=c(4,4,1,0))
ggplot(data = PCAvalues1, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 3) +
  geom_segment(data = PCAloadings1, aes(x = 0, y = 0, xend = (PC1 * 10), yend = (PC2 * 10)), 
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  annotate("text", x = (PCAloadings1$PC1 * 11), y = (PCAloadings1$PC2 * 11), 
           label = PCAloadings1$Variables) +
  scale_color_manual(values = c(
    "Cattle" = "orange", 
    "Human" = "red", 
    "Saithe" = "black", 
    "Cod" = "hotpink", 
    "Terr_Plants" = "darkgreen", 
    "Seagrasses" = "green", 
    "Bacteria" = "brown", 
    "Algae" = "blue"
  )) +
  labs(x = x_label, y = y_label) +  # Use the new labels here
  xlim(-10, 10) +
  ylim(-10, 10) +
  coord_fixed(ratio = 1)

dev.off()

png(file = "LarsenPCAEA.png", height = 2000, width = 2500, units = "px", res=300)
par(mar=c(4,4,1,0))
ggplot(data = PCAvalues1, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 3) +
  geom_segment(data = PCAloadings1, aes(x = 0, y = 0, xend = (PC1 * 10), yend = (PC2 * 10)), 
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  annotate("text", x = (PCAloadings1$PC1 * 11), y = (PCAloadings1$PC2 * 11), 
           label = PCAloadings1$Variables) +
  scale_color_manual(values = c(
    "Cattle" = "orange", 
    "Human" = "red", 
    "Saithe" = "black", 
    "Cod" = "hotpink", 
    "Terr_Plants" = "darkgreen", 
    "Seagrasses" = "green", 
    "Bacteria" = "brown", 
    "Algae" = "blue"
  )) +
  labs(x = x_label, y = y_label) +  # Use the new labels here
  xlim(-10, 10) +
  ylim(-10, 10) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

dev.off()


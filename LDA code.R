library(googlesheets4)
library(ggplot2)
library (devtools)
library(tidyverse)
library(dplyr)

#12. LDA analysis

Prepare Larsen data for LDA analysis
```{r}
#Select just amino acids for PCA
library(caret)
library(tidyverse)
#Create a normalised training set from Larsen's primary producters
gs4_deauth()
Larsen <<- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1Ugkt04-R1sxSgnc3U9tIR-MpJkI34AtQgu-rdAjlyEE/edit?gid=0#gid=0", sheet="Larsen")
#Subset non primary producers
Larsen_pp <- filter(Larsen, PP == "Y")
#Select group and amino acids
Larsen_pp <- Larsen_pp %>%  dplyr::select(Group, ile_d13C, leu_d13C, lys_d13C, val_d13C, phe_d13C)
Larsen_pp <- na.omit(Larsen_pp)
Larsen_pp$Group  <-as.factor (Larsen_pp$Group)
#omit Fungi
#Larsen_pp <-Larsen_pp %>% filter(Group!="Fungi")

```


Normalise the amino acids in Larsen reference dataset
```{r}
#Create a new coloumn containing the row mean of the key amino acid for each row. .
Larsen_pp <- Larsen_pp %>% dplyr::rowwise() %>% dplyr::mutate(row_mean = mean(c(ile_d13C, leu_d13C, lys_d13C, val_d13C, phe_d13C)))
#Subtract the row mean from each of the row values to normalize
Larsen_pp_n <-Larsen_pp %>% dplyr::rowwise() %>% dplyr::mutate(across(c( ile_d13C, leu_d13C,  lys_d13C, val_d13C, phe_d13C), ~ .x - row_mean))
Larsen_pp_n <- Larsen_pp_n %>%  dplyr::select(Group, ile_d13C, leu_d13C, lys_d13C, val_d13C, phe_d13C)#drop means
print (Larsen_pp_n)

```

Partition the data and create train and test sets

```{r}
ind <- sample(2, nrow(Larsen_pp_n),
              replace = TRUE,
              prob = c(0.9, 0.1))
train.data <- Larsen_pp_n[ind==1,]
test.data <- Larsen_pp_n[ind==2,]
```

Create LDA model
```{r}
library(MASS)
library(caret)
# Fit LDA model to the training data
lda_model <- lda(Group ~ ., data = train.data)
```
Evaluate the LDA  performance
```{r}
# Predict on the test data using the LDA model
lda_predictions <- predict(lda_model, test.data)

#lda_predictions$class
confusion_matrix <- confusionMatrix(lda_predictions$class, test.data$Group)
print(confusion_matrix)
```
all_data <<- read_sheet(ss="https://docs.google.com/spreadsheets/d/1Ugkt04-R1sxSgnc3U9tIR-MpJkI34AtQgu-rdAjlyEE/edit?gid=1734888199#gid=1734888199", sheet="all_data")
                        

Call archaeological data and normalise
```{r}
data <- all_data %>% dplyr::rowwise() %>% dplyr::mutate(row_mean = mean(c(ile_d13C, leu_d13C, lys_d13C, val_d13C, phe_d13C)))
#Subtract the row mean from each of the row values to normalize
data <- data %>% dplyr::rowwise() %>% dplyr::mutate(across(c( ile_d13C, leu_d13C,  lys_d13C, val_d13C, phe_d13C), ~ .x - row_mean))
data_n <- data %>%  dplyr::select(category, ile_d13C, leu_d13C, lys_d13C, val_d13C, phe_d13C)#drop means
data_n <- data_n %>% rename(Group=category)
data_n$Group <-as.factor (data_n$Group)
print (data_n)
```

```{r}
# Predict the classes for the unknown samples using the LDA model
lda_Danish_predictions <- predict(lda_model, newdata = data_n)
#Print the predicted classes for the unknown samples
predicted_Danish_classes <- lda_Danish_predictions$class
#print (predicted_Danish_classes)
```
Create two dataframes for plotting and combine

```{r}
Refs <- data.frame(LD1 = predict(lda_model, newdata = train.data)$x[, 1], LD2 = predict(lda_model, data = train.data)$x[, 2], Type = "Training", group = as.factor(train.data$Group))
Samples <- data.frame (LD1 = lda_Danish_predictions$x[, 1], LD2 = lda_Danish_predictions$x[, 2], Type = "Samples", group= data_n$Group)
Combined <- rbind (Refs, Samples)
#print (Combined)
```



#Find means for automatic labelling
``` {r}
#Find means for automatic labelling
means <- Refs %>% group_by(group) %>%
  summarise(dplyr::across(c(LD1, LD2), mean))
means
#means <- means %>%mutate(cats=c("Porcine", "Ruminant", "Dairy", "Freshwater", "Marine"))#Sh
```


Plot the LDA
```{r}
# Plot the combined data
p <- ggplot() +
  #geom_point(data=Refs, aes(x =LD1, y =LD2, fill = group), pch=22, colour = "black", size = 2) +
  stat_ellipse(geom= "polygon", data=Refs, alpha=0.1,
               aes (x = LD1, y = LD2, group = group),
               level = 0.68)+
  geom_point(data=Samples, aes(x =LD1, y =LD2, fill = group), pch=21, colour = "black", size = 3) +
  geom_text(data=means, aes(y=LD2, x=LD1,label=group), hjust=0, size=4, color="black")+
  labs(title = "LDA Plot", x = "Linear Discriminant 1", y = "Linear Discriminant 2", fill = "Group", shape = "Type") +
  theme_bw()

p
dev.off()
getwd()
setwd("C:/Users/annie/OneDrive/R")

png(file = "LDAanalysis.png", height = 1500, width = 2000, units = "px", res=300)
par(mar=c(4,4,1,0))
ggplot() +
  #geom_point(data=Refs, aes(x =LD1, y =LD2, fill = group), pch=22, colour = "black", size = 2) +
  stat_ellipse(geom= "polygon", data=Refs, alpha=0.1,
               aes (x = LD1, y = LD2, group = group),
               level = 0.68)+
  geom_point(data=Samples, aes(x =LD1, y =LD2, fill = group), pch=21, colour = "black", size = 3) +
  geom_text(data=means, aes(y=LD2, x=LD1,label=group), hjust=0, size=4, color="black")+
  labs(title = "LDA Plot", x = "Linear Discriminant 1", y = "Linear Discriminant 2", fill = "Group", shape = "Type") +
  theme_bw()
dev.off()

png(file = "LDAanalysisbigger.png", height = 2000, width = 2500, units = "px", res=300)
par(mar=c(4,4,1,0))
ggplot() +
  # Add ellipses for reference data
  stat_ellipse(geom= "polygon", data=Refs, alpha=0.1,
               aes(x = LD1, y = LD2, group = group),
               level = 0.68) +
  # Add sample points with custom colors
  geom_point(data=Samples, aes(x = LD1, y = LD2, fill = group), pch=21, colour = "black", size = 3) +
  # Add group labels
  geom_text(data=means, aes(y = LD2, x = LD1, label = group), hjust = 0, size = 4, color = "black") +
  # Custom color scale
  scale_fill_manual(values = c("Human" = "red", "Cattle" = "orange", "Saithe" = "blue3", "Cod" = "pink")) +
  # Add titles and labels
  labs(title = "LDA Plot", x = "Linear Discriminant 1", y = "Linear Discriminant 2", fill = "Group", shape = "Type") +
  theme_bw()
dev.off()

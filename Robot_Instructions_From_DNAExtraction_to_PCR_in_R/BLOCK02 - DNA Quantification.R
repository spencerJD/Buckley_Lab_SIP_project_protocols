library(plyr)
library(ggplot2)

## Input Data
x <- read.table(file = "pico.green.soil.health.allplates.output", sep = "\t", stringsAsFactors = F, header = T)

## User Input
dilution_factor = 100

# Calculate Standard Curve and Calculate Concentrations
x$Adj.Conc <- NA
for (plate in unique(x$Plate)){
  standards <- subset(x, Plate == plate & Type == "Standard")
  
  # Get regression coefficients
  b = coef(lm(Concentration ~ RFU, data=standards))[1]
  m = coef(lm(Concentration ~ RFU, data=standards))[2]
  
  # Correct Concentration in Main Datasheet (note: the "Concentration" is supposed to be the adjusted )
  x[which(x$Plate == plate),"Adj.Conc"] <- (m*x[which(x$Plate == plate),"RFU"]+b)*dilution_factor
}

# Remove Standards
x <- subset(x, Type != "Standard")

# Breakup Sample Factors
x$Combo <- x$Plate
for (n in 1:nrow(x)){
  x$Group[n] <- unlist(strsplit(x$Plate[n], ":"))[2]
  x$Plate[n] <- unlist(strsplit(x$Plate[n], ":"))[1]
}

# Fix Plate Names
x$Plate <- gsub("Plat3","Plate3", x$Plate)
x$Plate <- gsub("Pla3","Plate3", x$Plate)
x$Plate <- gsub("Plat4","Plate4", x$Plate)

# Average replicates
x <- ddply(x, ~ Plate + Sample + Group, summarise, Avg.Conc = mean(Adj.Conc))
head(x)

# Set negative values to zero
x$Avg.Conc[which(x$Avg.Conc < 0)] <- 0

# Get Rid of half-plate when one uses the second half of a plate for a separate set of smamples (i.e. keep only samples 1-24 on the last set of a plate, from 86-96)
x <- x[-which(x$Sample > 16 & x$Group == "81-96"),]

# Plot Histogram
qplot(x$Avg.Conc, geom="histogram")

# How Many Low DNA
length(which(x$Avg.Conc<=1))
length(which(x$Avg.Conc>1))

# Look for Plate Effects
# by plate
ddply(x, ~ Plate, summarise, average = mean(Avg.Conc))
TukeyHSD(aov(lm(Avg.Conc ~ Plate, x)))

# by Pico Green group
ddply(x, ~ Plate + Group, summarise, average = mean(Avg.Conc))
TukeyHSD(aov(lm(Avg.Conc ~ Plate*Group, x)))

# Mean Average Concentration
mean(x$Avg.Conc)

# Associate SampleIDs with DNA concentrations
access <- read.csv(file="sample2data.PicoGreen.csv",stringsAsFactors = F, header=T)
x <- merge(x, access, by =c("Plate","Group","Sample"))
colnames(x) <- c("DNA_Plate","Pico_Group","Pico_Sample","DNA","DNA_Well","DNA_Sample","SampleName","SampleID")
saveRDS(x, file="post.pico.plates1-4.rds")

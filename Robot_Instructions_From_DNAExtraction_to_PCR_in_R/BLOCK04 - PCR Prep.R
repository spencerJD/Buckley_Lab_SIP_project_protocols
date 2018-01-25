library(reshape2)
library("xlsx")  # Make sure you have Java 64 bit installed or 'RJava' will fail to load

make.plate <- function(){
  ids <- c(LETTERS, letters)
  n <- 96
  nrow <- 8
  samples <- character(n)
  samples[seq_along(ids)] <- ids
  samples <- matrix(samples, nrow=nrow)
  colnames(samples) <- seq_len(n/nrow)
  rownames(samples) <- LETTERS[seq_len(nrow)]
  samples <- melt(samples)
  samples$position <- paste0(samples$Var1, samples$Var2)
  samples <- samples[,c("position")]
  
  return(samples)
}

############################################################
## Prepare PCR Instructions for Buckley Lab Robot PCR Method

### This script was written to take a 96-well plate with 94 samples (wells H12 and G12 are positive and negative control DNA)
##  and generate PCR plates that include the following:
# i) ~32 samples in duplicate
# ii) qPCR standard curve in triplicate (0 (i.e. pure water), X, X, X, X, X, X, X dilutions).

## The qPCR Standard is added separately to each plate after the samples have been arrayed.
## ~20 uL of each qPCR dilution standard needs to be arrayed in a separate 96 well plate in A1 -> H1
## NOTE: Include pure water as your no-template standard

#### PART I: Arraying Samples
####

# Read in Data
x <- readRDS(file="diluted.template.1-4.rds")
x <- subset(x, FinalTemplateVol != "NA")

# Import PCR arraying database
array <- readRDS(file = "pcr.plate.array.duplicate.rds")

# Initialize other parameters
generic_plate <- make.plate()
x$PCRWell <- NA
x$PCR_Plate <- NA

##############################################################
# Make 'qPCR_plate_layout" by assigning duplicate target wells

for (plate in unique(x$DILUTION_Plate)){
  count = 1
  plate_counter = 1
  
  # Subset Plate
  foo <- subset(x, DILUTION_Plate == plate)
  
  # Run Through All Samples
  for (n in 1:nrow(foo)){
    if (count == 1){
      reps <- as.vector(array[n,1:2])
      PCR_sub_plate <- array[n,3]
      row <- rep(substr(reps[1],1,1),2)
      columns <- gsub(row[1], "", reps)
      PCR.df3 <- data.frame(Row = row, Column = columns, Target_Name = "", Sample_Name = foo[n,"SampleID"], stringsAsFactors = F)

    } else {
      reps <- as.vector(array[n,1:2])
      PCR_plate <- array[n,3]
      row <- rep(substr(reps[1],1,1),2)
      columns <- gsub(row[1], "", reps)
      PCR.df3 <- rbind(PCR.df3, data.frame(Row = row, Column = columns, Target_Name = "", Sample_Name = foo[n,"SampleID"], stringsAsFactors = F))
    }

    # Save PCR Location to Main Dataframe
    if (is.na(PCR_sub_plate) == TRUE){ # Some samples will have to be re-arrayed and done again (flaw upstream, I should only allow 90 samples per dilution plate)
      x[grep(paste("^",foo[n,"SampleID"],"$",sep=""), x[,"SampleID"]),"PCR_Plate"] <- NA
      x[grep(paste("^",foo[n,"SampleID"],"$",sep=""), x[,"SampleID"]),"PCRWell"] <- NA
    } else {
      x[grep(paste("^",foo[n,"SampleID"],"$",sep=""), x[,"SampleID"]),"PCR_Plate"] <- paste(plate,PCR_sub_plate,sep=".")
      x[grep(paste("^",foo[n,"SampleID"],"$",sep=""), x[,"SampleID"]),"PCRWell"] <- paste(row[1],columns[1],sep="")
    }
    
    # Punch Out Each PCR Subplate after 32 samples (ie. 1 x 94 sample plate) or when you run out of samples
    if (count == 32 | nrow(foo) == plate_counter){
      
      # Write Output
      #write.csv(PCR.df3, file = paste(plate,PCR_sub_plate,"PCR.plate_layout.csv",sep="."), row.names = F)
      write.xlsx(PCR.df3, file = paste(plate,PCR_sub_plate,"PCR.xlsx",sep="."), sheetName= "qPCR_plate_layout", col.names= T, row.names = F)
      
      # Add Pre-fab Output for Adding Standards to Plate
      #write.csv(readRDS(file = "qPCR.plate_layout.rds"), file = paste(plate,PCR_sub_plate,"qPCR.standards.plate_layout.csv",sep="."), row.names = F)
      write.xlsx(readRDS(file = "qPCR.plate_layout.rds"), file = paste(plate,PCR_sub_plate,"standards.xlsx",sep="."), sheetName= "qPCR_plate_layout", col.names= T, row.names = F)
      
      # Re-set Counter
      count = 1
      
    } else {
      # Update Counter
      count = count + 1
    }
    
    # Keep tally of total samples array in order to catch the
    plate_counter = plate_counter + 1
  }    
}

saveRDS(x, file="PCR.plates.1-4.rds")


###########################
## Prep 'sample_primer_map'

# Import Data
x <- readRDS(file = "PCR.plates.1-4.rds")

# Read in Primer Plate Well ID
primer_plate1 <- data.frame(PrimerWell = make.plate(), PrimerID = seq(1, 96, by = 1), stringsAsFactors = F)
primer_plate2 <- data.frame(PrimerWell = make.plate(), PrimerID = seq(97, 192, by = 1), stringsAsFactors = F)
primer_split1 <- data.frame(PrimerWell = make.plate(), PrimerID = c(seq(49, 96, by = 1),seq(97, 144, by = 1)), stringsAsFactors = F)
primer_split2 <- data.frame(PrimerWell = make.plate(), PrimerID = c(seq(145, 192, by = 1),seq(1, 48, by = 1)), stringsAsFactors = F)

# Set Primer Numbers
double_primers <- rep(seq(1,round(nrow(x)/2),by=1),2)
x$PrimerNumber <- double_primers[1:nrow(x)]

for (plate in unique(x$PCR_Plate)){
  count = 1
  
  # Subset Plate
  foo <- subset(x, PCR_Plate == plate)
  
  # Get Primer Plate Information in order to Determine Primer Well
  if (length(which(foo$PrimerNumber %in% primer_plate1$PrimerID)) == nrow(foo)){
    primer_plate <- primer_plate1
    which_plate <- "PrimerPlate1"
  } else if (length(which(foo$PrimerNumber %in% primer_split1$PrimerID)) == nrow(foo)){
    primer_plate <- primer_split1
    which_plate <- "SplitPlate1"
  } else if (length(which(foo$PrimerNumber %in% primer_split2$PrimerID)) == nrow(foo)){
    primer_plate <- primer_split2
    which_plate <- "SplitPlate2"
  } else if (length(which(foo$PrimerNumber %in% primer_plate2$PrimerID)) == nrow(foo)){
    primer_plate <- primer_plate2
    which_plate <- "PrimerPlate2"
  }
  
  # Run Through All Samples
  for (n in 1:nrow(foo)){
    PrimerWell <- primer_plate[which(primer_plate$PrimerID == foo[n,"PrimerNumber"]),"PrimerWell"]
    
    if (count == 1){
      PCR.df2 <- data.frame(Sample_Number = foo[n,"SampleID"], Sample_Well_ID = foo[n,"DILUTION_Well"], ul_to_add = foo[n,"FinalTemplateVol"], Primer_Number = foo[n,"PrimerNumber"], Primer_Well_ID = PrimerWell, ul_water_to_add = foo[n,"water_to_add"])
    } else {
      PCR.df2 <- rbind(PCR.df2, data.frame(Sample_Number = foo[n,"SampleID"], Sample_Well_ID = foo[n,"DILUTION_Well"], ul_to_add = foo[n,"FinalTemplateVol"], Primer_Number = foo[n,"PrimerNumber"], Primer_Well_ID = PrimerWell, ul_water_to_add = foo[n,"water_to_add"]))
    }
    
    # Update Counters
    count = count + 1
  }
  
  # Determine which primer plate to use
  PCR.df2$PrimerPlate <- which_plate
  
  # Write Output
  #write.csv(PCR.df2, file = paste(plate,"PCR.sample_primer_map.csv",sep="."), row.names = F)
  write.xlsx(PCR.df2, file = paste(plate,"PCR.xlsx",sep="."), sheetName= "sample_primer_map", col.names= T, row.names = F, append=TRUE)
  
  # Select 8 Primers from Current Plate at random
  random_primers <- c(sample(min(PCR.df2$Primer_Number):max(PCR.df2$Primer_Number), 1),sample(min(PCR.df2$Primer_Number):max(PCR.df2$Primer_Number), 1),sample(min(PCR.df2$Primer_Number):max(PCR.df2$Primer_Number), 1),sample(min(PCR.df2$Primer_Number):max(PCR.df2$Primer_Number), 1),sample(min(PCR.df2$Primer_Number):max(PCR.df2$Primer_Number), 1),sample(min(PCR.df2$Primer_Number):max(PCR.df2$Primer_Number), 1),sample(min(PCR.df2$Primer_Number):max(PCR.df2$Primer_Number), 1),sample(min(PCR.df2$Primer_Number):max(PCR.df2$Primer_Number), 1)) 
  random_primer_wells <- c(primer_plate[which(primer_plate$PrimerID == random_primers[1]),"PrimerWell"],primer_plate[which(primer_plate$PrimerID == random_primers[2]),"PrimerWell"],primer_plate[which(primer_plate$PrimerID == random_primers[3]),"PrimerWell"],primer_plate[which(primer_plate$PrimerID == random_primers[4]),"PrimerWell"],primer_plate[which(primer_plate$PrimerID == random_primers[5]),"PrimerWell"],primer_plate[which(primer_plate$PrimerID == random_primers[6]),"PrimerWell"],primer_plate[which(primer_plate$PrimerID == random_primers[7]),"PrimerWell"],primer_plate[which(primer_plate$PrimerID == random_primers[8]),"PrimerWell"])
  
  qPCR_standards <- data.frame(Sample_Number = seq(1,8,by=1), Sample_Well_ID = c("A1","B1","C1","D1","E1","F1","G1","H1"), ul_to_add = 5, Primer_Number = random_primers, Primer_Well_ID = random_primer_wells, ul_water_to_add = 4.4)
  
  # Output qPCR Standard Array Instructions
  #write.csv(qPCR_standards, file = paste(plate,"qPCR.standards.sample_primer_map.csv",sep="."), row.names = F)
  write.xlsx(qPCR_standards, file = paste(plate,"standards.xlsx",sep="."), sheetName= "sample_primer_map", col.names= T, row.names = F, append=TRUE)
}


##########################
## Prep 'DataForRobot' Tab

# Import Data
x <- readRDS(file = "PCR.plates.1-4.rds")

# Initialize other parameters
master_mix = 13.1  # uL

for (plate in unique(x$PCR_Plate)){
  count = 1
  
  # Subset Plate
  foo <- subset(x, PCR_Plate == plate)

  # Determine Final Sample Well
  final_well <- foo$PCRWell[length(foo$PCRWell)]
  final_well <- paste(substr(final_well[1],1,1), as.numeric(substr(final_well[1],2,2))+1,sep="") # Add a column to the last PCRWell (since the PCR well is only the first position of the duplicate)
  
  # Array Master Mix
  array_plate <- make.plate()
  array_plate <- array_plate[1:which(array_plate == final_well)]
  for (well in array_plate){
    if (count == 1){
      PCR.df1a <- data.frame(Name = "MasterMix_1", Type = "Mix", TargetWell = well, Volume = master_mix)
      count = count + 1 
    } else {
      PCR.df1a <- rbind(PCR.df1a, data.frame(Name = "MasterMix_1", Type = "Mix", TargetWell = well, Volume = master_mix))
    }
  }
  
  # Array DNA Template
  count = 1

  for (well in unique(foo$PCRWell)){
    well_row <- substr(well[1],1,1)
    well_column <- as.numeric(substr(well[1],2,2))
    
    # Add duplicate sample  
    sampleID <- foo[which(foo$PCRWell == well),"SampleID"]
    template_volume <- foo[which(foo$PCRWell == well),"FinalTemplateVol"]
    duplicate_well <- paste(well_row, well_column+1, sep="")

    if (count == 1){
      PCR.df1b <- data.frame(Name = sampleID, Type = "Sample", TargetWell = well, Volume = template_volume, WellRow = well_row, WellColumn = well_column, stringsAsFactors = F)
      PCR.df1b <- rbind(PCR.df1b, data.frame(Name = sampleID, Type = "Sample", TargetWell = duplicate_well, Volume = template_volume, WellRow = well_row, WellColumn = well_column+1))
      count = count + 1  
    } else {
      PCR.df1b <- rbind(PCR.df1b, data.frame(Name = sampleID, Type = "Sample", TargetWell = well, Volume = template_volume, WellRow = well_row, WellColumn = well_column))
      PCR.df1b <- rbind(PCR.df1b, data.frame(Name = sampleID, Type = "Sample", TargetWell = duplicate_well, Volume = template_volume, WellRow = well_row, WellColumn = well_column+1))
    }
  }
  
  # Order DNA Template Dataframe
  PCR.df1b <- PCR.df1b[order(PCR.df1b$WellColumn),]
  PCR.df1b <- PCR.df1b[,1:4]

  # Bind together two dataframes for the 'DataForRobot' tab
  PCR.df1 <- rbind(PCR.df1a, PCR.df1b)
  
  # Write Output for Samples
  #write.csv(PCR.df1, file = paste(plate,"PCR.DataForRobot.csv",sep="."), row.names = F)
  write.xlsx(PCR.df1, file = paste(plate,"PCR.xlsx",sep="."), sheetName= "DataForRobot", col.names= T, row.names = F, append=TRUE)
  
  # Add Pre-fab Output for Adding Standards to Plate
  #write.csv(readRDS(file = "qPCR.DataForRobot.rds"), file = paste(plate,"qPCR.standards.DataForRobot.csv",sep="."), row.names = F)
  write.xlsx(readRDS(file = "qPCR.DataForRobot.rds"), file = paste(plate,"standards.xlsx",sep="."), sheetName= "DataForRobot", col.names= T, row.names = F, append=TRUE)
}


##########################
## Append 'ComponentsSummary' Tab

for (plate in unique(x$PCR_Plate)){
  # Write Output for Samples
  #write.csv(readRDS(file = "ComponentsSummary.rds"), file = paste(plate,"PCR.ComponentsSummary.csv",sep="."), row.names = F)
  write.xlsx(readRDS(file = "ComponentsSummary.rds"), file = paste(plate,"PCR.xlsx",sep="."), sheetName= "ComponentsSummary", col.names= T, row.names = F, append=TRUE)
  
  # Add Pre-fab Output for Adding Standards to Plate
  #write.csv(readRDS(file = "ComponentsSummary.rds"), file = paste(plate,"qPCR.standards.ComponentsSummary.csv",sep="."), row.names = F)
  write.xlsx(readRDS(file = "ComponentsSummary.rds"), file = paste(plate,"standards.xlsx",sep="."), sheetName= "ComponentsSummary", col.names= T, row.names = F, append=TRUE)
}

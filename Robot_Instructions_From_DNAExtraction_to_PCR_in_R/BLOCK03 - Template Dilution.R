library(reshape2)

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

## Import Data
x <- readRDS(file="post.pico.plates1-4.rds")
x <- subset(x, SampleName != "Blank")  # Remove Blanks

##################################
## Prepare Dilution Plates for PCR

## Note: There will be a file produced for two stages of robot handling:
# i) the "*.Pooling.proto.array.DNA.for.dilution.csv" is instructions to re-array the DNA and is the SECOND STEP of the robot handling (though it is created first in this script)
# ii) the "*.Pooling.proto.WATER.for.dilution.csv" is instructions to array H2O for the dilution, and is the FIRST STEP

# Calculate Volume Needed for 2ng of Template
x$TemplateVol <- 2/x$DNA

# Designate Failed Extractions
x[which(x$DNA < 0.05),"TemplateVol"] <- NA

# Cap the Template Volume at 9.4 uL for low DNA samples (max. template for PCR reaction)
x[which(x$TemplateVol > 9.4),"TemplateVol"] <- 9.4

## Prepare Dilutions of High Concentration DNA to avoid Robot Pipetting Volumes < 5 uL
x$Diluted_DNA <- x$DNA
x$TemplateDiluteVol <- 0

########################
### CUSTOMIZABLE Section
########################

# Dilute Sample with concentrations <= 5 uL by 1:2  (i.e. 20 uL in 20 uL)
x[which(x$TemplateVol <= 4),"TemplateDiluteVol"] <- 20 
x[which(x$TemplateVol <= 4),"Diluted_DNA"] <- x[which(x$TemplateVol <= 4),"DNA"]/2

# Dilute Those <= 2 uL by 1:5 (i.e. 10 in 40 uL)
x[which(x$TemplateVol <= 1.75),"TemplateDiluteVol"] <- 40
x[which(x$TemplateVol <= 1.75),"Diluted_DNA"] <-  x[which(x$TemplateVol <= 1.75),"DNA"]/5

# Dilute Those <= 1 uL by 1:10 (i.e. 10 in 90 uL)
x[which(x$TemplateVol <= 0.8),"TemplateDiluteVol"] <- 90
x[which(x$TemplateVol <= 0.8),"Diluted_DNA"] <- x[which(x$TemplateVol <= 0.8),"DNA"]/10

# Dilute Those <= 0.5 uL by 1:20 (i.e. 5 in 95 uL)
x[which(x$TemplateVol <= 0.4),"TemplateDiluteVol"] <- 95
x[which(x$TemplateVol <= 0.4),"Diluted_DNA"] <- x[which(x$TemplateVol <= 0.4),"DNA"]/20

# Adjust Template Volume based on Dilutions
x$FinalTemplateVol <- x$TemplateVol
x$FinalTemplateVol[which(x$TemplateDiluteVol > 0)] <- 2/x$Diluted_DNA[which(x$TemplateDiluteVol > 0)]

# Add template for any sample less than 2 uL (because I don't trust the robot to pipette small volumes) 
x$FinalTemplateVol[which(x$FinalTemplateVol < 2)] <- 2

# Round volumes to one decimal place
x$FinalTemplateVol <- round(x$FinalTemplateVol, 1)

# Cap Final Template Volume at 9.4
x[which(x$FinalTemplateVol > 9.4),"FinalTemplateVol"] <- 9.4

# Calculate Water to Add
x$water_to_add <- 9.4-x$FinalTemplateVol

#####################
## END Of Customizing
#####################

###################################################################
## Prepare Instructions for Robot to Perform DNA Template Dilutions

# Order by DNA_Sample Order (i.e. 1-96 from A1->H1, B2->H2 etc.)
x <- x[order(x$DNA_Plate, x$DNA_Sample),]

# Remove counters as a precaution
rm("destination")
rm("destination_counter")
rm("destination_well")

# Set
x$DILUTION_Plate <- NA
x$DILUTION_Well <- NA

## Write Instructions to array DNA into dilution plates (Second Step for Robot after Adding water) (Plate Pooling Protocol)
## 
for (plate in unique(x$DNA_Plate)){
  # Re-set Counter
  count = 1
  
  # Subset Plate
  foo <- subset(x, DNA_Plate == plate)

  # Run Through All Samples
  for (n in 1:nrow(foo)){

    # Only Keep Samples with Template to Add (i.e. successful DNA extractions as defined in BLOCK02)
    if (is.na(foo$TemplateVol[n]) == FALSE){

      # Assign Destination Plate and Well
      if (exists("destination") == FALSE){
        destination_counter = 1
        destination_plate = 1
        destination <- make.plate()
        destination_well <- destination[destination_counter]
      } else if (destination_counter > 96) {  
        destination_counter = 1
        destination_plate = destination_plate + 1
        destination_well <- destination[destination_counter]
        count = 1
      } else {
        destination_well <- destination[destination_counter]
      }
      
      # Set Transfer Volume
      if (foo[n,"TemplateDiluteVol"] == 95){
        Transfer_volume = 5
      } else if (foo[n,"TemplateDiluteVol"] %in% c(40,90)){
        Transfer_volume = 10
      } else if (foo[n,"TemplateDiluteVol"] %in% c(20)){
        Transfer_volume = 20
      } else if (foo[n,"TemplateDiluteVol"] %in% c(0)){
        Transfer_volume = 30
      }
      
      # Write Instructions to re-array DNA with the 'plate pooling' protocol
      if (count == 1){
        transfer_csv <- data.frame(source_labware = "Plate2_HS", source_position = foo[n,"DNA_Well"], target_labware = "Plate10", target_position = destination_well, volume = Transfer_volume, sampleID = foo[n,"SampleID"], SourcePlate = plate)
      } else {
        transfer_csv <- rbind(transfer_csv, data.frame(source_labware = "Plate2_HS", source_position = foo[n,"DNA_Well"], target_labware = "Plate10", target_position = destination_well, volume = Transfer_volume, sampleID = foo[n,"SampleID"], SourcePlate = plate))
      }

      # Update Counters
      destination_counter = destination_counter + 1
      count = count + 1
      
      # Update Main Datasheet with TEMPLATE info
      source <- which(x$DNA_Plate == plate & x$DNA_Well == foo[n,"DNA_Well"])
      x$DILUTION_Plate[source] <- paste("Plate",destination_plate,sep="")
      x$DILUTION_Well[source] <- destination_well
      
      # Write Output
      write.csv(transfer_csv, file = paste(plate,"to",destination_plate,"Pooling.proto.array.DNA.for.dilution.csv",sep="."), row.names = F)
      
    } 
  }
}

# Save Information about location in PCR wells
saveRDS(x, file="diluted.template.1-4.rds")

## Write Instructions to array WATER into dilution plates (Plate Pooling Protocol)
##
x <- readRDS(file="diluted.template.1-4.rds")
x <- subset(x, FinalTemplateVol != "NA")

for (plate in unique(x$DILUTION_Plate)){
  # Set counter
  count = 1
  total_volume = 0  # Total Molecular Grade H2O Needed for Each Run
  
  # Subset to Plate
  foo <- subset(x, DILUTION_Plate == plate)
  
  # Run through all values
  for (n in 1:nrow(foo)){
    
    # Track to make sure the total volume of water for dilutions doesn't exceed 1500 uL
    if (total_volume < 1500){
      tube = 1
    } else if (total_volume > 1500 & total_volume < 3000){
      tube = 2
    } else if (total_volume > 3000 & total_volume < 4500){
      tube = 3
    } else if (total_volume > 4500 & total_volume < 6000){
      tube = 4
    } else if (total_volume > 6000 & total_volume < 7500){
      tube = 5
    }

    # Get Dilution Volume for Sample    
    Dilution_volume = foo[n,"TemplateDiluteVol"]
    
    # Write Instructions for Adding Water Dilution
    if (Dilution_volume > 0){  # Only include transfers of actual volumes of water (saves on tips)
      if (count == 1) {
        if (Dilution_volume > 50){   ## The Plate Pooling Protocol used 50uL tips, so volumes of 50 require duplicate entries
          dilution_csv <- data.frame(source_labware = "Tubes", source_position = 1, target_labware = "Plate10", target_position = foo[n,"DILUTION_Well"], volume = 50, sampleID = foo[n,"SampleID"])
          dilution_csv <- rbind(dilution_csv, data.frame(source_labware = "Tubes", source_position = 1, target_labware = "Plate10", target_position = foo[n,"DILUTION_Well"], volume = Dilution_volume-50, sampleID = foo[n,"SampleID"]))
        } else {
          dilution_csv <- data.frame(source_labware = "Tubes", source_position = 1, target_labware = "Plate10", target_position = foo[n,"DILUTION_Well"], volume = Dilution_volume, sampleID = foo[n,"SampleID"])
        }
        count = count + 1
      } else {
        if (Dilution_volume > 50){
          dilution_csv <- rbind(dilution_csv, data.frame(source_labware = "Tubes", source_position = tube, target_labware = "Plate10", target_position = foo[n,"DILUTION_Well"], volume = 50, sampleID = foo[n,"SampleID"]))
          dilution_csv <- rbind(dilution_csv, data.frame(source_labware = "Tubes", source_position = tube, target_labware = "Plate10", target_position = foo[n,"DILUTION_Well"], volume = Dilution_volume-50, sampleID = foo[n,"SampleID"]))
        } else {
          dilution_csv <- rbind(dilution_csv, data.frame(source_labware = "Tubes", source_position = tube, target_labware = "Plate10", target_position = foo[n,"DILUTION_Well"], volume = Dilution_volume, sampleID = foo[n,"SampleID"]))
        }
      }
    }

    # Update Counters
    total_volume = total_volume + Dilution_volume
  }

  # Write Output  
  write.csv(dilution_csv, file = paste(plate,"Pooling.proto.WATER.for.dilution.csv",sep=""), row.names = F)    
}
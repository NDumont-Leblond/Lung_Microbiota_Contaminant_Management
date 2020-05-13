#This script enable the removal of contaminating OTUs in microbiota dataset based of experimental controls

#As part of a validation protocol for multiple extraction techniques, the script takes in account the extraction methods used as the data set is in triplicat (one for each extraction method). This step can be removed when analysing microbiota samples in a real-life situation.
#The metadata file must contain a column for the patient ID (e.g.1,2,3,etc.) and the extraction method (e.g.:Blood, Microbial, Powersoil)

##Underlying principal
#The script first selects a patient and an extraction methods. It then extracts the related sample names. One OTU at a time (row), the script then check the OTU count for the cancerous tissue sample, the healthy tissue and the control sample for each patient and method pairs. 
#It only enters the loop if the controls a positive count for that particular OTU. If the OTU counts of both tissue is null, it escapes the loop. If OTU counts of either of the sample is positive, it checks if the relative abundance of reads in the control to that of the sample sample is lower than 0,1%. If so, the reads for that OTU are kept.
#Otherwise, the reads of both samples related to that patient for that OTU are removed. 
#The script loops for each combination of patient ID and method.

library(phyloseq)

#----Import data----
#Modify to following variables according to your data
#Set working directory (the rest of the script requires all the files to be located in that directory)
setwd("path")
#Name of mothur shared files
Shared <- "name.without.extension(.shared)"
#Name of the taxonomy information file
Taxonomy <- "name.of.tax.file.created.by.mothur.without.extension"
#The lowest level of identification (Genus = FALSE, Species = TRUE)
Species <- FALSE
#Name of the metadata (It must be a tabulated file. The first row must be column names (avoid spaces,numbers and special characters). The first column must be the matching samples names to the shared and taxonomy files)
#The metadata must include at least one column that discriminates the controls from samples (e.g.:colname = Type.sample, factors = c("Samples","Controls")) and an other one that discriminates the two experimental conditions (eg.:colname = Conditions, factors = c("Before","After")) 
Metadata <- "metadata.txt"

metadata <- read.delim(Metadata, row.names=1, sep = "\t")
OTUt <- import_mothur(mothur_shared_file = Shared ,mothur_constaxonomy_file = Taxonomy)
if (Species == FALSE){
  taxvec <- c("Domain","Phylum","Class","Order","Family","Genus")
} else if (Species == TRUE) {
  taxvec <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
}
colnames(OTUt@tax_table) <- taxvec
OTUt <- subset_taxa(OTUt, Domain == "Bacteria")
sample_data(object = OTUt) <- metadata

####Contaminating OTUs removal####
#Remove the otu from samples associated with the control (removal of otu base on a minimum proportion)

#Create otu table in dataframe format. This greatly increase the speed of the analysis, compare to referencing the logical tests directly in the phyloseq object. It also allows to keep the otu table as empirical counts
OTUt_table <- OTUt
#Normalize by sequencng depth
OTUt_tab <- transform_sample_counts(OTUt_table, function(x) x / sum(x) )
OTUt_tab <- as.data.frame(OTUt_tab@otu_table)

#Remove contaminants
#Patient ID
patient <- c(1,2,3,4,5)
#Extraction method ID
kit <- c("Powersoil","Microbial","Blood") 
nbotu <- rownames(OTUt@otu_table)

for (a in patient){
  for (b in kit){
    #Get names of samples (same patient & same kit)
    ctrl <- rownames(OTUt@sam_data[OTUt@sam_data$Patient==a & OTUt@sam_data$Extraction.kit==b & OTUt@sam_data$Type == "Reg" & OTUt@sam_data$Type.tissue=="Control",])
    echc <- rownames(OTUt@sam_data[OTUt@sam_data$Patient==a & OTUt@sam_data$Extraction.kit==b & OTUt@sam_data$Type == "Reg" & OTUt@sam_data$Type.tissue=="Cancerous",])
    echh <- rownames(OTUt@sam_data[OTUt@sam_data$Patient==a & OTUt@sam_data$Extraction.kit==b & OTUt@sam_data$Type == "Reg" & OTUt@sam_data$Type.tissue=="Healthy",])
    #If the counts of the control is higher than 0, the count of samples associated with that control are removed (==0)
    for (c in nbotu){
      if (OTUt_tab[c,ctrl] > 0){ #enter loop only if the controls has reads
        
        if (OTUt_tab[c,echc] == 0 & OTUt_tab[c,echh] == 0){ 
            next #If both samples have no reads of that otu, the loop stops
          } else if (OTUt_tab[c,echc] == 0 & OTUt_tab[c,echh] > 0){ #If only the healthy sample has reads for that otu
            if (!(OTUt_tab[c,ctrl]/OTUt_tab[c,echh] < 0.001)){   #If the proportion of reads ctrl to sample is lower than 0,1% the reads are kept
            OTUt@otu_table[c,echh] <- 0 #Otherwise the reads are removed for the sample
            }
          } else if (OTUt_tab[c,echc] > 0 & OTUt_tab[c,echh] == 0){ #If only the cancerous sample has reads for that otu
            if (!(OTUt_tab[c,ctrl]/OTUt_tab[c,echc] < 0.001)){ #If the proportion of reads ctrl to sample is lower than 0,1% the reads are kept
            OTUt@otu_table[c,echc] <- 0 #Otherwise the reads are removed for the sample
            }
          } else if (OTUt_tab[c,echc] > 0 & OTUt_tab[c,echh] > 0){ #If only the healthy sample has reads for that otu
            if (!(OTUt_tab[c,ctrl]/OTUt_tab[c,echc] < 0.001 | OTUt_tab[c,ctrl]/OTUt_tab[c,echh] < 0.001)){  #If the proportion of reads ctrl to sample is lower than 0,1% the reads are kept
            OTUt@otu_table[c,echh] <- 0 #Otherwise the reads are removed for the sample
            OTUt@otu_table[c,echc] <- 0 #Otherwise the reads are removed for the sample
            
            }
          }
        }
      }
    }
  }
}

#The script removes the OTUs identified as contaminants directly in the OTUt object

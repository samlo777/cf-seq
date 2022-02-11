#Notes: This code does not need to be executed in order for the app to run - it is simply for reference. The code here was initially used to
#gather KEGG pathway data for species present in the app, but the data objects were subsequently saved as external files, are present
#in the Github directory, and are loaded in from that directory in the 'Data Setup.R' script. Each call to the KEGG database, which retrieves
#KEGG pathway information for all genes of a single bacterial/fungal strain, can take up to several hours to run.

#1. Install necessary packages
{
  #BiocManager::install("edgeR")
  library(edgeR)
  
  #BiocManager::install("GEOquery")
  library(GEOquery)
  
  #install.packages("stringr")
  library(stringr)
  
  #BiocManager::install("KEGGREST")
  library(KEGGREST)
  
  #install.packages("stringr")
  library(stringr)
}

#2. KEGG_Data: Contains KEGG pathway info for compatible bacterial/fungal strains.
{
  
  #Overview: See all species with KEGG pathway annotations in the KEGG database
  {
    org <- keggList("organism")
    org <- as.data.frame(org)
    View(org)
  }
  
  #Aspergillus Fumigatus 
  {
    
    #1. Load in KEGG database information for selected strain
    List <- keggList("afm") 
    
    #2. Gather all specific pathway data of interest into a list object
    Pathway_Data <- vector(mode = "list") 
    for (i in 1:length(List)) {
      for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
        if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
          Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
          names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
          print(Pathway_Data)[length(Pathway_Data)]
        }
      }
    }
    
    #3. Convert large list into more concise data frame
    DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
    colnames(DF) <- c("Gene", "Pathway")
    DF$Gene <- names(Pathway_Data)
    for (i in 1:length(Pathway_Data)) {
      DF$Pathway[i] <- Pathway_Data[[i]]
    }
    
    #4. Export data frame and save outside of R
    write.csv(DF, "PathwayData/AFPaths.csv")
    
    #5. Read in data frame for use in the application
    AFPaths <- read.csv("PathwayData/AFPaths.csv")
    AFPaths <- AFPaths[,-1]
    AFPaths$Gene <- sub("afm:", "", AFPaths$Gene)
    
  }
  
  #Bacteroides Species
  {
    
    #B. thetaiotaomicron VPI-5482
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("bth") 
      length(List)
      
      #2. Gather all specific pathway data of interest into a list object
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/BTVPI-5482Paths.csv")
      
      #5. Read in data frame for use in the application
      BTVPI5482Paths <- read.csv("PathwayData/BTVPI-5482Paths.csv")
      BTVPI5482Paths <- BTVPI5482Paths[,-1]
      BTVPI5482Paths$Gene <- sub("bth:", "", BTVPI5482Paths$Gene)
    }
    
    #B. xylanixolvens
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("bxy") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/BXYPaths.csv")
      
      #5. Read in data frame for use in the application
      BXYPaths <- read.csv("PathwayData/BXYPaths.csv")
      BXYPaths <- BXYPaths[,-1]
      BXYPaths$Gene <- sub("bxy:", "", BXYPaths$Gene)
    }
    
  }
  
  #Burkholderia Species 
  {
    
    #B. cepacia AU1054 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("bcn") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/BCNPaths.csv")
      
      #5. Read in data frame for use in the application
      BCNPaths <- read.csv("PathwayData/BCNPaths.csv")
      BCNPaths <- BCNPaths[,-1]
      BCNPaths$Gene <- sub("bcn:", "", BCNPaths$Gene)
    }
    
    #B. Cepacia HI2424
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("bch") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/BCHPaths.csv")
      
      #5. Read in data frame for use in the application
      BCHPaths <- read.csv("PathwayData/BCHPaths.csv")
      BCHPaths <- BCHPaths[,-1]
      BCHPaths$Gene <- sub("bch:", "", BCHPaths$Gene)
    }
    
    #B. Pseudomallei K96243
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("bps") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/BPSPaths.csv")
      
      #5. Read in data frame for use in the application
      BPSPaths <- read.csv("PathwayData/BPSPaths.csv")
      BPSPaths <- BPSPaths[,-1]
      BPSPaths$Gene <- sub("bps:", "", BPSPaths$Gene)
    }
    
    #B. Thailandensis E264
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("bte") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/BTEPaths.csv")
      
      #5. Read in data frame for use in the application
      BTEPaths <- read.csv("PathwayData/BTEPaths.csv")
      BTEPaths <- BTEPaths[,-1]
      BTEPaths$Gene <- sub("bte:", "", BTEPaths$Gene)
    }
    
  }
  
  #Candida Albicans 
  {
    #1. Load in KEGG database information for selected strain
    List <- keggList("cal") 
    
    #2. Gather all specific pathway data of interest into a list object
    Pathway_Data <- vector(mode = "list") 
    for (i in 1:length(List)) {
      for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
        if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
          Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
          names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
          print(Pathway_Data)[length(Pathway_Data)]
        }
      }
    }
    
    #3. Convert large list into more concise data frame
    DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
    colnames(DF) <- c("Gene", "Pathway")
    DF$Gene <- names(Pathway_Data)
    for (i in 1:length(Pathway_Data)) {
      DF$Pathway[i] <- Pathway_Data[[i]]
    }
    
    #4. Export data frame and save outside of R
    write.csv(DF, "PathwayData/CALPaths.csv")
    
    #5. Read in data frame for use in the application
    CALPaths <- read.csv("PathwayData/CALPaths.csv")
    CALPaths <- CALPaths[,-1]
    CALPaths$Gene <- sub("cal:", "", CALPaths$Gene)
  }
  
  #Clostridioides difficile [NA]
  
  #Fusobacterium nucleatum [NA]
  
  #Haemophilus influenzae 
  {
    
    #H. influenza Rd KW20
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("hin") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/HINPaths.csv")
      
      #5. Read in data frame for use in the application
      HINPaths <- read.csv("PathwayData/HINPaths.csv")
      HINPaths <- HINPaths[,-1]
      HINPaths$Gene <- sub("hin:", "", HINPaths$Gene)
    }
    
    #H. influenza 86-028NP 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("hit") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/HITPaths.csv")
      
      #5. Read in data frame for use in the application
      HITPaths <- read.csv("PathwayData/HITPaths.csv")
      HITPaths <- HITPaths[,-1]
      HITPaths$Gene <- sub("hit:", "", HITPaths$Gene)
    }
    
    #H. influenza R2866
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("hiz") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/HIZPaths.csv")
      
      #5. Read in data frame for use in the application
      HIZPaths <- read.csv("PathwayData/HIZPaths.csv")
      HIZPaths <- HIZPaths[,-1]
      HIZPaths$Gene <- sub("hiz:", "", HIZPaths$Gene)
    }
    
    #H. influenza 723
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("hix") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/HIXPaths.csv")
      
      #5. Read in data frame for use in the application
      HIXPaths <- read.csv("PathwayData/HIXPaths.csv")
      HIXPaths <- HIXPaths[,-1]
      HIXPaths$Gene <- sub("hix:", "", HIXPaths$Gene)
    }
    
  }
  
  #Mycobacterium abscessus [NA]
  
  #Pseudomonas aeruginosa
  {
    #P. aeruginosa PA14 
    {
      List <- keggList("pau") #pau = UCBPP-PA14

      Pathway_Data <- vector(mode = "list") 
      for (i in 1: length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(PA14List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(PA14List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #Convert large vector into single data frame
      {
        DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
        colnames(DF) <- c("Gene", "Pathway")
        
        DF$Gene <- names(Pathway_Data)
        
        for (i in 1:length(Pathway_Data)) {
          DF$Pathway[i] <- Pathway_Data[[i]]
        }
      }
      
      #Export data to save outside of R (a hedge against accidental deletion) and re-upload as needed
      {
        write.csv(DF, "PathwayData/PA14Paths.csv") #In case the data frame created above gets deleted
        PA14Paths <- read.csv("PathwayData/PA14Paths.csv")
        PA14Paths <- PA14Paths[,-1]
        PA14Paths$Gene <- sub("pau:", "", PA14Paths$Gene)
      }
    }
    
    #P. aeruginosa PAO1
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("pae") 
      lng
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      } #Remove
      
      for (i in 1: length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/PAEPaths.csv")
      
      #5. Read in data frame for use in the application
      PAEPaths <- read.csv("PathwayData/PAEPaths.csv")
      PAEPaths <- PAEPaths[,-1]
      PAEPaths$Gene <- sub("pae:", "", PAEPaths$Gene)
    }
    
    #P. aeruginosa MTB-1
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("paem") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/PAEMPaths.csv")
      
      #5. Read in data frame for use in the application
      PAEMPaths <- read.csv("PathwayData/PAEMPaths.csv")
      PAEMPaths <- PAEMPaths[,-1]
      PAEMPaths$Gene <- sub("paem:", "", PAEMPaths$Gene)
    }
    
    #P. aeruginosa B136-33
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("psg") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/PSGPaths.csv")
      
      #5. Read in data frame for use in the application
      PSGPaths <- read.csv("PathwayData/PSGPaths.csv")
      PSGPaths <- PSGPaths[,-1]
      PSGPaths$Gene <- sub("psg:", "", PSGPaths$Gene)
    }
    
  }
  
  #Porphyromonas Species 
  {
    
    #P. gingivalis W83 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("pgi") 
      length(List)
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/PGIPaths.csv")
      
      #5. Read in data frame for use in the application
      PGIPaths <- read.csv("PathwayData/PGIPaths.csv")
      PGIPaths <- PGIPaths[,-1]
      PGIPaths$Gene <- sub("pgi:", "", PGIPaths$Gene)
    }
    
    #P. gingivalis ATCC 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("pgn") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/PGNPaths.csv")
      
      #5. Read in data frame for use in the application
      PGNPaths <- read.csv("PathwayData/PGNPaths.csv")
      PGNPaths <- PGNPaths[,-1]
      PGNPaths$Gene <- sub("pgn:", "", PGNPaths$Gene)
    }
    
  }
  
  #Staphylococcus aureus 
  {
    
    #S. aureus JDK6008
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("suk") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SUKPaths.csv")
      
      #5. Read in data frame for use in the application
      SUKPaths <- read.csv("PathwayData/SUKPaths.csv")
      SUKPaths <- SUKPaths[,-1]
      SUKPaths$Gene <- sub("suk:", "", SUKPaths$Gene)
    }
  
    #S. aureus RF122 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("sab") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SABPaths.csv")
      
      #5. Read in data frame for use in the application
      SABPaths <- read.csv("PathwayData/SABPaths.csv")
      SABPaths <- SABPaths[,-1]
      SABPaths$Gene <- sub("sab:", "", SABPaths$Gene)
    }
    
    #S. aureus NCTC8325 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("sao") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SAOPaths.csv")
      
      #5. Read in data frame for use in the application
      SAOPaths <- read.csv("PathwayData/SAOPaths.csv")
      SAOPaths <- SAOPaths[,-1]
      SAOPaths$Gene <- sub("sao:", "", SAOPaths$Gene)
    }
    
    #S. aureus Newman
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("sae")
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SAEPaths.csv")
      
      #5. Read in data frame for use in the application
      SAEPaths <- read.csv("PathwayData/SAEPaths.csv")
      SAEPaths <- SAEPaths[,-1]
      SAEPaths$Gene <- sub("sae:", "", SAEPaths$Gene)
    }
    
    #S. aureus USA300
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("sax") 
      length(List)
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 

      for (i in 1: length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SAXPaths.csv")
      
      #5. Read in data frame for use in the application
      SAXPaths <- read.csv("PathwayData/SAXPaths.csv")
      SAXPaths <- SAXPaths[,-1]
      SAXPaths$Gene <- sub("sax:", "", SAXPaths$Gene)
      
    }
    
  }
  
  #Stenotrophomonas maltophilia 
  {
    
    #S. maltophilia K279a 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("sml") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SMLPaths.csv")
      
      #5. Read in data frame for use in the application
      SMLPaths <- read.csv("PathwayData/SMLPaths.csv")
      SMLPaths <- SMLPaths[,-1]
      SMLPaths$Gene <- sub("sml:", "", SMLPaths$Gene)
    }
    
    #S. maltophilia D457 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("smz") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SMZPaths.csv")
      
      #5. Read in data frame for use in the application
      SMZPaths <- read.csv("PathwayData/SMZPaths.csv")
      SMZPaths <- SMZPaths[,-1]
      SMZPaths$Gene <- sub("smz:", "", SMZPaths$Gene)
    }
    
  }
  
  #Streptococcus Species 
  {
    
    #S. pyrogenes M1GAS
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("spy") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SPYPaths.csv")
      
      #5. Read in data frame for use in the application
      SPYPaths <- read.csv("PathwayData/SPYPaths.csv")
      SPYPaths <- SPYPaths[,-1]
      SPYPaths$Gene <- sub("spy:", "", SPYPaths$Gene)
    }
    
    #S. salivarius HSISS4 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("ssah") 
      length(List)
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SSAHPaths.csv")
      
      #5. Read in data frame for use in the application
      SSAHPaths <- read.csv("PathwayData/SSAHPaths.csv")
      SSAHPaths <- SSAHPaths[,-1]
      SSAHPaths$Gene <- sub("ssah:", "", SSAHPaths$Gene)
    }
    
    #S. mutans UA159 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("smu") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SMUPaths.csv")
      
      #5. Read in data frame for use in the application
      SMUPaths <- read.csv("PathwayData/SMUPaths.csv")
      SMUPaths <- SMUPaths[,-1]
      SMUPaths$Gene <- sub("smu:", "", SMUPaths$Gene)
    }
    
    #S. pneumoniae D39 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("spd") 
      
      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SPDPaths.csv")
      
      #5. Read in data frame for use in the application
      SPDPaths <- read.csv("PathwayData/SPDPaths.csv")
      SPDPaths <- SPDPaths[,-1]
      SPDPaths$Gene <- sub("spd:", "", SPDPaths$Gene)
    }
    
    #S. agalactiae A909 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("sak") 

      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SAKPaths.csv")
      
      #5. Read in data frame for use in the application
      SAKPaths <- read.csv("PathwayData/SAKPaths.csv")
      SAKPaths <- SAKPaths[,-1]
      SAKPaths$Gene <- sub("sak:", "", SAKPaths$Gene)
    }
    
    #S. gallollyticus UCN34 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("sga") 

      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SGAPaths.csv")
      
      #5. Read in data frame for use in the application
      SGAPaths <- read.csv("PathwayData/SGAPaths.csv")
      SGAPaths <- SGAPaths[,-1]
      SGAPaths$Gene <- sub("sga:", "", SGAPaths$Gene)
    }
    
    #S. suis 05ZYH33 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("ssu") 

      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SSUPaths.csv")
      
      #5. Read in data frame for use in the application
      SSUPaths <- read.csv("PathwayData/SSUPaths.csv")
      SSUPaths <- SSUPaths[,-1]
      SSUPaths$Gene <- sub("ssu:", "", SSUPaths$Gene)
    }
    
    #S. gordonii 
    {
      #1. Load in KEGG database information for selected strain
      List <- keggList("sgo") 

      #2. Gather all specific pathway data of interest into a list object
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            names(Pathway_Data)[length(Pathway_Data)] <- names(List)[i]
            print(Pathway_Data)[length(Pathway_Data)]
          }
        }
      }
      
      #3. Convert large list into more concise data frame
      DF <- data.frame(matrix(NA, nrow = length(Pathway_Data), ncol = 2))
      colnames(DF) <- c("Gene", "Pathway")
      DF$Gene <- names(Pathway_Data)
      for (i in 1:length(Pathway_Data)) {
        DF$Pathway[i] <- Pathway_Data[[i]]
      }
      
      #4. Export data frame and save outside of R
      write.csv(DF, "PathwayData/SGOPaths.csv")
      
      #5. Read in data frame for use in the application
      SGOPaths <- read.csv("PathwayData/SGOPaths.csv")
      SGOPaths <- SGOPaths[,-1]
      SGOPaths$Gene <- sub("sgo:", "", SGOPaths$Gene)
    }
    
  }
  
  
}
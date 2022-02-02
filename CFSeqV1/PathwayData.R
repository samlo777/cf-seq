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
  
  #Aspergillus Fumigatus [Done]
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
    write.csv(DF, "AFPaths.csv")
    
    #5. Read in data frame for use in the application
    AFPaths <- read.csv("AFPaths.csv")
    AFPaths <- AFPaths[,-1]
    AFPaths$Gene <- sub("afm:", "", AFPaths$Gene)
    
  }
  
  #Bacteroides Species [Done]
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
      write.csv(DF, "BTVPI-5482Paths.csv")
      
      #5. Read in data frame for use in the application
      BTVPI5482Paths <- read.csv("BTVPI-5482Paths.csv")
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
      write.csv(DF, "BXYPaths.csv")
      
      #5. Read in data frame for use in the application
      BXYPaths <- read.csv("BXYPaths.csv")
      BXYPaths <- BXYPaths[,-1]
      BXYPaths$Gene <- sub("bxy:", "", BXYPaths$Gene)
    }
    
  }
  
  #Burkholderia Species [Done]
  {
    
    #B. cepacia AU1054 [Done]
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
      write.csv(DF, "BCNPaths.csv")
      
      #5. Read in data frame for use in the application
      BCNPaths <- read.csv("BCNPaths.csv")
      BCNPaths <- BCNPaths[,-1]
      BCNPaths$Gene <- sub("bcn:", "", BCNPaths$Gene)
    }
    
    #B. Cepacia HI2424 [Done]
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
      write.csv(DF, "BCHPaths.csv")
      
      #5. Read in data frame for use in the application
      BCHPaths <- read.csv("BCHPaths.csv")
      BCHPaths <- BCHPaths[,-1]
      BCHPaths$Gene <- sub("bch:", "", BCHPaths$Gene)
    }
    
    #B. Pseudomallei K96243 [Done]
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
      write.csv(DF, "BPSPaths.csv")
      
      #5. Read in data frame for use in the application
      BPSPaths <- read.csv("BPSPaths.csv")
      BPSPaths <- BPSPaths[,-1]
      BPSPaths$Gene <- sub("bps:", "", BPSPaths$Gene)
    }
    
    #B. Thailandensis E264 [Done]
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
      write.csv(DF, "BTEPaths.csv")
      
      #5. Read in data frame for use in the application
      BTEPaths <- read.csv("BTEPaths.csv")
      BTEPaths <- BTEPaths[,-1]
      BTEPaths$Gene <- sub("bte:", "", BTEPaths$Gene)
    }
    
  }
  
  #Candida Albicans [Done]
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
    write.csv(DF, "CALPaths.csv")
    
    #5. Read in data frame for use in the application
    CALPaths <- read.csv("CALPaths.csv")
    CALPaths <- CALPaths[,-1]
    CALPaths$Gene <- sub("cal:", "", CALPaths$Gene)
  }
  
  #Clostridioides difficile [NA]
  
  #Fusobacterium nucleatum [NA]
  
  #Haemophilus influenzae [Done]
  {
    
    #H. influenza Rd KW20 [Done]
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
      write.csv(DF, "HINPaths.csv")
      
      #5. Read in data frame for use in the application
      HINPaths <- read.csv("HINPaths.csv")
      HINPaths <- HINPaths[,-1]
      HINPaths$Gene <- sub("hin:", "", HINPaths$Gene)
    }
    
    #H. influenza 86-028NP [Done]
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
      write.csv(DF, "HITPaths.csv")
      
      #5. Read in data frame for use in the application
      HITPaths <- read.csv("HITPaths.csv")
      HITPaths <- HITPaths[,-1]
      HITPaths$Gene <- sub("hit:", "", HITPaths$Gene)
    }
    
    #H. influenza R2866 [Done]
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
      write.csv(DF, "HIZPaths.csv")
      
      #5. Read in data frame for use in the application
      HIZPaths <- read.csv("HIZPaths.csv")
      HIZPaths <- HIZPaths[,-1]
      HIZPaths$Gene <- sub("hiz:", "", HIZPaths$Gene)
    }
    
    #H. influenza 723 [Done]
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
      write.csv(DF, "HIXPaths.csv")
      
      #5. Read in data frame for use in the application
      HIXPaths <- read.csv("HIXPaths.csv")
      HIXPaths <- HIXPaths[,-1]
      HIXPaths$Gene <- sub("hix:", "", HIXPaths$Gene)
    }
    
  }
  
  #Mycobacterium abscessus [NA]
  
  #Pseudomonas aeruginosa [Done]
  {
    #P. aeruginosa PA14 [Done]
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
        write.csv(DF, "PA14Paths.csv") #In case the data frame created above gets deleted
        PA14Paths <- read.csv("PA14Paths.csv")
        PA14Paths <- PA14Paths[,-1]
        PA14Paths$Gene <- sub("pau:", "", PA14Paths$Gene)
      }
    }
    
    #P. aeruginosa PAO1 [Done]
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
      write.csv(DF, "PAEPaths.csv")
      
      #5. Read in data frame for use in the application
      PAEPaths <- read.csv("PAEPaths.csv")
      PAEPaths <- PAEPaths[,-1]
      PAEPaths$Gene <- sub("pae:", "", PAEPaths$Gene)
    }
    
    #P. aeruginosa MTB-1 [Done]
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
      write.csv(DF, "PAEMPaths.csv")
      
      #5. Read in data frame for use in the application
      PAEMPaths <- read.csv("PAEMPaths.csv")
      PAEMPaths <- PAEMPaths[,-1]
      PAEMPaths$Gene <- sub("paem:", "", PAEMPaths$Gene)
    }
    
    #P. aeruginosa B136-33 [Done]
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
      write.csv(DF, "PSGPaths.csv")
      
      #5. Read in data frame for use in the application
      PSGPaths <- read.csv("PSGPaths.csv")
      PSGPaths <- PSGPaths[,-1]
      PSGPaths$Gene <- sub("psg:", "", PSGPaths$Gene)
    }
    
  }
  
  #Porphyromonas Species [Done]
  {
    
    #P. gingivalis W83 [Done]
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
      write.csv(DF, "PGIPaths.csv")
      
      #5. Read in data frame for use in the application
      PGIPaths <- read.csv("PGIPaths.csv")
      PGIPaths <- PGIPaths[,-1]
      PGIPaths$Gene <- sub("pgi:", "", PGIPaths$Gene)
    }
    
    #P. gingivalis ATCC [Done]
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
      write.csv(DF, "PGNPaths.csv")
      
      #5. Read in data frame for use in the application
      PGNPaths <- read.csv("PGNPaths.csv")
      PGNPaths <- PGNPaths[,-1]
      PGNPaths$Gene <- sub("pgn:", "", PGNPaths$Gene)
    }
    
  }
  
  #Staphylococcus aureus [Done]
  {
    
    #S. aureus JDK6008 [Done]
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
      write.csv(DF, "SUKPaths.csv")
      
      #5. Read in data frame for use in the application
      SUKPaths <- read.csv("SUKPaths.csv")
      SUKPaths <- SUKPaths[,-1]
      SUKPaths$Gene <- sub("suk:", "", SUKPaths$Gene)
    }
    
    #S. aureus RF122 [Done]
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
      write.csv(DF, "SABPaths.csv")
      
      #5. Read in data frame for use in the application
      SABPaths <- read.csv("SABPaths.csv")
      SABPaths <- SABPaths[,-1]
      SABPaths$Gene <- sub("sab:", "", SABPaths$Gene)
    }
    
    #S. aureus NCTC8325 [Done]
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
      write.csv(DF, "SAOPaths.csv")
      
      #5. Read in data frame for use in the application
      SAOPaths <- read.csv("SAOPaths.csv")
      SAOPaths <- SAOPaths[,-1]
      SAOPaths$Gene <- sub("sao:", "", SAOPaths$Gene)
    }
    
    #S. aureus Newman [Done]
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
      write.csv(DF, "SAEPaths.csv")
      
      #5. Read in data frame for use in the application
      SAEPaths <- read.csv("SAEPaths.csv")
      SAEPaths <- SAEPaths[,-1]
      SAEPaths$Gene <- sub("sae:", "", SAEPaths$Gene)
    }
    
    #S. aureus USA300 [Done]
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
      write.csv(DF, "SAXPaths.csv")
      
      #5. Read in data frame for use in the application
      SAXPaths <- read.csv("PathwayData/SAXPaths.csv")
      SAXPaths <- SAXPaths[,-1]
      SAXPaths$Gene <- sub("sax:", "", SAXPaths$Gene)
      
    }
    
  }
  
  #Stenotrophomonas maltophilia [Done]
  {
    
    #S. maltophilia K279a [Done]
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
      write.csv(DF, "SMLPaths.csv")
      
      #5. Read in data frame for use in the application
      SMLPaths <- read.csv("SMLPaths.csv")
      SMLPaths <- SMLPaths[,-1]
      SMLPaths$Gene <- sub("sml:", "", SMLPaths$Gene)
    }
    
    #S. maltophilia D457 [Done]
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
      write.csv(DF, "SMZPaths.csv")
      
      #5. Read in data frame for use in the application
      SMZPaths <- read.csv("SMZPaths.csv")
      SMZPaths <- SMZPaths[,-1]
      SMZPaths$Gene <- sub("smz:", "", SMZPaths$Gene)
    }
    
  }
  
  #Streptococcus Species [Done]
  {
    
    #S. pyrogenes M1GAS [Done]
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
      write.csv(DF, "SPYPaths.csv")
      
      #5. Read in data frame for use in the application
      SPYPaths <- read.csv("SPYPaths.csv")
      SPYPaths <- SPYPaths[,-1]
      SPYPaths$Gene <- sub("spy:", "", SPYPaths$Gene)
    }
    
    #S. salivarius HSISS4 [Done]
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
      write.csv(DF, "SSAHPaths.csv")
      
      #5. Read in data frame for use in the application
      SSAHPaths <- read.csv("SSAHPaths.csv")
      SSAHPaths <- SSAHPaths[,-1]
      SSAHPaths$Gene <- sub("ssah:", "", SSAHPaths$Gene)
    }
    
    #S. mutans UA159 [Done]
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
      write.csv(DF, "SMUPaths.csv")
      
      #5. Read in data frame for use in the application
      SMUPaths <- read.csv("SMUPaths.csv")
      SMUPaths <- SMUPaths[,-1]
      SMUPaths$Gene <- sub("smu:", "", SMUPaths$Gene)
    }
    
    #S. pneumoniae D39 [Done]
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
      write.csv(DF, "SPDPaths.csv")
      
      #5. Read in data frame for use in the application
      SPDPaths <- read.csv("SPDPaths.csv")
      SPDPaths <- SPDPaths[,-1]
      SPDPaths$Gene <- sub("spd:", "", SPDPaths$Gene)
    }
    
    #S. agalactiae A909 [Done]
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
      write.csv(DF, "SAKPaths.csv")
      
      #5. Read in data frame for use in the application
      SAKPaths <- read.csv("SAKPaths.csv")
      SAKPaths <- SAKPaths[,-1]
      SAKPaths$Gene <- sub("sak:", "", SAKPaths$Gene)
    }
    
    #S. gallollyticus UCN34 [Done]
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
      write.csv(DF, "SGAPaths.csv")
      
      #5. Read in data frame for use in the application
      SGAPaths <- read.csv("SGAPaths.csv")
      SGAPaths <- SGAPaths[,-1]
      SGAPaths$Gene <- sub("sga:", "", SGAPaths$Gene)
    }
    
    #S. suis 05ZYH33 [Done]
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
      write.csv(DF, "SSUPaths.csv")
      
      #5. Read in data frame for use in the application
      SSUPaths <- read.csv("SSUPaths.csv")
      SSUPaths <- SSUPaths[,-1]
      SSUPaths$Gene <- sub("ssu:", "", SSUPaths$Gene)
    }
    
    #S. gordonii [Done]
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
      write.csv(DF, "SGOPaths.csv")
      
      #5. Read in data frame for use in the application
      SGOPaths <- read.csv("SGOPaths.csv")
      SGOPaths <- SGOPaths[,-1]
      SGOPaths$Gene <- sub("sgo:", "", SGOPaths$Gene)
    }
    
  }
  
  
}
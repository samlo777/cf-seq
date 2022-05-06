#Notes: Here, the focus is on gathering annotation data for the strains where appropriate gene identifiers
#are available. The strain annotations have been fashioned into a table for each amenable strain
#The strain tables have different types of columns depending on what annotations are available

#Don't try to run this file all at once - it will literally take weeks of computational time, interfacing with multiple databases, 
#to gather all of the annotations. This code is best used as a reference. To access the tables that it produces, just grab them
#out of the GitHub directory from which you downloaded this code (in the Pathway_Data sub-directory)

#1. Install Necessary Packages
{
  #BiocManager::install("KEGGREST")
  library(KEGGREST)
  
  #install.packages("stringr")
  library(stringr)
  
  #install.packages("data.table")
  library(data.table)
  
  #BiocManager::install("UniProt.ws")
  library(UniProt.ws)
}

#2a. Overview: View all species with KEGG pathway annotations in the KEGG database
{
  org <- keggList("organism")
  org <- as.data.frame(org)
  View(org)
}

#2b. Gather KEGG data for accessible strains with locus IDs matching studies in app (16) 
{
  
  #Strain: B. thetaiotaomicron, VPI-5482 [Done]
  {
    List <- keggList("bth") 
    length(List) #4911
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "BTHPaths.csv") #In case the data frame created above gets deleted
    }
    
  }
  
  #Strain: B. cenocepacia, AU1054 [Done]
  {
    List <- keggList("bcn") 
    length(List) #6618
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "BCNPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: B. pseudomallei, K96243 [Done]
  {
    List <- keggList("bps") 
    length(List) #5935
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "BPSPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: C. albicans (not strain specific) [Done]
  {
    List <- keggList("cal") 
    length(List) #6162
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "CALPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: H. influenza, Strain 723 [Done]
  {
    List <- keggList("hix") 
    length(List) #1868
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "HIXPaths.csv") #In case the data frame created above gets deleted
    }
    
  }
  
  #Strain: P. aeruginosa, PAO1 [Done]
  {
    List <- keggList("pae") 
    length(List) #5697
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "PAEPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: P. aeruginosa, PA14 [Done]
  {
    List <- keggList("pau") 
    length(List) #5977
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "PAUPaths.csv") #In case the data frame created above gets deleted
    }
    
  }
  
  #Strain: P. gingivalis, W83 [Done]
  {
    List <- keggList("pgi") 
    length(List) #2015
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "PGIPaths.csv") #In case the data frame created above gets deleted
    }
    
  }
  
  #Strain: P. gingivalis, ATCC [Done]
  {
    List <- keggList("pgn") 
    length(List) #2155
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "PGNPaths.csv") #In case the data frame created above gets deleted
    }
    
  }
  
  #Strain: S. aureus, Newman [Done]
  {
    List <- keggList("sae") 
    length(List) #2697
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SAEPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: S. aureus, USA300 [Done]
  {
    List <- keggList("sax") 
    length(List) #2841
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SAXPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: S. aureus, NCTC8325 [Done]
  {
    List <- keggList("sao") 
    length(List) #2872
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SAOPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: S. aureus, JDK6008 [Done]
  {
    List <- keggList("suk") 
    length(List) #2836
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SUKPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: S. maltophilia, D457 
  {
    List <- keggList("smz") 
    length(List) #4239
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SMZPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: S. pneumonia, D39 [Done]
  {
    List <- keggList("spd") 
    length(List) #2068
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SPDPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: S. salivarius, HSISS4 [Done]
  {
    List <- keggList("ssah") 
    length(List) #1989
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SSAHPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
}

#2c. Gather KEGG data for all other strains in the app (lacking accessible locus IDs) (12) 
{
  #Strain: A. fumigatus (not strain-specific) [Done]
  {
    List <- keggList("afm") 
    length(List) #9630
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "AFMPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: B. Xylanixolvens (not strain-specific) [Done]
  {
    List <- keggList("bxy") 
    length(List) #4466
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "BXYPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: B. Thaliandensis E264 [Done]
  {
    List <- keggList("bte") 
    length(List) #5714
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "BTEPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: B. Cenocepacia HI2424 [Done]
  {
    List <- keggList("bch") 
    length(List) #7031
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "BCHPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain H. Influenza: R2866 [Done]
  {
    List <- keggList("hiz") 
    length(List) #1892
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "HIZPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain H. Influenza: 86-028NP [Done]
  {
    List <- keggList("hit") 
    length(List) #1899
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "HITPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain H. Influenza: Rd KW20 [Done]
  {
    List <- keggList("hin") 
    length(List) #1775
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "HINPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: P. Aeruginosa MTB-1 [Done]
  {
    List <- keggList("paem") 
    length(List) #6186
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "PAEMPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #Strain: P. Aeruginosa B136-33 [Done]
  {
    List <- keggList("psg") 
    length(List) #5904
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "PSGPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #S. Aureus: RF122 [Done]
  {
    List <- keggList("sab") 
    length(List) #2662
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SABPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #S. Maltophilia, K279a [Done]
  {
    List <- keggList("sml") 
    length(List) #4440
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SMLPaths.csv") #In case the data frame created above gets deleted
    }
  }
  
  #S. Pyogenes: M1 GAS [Done]
  {
    List <- keggList("spy") 
    length(List) #1801
    
    #Gather pathway data
    {
      
      Pathway_Data <- vector(mode = "list") 
      for (i in 1:length(List)) {
        for (j in 1:length(keggGet(names(List)[i])[[1]]$PATHWAY)) {
          if (is.null(keggGet(names(List)[i])[[1]]$PATHWAY) != TRUE) {
            Pathway_Data <- c(Pathway_Data, keggGet(names(List)[i])[[1]]$PATHWAY[[j]])
            if (is.null(keggGet(names(List)[i])[[1]]$SYMBOL) != TRUE) {
              GeneSymbol <- keggGet(names(List)[i])[[1]]$SYMBOL
            } else {
              GeneSymbol <- ""
            }
            if (is.null(keggGet(names(List)[i])[[1]]$NAME) != TRUE) {
              GeneName <- keggGet(names(List)[i])[[1]]$NAME
            } else {
              GeneName <- ""
            }
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(List)[i], sep = " ", GeneSymbol)
            names(Pathway_Data)[length(Pathway_Data)] <- paste0(names(Pathway_Data)[length(Pathway_Data)], sep = " ", GeneName)
            print(Pathway_Data)[length(Pathway_Data)]
          }
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
      write.csv(DF, "SPYPaths.csv") #In case the data frame created above gets deleted
    }
  }
}

#3. Add in additional functional annotations (COG) 
{
  
  COGAnnotations <- read.csv("Pathway_Data/COG_Annotations/cog-20.cog.csv") 
  COGStrains <- read.csv("Pathway_Data/COG_Annotations/cog-20.org.csv") 
  COGDescriptions <- read.table("Pathway_Data/COG_Annotations/cog-20.def.tab.txt", sep = "\t")
  
  #Note: If you have downloaded the code from Github, the cog-20.cog.csv and cog-20.org.csv are not present because they were too large
  #to upload. If you want to run this code, you can download those files yourself from this site: https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/
  #and drag them into the appropriate directory

  #1. Select out just those COG ID's that are relevant to the speices/strains we have in the app (use COGStrains file for guidance)
  
  #With KEGG Annotations
  #A. Bacteroides_thetaiotaomicron_VPI-5482: GCF_000011065.1
  #B. Pseudomonas_aeruginosa_PAO1: GCF_000006765.1
  #C. Porphyromonas_gingivalis_W83: GCF_000007585.1
  
  #Without KEGG Annotations
  #A. Clostridioides_difficile_630: GCF_000009205.2
  #B. Stenotrophomonas_maltophilia_K279a: GCF_000072485.1
  
  #2. See how well eligible COG ID's map to locus IDs from studies in app
  
  BT <- COGAnnotations[which(COGAnnotations$GCA_000007185.1 == "GCF_000011065.1"),] #e.g., BT_3759
  View(BT) #Matches KEGG ID
  PA <- COGAnnotations[which(COGAnnotations$GCA_000007185.1 == "GCF_000006765.1"),] #e.g., PA3039
  View(PA) #Matches KEGG ID
  PG <- COGAnnotations[which(COGAnnotations$GCA_000007185.1 == "GCF_000007585.1"),] #e.g., PG_RS02480
  View(PG) #Doesn't match KEGG ID or study locus ID (not useful)
  CD <- COGAnnotations[which(COGAnnotations$GCA_000007185.1 == "GCF_000009205.2"),] #e.g., CD630_20340
  View(CD) #Matches study locus ID (can include just COG annotations, not KEGG annotations)
  SM <- COGAnnotations[which(COGAnnotations$GCA_000007185.1 == "GCF_000072485.1"),] #e.g., SMLT_RS18455
  View(SM) #Matches study locus ID (can include just COG annotations, not KEGG annotations)
  
  #Adjust column names to reflect readme file for COG annotations
  colnames(BT) <- c("Locus_ID", "NCBI_Assembly_ID", "Protein_ID", "Protein_Length", "COG_Footprint_Coordinates", 
                    "Length_COG_Footprint", "COG_ID", "Reserved", "COG_Membership_Class","PSI_BLAST_Bit_Score",
                    "PSI-BLAST_e-value","COG_Profile_Length","Protein_Footprint_Coordinates")
  colnames(PA) <- c("Locus_ID", "NCBI_Assembly_ID", "Protein_ID", "Protein_Length", "COG_Footprint_Coordinates", 
                    "Length_COG_Footprint", "COG_ID", "Reserved", "COG_Membership_Class","PSI_BLAST_Bit_Score",
                    "PSI-BLAST_e-value","COG_Profile_Length","Protein_Footprint_Coordinates")
  colnames(PG) <- c("Locus_ID", "NCBI_Assembly_ID", "Protein_ID", "Protein_Length", "COG_Footprint_Coordinates", 
                    "Length_COG_Footprint", "COG_ID", "Reserved", "COG_Membership_Class","PSI_BLAST_Bit_Score",
                    "PSI-BLAST_e-value","COG_Profile_Length","Protein_Footprint_Coordinates")
  colnames(CD) <- c("Locus_ID", "NCBI_Assembly_ID", "Protein_ID", "Protein_Length", "COG_Footprint_Coordinates", 
                    "Length_COG_Footprint", "COG_ID", "Reserved", "COG_Membership_Class","PSI_BLAST_Bit_Score",
                    "PSI-BLAST_e-value","COG_Profile_Length","Protein_Footprint_Coordinates")
  colnames(SM) <- c("Locus_ID", "NCBI_Assembly_ID", "Protein_ID", "Protein_Length", "COG_Footprint_Coordinates", 
                    "Length_COG_Footprint", "COG_ID", "Reserved", "COG_Membership_Class","PSI_BLAST_Bit_Score",
                    "PSI-BLAST_e-value","COG_Profile_Length","Protein_Footprint_Coordinates")
  
  #Now simplify this table - pull out the COG info we want and merge with the COG descriptions
  BT1 <- BT[,c(1,7)]
  PA1 <- PA[,c(1,7)]
  PG1 <- PG[,c(1,7)]
  CD1 <- CD[,c(1,7)]
  SM1 <- SM[,c(1,7)]
  
  View(COGDescriptions)
  COGDescriptions1 <- COGDescriptions[,c(1,5)]
  colnames(COGDescriptions1) <- c("COG_ID", "Description")
  
  BT2 <- merge(BT1, COGDescriptions1, by = "COG_ID")
  PA2 <- merge(PA1, COGDescriptions1, by = "COG_ID")
  PG2 <- merge(PG1, COGDescriptions1, by = "COG_ID")
  CD2 <- merge(CD1, COGDescriptions1, by = "COG_ID")
  SM2 <- merge(SM1, COGDescriptions1, by = "COG_ID")
  
}

#4. Add in additional functional annotations (Uniprot) 
{
  S <- availableUniprotSpecies()
  
  SpSearch = c("Neosartorya fumigata", "Bacteroides thetaiotaomicron",
               "Bacteroides xylanisolvens",
               "Burkholderia thailandensis", "Burkholderia thailandensis",
               "Burkholderia pseudomallei", "Candida albicans", "Clostridioides difficile",
               "Fusobacterium nucleatum", "Haemophilus influenzae",
               "Mycobacteroides abscessus", "Pseudomonas aeruginosa",
               "Porphyromonas gingivalis", "Staphylococcus aureus",
               "Stenotrophomonas maltophilia", "Streptococcus pneumoniae",
               "Streptococcus salivarius","Streptococcus pyogenes", 
               "Streptococcus sanguinis"
  )
  
  SpRes <- lapply(SpSearch, function (x){
    S$`Species name`[grep(x, S$`Species name`, ignore.case = T)]
  } )
  
  names(SpRes) <- SpSearch
  
  StdCols <- c("UNIPROTKB_ID","UNIPROTKB", "ORGANISM", "ORGANISM-ID", "ID","EMBL/GENBANK/DDBJ", "EMBL/GENBANK/DDBJ_CDS",
               "GENENAME", "GENES", "GO", "GO-ID","KEGG", "KEGGID", 
               "EGGNOG", "EGGNOGO_ID", "PSEUDOCAP", "PSEUDOCAP_ID")
  
  myStrains <- c("Neosartorya fumigata (strain ATCC MYA-4609 / Af293 / CBS 101355 / FGSC A1100)",
                 "Neosartorya fumigata (strain CEA10 / CBS 144.89 / FGSC A1163)",
                 "Bacteroides thetaiotaomicron (strain ATCC 29148 / DSM 2079 / JCM 5827 / CCUG 10774 / NCTC 10582 / VPI-5482 / E50)",
                 "Burkholderia pseudomallei (strain K96243)",
                 "Burkholderia thailandensis (strain ATCC 700388 / DSM 13276 / CIP 106301 / E264)",
                 "Candida albicans (strain SC5314 / ATCC MYA-2876)",
                 "Clostridioides difficile (strain 630)",
                 "Haemophilus influenzae (strain 86-028NP)",
                 "Haemophilus influenzae (strain ATCC 51907 / DSM 11121 / KW20 / Rd)",
                 "Mycobacteroides abscessus (strain ATCC 19977 / DSM 44196 / CIP 104536 / JCM 13569 / NCTC 13031 / TMC 1543)",
                 "Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)",
                 "Pseudomonas aeruginosa (strain UCBPP-PA14)",
                 "Porphyromonas gingivalis (strain ATCC BAA-308 / W83)",
                 "Porphyromonas gingivalis (strain ATCC 33277 / DSM 20709 / CIP 103683 / JCM 12257 / NCTC 11834 / 2561)",
                 "Staphylococcus aureus (strain NCTC 8325 / PS 47)",
                 "Staphylococcus aureus (strain Newman)",
                 "Staphylococcus aureus (strain bovine RF122 / ET3-1)",
                 "Staphylococcus aureus (strain USA300 / TCH1516)",
                 "Stenotrophomonas maltophilia (strain K279a)",
                 "Streptococcus pneumoniae serotype 2 (strain D39 / NCTC 7466)",
                 "Streptococcus sanguinis (strain SK36)")
  
  createAnnotFileAll <- function(sname){
    taxon = S[S$`Species name` == sname,"taxon ID"]
    ws <- UniProt.ws(taxId=taxon)
    column_types <- columns(ws)
    k <- keys(ws, "UNIPROTKB_ID")
    a <- select(ws, k, StdCols[StdCols %in% column_types], "UNIPROTKB_ID")
    write.csv(a, file = paste(make.names(sname), "all", "csv", sep = '.'))
  }
  
  lapply(myStrains, createAnnotFileAll) #Create UniProt annotation files [Note: file names were altered outside of R before being merged with KEGG + COG annotations]
  
}

#5. Clean up + merge annotation data objects (KEGG + COG + Uniprot)
{
  
  #1. Strain: A. fumigatus, Af293 [KEGG + Uniprot]
  {
    #Create empty annotation data frame
    AF_Af293 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    AF_Af293_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/AFMPaths.csv")
    AF_Af293_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/AF_Af293.csv")
    AF_Af293_COG <- NULL
    
    #Merge all annotations into master table
    AF_Af293[1:nrow(AF_Af293_KEGG),] <- NA
    
    AF_Af293$KEGG_Pathway <- AF_Af293_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(AF_Af293)) { #Add in gene locus ID
      AF_Af293$Gene_Locus_ID[i] <- str_split(AF_Af293_KEGG$Gene, " ", 2)[[i]][1]
      AF_Af293$Gene_Locus_ID[i] <- substr(AF_Af293$Gene_Locus_ID[i], 5, nchar(AF_Af293$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(AF_Af293)) { #Add in gene symbol and gene function
      Gene <- str_split(AF_Af293_KEGG$Gene, " ", 2)[[i]][2]
      AF_Af293$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      AF_Af293$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    AF_Af293 <- as.data.table(AF_Af293)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(AF_Af293)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(AF_Af293_UNI)) {
      AF_Af293_UNI$KEGG[i] <- substr(AF_Af293_UNI$KEGG[i], 5, nchar(AF_Af293_UNI$KEGG[i]))
    }
    colnames(AF_Af293_UNI)[13] <- "Gene_Locus_ID"
    AF_Af293_UNI <- AF_Af293_UNI[,c(3,11,13,14)]
    AF_Af293 <- merge(AF_Af293, AF_Af293_UNI, by = "Gene_Locus_ID", all = T)
    AF_Af293$Uniprot_ID <- AF_Af293$UNIPROTKB
    AF_Af293$GO_Terms <- AF_Af293$GO
    AF_Af293$eggNOG <- AF_Af293$EGGNOG
    AF_Af293 <- AF_Af293[,-c(8:10)]
    
    #Remove duplicates
    AF_Af293 <- unique(AF_Af293)
  }
  
  #2. Strain: A. fumigatus, CEA10 + A1163 [Uniprot]
  {
    #Create empty annotation data frame
    AF_CEA10_A1163 <- data.frame(Uniprot_ID = character(),
                          GO_Terms = character(),
                          eggNOG = character(),
                          stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    AF_CEA10_A1163_KEGG <- NULL
    AF_CEA10_A1163_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/AF_CEA10_A1163.csv")
    AF_CEA10_A1163_COG <- NULL
    
    #Add in Uniprot annotations
    AF_CEA10_A1163[1:nrow(AF_CEA10_A1163_UNI),] <- NA
    
    AF_CEA10_A1163$Uniprot_ID <- AF_CEA10_A1163_UNI$UNIPROTKB
    AF_CEA10_A1163$GO_Terms <- AF_CEA10_A1163_UNI$GO
    AF_CEA10_A1163$eggNOG <- AF_CEA10_A1163_UNI$EGGNOG
    
  }
  
  #3. Strain: B. cenocepacia, AU1054 [KEGG]
  {
    #Create empty annotation data frame
    BC_AU1054 <- data.frame(Gene_Locus_ID = character(),
                          Gene_Symbol = character(),
                          Gene_Function = character(),
                          KEGG_Pathway = character(),
                          stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    BC_AU1054_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/BCNPaths.csv")
    BC_AU1054_UNI <- NULL
    BC_AU1054_COG <- NULL
    
    #Merge all annotations into master table
    BC_AU1054[1:nrow(BC_AU1054_KEGG),] <- NA
    
    BC_AU1054$KEGG_Pathway <- BC_AU1054_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(BC_AU1054)) { #Add in gene locus ID
      BC_AU1054$Gene_Locus_ID[i] <- str_split(BC_AU1054_KEGG$Gene, " ", 2)[[i]][1]
      BC_AU1054$Gene_Locus_ID[i] <- substr(BC_AU1054$Gene_Locus_ID[i], 5, nchar(BC_AU1054$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(BC_AU1054)) { #Add in gene symbol and gene function
      Gene <- str_split(BC_AU1054_KEGG$Gene, " ", 2)[[i]][2]
      BC_AU1054$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      BC_AU1054$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    BC_AU1054 <- as.data.table(BC_AU1054)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(BC_AU1054)[4] <- "KEGG_Pathway"
    
  }
  
  #4. Strain: B. cenocepacia, HI2424 [KEGG]
  {
    #Create empty annotation data frame
    BC_HI2424 <- data.frame(Gene_Locus_ID = character(),
                            Gene_Symbol = character(),
                            Gene_Function = character(),
                            KEGG_Pathway = character(),
                            stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    BC_HI2424_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/BCHPaths.csv")
    BC_HI2424_UNI <- NULL
    BC_HI2424_COG <- NULL
    
    #Merge all annotations into master table
    BC_HI2424[1:nrow(BC_HI2424_KEGG),] <- NA
    
    BC_HI2424$KEGG_Pathway <- BC_HI2424_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(BC_HI2424)) { #Add in gene locus ID
      BC_HI2424$Gene_Locus_ID[i] <- str_split(BC_HI2424_KEGG$Gene, " ", 2)[[i]][1]
      BC_HI2424$Gene_Locus_ID[i] <- substr(BC_HI2424$Gene_Locus_ID[i], 5, nchar(BC_HI2424$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(BC_HI2424)) { #Add in gene symbol and gene function
      Gene <- str_split(BC_HI2424_KEGG$Gene, " ", 2)[[i]][2]
      BC_HI2424$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      BC_HI2424$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    BC_HI2424 <- as.data.table(BC_HI2424)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(BC_HI2424)[4] <- "KEGG_Pathway"
    
  }
  
  #5. Strain: B. pseudomallei, K95243 [KEGG + Uniprot + COG] 
  {
    #Create empty annotation data frame
    BP_K95243 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    BP_K95243_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/BPSPaths.csv")
    BP_K95243_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/BP_K95243.csv")
    BP_K95243_COG <- NULL
    
    #Merge all annotations into master table
    BP_K95243[1:nrow(BP_K95243_KEGG),] <- NA
    
    BP_K95243$KEGG_Pathway <- BP_K95243_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(BP_K95243)) { #Add in gene locus ID
      BP_K95243$Gene_Locus_ID[i] <- str_split(BP_K95243_KEGG$Gene, " ", 2)[[i]][1]
      BP_K95243$Gene_Locus_ID[i] <- substr(BP_K95243$Gene_Locus_ID[i], 5, nchar(BP_K95243$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(BP_K95243)) { #Add in gene symbol and gene function
      Gene <- str_split(BP_K95243_KEGG$Gene, " ", 2)[[i]][2]
      BP_K95243$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      BP_K95243$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    BP_K95243 <- as.data.table(BP_K95243)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(BP_K95243)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(BP_K95243_UNI)) {
      BP_K95243_UNI$KEGG[i] <- substr(BP_K95243_UNI$KEGG[i], 5, nchar(BP_K95243_UNI$KEGG[i]))
    }
    colnames(BP_K95243_UNI)[13] <- "Gene_Locus_ID"
    BP_K95243_UNI <- BP_K95243_UNI[,c(3,11,13,14)]
    BP_K95243 <- merge(BP_K95243, BP_K95243_UNI, by = "Gene_Locus_ID", all = T)
    BP_K95243$Uniprot_ID <- BP_K95243$UNIPROTKB
    BP_K95243$GO_Terms <- BP_K95243$GO
    BP_K95243$eggNOG <- BP_K95243$EGGNOG
    BP_K95243 <- BP_K95243[,-c(8:10)]
    
    #Remove duplicates
    BP_K95243 <- unique(BP_K95243)

  }
  
  #6. Strain: B. thaliandensis, E264 [KEGG + Uniprot]
  {
    #Create empty annotation data frame
    BT_E264 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    BT_E264_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/BTEPaths.csv")
    BT_E264_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/BT_E264.csv")
    BT_E264_COG <- NULL
    
    #Merge all annotations into master table
    BT_E264[1:nrow(BT_E264_KEGG),] <- NA
    
    BT_E264$KEGG_Pathway <- BT_E264_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(BT_E264)) { #Add in gene locus ID
      BT_E264$Gene_Locus_ID[i] <- str_split(BT_E264_KEGG$Gene, " ", 2)[[i]][1]
      BT_E264$Gene_Locus_ID[i] <- substr(BT_E264$Gene_Locus_ID[i], 5, nchar(BT_E264$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(BT_E264)) { #Add in gene symbol and gene function
      Gene <- str_split(BT_E264_KEGG$Gene, " ", 2)[[i]][2]
      BT_E264$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      BT_E264$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    BT_E264 <- as.data.table(BT_E264)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(BT_E264)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(BT_E264_UNI)) {
      BT_E264_UNI$KEGG[i] <- substr(BT_E264_UNI$KEGG[i], 5, nchar(BT_E264_UNI$KEGG[i]))
    }
    colnames(BT_E264_UNI)[13] <- "Gene_Locus_ID"
    BT_E264_UNI <- BT_E264_UNI[,c(3,11,13,14)]
    BT_E264 <- merge(BT_E264, BT_E264_UNI, by = "Gene_Locus_ID", all = T)
    BT_E264$Uniprot_ID <- BT_E264$UNIPROTKB
    BT_E264$GO_Terms <- BT_E264$GO
    BT_E264$eggNOG <- BT_E264$EGGNOG
    BT_E264 <- BT_E264[,-c(8:10)]
    
    #Remove duplicates
    BT_E264 <- unique(BT_E264)
  }
  
  #7. Strain: B. thetaiotaomnicron, VPI-5482 [KEGG + Uniprot + COG] 
  {
    #Create empty annotation data frame
    BT_VPI5482 <- data.frame(Gene_Locus_ID = character(),
                          Gene_Symbol = character(),
                          Gene_Function = character(),
                          KEGG_Pathway = character(),
                          Uniprot_ID = character(),
                          GO_Terms = character(),
                          eggNOG = character(),
                          COG_ID = character(),
                          COG_Category = character(),
                          stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    BT_VPI5482_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/BTHPaths.csv")
    BT_VPI5482_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/BT_VPI-5482.csv")
    BT_VPI5482_COG <- BT2
    
    #Merge KEGG + COG + Uniprot annotations into master table
    BT_VPI5482[1:nrow(BT_VPI5482_KEGG),] <- NA
    
    BT_VPI5482$KEGG_Pathway <- BT_VPI5482_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(BT_VPI5482)) { #Add in gene locus ID
      BT_VPI5482$Gene_Locus_ID[i] <- str_split(BT_VPI5482_KEGG$Gene, " ", 2)[[i]][1]
      BT_VPI5482$Gene_Locus_ID[i] <- substr(BT_VPI5482$Gene_Locus_ID[i], 5, nchar(BT_VPI5482$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(BT_VPI5482)) { #Add in gene symbol and gene function
      Gene <- str_split(BT_VPI5482_KEGG$Gene, " ", 2)[[i]][2]
      BT_VPI5482$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      BT_VPI5482$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    colnames(BT_VPI5482_COG)[2] <- "Gene_Locus_ID" #Add in COG ID's and COG Categories
    BT_VPI5482 <- merge(BT_VPI5482, BT_VPI5482_COG, by = "Gene_Locus_ID", all = T)
    BT_VPI5482$COG_Category <- BT_VPI5482$Description
    BT_VPI5482$COG_ID.x <- BT_VPI5482$COG_ID.y
    BT_VPI5482 <- BT_VPI5482[,-c(10,11)]
    colnames(BT_VPI5482)[8] <- "COG_ID"
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    BT_VPI5482 <- as.data.table(BT_VPI5482)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG, COG_ID, COG_Category)]
    colnames(BT_VPI5482)[9] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(BT_VPI5482_UNI)) {
      BT_VPI5482_UNI$KEGG[i] <- substr(BT_VPI5482_UNI$KEGG[i], 5, nchar(BT_VPI5482_UNI$KEGG[i]))
    }
    colnames(BT_VPI5482_UNI)[13] <- "Gene_Locus_ID"
    BT_VPI5482_UNI <- BT_VPI5482_UNI[,c(3,11,13,14)]
    BT_VPI5482 <- merge(BT_VPI5482, BT_VPI5482_UNI, by = "Gene_Locus_ID", all = T)
    BT_VPI5482$Uniprot_ID <- BT_VPI5482$UNIPROTKB
    BT_VPI5482$GO_Terms <- BT_VPI5482$GO
    BT_VPI5482$eggNOG <- BT_VPI5482$EGGNOG
    BT_VPI5482 <- BT_VPI5482[,-c(10:12)]
    
    #Remove duplicates
    BT_VPI5482 <- unique(BT_VPI5482)
  }
  
  #8. Strain: B. xylanisolvens (General) [KEGG]
  {
    #Create empty annotation data frame
    BX_General <- data.frame(Gene_Locus_ID = character(),
                            Gene_Symbol = character(),
                            Gene_Function = character(),
                            KEGG_Pathway = character(),
                            stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    BX_General_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/BXYPaths.csv")
    BX_General_UNI <- NULL
    BX_General_COG <- NULL
    
    #Merge all annotations into master table
    BX_General[1:nrow(BX_General_KEGG),] <- NA
    
    BX_General$KEGG_Pathway <- BX_General_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(BX_General)) { #Add in gene locus ID
      BX_General$Gene_Locus_ID[i] <- str_split(BX_General_KEGG$Gene, " ", 2)[[i]][1]
      BX_General$Gene_Locus_ID[i] <- substr(BX_General$Gene_Locus_ID[i], 5, nchar(BX_General$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(BX_General)) { #Add in gene symbol and gene function
      Gene <- str_split(BX_General_KEGG$Gene, " ", 2)[[i]][2]
      BX_General$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      BX_General$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    BX_General <- as.data.table(BX_General)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(BX_General)[4] <- "KEGG_Pathway"
    
  }
  
  #9. Strain: C. albicans (General) [KEGG] 
  {
    #Create empty annotation data frame
    CA_General <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    CA_General_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/CALPaths.csv")
    CA_General_UNI <- NULL
    CA_General_COG <- NULL
    
    #Merge all annotations into master table
    CA_General[1:nrow(CA_General_KEGG),] <- NA
    
    CA_General$KEGG_Pathway <- CA_General_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(CA_General)) { #Add in gene locus ID
      CA_General$Gene_Locus_ID[i] <- str_split(CA_General_KEGG$Gene, " ", 2)[[i]][1]
      CA_General$Gene_Locus_ID[i] <- substr(CA_General$Gene_Locus_ID[i], 5, nchar(CA_General$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(CA_General)) { #Add in gene symbol and gene function
      Gene <- str_split(CA_General_KEGG$Gene, " ", 2)[[i]][2]
      CA_General$Gene_Symbol[i] <- str_split(Gene, "\\(", 2)[[1]][1]
      CA_General$Gene_Function[i] <- paste0("(", str_split(Gene, "\\(", 2)[[1]][2])
    }
    

    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    CA_General <- as.data.table(CA_General)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(CA_General)[4] <- "KEGG_Pathway"
    
  }
  
  #10. Strain: C. albicans, SC5314 [Uniprot]
  {
    #Create empty annotation data frame
    CA_SC5314 <- data.frame(Uniprot_ID = character(),
                                 GO_Terms = character(),
                                 eggNOG = character(),
                                 stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    CA_SC5314_KEGG <- NULL
    CA_SC5314_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/CA_SC5314.csv")
    CA_SC5314_COG <- NULL
    
    #Add in Uniprot annotations
    CA_SC5314[1:nrow(CA_SC5314_UNI),] <- NA
    
    CA_SC5314$Uniprot_ID <- CA_SC5314_UNI$UNIPROTKB
    CA_SC5314$GO_Terms <- CA_SC5314_UNI$GO
    CA_SC5314$eggNOG <- CA_SC5314_UNI$EGGNOG
    
  }
  
  #11. Strain: C. difficile, Strain 630 [COG]
  {
    #Create empty annotation data frame
    CD_630 <- data.frame(Gene_Locus_ID = character(),
                          COG_ID = character(),
                          COG_Category = character(),
                          stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    CD_630_KEGG <- NULL
    CD_630_UNI <- NULL
    CD_630_COG <- CD2
    
    #Merge all annotations into master table
    CD_630[1:nrow(CD_630_COG),] <- NA
    
    CD_630$Gene_Locus_ID <- CD_630_COG$Locus_ID
    CD_630$COG_ID <- CD_630_COG$COG_ID
    CD_630$COG_Category <- CD_630_COG$Description
    
    View(CD_630)
    
  }
  
  #12. Strain: H. influenza, Strain 723 [KEGG]
  {
    #Create empty annotation data frame
    HI_723 <- data.frame(Gene_Locus_ID = character(),
                            Gene_Symbol = character(),
                            Gene_Function = character(),
                            KEGG_Pathway = character(),
                            stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    HI_723_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/HIXPaths.csv")
    HI_723_UNI <- NULL
    HI_723_COG <- NULL
    
    #Merge all annotations into master table
    HI_723[1:nrow(HI_723_KEGG),] <- NA
    
    HI_723$KEGG_Pathway <- HI_723_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(HI_723)) { #Add in gene locus ID
      HI_723$Gene_Locus_ID[i] <- str_split(HI_723_KEGG$Gene, " ", 2)[[i]][1]
      HI_723$Gene_Locus_ID[i] <- substr(HI_723$Gene_Locus_ID[i], 5, nchar(HI_723$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(HI_723)) { #Add in gene symbol and gene function
      Gene <- str_split(HI_723_KEGG$Gene, " ", 2)[[i]][2]
      HI_723$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      HI_723$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    HI_723 <- as.data.table(HI_723)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(HI_723)[4] <- "KEGG_Pathway"
    
  }
  
  #13. Strain: H. influenza, R2866 [KEGG]
  {
    #Create empty annotation data frame
    HI_R2866 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    HI_R2866_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/HIZPaths.csv")
    HI_R2866_UNI <- NULL
    HI_R2866_COG <- NULL
    
    #Merge all annotations into master table
    HI_R2866[1:nrow(HI_R2866_KEGG),] <- NA
    
    HI_R2866$KEGG_Pathway <- HI_R2866_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(HI_R2866)) { #Add in gene locus ID
      HI_R2866$Gene_Locus_ID[i] <- str_split(HI_R2866_KEGG$Gene, " ", 2)[[i]][1]
      HI_R2866$Gene_Locus_ID[i] <- substr(HI_R2866$Gene_Locus_ID[i], 5, nchar(HI_R2866$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(HI_R2866)) { #Add in gene symbol and gene function
      Gene <- str_split(HI_R2866_KEGG$Gene, " ", 2)[[i]][2]
      HI_R2866$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      HI_R2866$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    HI_R2866 <- as.data.table(HI_R2866)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(HI_R2866)[4] <- "KEGG_Pathway"
    
  }
  
  #14. Strain: H. influenza, 86-028NP [KEGG + Uniprot] 
  {
    #Create empty annotation data frame
    HI_86028NP <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    HI_86028NP_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/HITPaths.csv")
    HI_86028NP_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/HI_86-028NP.csv")
    HI_86028NP_COG <- NULL
    
    #Merge all annotations into master table
    HI_86028NP[1:nrow(HI_86028NP_KEGG),] <- NA
    
    HI_86028NP$KEGG_Pathway <- HI_86028NP_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(HI_86028NP)) { #Add in gene locus ID
      HI_86028NP$Gene_Locus_ID[i] <- str_split(HI_86028NP_KEGG$Gene, " ", 2)[[i]][1]
      HI_86028NP$Gene_Locus_ID[i] <- substr(HI_86028NP$Gene_Locus_ID[i], 5, nchar(HI_86028NP$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(HI_86028NP)) { #Add in gene symbol and gene function
      Gene <- str_split(HI_86028NP_KEGG$Gene, " ", 2)[[i]][2]
      HI_86028NP$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      HI_86028NP$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    HI_86028NP <- as.data.table(HI_86028NP)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(HI_86028NP)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(HI_86028NP_UNI)) {
      HI_86028NP_UNI$KEGG[i] <- substr(HI_86028NP_UNI$KEGG[i], 5, nchar(HI_86028NP_UNI$KEGG[i]))
    }
    colnames(HI_86028NP_UNI)[13] <- "Gene_Locus_ID"
    HI_86028NP_UNI <- HI_86028NP_UNI[,c(3,11,13,14)]
    HI_86028NP <- merge(HI_86028NP, HI_86028NP_UNI, by = "Gene_Locus_ID", all = T)
    HI_86028NP$Uniprot_ID <- HI_86028NP$UNIPROTKB
    HI_86028NP$GO_Terms <- HI_86028NP$GO
    HI_86028NP$eggNOG <- HI_86028NP$EGGNOG
    HI_86028NP <- HI_86028NP[,-c(8:10)]
    
    #Remove duplicates
    HI_86028NP <- unique(HI_86028NP)
    
  }
  
  #15. Strain: H. influenza, Rd KW20 [KEGG + Uniprot + COG]
  {
    #Create empty annotation data frame
    HI_RdKW20 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    HI_RdKW20_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/HINPaths.csv")
    HI_RdKW20_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/HI_RdKW20.csv")
    HI_RdKW20_COG <- NULL
    
    #Merge all annotations into master table
    HI_RdKW20[1:nrow(HI_RdKW20_KEGG),] <- NA
    
    HI_RdKW20$KEGG_Pathway <- HI_RdKW20_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(HI_RdKW20)) { #Add in gene locus ID
      HI_RdKW20$Gene_Locus_ID[i] <- str_split(HI_RdKW20_KEGG$Gene, " ", 2)[[i]][1]
      HI_RdKW20$Gene_Locus_ID[i] <- substr(HI_RdKW20$Gene_Locus_ID[i], 5, nchar(HI_RdKW20$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(HI_RdKW20)) { #Add in gene symbol and gene function
      Gene <- str_split(HI_RdKW20_KEGG$Gene, " ", 2)[[i]][2]
      HI_RdKW20$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      HI_RdKW20$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    HI_RdKW20 <- as.data.table(HI_RdKW20)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(HI_RdKW20)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(HI_RdKW20_UNI)) {
      HI_RdKW20_UNI$KEGG[i] <- substr(HI_RdKW20_UNI$KEGG[i], 5, nchar(HI_RdKW20_UNI$KEGG[i]))
    }
    colnames(HI_RdKW20_UNI)[13] <- "Gene_Locus_ID"
    HI_RdKW20_UNI <- HI_RdKW20_UNI[,c(3,11,13,14)]
    HI_RdKW20 <- merge(HI_RdKW20, HI_RdKW20_UNI, by = "Gene_Locus_ID", all = T)
    HI_RdKW20$Uniprot_ID <- HI_RdKW20$UNIPROTKB
    HI_RdKW20$GO_Terms <- HI_RdKW20$GO
    HI_RdKW20$eggNOG <- HI_RdKW20$EGGNOG
    HI_RdKW20 <- HI_RdKW20[,-c(8:10)]

    #Remove duplicates
    HI_RdKW20 <- unique(HI_RdKW20)
  }
  
  #16. Strain: M. abscessus, ATCC 19977 [Uniprot] 
  {
    #Create empty annotation data frame
    MA_ATCC19977 <- data.frame(Uniprot_ID = character(),
                                 GO_Terms = character(),
                                 eggNOG = character(),
                                 stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    MA_ATCC19977_KEGG <- NULL
    MA_ATCC19977_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/MA_ATCC19977.csv")
    MA_ATCC19977_COG <- NULL
    
    #Add in Uniprot annotations
    MA_ATCC19977[1:nrow(MA_ATCC19977_UNI),] <- NA
    
    MA_ATCC19977$Uniprot_ID <- MA_ATCC19977_UNI$UNIPROTKB
    MA_ATCC19977$GO_Terms <- MA_ATCC19977_UNI$GO
    MA_ATCC19977$eggNOG <- MA_ATCC19977_UNI$EGGNOG
    
  }
  
  #17. Strain: P. aeruginosa, PAO1 [KEGG + Uniprot + COG]
  {
    #Create empty annotation data frame
    PA_PAO1 <- data.frame(Gene_Locus_ID = character(),
                          Gene_Symbol = character(),
                          Gene_Function = character(),
                          KEGG_Pathway = character(),
                          Uniprot_ID = character(),
                          GO_Terms = character(),
                          eggNOG = character(),
                          COG_ID = character(),
                          COG_Category = character(),
                          stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    PA_PAO1_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/PAEPaths.csv")
    PA_PAO1_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/PA_PAO1.csv")
    PA_PAO1_COG <- PA2
    
    #Merge KEGG + COG + Uniprot annotations into master table
    PA_PAO1[1:nrow(PA_PAO1_KEGG),] <- NA
    
    PA_PAO1$KEGG_Pathway <- PA_PAO1_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(PA_PAO1)) { #Add in gene locus ID
      PA_PAO1$Gene_Locus_ID[i] <- str_split(PA_PAO1_KEGG$Gene, " ", 2)[[i]][1]
      PA_PAO1$Gene_Locus_ID[i] <- substr(PA_PAO1$Gene_Locus_ID[i], 5, nchar(PA_PAO1$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(PA_PAO1)) { #Add in gene symbol and gene function
      Gene <- str_split(PA_PAO1_KEGG$Gene, " ", 2)[[i]][2]
      PA_PAO1$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      PA_PAO1$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    colnames(PA_PAO1_COG)[2] <- "Gene_Locus_ID" #Add in COG ID's and COG Categories
    PA_PAO1 <- merge(PA_PAO1, PA_PAO1_COG, by = "Gene_Locus_ID", all = T)
    PA_PAO1$COG_Category <- PA_PAO1$Description
    PA_PAO1$COG_ID.x <- PA_PAO1$COG_ID.y
    PA_PAO1 <- PA_PAO1[,-c(10,11)]
    colnames(PA_PAO1)[8] <- "COG_ID"
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    PA_PAO1 <- as.data.table(PA_PAO1)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG, COG_ID, COG_Category)]
    colnames(PA_PAO1)[9] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(PA_PAO1_UNI)) {
      PA_PAO1_UNI$KEGG[i] <- substr(PA_PAO1_UNI$KEGG[i], 5, nchar(PA_PAO1_UNI$KEGG[i]))
    }
    colnames(PA_PAO1_UNI)[13] <- "Gene_Locus_ID"
    PA_PAO1_UNI <- PA_PAO1_UNI[,c(3,11,13,14)]
    PA_PAO1 <- merge(PA_PAO1, PA_PAO1_UNI, by = "Gene_Locus_ID", all = T)
    PA_PAO1$Uniprot_ID <- PA_PAO1$UNIPROTKB
    PA_PAO1$GO_Terms <- PA_PAO1$GO
    PA_PAO1$eggNOG <- PA_PAO1$EGGNOG
    PA_PAO1 <- PA_PAO1[,-c(10:12)]
    
    #Remove duplicates
    PA_PAO1 <- unique(PA_PAO1)
  }
  
  #18. Strain: P. aeruginosa, PA14 [KEGG + Uniprot]
  {
    #Create empty annotation data frame
    PA_PA14 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    PA_PA14_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/PAUPaths.csv")
    PA_PA14_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/PA_PA14.csv")
    PA_PA14_COG <- NULL
    
    #Merge all annotations into master table
    PA_PA14[1:nrow(PA_PA14_KEGG),] <- NA
    
    PA_PA14$KEGG_Pathway <- PA_PA14_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(PA_PA14)) { #Add in gene locus ID
      PA_PA14$Gene_Locus_ID[i] <- str_split(PA_PA14_KEGG$Gene, " ", 2)[[i]][1]
      PA_PA14$Gene_Locus_ID[i] <- substr(PA_PA14$Gene_Locus_ID[i], 5, nchar(PA_PA14$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(PA_PA14)) { #Add in gene symbol and gene function
      Gene <- str_split(PA_PA14_KEGG$Gene, " ", 2)[[i]][2]
      PA_PA14$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      PA_PA14$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    PA_PA14 <- as.data.table(PA_PA14)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(PA_PA14)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(PA_PA14_UNI)) {
      PA_PA14_UNI$KEGG[i] <- substr(PA_PA14_UNI$KEGG[i], 5, nchar(PA_PA14_UNI$KEGG[i]))
    }
    colnames(PA_PA14_UNI)[13] <- "Gene_Locus_ID"
    PA_PA14_UNI <- PA_PA14_UNI[,c(3,11,13,14)]
    PA_PA14 <- merge(PA_PA14, PA_PA14_UNI, by = "Gene_Locus_ID", all = T)
    PA_PA14$Uniprot_ID <- PA_PA14$UNIPROTKB
    PA_PA14$GO_Terms <- PA_PA14$GO
    PA_PA14$eggNOG <- PA_PA14$EGGNOG
    PA_PA14 <- PA_PA14[,-c(8:10)]
    
    #Remove duplicates
    PA_PA14 <- unique(PA_PA14)
  }
  
  #19. Strain: P. aeruginosa, MTB-1 [KEGG]
  {
    #Create empty annotation data frame
    PA_MTB1 <- data.frame(Gene_Locus_ID = character(),
                            Gene_Symbol = character(),
                            Gene_Function = character(),
                            KEGG_Pathway = character(),
                            stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    PA_MTB1_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/PAEMPaths.csv")
    PA_MTB1_UNI <- NULL
    PA_MTB1_COG <- NULL
    
    #Merge all annotations into master table
    PA_MTB1[1:nrow(PA_MTB1_KEGG),] <- NA
    
    PA_MTB1$KEGG_Pathway <- PA_MTB1_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(PA_MTB1)) { #Add in gene locus ID
      PA_MTB1$Gene_Locus_ID[i] <- str_split(PA_MTB1_KEGG$Gene, " ", 2)[[i]][1]
      PA_MTB1$Gene_Locus_ID[i] <- substr(PA_MTB1$Gene_Locus_ID[i], 6, nchar(PA_MTB1$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(PA_MTB1)) { #Add in gene symbol and gene function
      Gene <- str_split(PA_MTB1_KEGG$Gene, " ", 2)[[i]][2]
      PA_MTB1$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      PA_MTB1$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    PA_MTB1 <- as.data.table(PA_MTB1)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(PA_MTB1)[4] <- "KEGG_Pathway"
    
  }
  
  #20. Strain: P. aeruginosa, B136-33 [KEGG]
  {
    #Create empty annotation data frame
    PA_B13633 <- data.frame(Gene_Locus_ID = character(),
                            Gene_Symbol = character(),
                            Gene_Function = character(),
                            KEGG_Pathway = character(),
                            stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    PA_B13633_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/PSGPaths.csv")
    PA_B13633_UNI <- NULL
    PA_B13633_COG <- NULL
    
    #Merge all annotations into master table
    PA_B13633[1:nrow(PA_B13633_KEGG),] <- NA
    
    PA_B13633$KEGG_Pathway <- PA_B13633_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(PA_B13633)) { #Add in gene locus ID
      PA_B13633$Gene_Locus_ID[i] <- str_split(PA_B13633_KEGG$Gene, " ", 2)[[i]][1]
      PA_B13633$Gene_Locus_ID[i] <- substr(PA_B13633$Gene_Locus_ID[i], 5, nchar(PA_B13633$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(PA_B13633)) { #Add in gene symbol and gene function
      Gene <- str_split(PA_B13633_KEGG$Gene, " ", 2)[[i]][2]
      PA_B13633$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      PA_B13633$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    PA_B13633 <- as.data.table(PA_B13633)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(PA_B13633)[4] <- "KEGG_Pathway"
    
  }
  
  #21. Strain: P. gingivalis, W83 [KEGG + Uniprot + COG]
  {
    #Create empty annotation data frame
    PG_W83 <- data.frame(Gene_Locus_ID = character(),
                          Gene_Symbol = character(),
                          Gene_Function = character(),
                          KEGG_Pathway = character(),
                          Uniprot_ID = character(),
                          GO_Terms = character(),
                          eggNOG = character(),
                          stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    PG_W83_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/PGIPaths.csv")
    PG_W83_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/PG_W83.csv")
    PG_W83_COG <- NULL
    
    #Merge KEGG + COG + Uniprot annotations into master table
    PG_W83[1:nrow(PG_W83_KEGG),] <- NA
    
    PG_W83$KEGG_Pathway <- PG_W83_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(PG_W83)) { #Add in gene locus ID
      PG_W83$Gene_Locus_ID[i] <- str_split(PG_W83_KEGG$Gene, " ", 2)[[i]][1]
      PG_W83$Gene_Locus_ID[i] <- substr(PG_W83$Gene_Locus_ID[i], 5, nchar(PG_W83$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(PG_W83)) { #Add in gene symbol and gene function
      Gene <- str_split(PG_W83_KEGG$Gene, " ", 2)[[i]][2]
      PG_W83$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      PG_W83$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    PG_W83 <- as.data.table(PG_W83)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(PG_W83)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(PG_W83_UNI)) {
      PG_W83_UNI$KEGG[i] <- substr(PG_W83_UNI$KEGG[i], 5, nchar(PG_W83_UNI$KEGG[i]))
    }
    colnames(PG_W83_UNI)[13] <- "Gene_Locus_ID"
    PG_W83_UNI <- PG_W83_UNI[,c(3,11,13,14)]
    PG_W83 <- merge(PG_W83, PG_W83_UNI, by = "Gene_Locus_ID", all = T)
    PG_W83$Uniprot_ID <- PG_W83$UNIPROTKB
    PG_W83$GO_Terms <- PG_W83$GO
    PG_W83$eggNOG <- PG_W83$EGGNOG
    PG_W83 <- PG_W83[,-c(8:10)]
    
    #Remove duplicates
    PG_W83 <- unique(PG_W83)
  }
  
  #22. Strain: P. gingivalis, ATCC 33277 [KEGG + Uniprot + COG]
  {
    #Create empty annotation data frame
    PG_ATCC33277 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    PG_ATCC33277_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/PGNPaths.csv")
    PG_ATCC33277_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/PG_ATCC33277.csv")
    PG_ATCC33277_COG <- NULL
    
    #Merge all annotations into master table
    PG_ATCC33277[1:nrow(PG_ATCC33277_KEGG),] <- NA
    
    PG_ATCC33277$KEGG_Pathway <- PG_ATCC33277_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(PG_ATCC33277)) { #Add in gene locus ID
      PG_ATCC33277$Gene_Locus_ID[i] <- str_split(PG_ATCC33277_KEGG$Gene, " ", 2)[[i]][1]
      PG_ATCC33277$Gene_Locus_ID[i] <- substr(PG_ATCC33277$Gene_Locus_ID[i], 5, nchar(PG_ATCC33277$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(PG_ATCC33277)) { #Add in gene symbol and gene function
      Gene <- str_split(PG_ATCC33277_KEGG$Gene, " ", 2)[[i]][2]
      PG_ATCC33277$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      PG_ATCC33277$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    PG_ATCC33277 <- as.data.table(PG_ATCC33277)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(PG_ATCC33277)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(PG_ATCC33277_UNI)) {
      PG_ATCC33277_UNI$KEGG[i] <- substr(PG_ATCC33277_UNI$KEGG[i], 5, nchar(PG_ATCC33277_UNI$KEGG[i]))
    }
    colnames(PG_ATCC33277_UNI)[13] <- "Gene_Locus_ID"
    PG_ATCC33277_UNI <- PG_ATCC33277_UNI[,c(3,11,13,14)]
    PG_ATCC33277 <- merge(PG_ATCC33277, PG_ATCC33277_UNI, by = "Gene_Locus_ID", all = T)
    PG_ATCC33277$Uniprot_ID <- PG_ATCC33277$UNIPROTKB
    PG_ATCC33277$GO_Terms <- PG_ATCC33277$GO
    PG_ATCC33277$eggNOG <- PG_ATCC33277$EGGNOG
    PG_ATCC33277 <- PG_ATCC33277[,-c(8:10)]
    
    #Remove duplicates
    PG_ATCC33277 <- unique(PG_ATCC33277)
  }
  
  #23. Strain: S. aureus, Newman [KEGG + Uniprot]
  {
    #Create empty annotation data frame
    SA_Newman <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SA_Newman_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SAEPaths.csv")
    SA_Newman_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/SA_Newman.csv")
    SA_Newman_COG <- NULL
    
    #Merge all annotations into master table
    SA_Newman[1:nrow(SA_Newman_KEGG),] <- NA
    
    SA_Newman$KEGG_Pathway <- SA_Newman_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SA_Newman)) { #Add in gene locus ID
      SA_Newman$Gene_Locus_ID[i] <- str_split(SA_Newman_KEGG$Gene, " ", 2)[[i]][1]
      SA_Newman$Gene_Locus_ID[i] <- substr(SA_Newman$Gene_Locus_ID[i], 5, nchar(SA_Newman$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SA_Newman)) { #Add in gene symbol and gene function
      Gene <- str_split(SA_Newman_KEGG$Gene, " ", 2)[[i]][2]
      SA_Newman$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SA_Newman$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    SA_Newman <- as.data.table(SA_Newman)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(SA_Newman)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(SA_Newman_UNI)) {
      SA_Newman_UNI$KEGG[i] <- substr(SA_Newman_UNI$KEGG[i], 5, nchar(SA_Newman_UNI$KEGG[i]))
    }
    colnames(SA_Newman_UNI)[13] <- "Gene_Locus_ID"
    SA_Newman_UNI <- SA_Newman_UNI[,c(3,11,13,14)]
    SA_Newman <- merge(SA_Newman, SA_Newman_UNI, by = "Gene_Locus_ID", all = T)
    SA_Newman$Uniprot_ID <- SA_Newman$UNIPROTKB
    SA_Newman$GO_Terms <- SA_Newman$GO
    SA_Newman$eggNOG <- SA_Newman$EGGNOG
    SA_Newman <- SA_Newman[,-c(8:10)]
    
    #Remove duplicates
    SA_Newman <- unique(SA_Newman)
  }
  
  #24. Strain: S. aureus, USA300 [KEGG + Uniprot]
  {
    #Create empty annotation data frame
    SA_USA300 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SA_USA300_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SAXPaths.csv")
    SA_USA300_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/SA_USA300.csv")
    SA_USA300_COG <- NULL
    
    #Merge all annotations into master table
    SA_USA300[1:nrow(SA_USA300_KEGG),] <- NA
    
    SA_USA300$KEGG_Pathway <- SA_USA300_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SA_USA300)) { #Add in gene locus ID
      SA_USA300$Gene_Locus_ID[i] <- str_split(SA_USA300_KEGG$Gene, " ", 2)[[i]][1]
      SA_USA300$Gene_Locus_ID[i] <- substr(SA_USA300$Gene_Locus_ID[i], 5, nchar(SA_USA300$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SA_USA300)) { #Add in gene symbol and gene function
      Gene <- str_split(SA_USA300_KEGG$Gene, " ", 2)[[i]][2]
      SA_USA300$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SA_USA300$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    SA_USA300 <- as.data.table(SA_USA300)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(SA_USA300)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(SA_USA300_UNI)) {
      SA_USA300_UNI$KEGG[i] <- substr(SA_USA300_UNI$KEGG[i], 5, nchar(SA_USA300_UNI$KEGG[i]))
    }
    colnames(SA_USA300_UNI)[13] <- "Gene_Locus_ID"
    SA_USA300_UNI <- SA_USA300_UNI[,c(3,11,13,14)]
    SA_USA300 <- merge(SA_USA300, SA_USA300_UNI, by = "Gene_Locus_ID", all = T)
    SA_USA300$Uniprot_ID <- SA_USA300$UNIPROTKB
    SA_USA300$GO_Terms <- SA_USA300$GO
    SA_USA300$eggNOG <- SA_USA300$EGGNOG
    SA_USA300 <- SA_USA300[,-c(8:10)]
    
    #Remove duplicates
    SA_USA300 <- unique(SA_USA300)
  }
  
  #25. Strain: S. aureus, NCTC 8325 [KEGG + Uniprot + COG]
  {
    #Create empty annotation data frame
    SA_NCTC8325 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SA_NCTC8325_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SAOPaths.csv")
    SA_NCTC8325_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/SA_NCTC8325.csv")
    SA_NCTC8325_COG <- NULL
    
    #Merge all annotations into master table
    SA_NCTC8325[1:nrow(SA_NCTC8325_KEGG),] <- NA
    
    SA_NCTC8325$KEGG_Pathway <- SA_NCTC8325_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SA_NCTC8325)) { #Add in gene locus ID
      SA_NCTC8325$Gene_Locus_ID[i] <- str_split(SA_NCTC8325_KEGG$Gene, " ", 2)[[i]][1]
      SA_NCTC8325$Gene_Locus_ID[i] <- substr(SA_NCTC8325$Gene_Locus_ID[i], 5, nchar(SA_NCTC8325$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SA_NCTC8325)) { #Add in gene symbol and gene function
      Gene <- str_split(SA_NCTC8325_KEGG$Gene, " ", 2)[[i]][2]
      SA_NCTC8325$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SA_NCTC8325$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    SA_NCTC8325 <- as.data.table(SA_NCTC8325)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(SA_NCTC8325)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(SA_NCTC8325_UNI)) {
      SA_NCTC8325_UNI$KEGG[i] <- substr(SA_NCTC8325_UNI$KEGG[i], 5, nchar(SA_NCTC8325_UNI$KEGG[i]))
    }
    colnames(SA_NCTC8325_UNI)[13] <- "Gene_Locus_ID"
    SA_NCTC8325_UNI <- SA_NCTC8325_UNI[,c(3,11,13,14)]
    SA_NCTC8325 <- merge(SA_NCTC8325, SA_NCTC8325_UNI, by = "Gene_Locus_ID", all = T)
    SA_NCTC8325$Uniprot_ID <- SA_NCTC8325$UNIPROTKB
    SA_NCTC8325$GO_Terms <- SA_NCTC8325$GO
    SA_NCTC8325$eggNOG <- SA_NCTC8325$EGGNOG
    SA_NCTC8325 <- SA_NCTC8325[,-c(8:10)]
    
    #Remove duplicates
    SA_NCTC8325 <- unique(SA_NCTC8325)
  }
  
  #26. Strain: S. aureus, JDK6008 [KEGG]
  {
    #Create empty annotation data frame
    SA_JDK6008 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SA_JDK6008_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SUKPaths.csv")
    SA_JDK6008_UNI <- NULL
    SA_JDK6008_COG <- NULL
    
    #Merge all annotations into master table
    SA_JDK6008[1:nrow(SA_JDK6008_KEGG),] <- NA
    
    SA_JDK6008$KEGG_Pathway <- SA_JDK6008_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SA_JDK6008)) { #Add in gene locus ID
      SA_JDK6008$Gene_Locus_ID[i] <- str_split(SA_JDK6008_KEGG$Gene, " ", 2)[[i]][1]
      SA_JDK6008$Gene_Locus_ID[i] <- substr(SA_JDK6008$Gene_Locus_ID[i], 5, nchar(SA_JDK6008$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SA_JDK6008)) { #Add in gene symbol and gene function
      Gene <- str_split(SA_JDK6008_KEGG$Gene, " ", 2)[[i]][2]
      SA_JDK6008$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SA_JDK6008$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    SA_JDK6008 <- as.data.table(SA_JDK6008)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(SA_JDK6008)[4] <- "KEGG_Pathway"
    
    
  }
  
  #27. Strain: S. aureus, RF122 [KEGG + Uniprot]
  {
    #Create empty annotation data frame
    SA_RF122 <- data.frame(Gene_Locus_ID = character(),
                         Gene_Symbol = character(),
                         Gene_Function = character(),
                         KEGG_Pathway = character(),
                         Uniprot_ID = character(),
                         GO_Terms = character(),
                         eggNOG = character(),
                         stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SA_RF122_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SABPaths.csv")
    SA_RF122_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/SA_RF122.csv")
    SA_RF122_COG <- NULL
    
    #Merge all annotations into master table
    SA_RF122[1:nrow(SA_RF122_KEGG),] <- NA
    
    SA_RF122$KEGG_Pathway <- SA_RF122_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SA_RF122)) { #Add in gene locus ID
      SA_RF122$Gene_Locus_ID[i] <- str_split(SA_RF122_KEGG$Gene, " ", 2)[[i]][1]
      SA_RF122$Gene_Locus_ID[i] <- substr(SA_RF122$Gene_Locus_ID[i], 5, nchar(SA_RF122$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SA_RF122)) { #Add in gene symbol and gene function
      Gene <- str_split(SA_RF122_KEGG$Gene, " ", 2)[[i]][2]
      SA_RF122$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SA_RF122$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    SA_RF122 <- as.data.table(SA_RF122)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(SA_RF122)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(SA_RF122_UNI)) {
      SA_RF122_UNI$KEGG[i] <- substr(SA_RF122_UNI$KEGG[i], 5, nchar(SA_RF122_UNI$KEGG[i]))
    }
    colnames(SA_RF122_UNI)[13] <- "Gene_Locus_ID"
    SA_RF122_UNI <- SA_RF122_UNI[,c(3,11,13,14)]
    SA_RF122 <- merge(SA_RF122, SA_RF122_UNI, by = "Gene_Locus_ID", all = T)
    SA_RF122$Uniprot_ID <- SA_RF122$UNIPROTKB
    SA_RF122$GO_Terms <- SA_RF122$GO
    SA_RF122$eggNOG <- SA_RF122$EGGNOG
    SA_RF122 <- SA_RF122[,-c(8:10)]
    
    #Remove duplicates
    SA_RF122 <- unique(SA_RF122)
  }
  
  #28. Strain: S. maltophilia, D457 [KEGG] 
  {
    #Create empty annotation data frame
    SM_D457 <- data.frame(Gene_Locus_ID = character(),
                            Gene_Symbol = character(),
                            Gene_Function = character(),
                            KEGG_Pathway = character(),
                            stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SM_D457_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SMZPaths.csv")
    SM_D457_UNI <- NULL
    SM_D457_COG <- NULL
    
    #Merge all annotations into master table
    SM_D457[1:nrow(SM_D457_KEGG),] <- NA
    
    SM_D457$KEGG_Pathway <- SM_D457_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SM_D457)) { #Add in gene locus ID
      SM_D457$Gene_Locus_ID[i] <- str_split(SM_D457_KEGG$Gene, " ", 2)[[i]][1]
      SM_D457$Gene_Locus_ID[i] <- substr(SM_D457$Gene_Locus_ID[i], 5, nchar(SM_D457$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SM_D457)) { #Add in gene symbol and gene function
      Gene <- str_split(SM_D457_KEGG$Gene, " ", 2)[[i]][2]
      SM_D457$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SM_D457$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    SM_D457 <- as.data.table(SM_D457)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(SM_D457)[4] <- "KEGG_Pathway"
    
  }
  
  #29. Strain: S. maltophilia, K279a [KEGG + Uniprot + COG]*
  {
    #Create empty annotation data frame
    SM_K279a <- data.frame(Gene_Locus_ID = character(),
                          Gene_Symbol = character(),
                          Gene_Function = character(),
                          KEGG_Pathway = character(),
                          Uniprot_ID = character(),
                          GO_Terms = character(),
                          eggNOG = character(),
                          COG_ID = character(),
                          COG_Category = character(),
                          stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SM_K279a_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SMLPaths.csv")
    SM_K279a_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/SM_K279a.csv")
    SM_K279a_COG <- SM2
    
    #Merge KEGG + COG + Uniprot annotations into master table
    SM_K279a[1:nrow(SM_K279a_KEGG),] <- NA
    
    SM_K279a$KEGG_Pathway <- SM_K279a_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SM_K279a)) { #Add in gene locus ID
      SM_K279a$Gene_Locus_ID[i] <- str_split(SM_K279a_KEGG$Gene, " ", 2)[[i]][1]
      SM_K279a$Gene_Locus_ID[i] <- substr(SM_K279a$Gene_Locus_ID[i], 5, nchar(SM_K279a$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SM_K279a)) { #Add in gene symbol and gene function
      Gene <- str_split(SM_K279a_KEGG$Gene, " ", 2)[[i]][2]
      SM_K279a$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SM_K279a$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    colnames(SM_K279a_COG)[2] <- "Gene_Locus_ID" #Add in COG ID's and COG Categories
    SM_K279a <- merge(SM_K279a, SM_K279a_COG, by = "Gene_Locus_ID", all = T)
    SM_K279a$COG_Category <- SM_K279a$Description
    SM_K279a$COG_ID.x <- SM_K279a$COG_ID.y
    SM_K279a <- SM_K279a[,-c(10,11)]
    colnames(SM_K279a)[8] <- "COG_ID"
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    SM_K279a <- as.data.table(SM_K279a)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG, COG_ID, COG_Category)]
    colnames(SM_K279a)[9] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(SM_K279a_UNI)) {
      SM_K279a_UNI$KEGG[i] <- substr(SM_K279a_UNI$KEGG[i], 5, nchar(SM_K279a_UNI$KEGG[i]))
    }
    colnames(SM_K279a_UNI)[13] <- "Gene_Locus_ID"
    SM_K279a_UNI <- SM_K279a_UNI[,c(3,11,13,14)]
    SM_K279a <- merge(SM_K279a, SM_K279a_UNI, by = "Gene_Locus_ID", all = T)
    SM_K279a$Uniprot_ID <- SM_K279a$UNIPROTKB
    SM_K279a$GO_Terms <- SM_K279a$GO
    SM_K279a$eggNOG <- SM_K279a$EGGNOG
    SM_K279a <- SM_K279a[,-c(10:12)]
    
    #Remove duplicates
    SM_K279a <- unique(SM_K279a)
    
  }
  
  #30. Strain: S. pneumonia D39 [KEGG + Uniprot + COG]
  {
    #Create empty annotation data frame
    SP_D39 <- data.frame(Gene_Locus_ID = character(),
                          Gene_Symbol = character(),
                          Gene_Function = character(),
                          KEGG_Pathway = character(),
                          Uniprot_ID = character(),
                          GO_Terms = character(),
                          eggNOG = character(),
                          stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SP_D39_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SPDPaths.csv")
    SP_D39_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/SP_D39.csv")
    SP_D39_COG <- NULL
    
    #Merge all annotations into master table
    SP_D39[1:nrow(SP_D39_KEGG),] <- NA
    
    SP_D39$KEGG_Pathway <- SP_D39_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SP_D39)) { #Add in gene locus ID
      SP_D39$Gene_Locus_ID[i] <- str_split(SP_D39_KEGG$Gene, " ", 2)[[i]][1]
      SP_D39$Gene_Locus_ID[i] <- substr(SP_D39$Gene_Locus_ID[i], 5, nchar(SP_D39$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SP_D39)) { #Add in gene symbol and gene function
      Gene <- str_split(SP_D39_KEGG$Gene, " ", 2)[[i]][2]
      SP_D39$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SP_D39$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line (to facilitate merge with Uniprot data)
    SP_D39 <- as.data.table(SP_D39)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function, Uniprot_ID, GO_Terms, eggNOG)]
    colnames(SP_D39)[7] <- "KEGG_Pathway"
    
    #Now merge in Uniprot Data
    for (i in 1:nrow(SP_D39_UNI)) {
      SP_D39_UNI$KEGG[i] <- substr(SP_D39_UNI$KEGG[i], 5, nchar(SP_D39_UNI$KEGG[i]))
    }
    colnames(SP_D39_UNI)[13] <- "Gene_Locus_ID"
    SP_D39_UNI <- SP_D39_UNI[,c(3,11,13,14)]
    SP_D39 <- merge(SP_D39, SP_D39_UNI, by = "Gene_Locus_ID", all = T)
    SP_D39$Uniprot_ID <- SP_D39$UNIPROTKB
    SP_D39$GO_Terms <- SP_D39$GO
    SP_D39$eggNOG <- SP_D39$EGGNOG
    SP_D39 <- SP_D39[,-c(8:10)]
    
    #Remove duplicates
    SP_D39 <- unique(SP_D39)
  }
  
  #31. Strain: S. pyogenes M1 GAS [KEGG]
  {
    #Create empty annotation data frame
    SP_M1GAS <- data.frame(Gene_Locus_ID = character(),
                            Gene_Symbol = character(),
                            Gene_Function = character(),
                            KEGG_Pathway = character(),
                            stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SP_M1GAS_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SPYPaths.csv")
    SP_M1GAS_UNI <- NULL
    SP_M1GAS_COG <- NULL
    
    #Merge all annotations into master table
    SP_M1GAS[1:nrow(SP_M1GAS_KEGG),] <- NA
    
    SP_M1GAS$KEGG_Pathway <- SP_M1GAS_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SP_M1GAS)) { #Add in gene locus ID
      SP_M1GAS$Gene_Locus_ID[i] <- str_split(SP_M1GAS_KEGG$Gene, " ", 2)[[i]][1]
      SP_M1GAS$Gene_Locus_ID[i] <- substr(SP_M1GAS$Gene_Locus_ID[i], 5, nchar(SP_M1GAS$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SP_M1GAS)) { #Add in gene symbol and gene function
      Gene <- str_split(SP_M1GAS_KEGG$Gene, " ", 2)[[i]][2]
      SP_M1GAS$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SP_M1GAS$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    SP_M1GAS <- as.data.table(SP_M1GAS)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(SP_M1GAS)[4] <- "KEGG_Pathway"
    
  }
  
  #32. Strain: S. salivarius, HSISS4 [KEGG]
  {
    #Create empty annotation data frame
    SS_HSISS4 <- data.frame(Gene_Locus_ID = character(),
                            Gene_Symbol = character(),
                            Gene_Function = character(),
                            KEGG_Pathway = character(),
                            stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SS_HSISS4_KEGG <- read.csv("./Pathway_Data/KEGG_Annotations/SSAHPaths.csv")
    SS_HSISS4_UNI <- NULL
    SS_HSISS4_COG <- NULL
    
    #Merge all annotations into master table
    SS_HSISS4[1:nrow(SS_HSISS4_KEGG),] <- NA
    
    SS_HSISS4$KEGG_Pathway <- SS_HSISS4_KEGG$Pathway #Add in KEGG Pathway Function
    
    for(i in 1:nrow(SS_HSISS4)) { #Add in gene locus ID
      SS_HSISS4$Gene_Locus_ID[i] <- str_split(SS_HSISS4_KEGG$Gene, " ", 2)[[i]][1]
      SS_HSISS4$Gene_Locus_ID[i] <- substr(SS_HSISS4$Gene_Locus_ID[i], 6, nchar(SS_HSISS4$Gene_Locus_ID[i]))
    }
    
    for (i in 1:nrow(SS_HSISS4)) { #Add in gene symbol and gene function
      Gene <- str_split(SS_HSISS4_KEGG$Gene, " ", 2)[[i]][2]
      SS_HSISS4$Gene_Symbol[i] <- str_split(Gene, " ", 2)[[1]][1]
      SS_HSISS4$Gene_Function[i] <- str_split(Gene, " ", 2)[[1]][2]
    }
    
    #Compress data frame so that KEGG pathway functions are all on one line 
    SS_HSISS4 <- as.data.table(SS_HSISS4)[, toString(KEGG_Pathway), by = list(Gene_Locus_ID, Gene_Symbol, Gene_Function)]
    colnames(SS_HSISS4)[4] <- "KEGG_Pathway"
    
  }
  
  #33. Strain: S. sanguinus, SK36 [Uniprot + COG]
  {
    #Create empty annotation data frame
    SS_SK36 <- data.frame(Uniprot_ID = character(),
                                 GO_Terms = character(),
                                 eggNOG = character(),
                                 stringsAsFactors = FALSE
    )
    
    #Load in annotation data
    SS_SK36_KEGG <- NULL
    SS_SK36_UNI <- read.csv("./Pathway_Data/UniProt_Annotations/SS_SK36.csv")
    SS_SK36_COG <- NULL
    
    #Add in Uniprot annotations
    SS_SK36[1:nrow(SS_SK36_UNI),] <- NA
    
    SS_SK36$Uniprot_ID <-SS_SK36_UNI$UNIPROTKB
    SS_SK36$GO_Terms <- SS_SK36_UNI$GO
    SS_SK36$eggNOG <- SS_SK36_UNI$EGGNOG
  
    #Remove duplicates
    SS_SK36 <- unique(SS_SK36)
  }
  
}

#6. Write annotation files into single directory [Data Setup.R loads these files into a data object that interfaces with the application]
{
  write.csv(AF_Af293, "./Pathway_Data/Full_Annotations/AF_Af293.csv") 
  write.csv(AF_CEA10_A1163, "./Pathway_Data/Full_Annotations/AF_CEA10_A1163.csv")
  write.csv(BC_AU1054, "./Pathway_Data/Full_Annotations/BC_AU1054.csv")
  write.csv(BC_HI2424, "./Pathway_Data/Full_Annotations/BC_HI2424.csv") 
  write.csv(BP_K95243, "./Pathway_Data/Full_Annotations/BP_K95243.csv")
  write.csv(BT_E264, "./Pathway_Data/Full_Annotations/BT_E264.csv") 
  write.csv(BT_VPI5482, "./Pathway_Data/Full_Annotations/BT_VPI5482.csv")
  write.csv(BX_General, "./Pathway_Data/Full_Annotations/BX_General.csv")
  write.csv(CA_General, "./Pathway_Data/Full_Annotations/CA_General.csv")
  write.csv(CA_SC5314, "./Pathway_Data/Full_Annotations/CA_SC5314.csv")
  write.csv(CD_630, "./Pathway_Data/Full_Annotations/CD_630.csv")
  write.csv(CD_R20291, "./Pathway_Data/Full_Annotations/CD_R20291.csv")
  write.csv(HI_723, "./Pathway_Data/Full_Annotations/HI_723.csv")
  write.csv(HI_R2866, "./Pathway_Data/Full_Annotations/HI_R2866.csv")
  write.csv(HI_86028NP, "./Pathway_Data/Full_Annotations/HI_86028NP.csv")
  write.csv(HI_RdKW20, "./Pathway_Data/Full_Annotations/HI_RdKW20.csv")
  write.csv(MA_ATCC19977, "./Pathway_Data/Full_Annotations/MA_ATCC19977.csv")
  write.csv(PA_PAO1, "./Pathway_Data/Full_Annotations/PA_PAO1.csv")
  write.csv(PA_PA14, "./Pathway_Data/Full_Annotations/PA_PA14.csv")
  write.csv(PA_MTB1, "./Pathway_Data/Full_Annotations/PA_MTB1.csv") #*
  write.csv(PA_B13633, "./Pathway_Data/Full_Annotations/PA_B13633.csv") 
  write.csv(PA_PAK, "./Pathway_Data/Full_Annotations/PA_PAK.csv")
  write.csv(PG_W83, "./Pathway_Data/Full_Annotations/PG_W83.csv")
  write.csv(PG_ATCC33277, "./Pathway_Data/Full_Annotations/PG_ATCC33277.csv")
  write.csv(SA_Newman, "./Pathway_Data/Full_Annotations/SA_Newman.csv")
  write.csv(SA_USA300, "./Pathway_Data/Full_Annotations/SA_USA300.csv")
  write.csv(SA_NCTC8325, "./Pathway_Data/Full_Annotations/SA_NCTC8325.csv")
  write.csv(SA_JDK6008, "./Pathway_Data/Full_Annotations/SA_JDK6008.csv")
  write.csv(SA_RF122, "./Pathway_Data/Full_Annotations/SA_RF122.csv")
  write.csv(SM_D457, "./Pathway_Data/Full_Annotations/SM_D457.csv")
  write.csv(SM_K279a, "./Pathway_Data/Full_Annotations/SM_K279a.csv")
  write.csv(SP_D39, "./Pathway_Data/Full_Annotations/SP_D39.csv")
  write.csv(SP_M1GAS, "./Pathway_Data/Full_Annotations/SP_M1GAS.csv")
  write.csv(SS_HSISS4, "./Pathway_Data/Full_Annotations/SS_HSISS4.csv")
  write.csv(SS_SK36, "./Pathway_Data/Full_Annotations/SS_SK36.csv")
  
  
}



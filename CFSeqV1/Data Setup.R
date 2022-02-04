#Note: If you would like to run the application from R after downloading the working directory from Github, you need to run through 
#all the steps on this file. These steps import data files that are included in the working directory from Github - so you need to 
#make sure they are downloaded and present in the working directory. These files constitute the whole compendium of RNA-seq data sets.

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

#2. Create CFSeq_Data object [Code will take less than 30 seconds to run]
{
  CFSeq_Data <- vector(mode = "list", 13)
  
  names(CFSeq_Data) <- c("Aspergillus_fumigatus", "Bacteroides_Species", "Burkholderia_Species","Candida_albicans","Clostridium_difficile", "Fusobacterium_nucleatum",
                         "Haemophilus_influenzae","Mycobacterium_abscessus", "Pseudomonas_aeruginosa",
                         "Porphyromonas_Species","Staphylococcus_aureus","Stenotrophomonas_maltophilia","Streptococcus_Species")
}

#3. Load in all count tables, design matrices, and metadata for each species [Code will take less than 30 seconds to run]
{
  
  #---
  #Aspergillus Fumigatus 
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/AFCounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/AFCounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/AFCounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[1]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[1]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[1]])) {
      CFSeq_Data[[1]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[1]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[1]])) {
      CFSeq_Data[[1]][[i]]$Count_Table <- read.csv(paste0("./CountTables/AFCounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[1]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/AFCounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[1]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/AFCounts/Additional_Metadata/", AMs[i]))
    
    }
    
  }
  
  #---
  #Bacteroides Species
  {
    
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/BacteroidesCounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/BacteroidesCounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/BacteroidesCounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[2]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[2]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[2]])) {
      CFSeq_Data[[2]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[2]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }

    #Load in data
    for (i in 1:length(CFSeq_Data[[2]])) {
      CFSeq_Data[[2]][[i]]$Count_Table <- read.csv(paste0("./CountTables/BacteroidesCounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[2]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/BacteroidesCounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[2]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/BacteroidesCounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
    
  #--- 
  #Burkholderia Species
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/BurkholderiaCounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/BurkholderiaCounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/BurkholderiaCounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[3]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[3]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[3]])) {
      CFSeq_Data[[3]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[3]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[3]])) {
      CFSeq_Data[[3]][[i]]$Count_Table <- read.csv(paste0("./CountTables/BurkholderiaCounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[3]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/BurkholderiaCounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[3]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/BurkholderiaCounts/Additional_Metadata/", AMs[i]))
    
    }
    
  }
  
  #--- 
  #Candida Albicans 
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/CACounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/CACounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/CACounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[4]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[4]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[4]])) {
      CFSeq_Data[[4]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[4]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[4]])) {
      CFSeq_Data[[4]][[i]]$Count_Table <- read.csv(paste0("./CountTables/CACounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[4]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/CACounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[4]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/CACounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Clostridioides difficile
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/CDCounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/CDCounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/CDCounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[5]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[5]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[5]])) {
      CFSeq_Data[[5]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[5]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[5]])) {
      CFSeq_Data[[5]][[i]]$Count_Table <- read.csv(paste0("./CountTables/CDCounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[5]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/CDCounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[5]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/CDCounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Fusobacterium nucleatum 
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/FNCounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/FNCounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/FNCounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[6]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[6]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[6]])) {
      CFSeq_Data[[6]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[6]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[6]])) {
      CFSeq_Data[[6]][[i]]$Count_Table <- read.csv(paste0("./CountTables/FNCounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[6]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/FNCounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[6]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/FNCounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Haemophilus influenzae
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/HICounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/HICounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/HICounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[7]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[7]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[7]])) {
      CFSeq_Data[[7]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[7]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[7]])) {
      CFSeq_Data[[7]][[i]]$Count_Table <- read.csv(paste0("./CountTables/HICounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[7]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/HICounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[7]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/HICounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Mycobacterium Abscessus
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/MACounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/MACounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/MACounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[8]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[8]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[8]])) {
      CFSeq_Data[[8]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[8]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[8]])) {
      CFSeq_Data[[8]][[i]]$Count_Table <- read.csv(paste0("./CountTables/MACounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[8]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/MACounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[8]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/MACounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Pseudomonas Aeruginosa
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/PACounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/PACounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/PACounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[9]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[9]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[9]])) {
      CFSeq_Data[[9]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[9]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[9]])) {
      CFSeq_Data[[9]][[i]]$Count_Table <- read.csv(paste0("./CountTables/PACounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[9]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/PACounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[9]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/PACounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Porphyromonas sp. 
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/PorphyromonasCounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/PorphyromonasCounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/PorphyromonasCounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[10]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[10]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[10]])) {
      CFSeq_Data[[10]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[10]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[10]])) {
      CFSeq_Data[[10]][[i]]$Count_Table <- read.csv(paste0("./CountTables/PorphyromonasCounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[10]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/PorphyromonasCounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[10]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/PorphyromonasCounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Staphylococcus aureus
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/SACounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/SACounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/SACounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[11]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[11]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[11]])) {
      CFSeq_Data[[11]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[11]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[11]])) {
      CFSeq_Data[[11]][[i]]$Count_Table <- read.csv(paste0("./CountTables/SACounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[11]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/SACounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[11]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/SACounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Stenotrophomonas maltophilia
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/SMCounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/SMCounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/SMCounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[12]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[12]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[12]])) {
      CFSeq_Data[[12]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[12]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[12]])) {
      CFSeq_Data[[12]][[i]]$Count_Table <- read.csv(paste0("./CountTables/SMCounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[12]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/SMCounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[12]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/SMCounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  #--- 
  #Streptococcus sp. 
  {
    
    #Add slots for studies
    Counts <- str_sort(list.files(path = "./CountTables/StrepCounts/Count_Tables", pattern = "*.csv"))
    DMs <- str_sort(list.files(path = "./CountTables/StrepCounts/Design_Matrices", pattern = "*.csv"))
    AMs <- str_sort(list.files(path = "./CountTables/StrepCounts/Additional_Metadata", pattern = "*.csv"))
    CFSeq_Data[[13]] <- vector(mode = "list", length(str_extract_all(string = Counts, pattern = "GSE[0-9]+")))
    names(CFSeq_Data[[13]]) <- str_extract_all(string = Counts, pattern = "GSE[0-9]+")
    
    for (i in 1:length(CFSeq_Data[[13]])) {
      CFSeq_Data[[13]][[i]] = vector(mode = "list", 3)
      names(CFSeq_Data[[13]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
    }
    
    #Load in data
    for (i in 1:length(CFSeq_Data[[13]])) {
      CFSeq_Data[[13]][[i]]$Count_Table <- read.csv(paste0("./CountTables/StrepCounts/Count_Tables/", Counts[i]))
      CFSeq_Data[[13]][[i]]$Design_Matrix <- read.csv(paste0("./CountTables/StrepCounts/Design_Matrices/", DMs[i]))
      CFSeq_Data[[13]][[i]]$Additional_Metadata <- read.csv(paste0("./CountTables/StrepCounts/Additional_Metadata/", AMs[i]))
      
    }
    
  }
  
  
}

#4. Perform DE Analysis By Species (Stored in CFSeq_Data) [Code will take approximately 5 minutes to run]
{
  #Part 1
  {
    
    #---
    #Aspergillus fumigatus
    {
      Species <- 1
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      } #Perform DE Analysis
    }
    
    #---
    #Bacteroides Species
    {
      Species <- 2
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      } 
    }
    
    #--- 
    #Burkholderia Species
    {
      Species <- 3
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      } 
    }
    
    #--- 
    #Candida Albicans
    {
      Species <- 4
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    }
    
    #--- 
    #Clostridioides difficile
    {
      Species <- 5
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    }
    
    #--- 
    #Fusobacterium nucleatum
    {
      Species <- 6
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      } 
    }
    
    #--- 
    #Haemophilus influenzae 
    {
      Species <- 7
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    }
    
    #--- 
    #Mycobacterium Abscessus
    {
      Species <- 8
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    }
    
  }
  
  #Part 2 
  {
    
    #--- 
    #Pseudomonas aeruginosa
    {
      Species <- 9
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {

        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    } 
    
    #--- 
    #Porphyromonas sp.
    {
      Species <- 10
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    }
    
    #--- 
    #Staphylococcus aureus
    {
      Species <- 11
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    }
    
    #--- 
    #Stenotrophomonas maltophilia 
    {
      Species <- 12
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    }
    
    #--- 
    #Streptococcus sp.
    {
      Species <- 13
      
      for (i in 1:(length(CFSeq_Data[[Species]]))) {
        
        if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate" & ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix) > 1) {
          length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[2]))
          #Create DGE Object
          CT <- CFSeq_Data[[Species]][[i]]$Count_Table
          rownames(CT) <- make.names(CT[,1], unique = TRUE)
          CT <- CT[,-1]
          for (a in 1:length(CT)) {
            if (typeof(CT[,a]) != "integer") {
              CT[,a] <- as.integer(CT[,a])
            }
          }
          CT <- na.omit(CT)
          DGE <- DGEList(counts = CT, remove.zeros = TRUE)
          #Set Factors
          Factors <- c()
          for (j in 1:length(which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate"))) {
            Factors[j] <- names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])[j]
          }
          
          #Determine class
          Class <- NULL
          if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2 |
              length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 3) {
            Class = 3
          } else if (length(names(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) > 3) {
            Class = 4
          } else if (length(unique(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) == 2) {
            Class = 1
          } else {
            Class = 2
          }
          
          #Set up Design Matrix
          if (is.null(Factors) != TRUE) {
            Contrasts <- vector(mode = "list", length = length(Factors))
          } else {
            Contrasts <- vector(mode = "list", length = 1)
          }
          for (k in 1:length(Contrasts)) {
            if (is.null(ncol(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])) != TRUE) {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")][,k])
            } else {
              Contrasts[[k]] <- factor(CFSeq_Data[[Species]][[i]]$Design_Matrix[, which(names(CFSeq_Data[[Species]][[i]]$Design_Matrix) != "Replicate")])
            }
          }
          if (length(Contrasts) == 1) {
            Var1 <- unlist(Contrasts[1])
            design <- model.matrix(~0+Var1)
          } else if (length(Contrasts) == 2) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            design <- model.matrix(~0+Var1 + Var2)
          } else if (length(Contrasts) == 3) {
            Var1 <- unlist(Contrasts[1])
            Var2 <- unlist(Contrasts[2])
            Var3 <- unlist(Contrasts[3])
            design <- model.matrix(~0+Var1 + Var2 + Var3)
          } else {
            design <- NULL
          }
          rownames(design) <- colnames(DGE)
          colnames(design) <- make.names(colnames(design), unique = FALSE, allow_ = TRUE)
          colnames(design) <- str_replace(colnames(design), "Var1", "")
          #Create outputs that differ by class
          if (Class == 1) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- colnames(fit$coefficients)
            qlf <- glmQLFTest(fit, contrast = c(-1,1))
            Results <- as.data.frame(qlf$table) 
            DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[2], "-", names[1], "Full Results")
            
            DE2 <- DE
            DE2$logFC <- DE$logFC*(-1)
            Results2 <- Results
            Results2$logFC <- Results$logFC*(-1)
            
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "DE Output")
            CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results2))
            names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[1], "-", names[2], "Full Results")
            
          } else if (Class == 2) {
            
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else if (Class == 3) {
            DM <- CFSeq_Data[[Species]][[i]]$Design_Matrix
            #Adjust design matrix to make new groups column (combining multiple factors)
            if (names(DM[1]) == "Replicate") {
              Group <- factor(paste(DM[,2],DM[,3], sep = "."))
            } else {
              Group <- factor(paste(DM[,1],DM[,2], sep = "."))
            }
            cbind(DM, Group = Group)
            design <- model.matrix(~0 + Group)
            colnames(design) <- levels(Group)
            #Now build analysis output
            keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
            DGE <- DGE[keep,]
            DGE <- calcNormFactors(DGE)
            DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 
            DGE <- estimateGLMTrendedDisp(DGE, design) 
            DGE <- estimateGLMTagwiseDisp(DGE, design) 
            fit <- glmQLFit(DGE, design)
            names <- make.names(colnames(fit$coefficients), unique = FALSE, allow_ = TRUE)
            for (l in 1:length(names)) {
              for (m in 1:length(names)) {
                if (l != m) {
                  Vec <- c(rep(0, ncol(design)))
                  Vec[l] = 1
                  Vec[m] = -1
                  qlf <- glmQLFTest(fit, contrast = Vec)
                  Results <- as.data.frame(qlf$table) 
                  DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(DE))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "DE Output")
                  CFSeq_Data[[Species]][[i]] <- c(CFSeq_Data[[Species]][[i]], list(Results))
                  names(CFSeq_Data[[Species]][[i]])[length(CFSeq_Data[[Species]][[i]])] <- paste(names[l], "-", names[m], "Full Results")
                }
              }
            }
            
          } else {
            
            CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
            CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
            
          }
          
        } else 
          
          CFSeq_Data[[Species]][[i]]$DE_Genes <- NULL
        CFSeq_Data[[Species]][[i]]$All_Genes <- NULL
      }
    }
    
  }
  
}

#5. Create gene-pathway info list object [Code will take less than 30 seconds to run]
{
  Files <- str_sort(list.files(path = "./PathwayData", pattern = "*.csv"))
  GenePathway_Data <- vector(mode = "list", length(Files))
  names(GenePathway_Data) <- Files
  
  for (i in 1:length(GenePathway_Data)) {
    GenePathway_Data[[i]] <- read.csv(paste0("./PathwayData/", Files[i]))
  }
  
  for (i in 1:length(GenePathway_Data)) {
    GenePathway_Data[[i]] <- GenePathway_Data[[i]][,-1]
    for (j in 1:nrow(GenePathway_Data[[i]])) {
      GenePathway_Data[[i]]$Gene[j] <- strsplit(GenePathway_Data[[i]]$Gene[j], ":")[[1]][2]
    }
  }
  
  Species_Numbers <- c(1, 3, 3, 3, 3, 2, 2, 4, 6, 7, 7, 7, 7, 9, 9, 9, 10, 10, 9, 11, 11, 13, 11, 11, 13, 13, 12, 13, 12, 13, 13, 13, 13, 11)
  
}

#The CFSeq_Data object stores all studies and associated analysis data, organized by species
View(CFSeq_Data)

#The GenePathway_Data object stores all KEGG pathway info, organized by strain
View(GenePathway_Data)



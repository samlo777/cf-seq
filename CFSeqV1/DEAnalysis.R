##Goal: Create differential expression analysis data objects for each species that can be displayed visually
## (as tables, volcano plots, MA plots) in app and used to filter studies on the basis of individual genes 
## being differentially expressed


#1. Perform DE Analysis By Species (Stored in CFSeq_Data)
{
#---
#Aspergillus fumigatus
{
Species <- 1
AFStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  AFStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
} #Create Study List  
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
#Bacteroides Species [No Studies Yet] 

#--- 
#Burkholderia Species
{
Species <- 2
BurkholderiaStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  BurkholderiaStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
} 
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
Species <- 3
CAStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  CAStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
Species <- 4
CDStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  CDStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
#Enterococcus sp. [No Studies Yet]

#--- 
#Fusobacterium nucleatum [No Studies Yet]

#--- 
#faecalibacterium prausnitzii 
{
Species <- 5
FPStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  FPStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
#haemophilus influenzae
{
Species <- 6
HIStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  HIStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
Species <- 7
MAStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  MAStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
#Pseudomonas Aeruginosa
{
Species <- 8
PAStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  PAStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
#Prevotella Intermedia [No Studies Yet]

#--- 
#Porphyromonas sp. 
{
Species <- 9
PorphyromonasStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  PorphyromonasStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
Species <- 10
SAStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  SAStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
Species <- 11
SMStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  SMStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
Species <- 12
StreptococcusStudyList <- vector()
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  StreptococcusStudyList[i] <- paste(CFSeq_Data[[Species]][[i]][[3]]$GEO.Accession, CFSeq_Data[[Species]][[i]]$Count_Table[,1][1], sep = " | Ex. ")
}
for (i in 1:(length(CFSeq_Data[[Species]]))) {
  
  if (colnames(CFSeq_Data[[Species]][[i]]$Design_Matrix)[1] == "Replicate") {
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
#Veillonella Parvula [No Studies Yet]
}
#Analysis outputs for each experimental comparison are loaded into the CFSeq_Data object at the level of individual studies
#and used as the basis for the DE analysis data tables, volcano plots, and MA plots displayed in app for each individual 
#study and comparison

#2.Create DE Table for each species (Stored in CFSeq_Genes)
{
  
  CFSeq_Genes <- vector(mode = "list")
  
  #Set up Data Frames for each species
  {
    AFDF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 1
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          AFDF <- rbind(AFDF, Rows)
        }
      }
    }
    
    BurkDF <- data.frame(Study = character(),
                         Comparison = character(),
                         Gene = character(),
                         logFC = double(),
                         logCPM = double(),
                         FValue = double(),
                         PValue = double(),
                         stringsAsFactors = FALSE)
    
    Species = 2
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          BurkDF <- rbind(BurkDF, Rows)
        }
      }
    }
    
    CADF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 3
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          CADF <- rbind(CADF, Rows)
        }
      }
    }
    
    CDDF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 4
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          CDDF <- rbind(CDDF, Rows)
        }
      }
    }
    
    FPDF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 5
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          FPDF <- rbind(FPDF, Rows)
        }
      }
    }
    
    HIDF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 6
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          HIDF <- rbind(HIDF, Rows)
        }
      }
    }
    
    MADF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 7
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          MADF <- rbind(MADF, Rows)
        }
      }
    }
    
    PADF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 8
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          PADF <- rbind(PADF, Rows)
        }
      }
    }
    
    PorphDF <- data.frame(Study = character(),
                          Comparison = character(),
                          Gene = character(),
                          logFC = double(),
                          logCPM = double(),
                          FValue = double(),
                          PValue = double(),
                          stringsAsFactors = FALSE)
    
    Species = 9
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          PorphDF <- rbind(PorphDF, Rows)
        }
      }
    }
    
    SADF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 10
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          SADF <- rbind(SADF, Rows)
        }
      }
    }
    
    SMDF <- data.frame(Study = character(),
                       Comparison = character(),
                       Gene = character(),
                       logFC = double(),
                       logCPM = double(),
                       FValue = double(),
                       PValue = double(),
                       stringsAsFactors = FALSE)
    
    Species = 11
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          SMDF <- rbind(SMDF, Rows)
        }
      }
    }
    
    StrepDF <- data.frame(Study = character(),
                          Comparison = character(),
                          Gene = character(),
                          logFC = double(),
                          logCPM = double(),
                          FValue = double(),
                          PValue = double(),
                          stringsAsFactors = FALSE)
    
    Species = 12
    for (i in 1:(length(CFSeq_Data[[Species]]))) {
      if (length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))]) != 0) {
        DEOutputs <- CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))] #Extract out all DE outputs
        for (j in 1:length(CFSeq_Data[[Species]][[i]][grep("DE Output", names(CFSeq_Data[[Species]][[i]]))])) {
          GEO <- rep(CFSeq_Data[[Species]][[i]]$Additional_Metadata$GEO.Accession, times = nrow(DEOutputs[[j]]))
          Comp <- rep(names(DEOutputs)[j], times = nrow(DEOutputs[[j]]))
          Rows <- data.frame(GEO, Comp, rownames(DEOutputs[[j]]), DEOutputs[[j]]$logFC, DEOutputs[[j]]$logCPM, DEOutputs[[j]]$F, DEOutputs[[j]]$PValue)
          colnames(Rows) <- c("Study","Comparison","Gene","logFC","logCPM","FValue","PValue")
          StrepDF <- rbind(StrepDF, Rows)
        }
      }
    }
    
  }
  
  
  #Create DE Gene Tables
  
  CFSeq_Genes[[1]] <- AFDF
  names(CFSeq_Genes)[1] <- "Aspergillus_fumigatus"
  CFSeq_Genes[[2]] <- BurkDF
  names(CFSeq_Genes)[2] <- "Burkholderia_species"
  CFSeq_Genes[[3]] <- CADF
  names(CFSeq_Genes)[3] <- "Candida_albicans"
  CFSeq_Genes[[4]] <- CDDF
  names(CFSeq_Genes)[4] <- "Clostridium_difficile"
  CFSeq_Genes[[5]] <- FPDF
  names(CFSeq_Genes)[5] <- "Faecalibacterium_prausnitzii"
  CFSeq_Genes[[6]] <- HIDF
  names(CFSeq_Genes)[6] <- "Haemophilus_influenzae"
  CFSeq_Genes[[7]] <- MADF
  names(CFSeq_Genes)[7] <- "Mycobacterium_abscessus"
  CFSeq_Genes[[8]] <- PADF
  names(CFSeq_Genes)[8] <- "Pseudomonas_aeruginosa"
  CFSeq_Genes[[9]] <- PorphDF
  names(CFSeq_Genes)[9] <- "Porphyromonas_species"
  CFSeq_Genes[[10]] <- SADF
  names(CFSeq_Genes)[10] <- "Staphylococcus_aureus"
  CFSeq_Genes[[11]] <- SMDF
  names(CFSeq_Genes)[11] <- "Stenotrophomonas_maltophilia"
  CFSeq_Genes[[12]] <- StrepDF
  names(CFSeq_Genes)[12] <- "Streptococcus_species"
  
}
View(CFSeq_Genes)


#The DE Tables for each species are large data tables contianing all DE genes for each comparison in one place. They are loaded
#into the CFSeq_Data object at the level of individual species. They are used for the purpose of filtering studies (in the study
#view tab) on the basis of individual genes being differentially expressed



#Goal ___
#Produce several large list-of-list objects that can be accessed by the app. These lists will contain all count tables, 
#study metadata, differential expression analysis data, and gene ID / KEGG pathway info for each study incorporated in CF-Seq. 

#Data Objects ___
#A. CFSeq_Data: Contains all study count tables, metadata, and differential expression analysis data
View(CFSeq_Data)

names(CFSeq_Data[1])
#1. Install necessary packages
{
BiocManager::install("edgeR")
library(edgeR)

#BiocManager::install("GEOquery")
library(GEOquery)
  
#install.packages("stringr")
library(stringr)
}

#2. Create CFSeq_Data object 
{
CFSeq_Data <- vector(mode = "list", 12)

names(CFSeq_Data) <- c("Aspergillus_fumigatus", "Burkholderia_Species","Candida_albicans","Clostridium_difficile", 
                       "Faecalibacterium_prausnitzii","Haemophilus_influenzae","Mycobacterium_abscessus", "Pseudomonas_aeruginosa",
                       "Porphyromonas_Species","Staphylococcus_aureus","Stenotrophomonas_maltophilia","Streptococcus_Species")
}
#Stores all studies, organized by species

#3. Load in all count tables, design matrices, and metadata for each species 
{

#---
#Achromobacter Xylosoxidans [No Studies Yet]
{#--- 

#Add slots for studies

#Check that all are formatted properly 

#Load in count tables

#Now load in design matrices 

#Now load in additional metadata
}

#---
#Aspergillus Fumigatus 
{#Add slots for studies
CFSeq_Data[[1]] = vector(mode = "list", 5)
names(CFSeq_Data[[1]]) <- c("GSE122391","GSE152682","GSE173349","GSE55648","GSE55943")

for (i in 1:length(CFSeq_Data[[1]])) {
  CFSeq_Data[[1]][[i]] = vector(mode = "list", 3)
  names(CFSeq_Data[[1]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
}

#Now load in count tables 
CFSeq_Data[[1]][[1]]$Count_Table <- read.csv("./CountTables/AFCounts/Count_Tables/GSE122391_counts.csv")
CFSeq_Data[[1]][[2]]$Count_Table <- read.csv("./CountTables/AFCounts/Count_Tables/GSE152682_counts.csv")
CFSeq_Data[[1]][[3]]$Count_Table <- read.csv("./CountTables/AFCounts/Count_Tables/GSE173349_counts.csv")
CFSeq_Data[[1]][[4]]$Count_Table <- read.csv("./CountTables/AFCounts/Count_Tables/GSE55648_counts.csv")
CFSeq_Data[[1]][[5]]$Count_Table <- read.csv("./CountTables/AFCounts/Count_Tables/GSE55943_counts.csv")

#Now load in design matrices 
CFSeq_Data[[1]][[1]]$Design_Matrix <- read.csv("./CountTables/AFCounts/Design_Matrices/GSE122391_design.csv")
CFSeq_Data[[1]][[2]]$Design_Matrix <- read.csv("./CountTables/AFCounts/Design_Matrices/GSE152682_design.csv")
CFSeq_Data[[1]][[3]]$Design_Matrix <- read.csv("./CountTables/AFCounts/Design_Matrices/GSE173349_design.csv")
CFSeq_Data[[1]][[4]]$Design_Matrix <- read.csv("./CountTables/AFCounts/Design_Matrices/GSE55648_design.csv")
CFSeq_Data[[1]][[5]]$Design_Matrix <- read.csv("./CountTables/AFCounts/Design_Matrices/GSE55943_design.csv")

#Now load in additional metadata
CFSeq_Data[[1]][[1]]$Additional_Metadata <- read.csv("./CountTables/AFCounts/Additional_Metadata/GSE122391_metadata.csv")
CFSeq_Data[[1]][[2]]$Additional_Metadata <- read.csv("./CountTables/AFCounts/Additional_Metadata/GSE152682_metadata.csv")
CFSeq_Data[[1]][[3]]$Additional_Metadata <- read.csv("./CountTables/AFCounts/Additional_Metadata/GSE173349_metadata.csv")
CFSeq_Data[[1]][[4]]$Additional_Metadata <- read.csv("./CountTables/AFCounts/Additional_Metadata/GSE55648_metadata.csv")
CFSeq_Data[[1]][[5]]$Additional_Metadata <- read.csv("./CountTables/AFCounts/Additional_Metadata/GSE55943_metadata.csv")

}

#---
#Bacteroides Species [No Studies Yet] 
{#--- 
#Bacteroides Species [No Studies Yet]

#Add slots for studies

#Check that all are formatted properly 

#Load in count tables

#Now load in design matrices 

#Now load in additional metadata
}

#--- 
#Burkholderia Species
{#Add slots for studies
CFSeq_Data[[2]] = vector(mode = "list", 2)
names(CFSeq_Data[[2]]) <- c("GSE19115","GSE77970")

for (i in 1:length(CFSeq_Data[[2]])) {
  CFSeq_Data[[2]][[i]] = vector(mode = "list", 3)
  names(CFSeq_Data[[2]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
}

#Load in count tables
CFSeq_Data[[2]][[1]]$Count_Table <- read.csv("./CountTables/BurkholderiaCounts/Count_Tables/GSE19115_counts.csv")
CFSeq_Data[[2]][[2]]$Count_Table <- read.csv("./CountTables/BurkholderiaCounts/Count_Tables/GSE77970_counts.csv")

#Now load in design matrices 
CFSeq_Data[[2]][[1]]$Design_Matrix <- read.csv("./CountTables/BurkholderiaCounts/Design_Matrices/GSE19115_design.csv")
CFSeq_Data[[2]][[2]]$Design_Matrix <- read.csv("./CountTables/BurkholderiaCounts/Design_Matrices/GSE77970_design.csv")

#Now load in additional metadata
CFSeq_Data[[2]][[1]]$Additional_Metadata <- read.csv("./CountTables/BurkholderiaCounts/Additional_Metadata/GSE19115_metadata.csv")
CFSeq_Data[[2]][[2]]$Additional_Metadata <- read.csv("./CountTables/BurkholderiaCounts/Additional_Metadata/GSE77970_metadata.csv")

}

#--- 
#Candida Albicans 
{#Add slots for studies
CFSeq_Data[[3]] = vector(mode = "list", 15)
names(CFSeq_Data[[3]]) <- c("GSE100737", "GSE102039", "GSE116533", "GSE123122", "GSE124137", "GSE125636", "GSE130948",
                            "GSE133611", "GSE138069", "GSE154488", "GSE158472",
                            "GSE173668", "GSE56174", "GSE86540", "GSE96965")

for (i in 1:length(CFSeq_Data[[3]])) {
  CFSeq_Data[[3]][[i]] = vector(mode = "list", 3)
  names(CFSeq_Data[[3]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
}

#Load in count tables
CFSeq_Data[[3]][[1]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE100737_counts.csv")
CFSeq_Data[[3]][[2]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE102039_counts.csv")
CFSeq_Data[[3]][[3]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE116533_counts.csv")
CFSeq_Data[[3]][[4]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE123122_counts.csv")
CFSeq_Data[[3]][[5]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE124137_counts.csv")
CFSeq_Data[[3]][[6]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE125636_counts.csv")
CFSeq_Data[[3]][[7]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE130948_counts.csv")
CFSeq_Data[[3]][[8]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE133611_counts.csv")
CFSeq_Data[[3]][[9]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE138069_counts.csv")
CFSeq_Data[[3]][[10]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE154488_counts.csv")
CFSeq_Data[[3]][[11]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE158472_counts.csv")
CFSeq_Data[[3]][[12]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE173668_counts.csv")
CFSeq_Data[[3]][[13]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE56174_counts.csv")
CFSeq_Data[[3]][[14]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE86540_counts.csv")
CFSeq_Data[[3]][[15]]$Count_Table <- read.csv("./CountTables/CACounts/Count_Tables/GSE96965_counts.csv")

#Now load in design matrices 
CFSeq_Data[[3]][[1]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE100737_design.csv")
CFSeq_Data[[3]][[2]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE102039_design.csv")
CFSeq_Data[[3]][[3]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE116533_design.csv")
CFSeq_Data[[3]][[4]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE123122_design.csv")
CFSeq_Data[[3]][[5]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE124137_design.csv")
CFSeq_Data[[3]][[6]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE125636_design.csv")
CFSeq_Data[[3]][[7]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE130948_design.csv")
CFSeq_Data[[3]][[8]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE133611_design.csv")
CFSeq_Data[[3]][[9]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE138069_design.csv")
CFSeq_Data[[3]][[10]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE154488_design.csv")
CFSeq_Data[[3]][[11]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE158472_design.csv")
CFSeq_Data[[3]][[12]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE173668_design.csv")
CFSeq_Data[[3]][[13]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE56174_design.csv")
CFSeq_Data[[3]][[14]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE86540_design.csv")
CFSeq_Data[[3]][[15]]$Design_Matrix <- read.csv("./CountTables/CACounts/Design_Matrices/GSE96965_design.csv")

#Now load in additional metadata
CFSeq_Data[[3]][[1]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE100737_metadata.csv")
CFSeq_Data[[3]][[2]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE102039_metadata.csv")
CFSeq_Data[[3]][[3]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE116533_metadata.csv")
CFSeq_Data[[3]][[4]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE123122_metadata.csv")
CFSeq_Data[[3]][[5]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE124137_metadata.csv")
CFSeq_Data[[3]][[6]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE125636_metadata.csv")
CFSeq_Data[[3]][[7]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE130948_metadata.csv")
CFSeq_Data[[3]][[8]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE133611_metadata.csv")
CFSeq_Data[[3]][[9]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE138069_metadata.csv")
CFSeq_Data[[3]][[10]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE154488_metadata.csv")
CFSeq_Data[[3]][[11]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE158472_metadata.csv")
CFSeq_Data[[3]][[12]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE173668_metadata.csv")
CFSeq_Data[[3]][[13]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE56174_metadata.csv")
CFSeq_Data[[3]][[14]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE86540_metadata.csv")
CFSeq_Data[[3]][[15]]$Additional_Metadata <- read.csv("./CountTables/CACounts/Additional_Metadata/GSE96965_metadata.csv")

}

#--- 
#Clostridioides difficile
{#Add slots for studies
CFSeq_Data[[4]] = vector(mode = "list", 6)
names(CFSeq_Data[[4]]) <- c("GSE86152","GSE86612","GSE103952","GSE107961","GSE120198","GSE135912")

for (i in 1:length(CFSeq_Data[[4]])) {
  CFSeq_Data[[4]][[i]] = vector(mode = "list", 3)
  names(CFSeq_Data[[4]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
}

#Load in count tables
CFSeq_Data[[4]][[1]]$Count_Table <- read.csv("./CountTables/CDCounts/Count_Tables/GSE86152_counts.csv")
CFSeq_Data[[4]][[2]]$Count_Table <- read.csv("./CountTables/CDCounts/Count_Tables/GSE86612_counts.csv")
CFSeq_Data[[4]][[3]]$Count_Table <- read.csv("./CountTables/CDCounts/Count_Tables/GSE103952_counts.csv")
CFSeq_Data[[4]][[4]]$Count_Table <- read.csv("./CountTables/CDCounts/Count_Tables/GSE107961_counts.csv")
CFSeq_Data[[4]][[5]]$Count_Table <- read.csv("./CountTables/CDCounts/Count_Tables/GSE120198_counts.csv")
CFSeq_Data[[4]][[6]]$Count_Table <- read.csv("./CountTables/CDCounts/Count_Tables/GSE135912_counts.csv")

#Now load in design matrices 
CFSeq_Data[[4]][[1]]$Design_Matrix <- read.csv("./CountTables/CDCounts/Design_Matrices/GSE86152_design.csv")
CFSeq_Data[[4]][[2]]$Design_Matrix <- read.csv("./CountTables/CDCounts/Design_Matrices/GSE86612_design.csv")
CFSeq_Data[[4]][[3]]$Design_Matrix <- read.csv("./CountTables/CDCounts/Design_Matrices/GSE103952_design.csv")
CFSeq_Data[[4]][[4]]$Design_Matrix <- read.csv("./CountTables/CDCounts/Design_Matrices/GSE107961_design.csv")
CFSeq_Data[[4]][[5]]$Design_Matrix <- read.csv("./CountTables/CDCounts/Design_Matrices/GSE120198_design.csv")
CFSeq_Data[[4]][[6]]$Design_Matrix <- read.csv("./CountTables/CDCounts/Design_Matrices/GSE135912_design.csv")

#Now load in additional metadata
CFSeq_Data[[4]][[1]]$Additional_Metadata <- read.csv("./CountTables/CDCounts/Additional_Metadata/GSE86152_metadata.csv")
CFSeq_Data[[4]][[2]]$Additional_Metadata <- read.csv("./CountTables/CDCounts/Additional_Metadata/GSE86612_metadata.csv")
CFSeq_Data[[4]][[3]]$Additional_Metadata <- read.csv("./CountTables/CDCounts/Additional_Metadata/GSE103952_metadata.csv")
CFSeq_Data[[4]][[4]]$Additional_Metadata <- read.csv("./CountTables/CDCounts/Additional_Metadata/GSE107961_metadata.csv")
CFSeq_Data[[4]][[5]]$Additional_Metadata <- read.csv("./CountTables/CDCounts/Additional_Metadata/GSE120198_metadata.csv")
CFSeq_Data[[4]][[6]]$Additional_Metadata <- read.csv("./CountTables/CDCounts/Additional_Metadata/GSE135912_metadata.csv")

}

#--- 
#Enterococcus sp. [No Studies Yet]
{#Add slots for studies

#Check that all are formatted properly 

#Load in count tables

#Now load in design matrices 

#Now load in additional metadata
}

#--- 
#Fusobacterium nucleatum [No Studies Yet]
{#Add slots for studies

#Check that all are formatted properly 

#Load in count tables

#Now load in design matrices 

#Now load in additional metadata
}

#--- 
#faecalibacterium prausnitzii 
{#Add slots for studies
  CFSeq_Data[[5]] = vector(mode = "list", 1)
  names(CFSeq_Data[[5]]) <- c("GSE168352")
  
  for (i in 1:length(CFSeq_Data[[5]])) {
    CFSeq_Data[[5]][[i]] = vector(mode = "list", 3)
    names(CFSeq_Data[[5]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
  }
  
#Load in count tables
CFSeq_Data[[5]][[1]]$Count_Table <- read.csv("./CountTables/FPCounts/Count_Tables/GSE168352_counts.csv")

#Now load in design matrices 
CFSeq_Data[[5]][[1]]$Design_Matrix <- read.csv("./CountTables/FPCounts/Design_Matrices/GSE168352_design.csv")

#Now load in additional metadata
CFSeq_Data[[5]][[1]]$Additional_Metadata <- read.csv("./CountTables/FPCounts/Additional_Metadata/GSE168352_metadata.csv")

}
  
#--- 
#haemophilus influenzae
{#Add slots for studies
  CFSeq_Data[[6]] = vector(mode = "list", 1)
  names(CFSeq_Data[[6]]) <- c("GSE125415")
  
  for (i in 1:length(CFSeq_Data[[6]])) {
    CFSeq_Data[[6]][[i]] = vector(mode = "list", 3)
    names(CFSeq_Data[[6]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
  }
  
#Load in count tables
  CFSeq_Data[[6]][[1]]$Count_Table <- read.csv("./CountTables/HICounts/Count_Tables/GSE125415_counts.csv")
  
#Now load in design matrices 
  CFSeq_Data[[6]][[1]]$Design_Matrix <- read.csv("./CountTables/HICounts/Design_Matrices/GSE125415_design.csv")

#Now load in additional metadata
  CFSeq_Data[[6]][[1]]$Additional_Metadata <- read.csv("./CountTables/HICounts/Additional_Metadata/GSE125415_metadata.csv")
  
}

#--- 
#Mycobacterium Abscessus
{#Add slots for studies
  CFSeq_Data[[7]] = vector(mode = "list", 1)
  names(CFSeq_Data[[7]]) <- c("GSE78787")
  
  for (i in 1:length(CFSeq_Data[[7]])) {
    CFSeq_Data[[7]][[i]] = vector(mode = "list", 3)
    names(CFSeq_Data[[7]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
  }
  
#Load in count tables
  CFSeq_Data[[7]][[1]]$Count_Table <- read.csv("./CountTables/MACounts/Count_Tables/GSE78787_counts.csv")
  
#Now load in design matrices 
  CFSeq_Data[[7]][[1]]$Design_Matrix <- read.csv("./CountTables/MACounts/Design_Matrices/GSE78787_design.csv")
  
#Now load in additional metadata
  CFSeq_Data[[7]][[1]]$Additional_Metadata <- read.csv("./CountTables/MACounts/Additional_Metadata/GSE78787_metadata.csv")

}

#--- 
#Pseudomonas Aeruginosa
{#Add slots for studies
  CFSeq_Data[[8]] = vector(mode = "list", 29)
  names(CFSeq_Data[[8]]) <- c("GSE81065","GSE86211","GSE87213","GSE99729","GSE99981","GSE103620","GSE118801",
                              "GSE121243","GSE122048","GSE123356","GSE124385",
                              "GSE125646","GSE130190","GSE136111","GSE138731","GSE139104",
                              "GSE142464","GSE144365","GSE148116","GSE148597","GSE148955",
                              "GSE152480","GSE153067","GSE156995","GSE163234","GSE163248","GSE163555",
                              "GSE166602","GSE181354")
  
  for (i in 1:length(CFSeq_Data[[8]])) {
    CFSeq_Data[[8]][[i]] = vector(mode = "list", 3)
    names(CFSeq_Data[[8]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
  }
  
#Load in count tables
  CFSeq_Data[[8]][[1]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE81065_counts.csv")
  CFSeq_Data[[8]][[2]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE86211_counts.csv")
  CFSeq_Data[[8]][[3]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE87213_counts.csv")
  CFSeq_Data[[8]][[4]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE99729_counts.csv") 
  CFSeq_Data[[8]][[5]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE99981_counts.csv")
  CFSeq_Data[[8]][[6]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE103620_counts.csv")
  CFSeq_Data[[8]][[7]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE118801_counts.csv")
  CFSeq_Data[[8]][[8]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE121243_counts.csv")
  CFSeq_Data[[8]][[9]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE122048_counts.csv")
  CFSeq_Data[[8]][[10]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE123356_counts.csv")
  CFSeq_Data[[8]][[11]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE124385_counts.csv")
  CFSeq_Data[[8]][[12]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE125646_counts.csv")
  CFSeq_Data[[8]][[13]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE130190_counts.csv")
  CFSeq_Data[[8]][[14]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE136111_counts.csv")
  CFSeq_Data[[8]][[15]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE138731_counts.csv")
  CFSeq_Data[[8]][[16]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE139104_counts.csv")
  CFSeq_Data[[8]][[17]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE142464_counts.csv")
  CFSeq_Data[[8]][[18]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE144365_counts.csv")
  CFSeq_Data[[8]][[19]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE148116_counts.csv")
  CFSeq_Data[[8]][[20]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE148597_counts.csv")
  CFSeq_Data[[8]][[21]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE148955_counts.csv")
  CFSeq_Data[[8]][[22]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE152480_counts.csv")
  CFSeq_Data[[8]][[23]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE153067_counts.csv")
  CFSeq_Data[[8]][[24]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE156995_counts.csv")
  CFSeq_Data[[8]][[25]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE163234_counts.csv")
  CFSeq_Data[[8]][[26]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE163248_counts.csv")
  CFSeq_Data[[8]][[27]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE163555_counts.csv")
  CFSeq_Data[[8]][[28]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE166602_counts.csv")
  CFSeq_Data[[8]][[29]]$Count_Table <- read.csv("./CountTables/PACounts/Count_Tables/GSE181354_counts.csv")
  
#Now load in design matrices 
  CFSeq_Data[[8]][[1]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE81065_design.csv")
  CFSeq_Data[[8]][[2]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE86211_design.csv")
  CFSeq_Data[[8]][[3]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE87213_design.csv")
  CFSeq_Data[[8]][[4]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE99729_design.csv")
  CFSeq_Data[[8]][[5]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE99981_design.csv")
  CFSeq_Data[[8]][[6]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE103620_design.csv")
  CFSeq_Data[[8]][[7]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE118801_design.csv")
  CFSeq_Data[[8]][[8]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE121243_design.csv")
  CFSeq_Data[[8]][[9]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE122048_design.csv")
  CFSeq_Data[[8]][[10]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE123356_design.csv")
  CFSeq_Data[[8]][[11]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE124385_design.csv")
  CFSeq_Data[[8]][[12]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE125646_design.csv")
  CFSeq_Data[[8]][[13]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE130190_design.csv")
  CFSeq_Data[[8]][[14]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE136111_design.csv")
  CFSeq_Data[[8]][[15]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE138731_design.csv")
  CFSeq_Data[[8]][[16]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE139104_design.csv")
  CFSeq_Data[[8]][[17]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE142464_design.csv")
  CFSeq_Data[[8]][[18]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE144365_design.csv")
  CFSeq_Data[[8]][[19]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE148116_design.csv")
  CFSeq_Data[[8]][[20]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE148597_design.csv")
  CFSeq_Data[[8]][[21]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE148955_design.csv")
  CFSeq_Data[[8]][[22]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE152480_design.csv")
  CFSeq_Data[[8]][[23]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE153067_design.csv")
  CFSeq_Data[[8]][[24]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE156995_design.csv")
  CFSeq_Data[[8]][[25]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE163234_design.csv")
  CFSeq_Data[[8]][[26]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE163248_design.csv")
  CFSeq_Data[[8]][[27]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE163555_design.csv")
  CFSeq_Data[[8]][[28]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE166602_design.csv")
  CFSeq_Data[[8]][[29]]$Design_Matrix <- read.csv("./CountTables/PACounts/Design_Matrices/GSE181354_design.csv")
  
#Now load in additional metadata
  CFSeq_Data[[8]][[1]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE81065_metadata.csv")
  CFSeq_Data[[8]][[2]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE86211_metadata.csv")
  CFSeq_Data[[8]][[3]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE87213_metadata.csv")
  CFSeq_Data[[8]][[4]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE99729_metadata.csv")
  CFSeq_Data[[8]][[5]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE99981_metadata.csv")
  CFSeq_Data[[8]][[6]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE103620_metadata.csv")
  CFSeq_Data[[8]][[7]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE118801_metadata.csv")
  CFSeq_Data[[8]][[8]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE121243_metadata.csv")
  CFSeq_Data[[8]][[9]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE122048_metadata.csv")
  CFSeq_Data[[8]][[10]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE123356_metadata.csv")
  CFSeq_Data[[8]][[11]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE124385_metadata.csv")
  CFSeq_Data[[8]][[12]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE125646_metadata.csv")
  CFSeq_Data[[8]][[13]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE130190_metadata.csv")
  CFSeq_Data[[8]][[14]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE136111_metadata.csv")
  CFSeq_Data[[8]][[15]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE138731_metadata.csv")
  CFSeq_Data[[8]][[16]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE139104_metadata.csv")
  CFSeq_Data[[8]][[17]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE142464_metadata.csv")
  CFSeq_Data[[8]][[18]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE144365_metadata.csv")
  CFSeq_Data[[8]][[19]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE148116_metadata.csv")
  CFSeq_Data[[8]][[20]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE148597_metadata.csv")
  CFSeq_Data[[8]][[21]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE148955_metadata.csv")
  CFSeq_Data[[8]][[22]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE152480_metadata.csv")
  CFSeq_Data[[8]][[23]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE153067_metadata.csv")
  CFSeq_Data[[8]][[24]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE156995_metadata.csv")
  CFSeq_Data[[8]][[25]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE163234_metadata.csv")
  CFSeq_Data[[8]][[26]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE163248_metadata.csv")
  CFSeq_Data[[8]][[27]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE163555_metadata.csv")
  CFSeq_Data[[8]][[28]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE166602_metadata.csv")
  CFSeq_Data[[8]][[29]]$Additional_Metadata <- read.csv("./CountTables/PACounts/Additional_Metadata/GSE181354_metadata.csv")
  
}
  
#--- 
#Prevotella Intermedia [No Studies Yet]
{#Add slots for studies

#Check that all are formatted properly 

#Load in count tables

#Now load in design matrices 

#Now load in additional metadata
}
  
#--- 
#Porphyromonas sp. 
{#Add slots for studies
  CFSeq_Data[[9]] = vector(mode = "list", 1)
  names(CFSeq_Data[[9]]) <- c("GSE124206")
  
  for (i in 1:length(CFSeq_Data[[9]])) {
    CFSeq_Data[[9]][[i]] = vector(mode = "list", 3)
    names(CFSeq_Data[[9]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
  }

#Load in count tables
  CFSeq_Data[[9]][[1]]$Count_Table <- read.csv("./CountTables/PorphyromonasCounts/Count_Tables/GSE124206_counts.csv")

#Now load in design matrices 
  CFSeq_Data[[9]][[1]]$Design_Matrix <- read.csv("./CountTables/PorphyromonasCounts/Design_Matrices/GSE124206_design.csv")
  
#Now load in additional metadata
  CFSeq_Data[[9]][[1]]$Additional_Metadata <- read.csv("./CountTables/PorphyromonasCounts/Additional_Metadata/GSE124206_metadata.csv")

}
  
#--- 
#Staphylococcus aureus
{#Add slots for studies
  CFSeq_Data[[10]] = vector(mode = "list", 8)
  names(CFSeq_Data[[10]]) <- c("GSE40864","GSE77473","GSE79407",
                               "GSE122048","GSE122065","GSE125741","GSE130777",
                               "GSE139659")
  
  for (i in 1:length(CFSeq_Data[[10]])) {
    CFSeq_Data[[10]][[i]] = vector(mode = "list", 3)
    names(CFSeq_Data[[10]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
  }
  
#Load in count tables
  CFSeq_Data[[10]][[1]]$Count_Table <- read.csv("./CountTables/SACounts/Count_Tables/GSE40864_counts.csv")
  CFSeq_Data[[10]][[2]]$Count_Table <- read.csv("./CountTables/SACounts/Count_Tables/GSE77473_counts.csv")
  CFSeq_Data[[10]][[3]]$Count_Table <- read.csv("./CountTables/SACounts/Count_Tables/GSE79407_counts.csv")
  CFSeq_Data[[10]][[4]]$Count_Table <- read.csv("./CountTables/SACounts/Count_Tables/GSE122048_counts.csv")
  CFSeq_Data[[10]][[5]]$Count_Table <- read.csv("./CountTables/SACounts/Count_Tables/GSE122065_counts.csv")
  CFSeq_Data[[10]][[6]]$Count_Table <- read.csv("./CountTables/SACounts/Count_Tables/GSE125741_counts.csv")
  CFSeq_Data[[10]][[7]]$Count_Table <- read.csv("./CountTables/SACounts/Count_Tables/GSE130777_counts.csv")
  CFSeq_Data[[10]][[8]]$Count_Table <- read.csv("./CountTables/SACounts/Count_Tables/GSE139659_counts.csv")

#Now load in design matrices 
  CFSeq_Data[[10]][[1]]$Design_Matrix <- read.csv("./CountTables/SACounts/Design_Matrices/GSE40864_design.csv")
  CFSeq_Data[[10]][[2]]$Design_Matrix <- read.csv("./CountTables/SACounts/Design_Matrices/GSE77473_design.csv")
  CFSeq_Data[[10]][[3]]$Design_Matrix <- read.csv("./CountTables/SACounts/Design_Matrices/GSE79407_design.csv")
  CFSeq_Data[[10]][[4]]$Design_Matrix <- read.csv("./CountTables/SACounts/Design_Matrices/GSE122048_design.csv")
  CFSeq_Data[[10]][[5]]$Design_Matrix <- read.csv("./CountTables/SACounts/Design_Matrices/GSE122065_design.csv")
  CFSeq_Data[[10]][[6]]$Design_Matrix <- read.csv("./CountTables/SACounts/Design_Matrices/GSE125741_design.csv")
  CFSeq_Data[[10]][[7]]$Design_Matrix <- read.csv("./CountTables/SACounts/Design_Matrices/GSE130777_design.csv")
  CFSeq_Data[[10]][[8]]$Design_Matrix <- read.csv("./CountTables/SACounts/Design_Matrices/GSE139659_design.csv")
  
#Now load in additional metadata
  CFSeq_Data[[10]][[1]]$Additional_Metadata <- read.csv("./CountTables/SACounts/Additional_Metadata/GSE40864_metadata.csv")
  CFSeq_Data[[10]][[2]]$Additional_Metadata <- read.csv("./CountTables/SACounts/Additional_Metadata/GSE77473_metadata.csv")
  CFSeq_Data[[10]][[3]]$Additional_Metadata <- read.csv("./CountTables/SACounts/Additional_Metadata/GSE79407_metadata.csv")
  CFSeq_Data[[10]][[4]]$Additional_Metadata <- read.csv("./CountTables/SACounts/Additional_Metadata/GSE122048_metadata.csv")
  CFSeq_Data[[10]][[5]]$Additional_Metadata <- read.csv("./CountTables/SACounts/Additional_Metadata/GSE122065_metadata.csv")
  CFSeq_Data[[10]][[6]]$Additional_Metadata <- read.csv("./CountTables/SACounts/Additional_Metadata/GSE125741_metadata.csv")
  CFSeq_Data[[10]][[7]]$Additional_Metadata <- read.csv("./CountTables/SACounts/Additional_Metadata/GSE130777_metadata.csv")
  CFSeq_Data[[10]][[8]]$Additional_Metadata <- read.csv("./CountTables/SACounts/Additional_Metadata/GSE139659_metadata.csv")
  
}
  
#--- 
#Stenotrophomonas maltophilia
{#Add slots for studies
  CFSeq_Data[[11]] = vector(mode = "list", 3)
  names(CFSeq_Data[[11]]) <- c("GSE121347","GSE125704","GSE141276")
  
  for (i in 1:length(CFSeq_Data[[11]])) {
    CFSeq_Data[[11]][[i]] = vector(mode = "list", 3)
    names(CFSeq_Data[[11]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
  }
  
#Load in count tables
  CFSeq_Data[[11]][[1]]$Count_Table <- read.csv("./CountTables/SMCounts/Count_Tables/GSE121347_counts.csv")
  CFSeq_Data[[11]][[2]]$Count_Table <- read.csv("./CountTables/SMCounts/Count_Tables/GSE125704_counts.csv")
  CFSeq_Data[[11]][[3]]$Count_Table <- read.csv("./CountTables/SMCounts/Count_Tables/GSE141276_counts.csv")
  
#Now load in design matrices 
  CFSeq_Data[[11]][[1]]$Design_Matrix <- read.csv("./CountTables/SMCounts/Design_Matrices/GSE121347_design.csv")
  CFSeq_Data[[11]][[2]]$Design_Matrix <- read.csv("./CountTables/SMCounts/Design_Matrices/GSE125704_design.csv")
  CFSeq_Data[[11]][[3]]$Design_Matrix <- read.csv("./CountTables/SMCounts/Design_Matrices/GSE141276_design.csv")
  
#Now load in additional metadata
  CFSeq_Data[[11]][[1]]$Additional_Metadata <- read.csv("./CountTables/SMCounts/Additional_Metadata/GSE121347_metadata.csv")
  CFSeq_Data[[11]][[2]]$Additional_Metadata <- read.csv("./CountTables/SMCounts/Additional_Metadata/GSE125704_metadata.csv")
  CFSeq_Data[[11]][[3]]$Additional_Metadata <- read.csv("./CountTables/SMCounts/Additional_Metadata/GSE141276_metadata.csv")
  
}
  
#--- 
#Streptococcus sp. 
{#Add slots for studies
  CFSeq_Data[[12]] = vector(mode = "list", 23)
  names(CFSeq_Data[[12]]) <- c("GSE67533","GSE86854","GSE120640","GSE128534",
                               "GSE139093","GSE142362","GSE142458","GSE153766",
                               "GSE158512","GSE160737","GSE164206","GSE167217","GSE167895","GSE167896",
                               "GSE167897","GSE167898","GSE167899","GSE167900","GSE167901","GSE167902",
                               "GSE167903","GSE173362","GSE181516")
  
  for (i in 1:length(CFSeq_Data[[12]])) {
    CFSeq_Data[[12]][[i]] = vector(mode = "list", 3)
    names(CFSeq_Data[[12]][[i]]) = c("Count_Table","Design_Matrix","Additional_Metadata")
  }
  
#Load in count tables
  CFSeq_Data[[12]][[1]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE67533_counts.csv")
  CFSeq_Data[[12]][[2]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE86854_counts.csv")
  CFSeq_Data[[12]][[3]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE120640_counts.csv")
  CFSeq_Data[[12]][[4]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE128534_counts.csv")
  CFSeq_Data[[12]][[5]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE139093_counts.csv")
  CFSeq_Data[[12]][[6]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE142362_counts.csv")
  CFSeq_Data[[12]][[7]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE142458_counts.csv")
  CFSeq_Data[[12]][[8]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE153766_counts.csv")
  CFSeq_Data[[12]][[9]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE158512_counts.csv")
  CFSeq_Data[[12]][[10]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE160737_counts.csv")
  CFSeq_Data[[12]][[11]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE164206_counts.csv")
  CFSeq_Data[[12]][[12]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167217_counts.csv")
  CFSeq_Data[[12]][[13]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167895_counts.csv")
  CFSeq_Data[[12]][[14]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167896_counts.csv")
  CFSeq_Data[[12]][[15]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167897_counts.csv")
  CFSeq_Data[[12]][[16]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167898_counts.csv")
  CFSeq_Data[[12]][[17]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167899_counts.csv")
  CFSeq_Data[[12]][[18]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167900_counts.csv")
  CFSeq_Data[[12]][[19]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167901_counts.csv")
  CFSeq_Data[[12]][[20]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167902_counts.csv")
  CFSeq_Data[[12]][[21]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE167903_counts.csv")
  CFSeq_Data[[12]][[22]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE173362_counts.csv")
  CFSeq_Data[[12]][[23]]$Count_Table <- read.csv("./CountTables/StrepCounts/Count_Tables/GSE181516_counts.csv")
  
#Now load in design matrices 
  CFSeq_Data[[12]][[1]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE67533_design.csv")
  CFSeq_Data[[12]][[2]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE86854_design.csv")
  CFSeq_Data[[12]][[3]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE120640_design.csv")
  CFSeq_Data[[12]][[4]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE128534_design.csv")
  CFSeq_Data[[12]][[5]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE139093_design.csv")
  CFSeq_Data[[12]][[6]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE142362_design.csv")
  CFSeq_Data[[12]][[7]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE142458_design.csv")
  CFSeq_Data[[12]][[8]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE153766_design.csv")
  CFSeq_Data[[12]][[9]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE158512_design.csv")
  CFSeq_Data[[12]][[10]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE160737_design.csv")
  CFSeq_Data[[12]][[11]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE164206_design.csv")
  CFSeq_Data[[12]][[12]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167217_design.csv")
  CFSeq_Data[[12]][[13]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167895_design.csv")
  CFSeq_Data[[12]][[14]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167896_design.csv")
  CFSeq_Data[[12]][[15]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167897_design.csv")
  CFSeq_Data[[12]][[16]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167898_design.csv")
  CFSeq_Data[[12]][[17]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167899_design.csv")
  CFSeq_Data[[12]][[18]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167900_design.csv")
  CFSeq_Data[[12]][[19]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167901_design.csv")
  CFSeq_Data[[12]][[20]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167902_design.csv")
  CFSeq_Data[[12]][[21]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE167903_design.csv")
  CFSeq_Data[[12]][[22]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE173362_design.csv")
  CFSeq_Data[[12]][[23]]$Design_Matrix <- read.csv("./CountTables/StrepCounts/Design_Matrices/GSE181516_design.csv")
  
#Now load in additional metadata
  CFSeq_Data[[12]][[1]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE67533_metadata.csv")
  CFSeq_Data[[12]][[2]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE86854_metadata.csv")
  CFSeq_Data[[12]][[3]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE120640_metadata.csv")
  CFSeq_Data[[12]][[4]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE128534_metadata.csv")
  CFSeq_Data[[12]][[5]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE139093_metadata.csv")
  CFSeq_Data[[12]][[6]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE142362_metadata.csv")
  CFSeq_Data[[12]][[7]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE142458_metadata.csv")
  CFSeq_Data[[12]][[8]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE153766_metadata.csv")
  CFSeq_Data[[12]][[9]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE158512_metadata.csv")
  CFSeq_Data[[12]][[10]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE160737_metadata.csv")
  CFSeq_Data[[12]][[11]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE164206_metadata.csv")
  CFSeq_Data[[12]][[12]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167217_metadata.csv")
  CFSeq_Data[[12]][[13]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167895_metadata.csv")
  CFSeq_Data[[12]][[14]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167896_metadata.csv")
  CFSeq_Data[[12]][[15]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167897_metadata.csv")
  CFSeq_Data[[12]][[16]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167898_metadata.csv")
  CFSeq_Data[[12]][[17]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167899_metadata.csv")
  CFSeq_Data[[12]][[18]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167900_metadata.csv")
  CFSeq_Data[[12]][[19]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167901_metadata.csv")
  CFSeq_Data[[12]][[20]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167902_metadata.csv")
  CFSeq_Data[[12]][[21]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE167903_metadata.csv")
  CFSeq_Data[[12]][[22]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE173362_metadata.csv")
  CFSeq_Data[[12]][[23]]$Additional_Metadata <- read.csv("./CountTables/StrepCounts/Additional_Metadata/GSE181516_metadata.csv")
  
}
  
#--- 
#Veillonella Parvula [No Studies Yet]
{#Add slots for studies

#Check that all are formatted properly 

#Load in count tables

#Now load in design matrices 

#Now load in additional metadata
}

}

#B. KEGG_Data: Contains KEGG pathway info for compatible bacterial strains. Currently, this includes: PA14, ...

#KEGG Data Extraction (with the R/Bioconductor package KEGGREST)
{
 
  #PA14
  {
  PA14List <- keggList("pau") #pau = UCBPP-PA14

  #Download all KEGG pathway info by gene ID
  {
  PA14_1000 <- vector(mode = "list")
  for (i in 1:1000) {
    for (j in 1:length(keggGet(names(PA14List)[i])[[1]]$PATHWAY)) {
      if (is.null(keggGet(names(PA14List)[i])[[1]]$PATHWAY) != TRUE) {
        PA14_1000 <- c(PA14_1000, keggGet(names(PA14List)[i])[[1]]$PATHWAY[[j]])
        names(PA14_1000)[length(PA14_1000)] <- names(PA14List)[i]
        print(PA14_1000)[length(PA14_1000)]
      }
    }
  } #Done
  
  for (i in 1001:2000) {
    for (j in 1:length(keggGet(names(PA14List)[i])[[1]]$PATHWAY)) {
      if (is.null(keggGet(names(PA14List)[i])[[1]]$PATHWAY) != TRUE) {
        PA14_1000 <- c(PA14_1000, keggGet(names(PA14List)[i])[[1]]$PATHWAY[[j]])
        names(PA14_1000)[length(PA14_1000)] <- names(PA14List)[i]
        print(PA14_1000)[length(PA14_1000)]
      }
    }
  } #Done
  
  for (i in 2001:3000) {
    for (j in 1:length(keggGet(names(PA14List)[i])[[1]]$PATHWAY)) {
      if (is.null(keggGet(names(PA14List)[i])[[1]]$PATHWAY) != TRUE) {
        PA14_1000 <- c(PA14_1000, keggGet(names(PA14List)[i])[[1]]$PATHWAY[[j]])
        names(PA14_1000)[length(PA14_1000)] <- names(PA14List)[i]
        print(PA14_1000)[length(PA14_1000)]
      }
    }
  } #Done
  
  for (i in 3001:4000) {
    for (j in 1:length(keggGet(names(PA14List)[i])[[1]]$PATHWAY)) {
      if (is.null(keggGet(names(PA14List)[i])[[1]]$PATHWAY) != TRUE) {
        PA14_1000 <- c(PA14_1000, keggGet(names(PA14List)[i])[[1]]$PATHWAY[[j]])
        names(PA14_1000)[length(PA14_1000)] <- names(PA14List)[i]
        print(PA14_1000)[length(PA14_1000)]
      }
    }
  } #Done
  
  for (i in 4001:5000) {
    for (j in 1:length(keggGet(names(PA14List)[i])[[1]]$PATHWAY)) {
      if (is.null(keggGet(names(PA14List)[i])[[1]]$PATHWAY) != TRUE) {
        PA14_1000 <- c(PA14_1000, keggGet(names(PA14List)[i])[[1]]$PATHWAY[[j]])
        names(PA14_1000)[length(PA14_1000)] <- names(PA14List)[i]
        print(PA14_1000)[length(PA14_1000)]
      }
    }
  } #Done
  
  for (i in 5001:5977) {
    for (j in 1:length(keggGet(names(PA14List)[i])[[1]]$PATHWAY)) {
      if (is.null(keggGet(names(PA14List)[i])[[1]]$PATHWAY) != TRUE) {
        PA14_1000 <- c(PA14_1000, keggGet(names(PA14List)[i])[[1]]$PATHWAY[[j]])
        names(PA14_1000)[length(PA14_1000)] <- names(PA14List)[i]
        print(PA14_1000)[length(PA14_1000)]
      }
    }
  }
  }
  
  #Convert large vector into single data frame
  {
  DF <- data.frame(matrix(NA, nrow = length(PA14_1000), ncol = 2))
  colnames(DF) <- c("Gene", "Pathway")
  
  DF$Gene <- names(PA14_1000)
  
  for (i in 1:length(PA14_1000)) {
    DF$Pathway[i] <- PA14_1000[[i]]
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
  
  #PAO1
  {
  PAO1List <- keggList("pae") # pae = PA01
  }
  
}


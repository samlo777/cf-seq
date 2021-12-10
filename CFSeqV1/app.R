#####

#Load Libraries
{
library(shiny)
library(shinydashboard)
library(DT) 
library(tidyverse)
#library(datamods) #Don't think this one is necessary
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(manipulateWidget)
library(scattermore)
}

#####

#UI Code 
{
ui <- dashboardPage(
  
  #Initialize Header
  dashboardHeader(title = "Welcome to CF-Seq!"),
  skin = "green",
  
  #Initialize Sidebar
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("How to Use", tabName = "How_to_Use", icon = icon("th"),  badgeLabel = "Start Here!", badgeColor = "green"),
      menuItem("Choose Studies", tabName = "Choose_Studies", icon = icon("dashboard")),
      menuItem("Run Analysis", tabName = "Run_Analysis", icon = icon("th"))
    )
  ),
  
  #Initialize Body
  dashboardBody(
    
    #Setup JS + CSS
    useShinyjs(), 
    
    #Initialize Tab Items
    tabItems(
    
    #Create Study Selection Window
    tabItem(tabName = "Choose_Studies",  
    includeScript(path = "./www/ButtonFunction.js"),
    
    fluidRow(
      column(12,
             strong(p(class = "strong",
                      "1. Choose your species"
             )),
             box(width = NULL, status = "warning",
                 selectizeInput(inputId = "Select_Species", label =  "Choose Your Species of Interest",
                             choices = c("[Pick a Species]", names(CFSeq_Data)),
                             selected = "[Pick a Species]"
                 ),
                 actionButton("RESET", "Reset Species")
             )
      )
    ), 
    
    fluidRow(
      column(8, 
             img(id = "Logo", src = "CFSeq_Logo.png", width = "100%", height = 520),
             ), 
      
      column(4, 
             strong(p(class = "strong",
                      "2. Choose your study filters"
             )),
             box(width = "100%", height = "100%", status = "warning",
                 
                 selectizeInput(inputId = "Select_SubSpecies", "Choose Your Sub-Species of Interest",
                                choices = c("[Select A Sub-Species]"),
                                selected = "[Select A Sub-Species]"
                 ),
                 selectizeInput(inputId = "Select_Strain", "Choose Your Strain of Interest",
                             choices = c("[Select A Strain]"),
                             selected = "[Select A Strain]"
                 ),
                 selectizeInput(inputId = "Select_Treatment", "Choose Your Treatment of Interest",
                               choices = c("[Select A Treatment]"),
                               selected = "[Select A Treatment]"
                 ),
                 selectizeInput(inputId = "Select_Medium", "Choose Your Medium of Interest",
                               choices = c("[Select A Medium]"),
                               selected = "[Select A Medium]"
                 ),
                 selectizeInput(inputId = "Select_Genotype", "Choose Your Gene Perturbation of Interest",
                                choices = c("[Select A Gene Perturbation]"),
                                selected = "[Select A Gene Perturbation]"
                 ),
                 p(class = "text-muted",
                   "Select from any of the conditions above",
                   "to refine your search."
                 ),
                 actionButton("button", "Finalize Inputs"),
                 actionButton("RESET_FILTER", "Reset Filters")
             )
    
             )
    ),
    
    br(),
    
    fluidRow(
      column(12,
        DT::dataTableOutput("table", width = "100%")
      )
    )
  ),
      
  #Create Run Analysis Window
  tabItem(tabName = "Run_Analysis", 
     
     fluidRow(
       column(12,
              DT::dataTableOutput("ResultTable", width = "100%")
              )
     ),
     
     fluidRow(
       column(8,
              tabsetPanel(id = "Analysis_Plots", type = "tabs",
                          tabPanel("Volcano Plot", plotlyOutput("VolcanoPlot")),
                          tabPanel("MA Plot", plotlyOutput("MAPlot"))
              )
              
       ),
       column(4,
              br(),
              br(),
              box(width = NULL, status = "warning",
                  selectizeInput(inputId = "Select_Comparison", label =  "Select A Comparison",
                              choices = c("[Select A Comparison]"),
                              selected = "[Select A Comparison]"
                  ),
                  selectizeInput(inputId = "Find_Gene", label =  "Find a Gene",
                              choices = c("[Select a Gene]"),
                              selected = "[Select a Gene]"
                  ),
                  selectizeInput(inputId = "Find_Pathway", label =  "Find a Pathway",
                              choices = c("[Select a Pathway]"),
                              selected = "[Select a Pathway]"
                  ),
                  downloadButton("downloadCountDT", "Download Study Count Table"),
                  downloadButton("fullAnalysisOutput", "Download All Analysis Data"),
                  actionButton("CLEAR_DATA2", "Reset Analysis")
              ))
     ),
     
     fluidRow(
       column(12,
              br(),
              box(width = NULL, status = "warning",
                  actionButton("CLEAR_DATA1", "Reset Analysis")
              ))
     )
     
  )

)
)
)
}

#####

#Additional Functions
{
create_btns <- function(x) {
  x %>%
    purrr::map_chr(~
                     paste0(
                   '<div class = "btn-group">
                    <button class="btn btn-default action-button btn-info action_button" id="StudyInfo_',
                       .x, '" type="button" onclick=get_id(this.id)><i class="fas fa-arrow-alt-circle-right" style="color: #19AA53; background-color: #EDEDED"></i></button></div>'
                     ))
}
}

#####

#Server Code
{
server <- function(input, output, session) {
  
  #1. Select a species and generate species-specific data
  {
  observeEvent(input$Select_Species, {
    
    #A. Initialize variables
    {
      VALUES <- reactiveValues(CFSeqStrain = NULL, CFSeqTreatment  = NULL, CFSeqMedium = NULL, CFSeqGenotype = NULL, CFSeqSubSpecies = NULL, CFSeq_Data = CFSeq_Data, S = NULL) 
      output$table <- NULL
      Strains <- vector()
      Treatments <- vector()
      Media <- vector()
      Genotype <- vector()
      Species2 <- vector()
      GeneIDList <- vector()
      GeneList <- vector()
      
    }
    
    if(input$Select_Species != "[Pick a Species]") {
      observe({
        updateSelectizeInput(session, inputId = "Select_Species", "Choose Your Species of Interest", choices = input$Select_Species,
                             selected = input$Select_Species)
      })
    } else {
      observe({
        updateSelectizeInput(session, inputId = "Select_Species", "Choose Your Species of Interest", choices = c("[Pick a Species]", names(CFSeq_Data)),
                             selected = "[Pick a Species]")
      })
    }
    
    #B. Update drop-down menus
    {
    if (input$Select_Species == "Aspergillus_fumigatus") {
      VALUES$S = 1
    } else if (input$Select_Species == "Burkholderia_Species") {
      VALUES$S = 2
    } else if (input$Select_Species == "Candida_albicans") {
      VALUES$S = 3
    } else if (input$Select_Species == "Clostridium_difficile") {
      VALUES$S = 4
    } else if (input$Select_Species == "Faecalibacterium_prausnitzii") {
      VALUES$S = 5
    } else if (input$Select_Species == "Haemophilus_influenzae") {
      VALUES$S = 6
    } else if (input$Select_Species == "Mycobacterium_abscessus") {
      VALUES$S = 7
    } else if (input$Select_Species == "Pseudomonas_aeruginosa") {
      VALUES$S = 8
    } else if (input$Select_Species == "Porphyromonas_Species") {
      VALUES$S = 9
    } else if (input$Select_Species == "Staphylococcus_aureus") {
      VALUES$S = 10
    } else if (input$Select_Species == "Stenotrophomonas_maltophilia") {
      VALUES$S = 11
    } else if (input$Select_Species == "Streptococcus_Species") {
      VALUES$S = 12
    } else {
      VALUES$S = NULL
    }
    
    
    if (is.null(VALUES$S) != TRUE) {
    
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
      Strains[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Strain
    }
    Strains <- na.omit(unique(Strains))
    if (length(Strains) != 0) {
      Strains <- unlist(strsplit(Strains, ","))
    }
    
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
      Treatments[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Treatment
    }
    Treatments <- na.omit(unique(Treatments))
    if (length(Treatments) != 0) {
      Treatments <- unlist(strsplit(Treatments, ","))
    }
    
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
      Media[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Medium
    }
    Media <- na.omit(unique(Media))
    if (length(Media) != 0) {
      Media <- unlist(strsplit(Media, ","))
    }
    
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
      Genotype[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Genotype
    }
    Genotype <- na.omit(unique(Genotype))
    if (length(Genotype) != 0) {
      Genotype <- unlist(strsplit(Genotype, ","))
    }
    
    if (is.null(VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Species) != TRUE) {
      for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
        Species2[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Species
      }
      Species2 <- na.omit(unique(Species2))
      if (length(Species2) != 0) {
        Species2 <- unlist(strsplit(Species2, ","))
      }
    } else {
      Species2 <- NULL
    }
    
    observe({
      updateSelectizeInput(session, inputId = "Select_SubSpecies", "Choose Your Sub-Species of Interest", choices = c("[Select A Sub-Species]", Species2),
                        selected = "[Select A Sub-Species]")
      updateSelectizeInput(session, inputId = "Select_Strain", "Choose Your Strain of Interest", choices = c("[Select A Strain]", Strains),
                        selected = "[Select A Strain]")
      updateSelectizeInput(session, inputId = "Select_Treatment", "Choose Your Treatment of Interest", choices = c("[Select A Treatment]", Treatments),
                        selected = "[Select A Treatment]")
      updateSelectizeInput(session, inputId = "Select_Medium", "Choose Your Medium of Interest", choices = c("[Select A Medium]", Media),
                        selected = "[Select A Medium]")
      updateSelectizeInput(session, inputId = "Select_Genotype", "Choose Your Gene Perturbation of Interest", choices = c("[Select A Gene Perturbation]", Genotype),
                        selected = "[Select A Gene Perturbation]")
    })
    
    }
      
    if (is.null(Species2) == TRUE) {
      shinyjs::hide(id = "Select_SubSpecies")
    } else {
      shinyjs::show(id = "Select_SubSpecies")
    }
      
    }
    
    if (input$Select_Species != "[Pick a Species]") {

    #C. Update studies shown
    {
        #Output starting data table 
        DF = data.frame(Date=character(),
                        GEO_Accession = character(),
                        Title=character(),
                        stringsAsFactors=FALSE)
        
        for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
          DF[i,] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]][,c(1:3)]
        }
        
        
        btn2 <- create_btns(1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) 
        DF = DF %>% dplyr::bind_cols(tibble("Get More Info + Run Analysis" = btn2))
        
        values <- reactiveValues(
          df = DF,
          dt_row = NULL,
          keep_track_id = nrow(df) + 1
        )
        
        output$table <- DT::renderDataTable(
          {
            shiny::isolate(values$df)
          },
          escape = F,
          rownames = FALSE,
          options = list(processing = FALSE)
        )
        
    }
      
    #D. Drop-down Inputs: Select Filtering inputs 
    {
      
      
      #Set Variables to Lock In Filters
      {
      CompleteS = FALSE
      CompleteT = FALSE
      CompleteM = FALSE
      CompleteG = FALSE
      CompleteSub = FALSE
      
      makeReactiveBinding("CompleteS")
      makeReactiveBinding("CompleteT")
      makeReactiveBinding("CompleteM")
      makeReactiveBinding("CompleteG")
      makeReactiveBinding("CompleteSub")
      }
      
      #Select a Strain
      {
      observeEvent(input$Select_Strain, {

        #Create Strains Data Object
        {
        VALUES$CFSeqStrain <- VALUES$CFSeq_Data[[VALUES$S]]
        
        if (input$Select_Strain != "[Select A Strain]") {
          
          if(names(CFSeq_Data[VALUES$S]) == input$Select_Species) {
          for (i in 1:(length(VALUES$CFSeqStrain))) {
            if (length(grep(input$Select_Strain, VALUES$CFSeqStrain[[i]][[3]]$Strain)) == 0) {
              names(VALUES$CFSeqStrain)[i] <- "remove"
            }
          }
          VALUES$CFSeqStrain <- VALUES$CFSeqStrain[names(VALUES$CFSeqStrain) %in% "remove" == FALSE]  
          
          } 
        }
        }
        
        #Reset Drop-Down menus
        {
        UTreatments <- c()
        NTreatments <- c()
        for (i in 1:length(VALUES$CFSeqStrain)) {
          UTreatments[i] <- VALUES$CFSeqStrain[[i]][[3]]$Treatment
        }
        if(length(na.omit(UTreatments)) != 0) {
          UTreatments <- unlist(strsplit(na.omit(UTreatments), ","))
        }
        for (i in 1:length(UTreatments)) {
          NTreatments[i] <- Treatments[grep(UTreatments[i], Treatments)]
        }

        UMedia <- c()
        NMedia <- c()
        for (i in 1:length(VALUES$CFSeqStrain)) {
          UMedia[i] <- VALUES$CFSeqStrain[[i]][[3]]$Medium
        }
        if(length(na.omit(UMedia)) != 0) {
          UMedia <- unlist(strsplit(na.omit(UMedia), ","))
        }
        for (i in 1:length(UMedia)) {
          NMedia[i] <- Media[grep(UMedia[i], Media)]
        }

        UGenotype <- c()
        NGenotype <- c()
        for (i in 1:length(VALUES$CFSeqStrain)) {
          UGenotype[i] <- VALUES$CFSeqStrain[[i]][[3]]$Genotype
        }
        if(length(na.omit(UGenotype)) != 0) {
          UGenotype <- unlist(strsplit(na.omit(UGenotype), ","))
        }
        for (i in 1:length(UGenotype)) {
          NGenotype[i] <- Genotype[grep(UGenotype[i], Genotype)]
        }
        
        if (is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
          USubSpecies <- c()
          NSubSpecies <- c()
          for (i in 1:length(VALUES$CFSeqStrain)) {
            USubSpecies[i] <- VALUES$CFSeqStrain[[i]][[3]]$Species
          }
          if(length(na.omit(USubSpecies)) != 0) {
            USubSpecies <- unlist(strsplit(na.omit(USubSpecies), ","))
          }
          for (i in 1:length(USubSpecies)) {
            NSubSpecies[i] <- Species2[grep(USubSpecies[i], Species2)]
          }
        }
        
          if (CompleteT == FALSE) {
            observe({
              updateSelectizeInput(session, inputId = "Select_Treatment", "Choose Your Treatment of Interest", choices = c("[Select A Treatment]", na.omit(NTreatments)),
                                   selected = "[Select A Treatment]")
            })
          }
          if (CompleteM == FALSE) {
            observe({
              updateSelectizeInput(session, inputId = "Select_Medium", "Choose Your Medium of Interest", choices = c("[Select A Medium]", na.omit(NMedia)),
                                   selected = "[Select A Medium]")
            })
          }
          if (CompleteG == FALSE) {
            observe({
              updateSelectizeInput(session, inputId = "Select_Genotype", "Choose Your Gene Perturbation of Interest", choices = c("[Select A Gene Perturbation]", na.omit(NGenotype)),
                                   selected = "[Select A Gene Perturbation]")
            })
          }
          
          if (CompleteSub == FALSE & is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
            observe({
              updateSelectizeInput(session, inputId = "Select_SubSpecies", "Choose Your Sub-Species of Interest", choices = c("[Select A Sub-Species]", na.omit(NSubSpecies)),
                                   selected = "[Select A Sub-Species]")
            })
          }
          
        
        }
        
        if (input$Select_Strain != "[Select A Strain]") {
          CompleteS <<- TRUE
        }
      
       #print(names(VALUES$CFSeqStrain))
        
        
        StrainSelect <<- input$Select_Strain
        
      })
      }
      
      #Select a Treatment 
      {
        observeEvent(input$Select_Treatment, {
          
          #Create Treatment Data Object
          {
          VALUES$CFSeqTreatment <- VALUES$CFSeq_Data[[VALUES$S]]
          if (input$Select_Treatment != "[Select A Treatment]") {
            
            if(names(CFSeq_Data[VALUES$S]) == input$Select_Species) {
            for (i in 1:(length(VALUES$CFSeqTreatment))) {
              if (length(grep(input$Select_Treatment, VALUES$CFSeqTreatment[[i]][[3]]$Treatment)) == 0) {
                names(VALUES$CFSeqTreatment)[i] <- "remove"
              }
            }
            VALUES$CFSeqTreatment <- VALUES$CFSeqTreatment[names(VALUES$CFSeqTreatment) %in% "remove" == FALSE]    
            }
          }
          }
          
          #Reset Drop-Down Menus
          {
              UStrains <- c()
              NStrains <- c()
              for (i in 1:length(VALUES$CFSeqTreatment)) {
                UStrains[i] <- VALUES$CFSeqTreatment[[i]][[3]]$Strain
              }
              if(length(na.omit(UStrains)) != 0) {
                UStrains <- unlist(strsplit(na.omit(UStrains), ","))
              }
              for (i in 1:length(UStrains)) {
                NStrains[i] <- Strains[grep(UStrains[i],Strains)]
              }

              UMedia <- c()
              NMedia <- c()
              for (i in 1:length(VALUES$CFSeqTreatment)) {
                UMedia[i] <- VALUES$CFSeqTreatment[[i]][[3]]$Medium
              }
              if(length(na.omit(UMedia)) != 0) {
                UMedia <- unlist(strsplit(na.omit(UMedia), ","))
              }
              for (i in 1:length(UMedia)) {
                NMedia[i] <- Media[grep(UMedia[i], Media)]
              }

              UGenotype <- c()
              NGenotype <- c()
              for (i in 1:length(VALUES$CFSeqTreatment)) {
                UGenotype[i] <- VALUES$CFSeqTreatment[[i]][[3]]$Genotype
              }
              if(length(na.omit(UGenotype)) != 0) {
                UGenotype <- unlist(strsplit(na.omit(UGenotype), ","))
              }
              for (i in 1:length(UGenotype)) {
                NGenotype[i] <- Genotype[grep(UGenotype[i], Genotype)]
              }
              
              if (is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
              USubSpecies <- c()
              NSubSpecies <- c()
              for (i in 1:length(VALUES$CFSeqTreatment)) {
                USubSpecies[i] <- VALUES$CFSeqTreatment[[i]][[3]]$Species
              }
              if(length(na.omit(USubSpecies)) != 0) {
                USubSpecies <- unlist(strsplit(na.omit(USubSpecies), ","))
              }
              for (i in 1:length(USubSpecies)) {
                NSubSpecies[i] <- Species2[grep(USubSpecies[i], Species2)]
              }
              }
              

                if (CompleteS == FALSE) {
                  observe({
                    updateSelectizeInput(session, inputId = "Select_Strain", "Choose Your Strain of Interest", choices = c("[Select A Strain]", na.omit(NStrains)),
                                         selected = "[Select A Strain]")
                  })
                }
                if (CompleteM == FALSE) {
                  observe({
                    updateSelectizeInput(session, inputId = "Select_Medium", "Choose Your Medium of Interest", choices = c("[Select A Medium]", na.omit(NMedia)),
                                         selected = "[Select A Medium]")
                  })
                }
                if (CompleteG == FALSE) {
                  observe({
                    updateSelectizeInput(session, inputId = "Select_Genotype", "Choose Your Gene Perturbation of Interest", choices = c("[Select A Gene Perturbation]", na.omit(NGenotype)),
                                         selected = "[Select A Gene Perturbation]")
                  })
                }
                
                if (CompleteSub == FALSE & is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
                  observe({
                    updateSelectizeInput(session, inputId = "Select_SubSpecies", "Choose Your Sub-Species of Interest", choices = c("[Select A Sub-Species]", na.omit(NSubSpecies)),
                                         selected = "[Select A Sub-Species]")
                  })
                }
              
              
                
          }
          
          if (input$Select_Treatment != "[Select A Treatment]") {
            CompleteT <<- TRUE
          }
          
          TreatmentSelect <<- input$Select_Treatment
          
        })  
      }
      
      
      #Select a Medium
      {
        observeEvent(input$Select_Medium, {
          
          #Create Medium Data Object
          {
          VALUES$CFSeqMedium <- VALUES$CFSeq_Data[[VALUES$S]]
          if (input$Select_Medium != "[Select A Medium]") {
            
            if(names(CFSeq_Data[VALUES$S]) == input$Select_Species) {
            for (i in 1:(length(VALUES$CFSeqMedium))) {
              if (length(grep(input$Select_Medium, VALUES$CFSeqMedium[[i]][[3]]$Medium)) == 0) {
                names(VALUES$CFSeqMedium)[i] <- "remove"
              }
            }
            VALUES$CFSeqMedium <- VALUES$CFSeqMedium[names(VALUES$CFSeqMedium) %in% "remove" == FALSE]    
            } 
          }
          }
          
          #Reset Drop-Down Menus
          {
            UStrains <- c()
            NStrains <- c()
            for (i in 1:length(VALUES$CFSeqMedium)) {
              UStrains[i] <- VALUES$CFSeqMedium[[i]][[3]]$Strain
            }
            if(length(na.omit(UStrains)) != 0) {
              UStrains <- unlist(strsplit(na.omit(UStrains), ","))
            }
            for (i in 1:length(UStrains)) {
              NStrains[i] <- Strains[grep(UStrains[i], Strains)]
            }
            
            UTreatments <- c()
            NTreatments <- c()
            for (i in 1:length(VALUES$CFSeqMedium)) {
              UTreatments[i] <- VALUES$CFSeqMedium[[i]][[3]]$Treatment
            }
            if(length(na.omit(UTreatments)) != 0) {
              UTreatments <- unlist(strsplit(na.omit(UTreatments), ","))
            }
            for (i in 1:length(UTreatments)) {
              NTreatments[i] <- Treatments[grep(UTreatments[i], Treatments)]
            }
            
            UGenotype <- c()
            NGenotype <- c()
            for (i in 1:length(VALUES$CFSeqMedium)) {
              UGenotype[i] <- VALUES$CFSeqMedium[[i]][[3]]$Genotype
            }
            if(length(na.omit(UGenotype)) != 0) {
              UGenotype <- unlist(strsplit(na.omit(UGenotype), ","))
            }
            for (i in 1:length(UGenotype)) {
              NGenotype[i] <- Genotype[grep(UGenotype[i], Genotype)]
            }
            
            if (is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
            USubSpecies <- c()
            NSubSpecies <- c()
            for (i in 1:length(VALUES$CFSeqMedium)) {
              USubSpecies[i] <- VALUES$CFSeqMedium[[i]][[3]]$Species
            }
            if(length(na.omit(USubSpecies)) != 0) {
              USubSpecies <- unlist(strsplit(na.omit(USubSpecies), ","))
            }
            for (i in 1:length(USubSpecies)) {
              NSubSpecies[i] <- Species2[grep(USubSpecies[i], Species2)]
            }
            }
            
      
              

              if (CompleteS == FALSE) {
                observe({
                  updateSelectizeInput(session, inputId = "Select_Strain", "Choose Your Strain of Interest", choices = c("[Select A Strain]", na.omit(NStrains)),
                                       selected = "[Select A Strain]")
                })
                
              }
              if (CompleteT == FALSE) {
                observe({
                  updateSelectizeInput(session, inputId = "Select_Treatment", "Choose Your Treatment of Interest", choices = c("[Select A Treatment]", na.omit(NTreatments)),
                                       selected = "[Select A Treatment]")
                })
              }
              if (CompleteG == FALSE) {
                observe({
                  updateSelectizeInput(session, inputId = "Select_Genotype", "Choose Your Gene Perturbation of Interest", choices = c("[Select A Gene Perturbation]", na.omit(NGenotype)),
                                       selected = "[Select A Gene Perturbation]")
                })
                
              }
              
              if (CompleteSub == FALSE & is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
                observe({
                  updateSelectizeInput(session, inputId = "Select_SubSpecies", "Choose Your Sub-Species of Interest", choices = c("[Select A Sub-Species]", na.omit(NSubSpecies)),
                                       selected = "[Select A Sub-Species]")
                })
              }
            
            
          }
          
          if (input$Select_Medium != "[Select A Medium]") {
            CompleteM <<- TRUE
            MTest <- TRUE
          }

          MediumSelect <<- input$Select_Medium
        })
      }
      
      
      #Select a Gene Perturbation
      {
      observeEvent(input$Select_Genotype, {
        
        #Create Genotype Data Object
        {
        VALUES$CFSeqGenotype <- VALUES$CFSeq_Data[[VALUES$S]]
        if (input$Select_Genotype != "[Select A Gene Perturbation]") {
          
          if(names(CFSeq_Data[VALUES$S]) == input$Select_Species) {
          for (i in 1:(length(VALUES$CFSeqGenotype))) {
            if (length(grep(input$Select_Genotype, VALUES$CFSeqGenotype[[i]][[3]]$Genotype)) == 0) {
              names(VALUES$CFSeqGenotype)[i] <- "remove"
            }
          }
          VALUES$CFSeqGenotype <- VALUES$CFSeqGenotype[names(VALUES$CFSeqGenotype) %in% "remove" == FALSE]    
          }
        }
        }
      
        #Reset Drop-Down Menus
        {
          UStrains <- c()
          NStrains <- c()
          for (i in 1:length(VALUES$CFSeqGenotype)) {
            UStrains[i] <- VALUES$CFSeqGenotype[[i]][[3]]$Strain
          }
          if(length(na.omit(UStrains)) != 0) {
            UStrains <- unlist(strsplit(na.omit(UStrains), ","))
          }
          for (i in 1:length(UStrains)) {
            NStrains[i] <- Strains[grep(UStrains[i], Strains)]
          }          
          UTreatments <- c()
          NTreatments <- c()
          for (i in 1:length(VALUES$CFSeqGenotype)) {
            UTreatments[i] <- VALUES$CFSeqGenotype[[i]][[3]]$Treatment
          }
          if(length(na.omit(UTreatments)) != 0) {
            UTreatments <- unlist(strsplit(na.omit(UTreatments), ","))
          }
          for (i in 1:length(UTreatments)) {
            NTreatments[i] <- Treatments[grep(UTreatments[i], Treatments)]
          }       
          
          UMedia <- c()
          NMedia <- c()
          for (i in 1:length(VALUES$CFSeqGenotype)) {
            UMedia[i] <- VALUES$CFSeqGenotype[[i]][[3]]$Medium
          }
          if(length(na.omit(UMedia)) != 0) {
            UMedia <- unlist(strsplit(na.omit(UMedia), ","))
          }
          for (i in 1:length(UMedia)) {
            NMedia[i] <- Media[grep(UMedia[i], Media)]
          }
          
          if (is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
          USubSpecies <- c()
          NSubSpecies <- c()
          for (i in 1:length(VALUES$CFSeqGenotype)) {
            USubSpecies[i] <- VALUES$CFSeqGenotype[[i]][[3]]$Species
          }
          if(length(na.omit(USubSpecies)) != 0) {
            USubSpecies <- unlist(strsplit(na.omit(USubSpecies), ","))
          }
          for (i in 1:length(USubSpecies)) {
            NSubSpecies[i] <- Species2[grep(USubSpecies[i], Species2)]
          }
          }
          
            
            if (CompleteS == FALSE) {
              observe({
                updateSelectizeInput(session, inputId = "Select_Strain", "Choose Your Strain of Interest", choices = c("[Select A Strain]", na.omit(NStrains)),
                                     selected = "[Select A Strain]")
              })
             
            }
            if (CompleteT == FALSE) {
              observe({
                updateSelectizeInput(session, inputId = "Select_Treatment", "Choose Your Treatment of Interest", choices = c("[Select A Treatment]", na.omit(NTreatments)),
                                     selected = "[Select A Treatment]")
              })
            }
            if (CompleteM == FALSE) {
              observe({
                updateSelectizeInput(session, inputId = "Select_Medium", "Choose Your Medium of Interest", choices = c("[Select A Medium]", na.omit(NMedia)),
                                     selected = "[Select A Medium]")
              })
              
            }
            
            if (CompleteSub == FALSE & is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
              observe({
                updateSelectizeInput(session, inputId = "Select_SubSpecies", "Choose Your Sub-Species of Interest", choices = c("[Select A Sub-Species]", na.omit(NSubSpecies)),
                                     selected = "[Select A Sub-Species]")
              })
            }
        
          
        }
        
        if (input$Select_Genotype != "[Select A Gene Perturbation]") {
          CompleteG <<- TRUE
          GTest <- TRUE
        }
        
        GeneSelect <<- input$Select_Genotype
        
      })
      }
      
      
      #Select a Sub-Species
      {
        observeEvent(input$Select_SubSpecies, {
          
          #Create Sub-Species Data Object
          {
            VALUES$CFSeqSubSpecies <- VALUES$CFSeq_Data[[VALUES$S]]
            if (input$Select_SubSpecies != "[Select A Sub-Species]" & is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
              
              if(names(CFSeq_Data[VALUES$S]) == input$Select_Species) {
              for (i in 1:(length(VALUES$CFSeqSubSpecies))) {
                if (length(grep(input$Select_SubSpecies, VALUES$CFSeqSubSpecies[[i]][[3]]$Species)) == 0) {
                  names(VALUES$CFSeqSubSpecies)[i] <- "remove"
                }
              }
              VALUES$CFSeqSubSpecies <- VALUES$CFSeqSubSpecies[names(VALUES$CFSeqSubSpecies) %in% "remove" == FALSE]    
              }
            }
          }
          
          
          if (is.null(VALUES$CFSeq_Data[[VALUES$S]][[1]][[3]]$Species) != TRUE) {
            #Reset Drop-Down Menus
            {
              UStrains <- c()
              NStrains <- c()
              for (i in 1:length(VALUES$CFSeqSubSpecies)) {
                UStrains[i] <- VALUES$CFSeqSubSpecies[[i]][[3]]$Strain
              }
              if(length(na.omit(UStrains)) != 0) {
                UStrains <- unlist(strsplit(na.omit(UStrains), ","))
              }
              for (i in 1:length(UStrains)) {
                NStrains[i] <- Strains[grep(UStrains[i], Strains)]
              }          
              UTreatments <- c()
              NTreatments <- c()
              for (i in 1:length(VALUES$CFSeqSubSpecies)) {
                UTreatments[i] <- VALUES$CFSeqSubSpecies[[i]][[3]]$Treatment
              }
              if(length(na.omit(UTreatments)) != 0) {
                UTreatments <- unlist(strsplit(na.omit(UTreatments), ","))
              }
              for (i in 1:length(UTreatments)) {
                NTreatments[i] <- Treatments[grep(UTreatments[i], Treatments)]
              }       
              
              UMedia <- c()
              NMedia <- c()
              for (i in 1:length(VALUES$CFSeqSubSpecies)) {
                UMedia[i] <- VALUES$CFSeqSubSpecies[[i]][[3]]$Medium
              }
              if(length(na.omit(UMedia)) != 0) {
                UMedia <- unlist(strsplit(na.omit(UMedia), ","))
              }
              for (i in 1:length(UMedia)) {
                NMedia[i] <- Media[grep(UMedia[i], Media)]
              }
              
              UGenotype <- c()
              NGenotype <- c()
              for (i in 1:length(VALUES$CFSeqSubSpecies)) {
                UGenotype[i] <- VALUES$CFSeqSubSpecies[[i]][[3]]$Genotype
              }
              if(length(na.omit(UGenotype)) != 0) {
                UGenotype <- unlist(strsplit(na.omit(UGenotype), ","))
              }
              for (i in 1:length(UGenotype)) {
                NGenotype[i] <- Genotype[grep(UGenotype[i], Genotype)]
              }
              

                
                if (CompleteS == FALSE) {
                  observe({
                    updateSelectizeInput(session, inputId = "Select_Strain", "Choose Your Strain of Interest", choices = c("[Select A Strain]", na.omit(NStrains)),
                                         selected = "[Select A Strain]")
                  })
                  
                }
                if (CompleteT == FALSE) {
                  observe({
                    updateSelectizeInput(session, inputId = "Select_Treatment", "Choose Your Treatment of Interest", choices = c("[Select A Treatment]", na.omit(NTreatments)),
                                         selected = "[Select A Treatment]")
                  })
                }
                if (CompleteM == FALSE) {
                  observe({
                    updateSelectizeInput(session, inputId = "Select_Medium", "Choose Your Medium of Interest", choices = c("[Select A Medium]", na.omit(NMedia)),
                                         selected = "[Select A Medium]")
                  })
                  
                }
                if (CompleteG == FALSE) {
                  observe({
                    updateSelectizeInput(session, inputId = "Select_Genotype", "Choose Your Gene Perturbation of Interest", choices = c("[Select A Gene Perturbation]", na.omit(NGenotype)),
                                         selected = "[Select A Gene Perturbation]")
                  })
                  
                }
                
              
            } 
          }
          
          
          if (input$Select_SubSpecies != "[Select A Sub-Species]") {
            CompleteSub <<- TRUE
            GTest <- TRUE
          }
          
          
          #print(names(VALUES$CFSeqGenotype))
        })
      }
      
      
      }
    
    #E. Action Button: Confirm filter selections
    {
      observeEvent(input$button, {
        
        if(names(CFSeq_Data[VALUES$S]) == input$Select_Species) {
        
        if (input$Select_Strain == "[Select A Strain]" & input$Select_Treatment == "[Select A Treatment]" & input$Select_Medium == "[Select A Medium]" & input$Select_Genotype == "[Select A Gene Perturbation]" & input$Select_SubSpecies == "[Select A Sub-Species]") {
          
          output$table <- DT::renderDataTable(
            {
              shiny::isolate(values$df)
            },
            escape = F,
            rownames = FALSE,
            options = list(processing = FALSE)
          )
          
          DF1 = data.frame(Date=character(),
                           GEO_Accession = character(),
                           Title=character(),
                           stringsAsFactors=FALSE)
          
          for (i in 1:length(VALUES$CFSeq_Data[[VALUES$S]])) {
            DF1[i,] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]][,c(1:3)]
          }
        
          btn1 <- create_btns(1:length(VALUES$CFSeq_Data[[VALUES$S]]))
          DF1 = DF1 %>% dplyr::bind_cols(tibble("Get More Info + Run Analysis" = btn1))
          
        } else {
        
        DF1 = data.frame(Date=character(),
                         GEO_Accession = character(),
                         Title=character(),
                         stringsAsFactors=FALSE)
        

        #Fill Data Frame with filter-selected studies
        {
        Approve = 0

        for (i in 1:length(CFSeq_Data[[VALUES$S]])) {
          
          #Reset Counter Variable
          Approve = 0
          
          #Sub-Species Check
          if (input$Select_SubSpecies != "[Select A Sub-Species]") {
            if (length(grep(input$Select_SubSpecies, VALUES$CFSeq_Data[[VALUES$S]][[i]]$Additional_Metadata$Species)) != 0) { 
              Approve = Approve
            } else {
              Approve = Approve + 1
            }
          } else {
            Approve = Approve
          }
          
          #Strain Check
          if (input$Select_Strain != "[Select A Strain]") {
            if (length(grep(input$Select_Strain, VALUES$CFSeq_Data[[VALUES$S]][[i]]$Additional_Metadata$Strain)) != 0) { 
              Approve = Approve
            } else {
              Approve = Approve + 1
            }
          } else {
            Approve = Approve
          }
          
          #Treatment Check
          if (input$Select_Treatment != "[Select A Treatment]") {
            if (length(grep(input$Select_Treatment, VALUES$CFSeq_Data[[VALUES$S]][[i]]$Additional_Metadata$Treatment)) != 0) { 
              Approve = Approve
            } else {
              Approve = Approve + 1
            }
          } else {
            Approve = Approve
          }

          #Medium Check
          if (input$Select_Medium != "[Select A Medium]") {
            if (length(grep(input$Select_Medium, VALUES$CFSeq_Data[[VALUES$S]][[i]]$Additional_Metadata$Medium)) != 0) { 
              Approve = Approve
            } else {
              Approve = Approve + 1
            }
          } else {
            Approve = Approve
          }
          
          #Genotype Check
          if (input$Select_Genotype != "[Select A Gene Perturbation]") {
            if (length(grep(input$Select_Genotype, VALUES$CFSeq_Data[[VALUES$S]][[i]]$Additional_Metadata$Genotype)) != 0) { 
              Approve = Approve
            } else {
              Approve = Approve + 1
            }
          } else {
            Approve = Approve
          }
          
          if (Approve == 0) {
            DF1[nrow(DF1)+1,] <- VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]][,c(1:3)]
          }
        }
        
          
        #Create additional columns for data frame 
        btn1 <- create_btns(1:nrow(DF1))
        DF1 = DF1 %>% dplyr::bind_cols(tibble("Get More Info + Run Analysis" = btn1))

        #Render data frame visually in app
        values1 <- reactiveValues(
          df = DF1,
          dt_row = NULL,
          keep_track_id = nrow(df) + 1
        )
        
        output$table <- DT::renderDataTable(
          {
            shiny::isolate(values1$df)
          },
          escape = F,
          rownames = FALSE,
          options = list(processing = FALSE)
        )
        }
        
        } 
        }
        
        DF1 <<- DF1
      })
    }
      
    #F. Reset Filter Inputs
    {
        observeEvent(input$RESET_FILTER, {
          observe({
            updateSelectizeInput(session, inputId = "Select_SubSpecies", choices = c("[Select A Sub-Species]", Species2),
                                 selected = "[Select A Sub-Species]")
          })
          observe({
            updateSelectizeInput(session, inputId = "Select_Strain", choices = c("[Select A Strain]", Strains),
                                 selected = "[Select A Strain]")
          })
          observe({
            updateSelectizeInput(session, inputId = "Select_Treatment", choices = c("[Select A Treatment]", Treatments),
                                 selected = "[Select A Treatment]")
          })
          observe({
            updateSelectizeInput(session, inputId = "Select_Medium", choices = c("[Select A Medium]", Media),
                                 selected = "[Select A Medium]")
          })
          observe({
            updateSelectizeInput(session, inputId = "Select_Genotype", choices = c("[Select A Gene Perturbation]", Genotype),
                                 selected = "[Select A Gene Perturbation]")
          })
          observe({
            updateSelectizeInput(session, inputId = "Select_SubSpecies", choices = c("[Select A Sub-Species]", Species2),
                                 selected = "[Select A Sub-Species]")
          })
          UStrains <- c() 
          NStrains <- c()
          UTreatments <- c()
          NTreatments <- c()
          UMedia <- c()
          NMedia <- c()
          UGenotype <- c()
          NGenotype <- c()
          CompleteS <<- FALSE
          CompleteT <<- FALSE
          CompleteM <<- FALSE
          CompleteG <<- FALSE
          CompleteSub <<- FALSE
          
          output$table <- DT::renderDataTable(
            {
              shiny::isolate(values$df)
            },
            escape = F,
            rownames = FALSE,
            options = list(processing = FALSE)
          )
        })
      }
      
    SpeciesID <<- VALUES$S
    Strains <<- Strains
    Treatments <<- Treatments
    Media <<- Media
    Genotype <<- Genotype
    Species2 <<- Species2
    SpeciesSelect1 <<- names(CFSeq_Data[VALUES$S])
    SpeciesSelect2 <<- input$Select_Species
    VALUES$CFSeqStrain <<- VALUES$CFSeqStrain
    VALUES$CFSeqTreatment <<- VALUES$CFSeqTreatment
    VALUES$CFSeqMedium <<- VALUES$CFSeqMedium
    VALUES$CFSeqGenotype <<- VALUES$CFSeqGenotype
    VALUES$CFSeqSubSpecies <<- VALUES$CFSeqSubSpecies
    
    } else {
      output$table <- NULL
    }
    

    
  })
 }
    
  #2. Action Button: More Details Dialogue + Run Analysis Window
  {
      observeEvent(input$current_id, {
       
        if(SpeciesSelect1 == SpeciesSelect2) {
          
          #Show Dialog Window with study info + Run Analysis Button 
          {
            shiny::req(!is.null(input$current_id) & stringr::str_detect(input$current_id, pattern = "StudyInfo"))
            
            if (StrainSelect == "[Select A Strain]" & TreatmentSelect == "[Select A Treatment]" & MediumSelect == "[Select A Medium]" & GeneSelect == "[Select A Gene Perturbation]") {
              FilterInput <- CFSeq_Data[[SpeciesID]]
            } else {
              FilterInput <- vector(mode = "list")
              for (i in 1:length(CFSeq_Data[[SpeciesID]])) {
                if (CFSeq_Data[[SpeciesID]][[i]][[3]]$GEO.Accession %in% DF1$GEO_Accession) {
                  FilterInput[[length(FilterInput)+1]] <- CFSeq_Data[[SpeciesID]][[i]]
                }
              }
            }
            FilterInput <<- FilterInput
            
            #Comparisons
            Comparison <<- names(FilterInput[[as.integer(substr(input$current_id,11,15))]][4:length(CFSeq_Data[[SpeciesID]][[as.integer(substr(input$current_id,11,15))]])])
            Comparison <<- Comparison[grep("Full Results", Comparison)]
            
            showModal(modalDialog(
              title = "Study Details",
              div(
                class = "text-center",
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Strain: ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Strain)
                  })
                ),
                br(),
                br(),
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Medium: ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Medium)
                  })
                ),
                br(),
                br(),
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Treatment: ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Treatment)
                  })
                ),
                br(),
                br(),
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Description: ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Description)
                  })
                ),
                br(),
                br(),
                div(
                  style = "display: inline-block;",
                  style = "color:red",
                  shiny::renderText({
                    if (length(Comparison) != 0) {
                      paste("Study design allows for differential expression analysis. Data visualization is available to view and download")
                    } else {
                      paste("Study design does not allow for differential expression analysis. Only count table is available for download")
                    }
                  })
                )
              ),
              footer = div(
                class = "pull-right container",
                shiny::actionButton(
                  inputId = "ShowAnalysis",
                  label = "Show Differential Expression Analysis",
                  icon = shiny::icon("chart-bar"),
                  class = "btn-info"
                ),
                shiny::actionButton(
                  inputId = "dismiss_modal",
                  label = "Close",
                  class = "btn-danger"
                )
              )
            ))
          }
          
        }
      })
      
  }
    
  #3. Press Show Analysis Button --> Show Analysis Screen with default comparison 
  {
      shiny::observeEvent(input$ShowAnalysis, {
        
        output$VolcanoPlot <- NULL
        
        #Show / Hide outputs based on study type
        {
          if (length(Comparison) == 0) {
            shinyjs::hide(id = "Select_Comparison")
            shinyjs::hide(id = "Find_Gene")
            shinyjs::hide(id = "Find_Pathway")
            shinyjs::hide(id = "downloadDataDT")
            shinyjs::hide(id = "downloadDataVP")
            shinyjs::hide(id = "downloadDataMA")
            shinyjs::hide(id = "Analysis_Plots")
            shinyjs::hide(id = "CLEAR_DATA1")
            shinyjs::show(id = "CLEAR_DATA2")
          } else {
            shinyjs::show(id = "Select_Comparison")
            shinyjs::show(id = "Find_Gene")
            shinyjs::show(id = "Find_Pathway")
            shinyjs::show(id = "downloadDataDT")
            shinyjs::show(id = "downloadDataVP")
            shinyjs::show(id = "downloadDataMA")
            shinyjs::show(id = "Analysis_Plots")
            shinyjs::show(id = "CLEAR_DATA1")
            shinyjs::hide(id = "CLEAR_DATA2")
          }
          
          if(SpeciesID != 8) {
            shinyjs::hide(id = "Find_Pathway")
          } else {
            shinyjs::show(id = "Find_Pathway")
          }
        }
        
        #Switch to Analysis Tab
        {
          newtab <- switch(input$tabs,
                           "Choose_Studies" = "Run_Analysis"
                           
          )
          updateTabItems(session, "tabs", newtab) 
        }
        
        #Update Comparison list and pathways for specific study 
        {
          #Comparisons (Again)
          Comparison <- names(FilterInput[[as.integer(substr(input$current_id,11,15))]][4:length(CFSeq_Data[[SpeciesID]][[as.integer(substr(input$current_id,11,15))]])])
          Comparison <- Comparison[grep("Full Results", Comparison)]
          
          #Pathways
          Pathways <- unique(PA14Paths$Pathway)
          
          observe({
            updateSelectInput(session, inputId = "Select_Comparison", "Choose Your Comparison of Interest", choices = c("[Select A Comparison]", Comparison),
                              selected = "[Select A Comparison]")
            updateSelectInput(session, inputId = "Find_Pathway", "Highlight a Pathway", choices = c("[Select a Pathway]", Pathways),
                              selected = "[Select a Pathway]")
          })
        }
        
        
        if (length(Comparison) != 0) {
          
          #User Selects a Comparison
          {
            observeEvent(input$Select_Comparison, {
              
              if (input$Select_Comparison != "[Select A Comparison]") {
                
                #Create data table, volcano plot, and MA plot
                {  
                  DF2 = data.frame(logFC=double(),
                                   logCPM = double(),
                                   FValue=double(),
                                   PValue=double(),
                                   stringsAsFactors=FALSE)
                  
                  DF2 <- FilterInput[[as.integer(substr(input$current_id,11,15))]][[which(names(FilterInput[[as.integer(substr(input$current_id,11,15))]]) == input$Select_Comparison)]]
                  
                  values3 <- reactiveValues(
                    df = DF2,
                    dt_row = NULL,
                    keep_track_id = nrow(df) + 1
                  )
                  
                  
                  #Output Data Table
                  output$ResultTable <- DT::renderDataTable(
                    {
                      shiny::isolate(values3$df)
                    },
                    escape = F,
                    rownames = TRUE,
                    options = list(processing = FALSE)
                  )
                  
                  #Create + Output Volcano Plot in App
                  {
                    
                    output$VolcanoPlot <- renderPlotly({
                      
                      plot_ly(values3$df, x = values3$df$logFC, y = -log10(values3$df$PValue),
                              text = ~paste("logFC: ", values3$df$logFC, '$<br>-log10(PValue):',
                                            -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                        layout(xaxis = list(title="logFC"),
                               yaxis = list(title="-log10(PValue)"))
                      
                      #GGTEST1 <-
                      #  ggplot(values3$df, aes(logFC, -log10(PValue), text=sprintf("Gene Name: %s", row.names(values3$df)))) + 
                      #  geom_point(alpha=0.3, size = 1)
                      
                      #GGTEST1
                      
                    })
                  }
                  
                  #Create + Output MA Plot in App
                  {
                    
                    output$MAPlot <- renderPlotly({
                      
                      plot_ly(values3$df, x = values3$df$logCPM, y = values3$df$logFC, 
                              text = ~paste("logCPM: ", values3$df$logCPM, '$<br>logFC:', 
                                            values3$df$logFC, '$<br>Gene Name:', row.names(values3$df))) %>%
                        layout(xaxis = list(title="logCPM"),
                               yaxis = list(title="logFC"))
                      
                      #GGTEST2 <-  
                      #  ggplot(values3$df, aes(logCPM, logFC, text=sprintf("Gene Name: %s", row.names(values3$df)))) + 
                      #  geom_point(alpha=0.3, size = 1)
                      
                      #GGTEST2
                      
                    })
                  }
                  
                  #Download Plots
                  {
                    GGVolcano <- values3$df %>%
                      ggplot(aes(logFC, y=-log10(PValue))) + 
                      geom_point(alpha=0.3)
                    
                    GGMA <- values3$df %>%
                      ggplot(aes(logCPM, y=logFC)) + 
                      geom_point(alpha=0.3)
                    
                    output$fullAnalysisOutput <- downloadHandler(
                      filename = function() {
                        paste0(input$Select_Comparison, ".zip")
                      },
                      content <- function(file) {
                        UserDir = tempdir()
                        write.csv(DF2, paste0(UserDir, "/Analysis_Table.csv"))
                        MyWidget = combineWidgets(ggplotly(GGVolcano), ggplotly(GGMA), ncol = 2)
                        saveWidget(MyWidget, paste0(UserDir, "/Analysis_Plots.html"))
                        zip(file, paste0(UserDir, "/Analysis_Table.csv"), paste0(UserDir, "/Analysis_Plots.html"), flags = "-r9Xj")
                      },
                      contentType = "application/zip")
                  }
                  
                  
                }
                
                #Update the gene list and pathways
                {
                  GeneList <- rownames(values3$df)
                  
                  observe({
                    updateSelectizeInput(session, inputId = "Find_Gene", "Find a Gene", choices = c(Choose = "", GeneList),
                                         options = list(maxOptions = 15000))
                    updateSelectizeInput(session, inputId = "Find_Pathway", "Find a Pathway", choices = c(Choose = "", Pathways),
                                         options = list(maxOptions = 15000))
                  })
                  
                }
                
                #Find a Gene on the Volcano or MA Plot
                {
                  observeEvent(input$Find_Gene, {
                    
                    if (input$Find_Gene != "" & input$Find_Gene != "[Select a Gene]") {
                      
                      Found_Gene <- values3$df[which(rownames(values3$df) == input$Find_Gene),]
                      
                      #Volcano Plot
                      {
                        
                        output$VolcanoPlot <- renderPlotly({
                          
                          print(Found_Gene)
                          VP1 <- plot_ly(type = "scatter", values3$df, x = values3$df$logFC, y = -log10(values3$df$PValue),
                                         text = ~paste("logFC: ", values3$df$logFC, '$<br>-log10(PValue):',
                                                       -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(showlegend = FALSE, xaxis = list(title="logFC"),
                                   yaxis = list(title="-log10(PValue)"))
                          
                          VP1 <- add_trace(VP1, x = Found_Gene$logFC, y = -log10(Found_Gene$PValue), type = "scatter", mode = "markers", color = I("red"), inherit = FALSE)
                          VP1
                          
                          #GGTEST3 <- values3$df %>%
                          #  ggplot(aes(logFC, -log10(PValue), text=sprintf("Gene Name: %s", row.names(values3$df)))) + 
                          #  geom_point(alpha=0.3) +
                          #  geom_point(data=Found_Gene, 
                          #             aes(logFC, -log10(PValue), text=sprintf("Gene Name: %s", row.names(Found_Gene))), 
                          #             color='blue',
                          #             size=2) 
                          
                          #GGTEST3
                          
                        })
                      }
                      
                      #MAPlot
                      {
                        
                        output$MAPlot <- renderPlotly({
                          
                          MA1 <- plot_ly(type = "scatter", values3$df, x = values3$df$logCPM, y = values3$df$logFC,
                                         text = ~paste("logCPM: ", values3$df$logCPM, '$<br>logFC:',
                                                       values3$df$logFC, '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="logCPM"),
                                   yaxis = list(title="logFC"))
                          
                          MA1 <- add_trace(MA1, x = Found_Gene$logCPM, y = Found_Gene$logFC, type = "scatter", mode = "markers", color = I("red"), inherit = FALSE)
                          MA1 
                          
                          #GGTEST4 <- values3$df %>%
                          #  ggplot(aes(logCPM, logFC, text=sprintf("Gene Name: %s", row.names(values3$df)))) + 
                          #  geom_point(alpha=0.3) +
                          #  geom_point(data=Found_Gene, 
                          #             aes(logCPM, logFC, text=sprintf("Gene Name: %s", row.names(Found_Gene))), 
                          #             color='blue',
                          #             size=2)
                          
                          #GGTEST4
                          
                        })
                      }
                      
                      #Download Plots
                      {
                        GGVolcano <- values3$df %>%
                          ggplot(aes(logFC, -log10(PValue))) + 
                          geom_point(alpha=0.3) +
                          geom_point(data=Found_Gene, 
                                     aes(logFC, -log10(PValue)), 
                                     color='blue',
                                     size=2)
                        
                        GGMA <- values3$df %>%
                          ggplot(aes(logCPM, logFC)) + 
                          geom_point(alpha=0.3) +
                          geom_point(data=Found_Gene, 
                                     aes(logCPM, logFC), 
                                     color='blue',
                                     size=2)
                        
                        output$fullAnalysisOutput <- downloadHandler(
                          filename = function() {
                            paste0(input$Select_Comparison, ".zip")
                          },
                          content <- function(file) {
                            UserDir = tempdir()
                            write.csv(DF2, paste0(UserDir, "/Analysis_Table.csv"))
                            MyWidget = combineWidgets(ggplotly(GGVolcano), ggplotly(GGMA), ncol = 2)
                            saveWidget(MyWidget, paste0(UserDir, "/Analysis_Plots.html"))
                            zip(file, paste0(UserDir, "/Analysis_Table.csv"), paste0(UserDir, "/Analysis_Plots.html"), flags = "-r9Xj")
                          },
                          contentType = "application/zip")
                      }
                      
                    }
                    
                  })
                }
                
                #Highlight a pathway on the Volcano or MA Plot
                {
                  observeEvent(input$Find_Pathway, {
                    
                    if (input$Find_Pathway != "" & input$Find_Pathway != "[Select a Pathway]") {
                      
                      Selected_Pathway <- PA14Paths[which(PA14Paths$Pathway == input$Find_Pathway),]
                      
                      Selected_Genes <- subset(values3$df, rownames(values3$df) %in% Selected_Pathway$Gene)
                      
                      #Volcano Plot
                      {
                        
                        output$VolcanoPlot <- renderPlotly({
                          
                          VP2 <- plot_ly(type = "scatter", values3$df, x = values3$df$logFC, y = -log10(values3$df$PValue),
                                         text = ~paste("logFC: ", values3$df$logFC, '$<br>-log10(PValue):',
                                                       -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="logFC"),
                                   yaxis = list(title="-log10(PValue)"))
                          
                          VP2 <- add_trace(VP2, x = Selected_Genes$logFC, y = -log10(Selected_Genes$PValue), type = "scatter", mode = "markers", color = I("green"), inherit = FALSE)
                          VP2
                          
                          #GGTEST5 <- values3$df %>%
                          #  ggplot(aes(logFC, -log10(PValue), text=sprintf("Gene Name: %s", row.names(values3$df)))) + 
                          #  geom_point(alpha=0.3) +
                          #  geom_point(data=Selected_Genes, 
                          #             aes(logFC, -log10(PValue), text=sprintf("Gene Name: %s", row.names(Selected_Genes))), 
                          #             color='green',
                          #             size=2) 
                          
                          #GGTEST5
                          
                        })
                      }
                      
                      #MA Plot
                      {
                        
                        output$MAPlot <- renderPlotly({
                          
                          MA2 <- plot_ly(type = "scatter", values3$df, x = values3$df$logCPM, y = values3$df$logFC,
                                         text = ~paste("logCPM: ", values3$df$logCPM, '$<br>logFC:',
                                                       values3$df$logFC, '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="logCPM"),
                                   yaxis = list(title="logFC"))
                          
                          MA2 <- add_trace(MA2, x = Selected_Genes$logCPM, y = Selected_Genes$logFC, type = "scatter", mode = "markers", color = I("green"), inherit = FALSE)
                          MA2
                          
                          #GGTEST6 <- values3$df %>%
                          #  ggplot(aes(logCPM, logFC, text=sprintf("Gene Name: %s", row.names(values3$df)))) + 
                          #  geom_point(alpha=0.3) +
                          #  geom_point(data=Selected_Genes, 
                          #             aes(logCPM, logFC, text=sprintf("Gene Name: %s", row.names(Selected_Genes))), 
                          #             color='green',
                          #             size=2)
                          
                          #GGTEST6
                          
                        })
                      }
                      
                      #Download Plots
                      {
                        GGVolcano <- values3$df %>%
                          ggplot(aes(logFC, -log10(PValue))) + 
                          geom_point(alpha=0.3) +
                          geom_point(data=Selected_Genes, 
                                     aes(logFC, -log10(PValue)), 
                                     color='green',
                                     size=2)
                        
                        GGMA <- values3$df %>%
                          ggplot(aes(logCPM, logFC)) + 
                          geom_point(alpha=0.3) +
                          geom_point(data=Selected_Genes, 
                                     aes(logCPM, logFC), 
                                     color='green',
                                     size=2)
                        
                        output$fullAnalysisOutput <- downloadHandler(
                          filename = function() {
                            paste0(input$Select_Comparison, ".zip")
                          },
                          content <- function(file) {
                            UserDir = tempdir()
                            write.csv(DF2, paste0(UserDir, "/Analysis_Table.csv"))
                            MyWidget = combineWidgets(ggplotly(GGVolcano), ggplotly(GGMA), ncol = 2)
                            saveWidget(MyWidget, paste0(UserDir, "/Analysis_Plots.html"))
                            zip(file, paste0(UserDir, "/Analysis_Table.csv"), paste0(UserDir, "/Analysis_Plots.html"), flags = "-r9Xj")
                          },
                          contentType = "application/zip")
                      }
                      
                    }
                  })
                }
                
                
              }
              
            })
          }
          
        } else {
          
          shiny::observeEvent(input$ShowAnalysis, {
            output$ResultTable <- DT::renderDataTable(
              {
                shiny::isolate(VALUES$CFSeq_Data[[SpeciesID]][[as.integer(substr(input$current_id,11,15))]]$Count_Table)
              },
              escape = F,
              rownames = TRUE,
              options = list(processing = FALSE)
            )
          })
          
        }
        
      })
  }
    
  #4. Reset Analysis Output
  {
      
      #When bottom reset button is present
      shiny::observeEvent(input$CLEAR_DATA1, {
        
        #Switch to Study View Tab
        {
          newtab2 <- switch(input$tabs,
                            "Run_Analysis" = "Choose_Studies"
          )
          updateTabItems(session, "tabs", newtab2) 
        }
        
        #Clear out data
        {
          observe({
            updateSelectInput(session, inputId = "Select_Comparison", "Choose Your Comparison of Interest", choices = c("[Select A Comparison]"),
                              selected = "[Select A Comparison]")
          })
          
          observe({
            updateSelectInput(session, inputId = "Find_Pathway", "Highlight a Pathway", choices = c("[Select a Pathway]"),
                              selected = "[Select a Pathway]")
          })
          
          observe({
            updateSelectizeInput(session, inputId = "Find_Gene", "Find a Gene", choices = c(Choose = ""),
                                 options = list(maxOptions = 15000))
          })
          
          output$ResultTable <- NULL
          output$VolcanoPlot <- NULL
          output$MAPlot <- NULL
          output$Comparison <-  NULL
          DF2 <- NULL
          Values3 <- NULL
          GeneList <- NULL
        }
        
      })
      
      #When side menu reset button is present
      shiny::observeEvent(input$CLEAR_DATA2, {
        
        #Switch to Study View Tab
        {
          newtab2 <- switch(input$tabs,
                            "Run_Analysis" = "Choose_Studies"
          )
          updateTabItems(session, "tabs", newtab2) 
        }
        
        #Clear out data
        {
          observe({
            updateSelectInput(session, inputId = "Select_Comparison", "Choose Your Comparison of Interest", choices = c("[Select A Comparison]"),
                              selected = "[Select A Comparison]")
          })
          
          observe({
            updateSelectInput(session, inputId = "Find_Pathway", "Highlight a Pathway", choices = c("[Select a Pathway]"),
                              selected = "[Select a Pathway]")
          })
          
          observe({
            updateSelectizeInput(session, inputId = "Find_Gene", "Find a Gene", choices = c(Choose = ""),
                                 options = list(maxOptions = 15000))
          })
          
          output$ResultTable <- NULL
          output$VolcanoPlot <- NULL
          output$MAPlot <- NULL
          output$
            Comparison <-  NULL
          DF2 <- NULL
          Values3 <- NULL
          GeneList <- NULL
        }
        
      })
      
    }
    
  #5. Dismiss Modal Dialog
  {
      shiny::observeEvent(input$dismiss_modal, {
        shiny::removeModal()
      })
    }
    
  #6. Download Data
  {
      #Download Count Table
      output$downloadCountDT <- downloadHandler(
        filename = function() {
          paste(VALUES$CFSeq_Data[[SpeciesID]][[as.integer(substr(input$current_id,11,15))]][[3]]$GEO.Accession, "_Count_Table", ".csv", sep = "")
        },
        content = function(file) {
          write.csv(CFSeq_Data[[SpeciesID]][[as.integer(substr(input$current_id,11,15))]][[1]], file, row.names = TRUE)
        }
      )
      
      #Download DE Analysis Table
      output$downloadDataDT <- downloadHandler(
        filename = function() {
          paste(input$Select_Comparison, "_Data_Table", ".csv", sep = "")
        },
        content = function(file) {
          write.csv(DF2, file, row.names = TRUE)
        }
      )  
      
      
    }
  
  #7. Reset Species
  {
    observeEvent(input$RESET, {
      observe({
        updateSelectizeInput(session, inputId = "Select_SubSpecies", "Choose Your Sub-Species of Interest", choices = c("[Select A Sub-Species]"),
                             selected = "[Select A Sub-Species]")
      })
      observe({
        updateSelectizeInput(session, inputId = "Select_Species", "Choose Your Species of Interest", choices = c("[Pick a Species]", names(CFSeq_Data)),
                             selected = "[Pick a Species]")
      })
      observe({
        updateSelectizeInput(session, inputId = "Select_Strain", "Choose Your Strain of Interest",
                             choices = c("[Select A Strain]"),
                             selected = "[Select A Strain]")
      })
      observe({
        updateSelectizeInput(session, inputId = "Select_Treatment", "Choose Your Treatment of Interest",
                             choices = c("[Select A Treatment]"),
                             selected = "[Select A Treatment]")
      })
      observe({
        updateSelectizeInput(session, inputId = "Select_Medium", "Choose Your Medium of Interest",
                             choices = c("[Select A Medium]"),
                             selected = "[Select A Medium]")
      })
      observe({
        updateSelectizeInput(session, inputId = "Select_Genotype", "Choose Your Gene Perturbation of Interest",
                             choices = c("[Select A Gene Perturbation]"),
                             selected = "[Select A Gene Perturbation]")
      })
      observe({
        updateSelectizeInput(session, inputId = "Find_Gene", "Find a Gene",
                             choices = c("[Select A Gene]"),
                             selected = "[Select A Gene]")
      })
      observe({
        updateSelectizeInput(session, inputId = "Find_Pathway", "Find a Pathway",
                             choices = c("[Select A Pathway]"),
                             selected = "[Select A Pathway]")
      })
      
      reset("current_id")
      reset("ShowAnalysis")
      reset("button")
      reset("dismiss_modal")

      #Re-initialize variables
      {
        VALUES <- reactiveValues(CFSeqStrain = NULL, CFSeqTreatment  = NULL, CFSeqMedium = NULL, CFSeqSubSpecies = NULL, DFSeqIntersect = NULL, s = NULL)
        output$table <- NULL
        output$ResultTable <- NULL
        output$VolcanoPlot <- NULL
        output$MAPlot <- NULL
        Strains <- vector()
        Treatments <- vector()
        Media <- vector()
        GeneList <- vector()
        GeneIDList <- vector()
        GeneTable <- NULL
        GeneFrame <- NULL
        Comparison <- NULL
        values = NULL
        values1 = NULL
        values2 = NULL
        values3 = NULL
        DF = NULL
        DF1 = NULL
        DF2 = NULL
        DF3 = NULL
        NMedia <- vector()
        UMedia <- vector()
        NStrains <- vector()
        UStrains <- vector()
        NTreatments <- vector()
        UTreatments <- vector()
        NGenotype <- vector()
        UGenotype <- vector()
        CompleteS <<- FALSE
        CompleteM <<- FALSE
        CompleteT <<- FALSE
        CompleteG <<- FALSE
        Approve = 0
        newtab = NULL
      }
      
    })
  }
  
}
}

#####

shinyApp(ui = ui, server = server)


#####

#Load Libraries
{
library(shiny)
library(shinydashboard)
library(DT) 
library(tidyverse)
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
    
    #Create manual window
    tabItem(tabName = "How_to_Use", 
            
      h1("CF-Seq User Guide"),
      br(),
      
      fluidRow(
        column(12,
            box(width = NULL, status = "warning",
               h2("Table of Contents"),
               p(" 1. Why Use CF-Seq?"),
               p(" 2. How CF-Seq Works"),
               p(" 3. Citing CF-Seq"),
               p(" 4. Launch the App")
        ))
      ), 
      
      fluidRow(
        column(12,
               box(width = NULL, status = "warning",
                  
                   h2("Why Use CF-Seq?"),
                   p(em("Note: Short on time? Scroll down to the end of the user manual to launch the app")),
                   p("The CF-Seq application was launched with a central goal in mind: to illuminate the landscape of prior published research on CF pathogens. Our aim is to take the motley collection of individual RNA-seq datasets stored in the Gene Expression Omnibus (GEO) – the go-to database for transcriptomics experiments – and show how they are connected. With this approach, we hope to make it more possible for any CF researcher to build on prior findings and generate new experimental hypotheses."),
                   p("With CF-Seq, you can identify how different bacterial strains, media, treatments, and gene perturbations have been employed across experiments on 14 CF pathogens. After choosing a study of interest, you can explore how each gene has been affected by the conditions of the experiment. And if there are multiple conditions – say multiple antibiotics, growth conditions, or media tested – for a single experiment, you can explore the effects of each condition independently. For many studies, you can also see how the genes on various biological pathways are up or downregulated to get a global picture of bacterial response to any given experiment. All of this can be done with zero computational experience!"),
                   p("CF-Seq currently includes a set of more than 150 CF pathogen studies. Over time we will expand not only the number of studies we include in the app, but the species we include (more CF pathogens and CF studies involving human cells and other model organism).")
                    
               ))
      ), 
      
      fluidRow(
        column(12,
               box(width = NULL, status = "warning",
                   
                   h2("How CF-Seq Works"),
                   p("Before you start using CF-Seq, make sure to read through the basic instructions below so you understand how to use it most effectively."),
                   p(em("[Viewing and Filtering Studies]")),
                   img(src = "StudyView.png", height = 600, width = "100%"),
                   br(),
                   br(),
                   p("Once you press the ‘Launch App’ button below, you will be transported to the ‘Study View’ panel of the app. This is where you can select a species and view all of its studies. Simply select your species of interest in the ‘Select a Species’ drop-down menu"),
                   p("Once you select a species, all of the studies for that CF pathogen will appear at the bottom of the screen. At this point, you can click the blue button next to any study, view more detailed metadata, and run analysis."),
                   p("But if you would like to filter available studies by experimental characteristics first, you can adjust any of the drop-down menus at the right-hand side of the app window. Once you finalize selections, you will see just those studies that meet your filtering criteria. At this point, you can change the filters, or click the blue button and run analysis on a study of interest"),
                   br(),
                   p(em("[Viewing Study Analysis]")),
                   img(src = "StudyAnalysis.png", height = 600, width = "100%"),
                   br(),
                   br(),
                   p("After pressing the run analysis button for any given study, you will be transported to the ‘Study Analysis’ window of the application. Here, you can see all of the experimental conditions that the CF pathogen was subjected to in the selected study, and pick two for which to compare gene expression and view differential expression analysis results. When you select your chosen comparison, and finalize that selection, you will be presented with a volcano plot and MA plot which show how the expression of each gene different between conditions (whether expression went up or down) and by how much. Furthermore, you will see a table of analysis with the P values, fold changes, and other metrics for each gene."),
                   p("Once a comparison is shown, you also have the option to highlight individual genes, genes that surpass a chosen statistical significance (P Value) or fold change, and genes on a selected KEGG biological pathway. Any of the plots and tables that you generate by filtering the data can be downloaded exactly as shown in the application (With any selected genes highlighted), as can the original study count table if you’d like to perform computational analysis yourself. When finished performing analysis for a particular study, you are welcome to check out another by pressing the ‘reset analysis’ button at the bottom of the ‘Study Analysis’ panel to be taken back to ‘Study View’ panel"),
                   br()
                   
               ))
      ), 
      
      

      fluidRow(
        column(12,
               box(width = NULL, status = "warning",
                   
                   h2("Citing CF-Seq"),
                   p("If you make use of the CF-Seq application to further your own experimental efforts any way – to explore prior findings, generate new hypotheses, or as inspiration for future data analysis applications, we kindly ask that you cite our associated paper:")
                   #PAPER CITATION HERE
                   
               ))
      ), 
      
      

      fluidRow(
        column(12,
               box(width = NULL, status = "warning",
                   
                   h2("Launch the App"),
                   p("Thanks for reading! Now go ahead and press the button below to start the application:"),
                   br(),
                   actionButton("Launch","Launch App", class = "btn-success btn-lg btn-block"),
                   
               ))
      )
          
    ),
            
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
                 actionButton("button", "Search"),
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
     
     
     fluidRow(id = "FA1",
       column(12,

              box(width = NULL, status = "warning",
              code(id = "No_Analysis", "Study design does not allow for differential expression analysis: 
              either there is only 1 replicate per condition or conditions have an unequal number of replicates.
              The study count table is still available to download if you wish to check it out for yourself"),
              DT::dataTableOutput("ComparisonTable", width = "100%"),
              br(),
              selectizeInput(inputId = "Select_Comparison", label =  "Select A Comparison",
                             choices = c("[Select A Comparison]"),
                             selected = "[Select A Comparison]"
              )),
              p(class = "text-muted", width = "100%", 
                "*Interpret comparison, [Comparison A] - [Comparison B], as follows: [Selected Gene] is upregulated/downregulated in [Comparison A]  vs. [Comparison B]",
                br(),
                "In other words, if [Selected Gene] has a fold change greater than zero, [Selected Gene] is upregulated in [Comparison A] vs. [Comparison B]",
                br(),
                "If [Selected Gene] has a fold change less than zero, [Selected Gene] is downregulated in [Comparison A] vs. [Comparison B]"
              ),
              br()
              )
     ),
     
     fluidRow(id = "FA2",
       column(8,
              tabsetPanel(id = "Analysis_Plots", type = "tabs",
                          tabPanel("Volcano Plot", plotlyOutput("VolcanoPlot")),
                          tabPanel("MA Plot", plotlyOutput("MAPlot"))
       )),
       
       column(4,
              br(),
              br(),
              box(width = NULL, status = "warning",
                  
                  selectizeInput(inputId = "Choose_P", label =  "Choose a P Value Cutoff",
                                 choices = c("[Select a Cutoff]", 0.001, 0.01, 0.05, 0.10, 0.15, 0.20),
                                 selected = "[Select a Cutoff]"
                  ),
                  selectizeInput(inputId = "Choose_FC", label =  "Choose a Fold Change Cutoff",
                                 choices = c("[Select a Cutoff]", 2, 3, 4, 5, 10, 50),
                                 selected = "[Select a Cutoff]"
                  ),
                  selectizeInput(inputId = "Find_Gene", label =  "Find a Gene",
                                 choices = c("[Select a Gene]"),
                                 selected = "[Select a Gene]"
                  ),
                  selectizeInput(inputId = "Find_Pathway", label =  "Find a Pathway",
                                 choices = c("[Select a Pathway]"),
                                 selected = "[Select a Pathway]"
                  ),
                  actionButton("RESET_ANALYSIS_OUTPUTS", "Reset Analysis Outputs"),
                  br(),
                  actionButton("View_Metadata", "View Study Metadata"),
                  downloadButton("downloadCountDT", "Download Study Count Table"),
                  downloadButton("fullAnalysisOutput", "Download All Analysis Data"),
                  actionButton("CLEAR_DATA2", "Reset Analysis")
              ))
     ),
     

     
     
     fluidRow(id = "FA3",
       column(12,
              DT::dataTableOutput("ResultTable", width = "100%")
       )
     ),
     
     fluidRow(id = "FA4",
       column(12,
              br(),
              box(width = NULL, status = "warning",
                  actionButton("View_Metadata2", "View Study Metadata"),
                  downloadButton("downloadCountDT2", "Download Study Count Table"),
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
  
  #1. Switch from user guide to study view 
  {
    
    shiny::observeEvent(input$Launch, {
      
      Newtab <- switch(input$tabs,
                       "How_to_Use" = "Choose_Studies"
                       
      )
      updateTabItems(session, "tabs", Newtab) 
      
    })
    
  }
  
  #2. Select a species and generate species-specific data
  {
  observeEvent(input$Select_Species, {
    
    #A. Initialize variables
    {
      VALUES <<- reactiveValues(CFSeqStrain = NULL, CFSeqTreatment  = NULL, CFSeqMedium = NULL, CFSeqGenotype = NULL, CFSeqSubSpecies = NULL, CFSeq_Data = CFSeq_Data, S = NULL) 
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
    } else if (input$Select_Species == "Bacteroides_Species") {
      VALUES$S = 2
    } else if (input$Select_Species == "Burkholderia_Species") {
      VALUES$S = 3
    } else if (input$Select_Species == "Candida_albicans") {
      VALUES$S = 4
    } else if (input$Select_Species == "Clostridium_difficile") {
      VALUES$S = 5
    } else if (input$Select_Species == "Fusobacterium_nucleatum") {
      VALUES$S = 6
    } else if (input$Select_Species == "Haemophilus_influenzae") {
      VALUES$S = 7
    } else if (input$Select_Species == "Mycobacterium_abscessus") {
      VALUES$S = 8
    } else if (input$Select_Species == "Pseudomonas_aeruginosa") {
      VALUES$S = 9
    } else if (input$Select_Species == "Porphyromonas_Species") {
      VALUES$S = 10
    } else if (input$Select_Species == "Staphylococcus_aureus") {
      VALUES$S = 11
    } else if (input$Select_Species == "Stenotrophomonas_maltophilia") {
      VALUES$S = 12
    } else if (input$Select_Species == "Streptococcus_Species") {
      VALUES$S = 13
    } else {
      VALUES$S = NULL
    }
    
    if (is.null(VALUES$S) != TRUE) {
    
    #Set Strains values
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
      Strains[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Strain
    }
    Strains <- na.omit(unique(Strains))
    if (length(Strains) != 0) {
      Strains <- unlist(strsplit(Strains, ","))
    }

    if (length(Strains) != 0) {
      for (i in 1:length(Strains)) {
        if (substr(Strains[i], 1, 1) == " ") {
          Strains[i] <- substr(Strains[i], 2, nchar(Strains[i])) 
        } 
        if (substr(Strains[i], nchar(Strains[i]), nchar(Strains[i])) == " ") {
          Strains[i] <- substr(Strains[i], 1, (nchar(Strains[i]) - 1)) 
        } 
      }
    }
    
    #Set Treatments values
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
      Treatments[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Treatment
    }
    Treatments <- na.omit(unique(Treatments))
    if (length(Treatments) != 0) {
      Treatments <- unlist(strsplit(Treatments, ","))
    }
    
    if (length(Treatments) != 0) {
      for (i in 1:length(Treatments)) {
        if (substr(Treatments[i], 1, 1) == " ") {
          Treatments[i] <- substr(Treatments[i], 2, nchar(Treatments[i])) 
        } 
        if (substr(Treatments[i], nchar(Treatments[i]), nchar(Treatments[i])) == " ") {
          Treatments[i] <- substr(Treatments[i], 1, (nchar(Treatments[i]) - 1)) 
        } 
      }
    }
    
    #Set Media values
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
      Media[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Medium
    }
    Media <- na.omit(unique(Media))
    if (length(Media) != 0) {
      Media <- unlist(strsplit(Media, ","))
    }
    
    if (length(Media) != 0) {
    for (i in 1:length(Media)) {
      if (substr(Media[i], 1, 1) == " ") {
        Media[i] <- substr(Media[i], 2, nchar(Media[i])) 
      } 
      if (substr(Media[i], nchar(Media[i]), nchar(Media[i])) == " ") {
        Media[i] <- substr(Media[i], 1, (nchar(Media[i]) - 1)) 
      } 
    }
    }
    
    #Set Gene perturbation values
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
      Genotype[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Genotype
    }
    Genotype <- na.omit(unique(Genotype))
    if (length(Genotype) != 0) {
      Genotype <- unlist(strsplit(Genotype, ","))
    }
    
    if (length(Genotype) != 0) {
    for (i in 1:length(Genotype)) {
      if (substr(Genotype[i], 1, 1) == " ") {
        Genotype[i] <- substr(Genotype[i], 2, nchar(Genotype[i])) 
      } 
      if (substr(Genotype[i], nchar(Genotype[i]), nchar(Genotype[i])) == " ") {
        Genotype[i] <- substr(Genotype[i], 1, (nchar(Genotype[i]) - 1)) 
      } 
    }
    }
    
    #Set Sub-Species values
    
    for (i in 1:(length(VALUES$CFSeq_Data[[VALUES$S]]))) {
    if (is.null(VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Species) != TRUE) {
      Species2[i] = VALUES$CFSeq_Data[[VALUES$S]][[i]][[3]]$Species
      Species2 <- na.omit(unique(Species2))
      if (length(Species2) != 0) {
        Species2 <- unlist(strsplit(Species2, ","))
        for (i in 1:length(Species2)) {
          if (substr(Species2[i], 1, 1) == " ") {
            Species2[i] <- substr(Species2[i], 2, nchar(Species2[i])) 
          } 
          if (substr(Species2[i], nchar(Species2[i]), nchar(Species2[i])) == " ") {
            Species2[i] <- substr(Species2[i], 1, (nchar(Species2[i]) - 1)) 
          } 
        }
      }
    } else {
      Species2 <- NULL
    }
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
      
      SetStrain <- "Fresh"
      SetMedium <- "Fresh"
      SetGene <- "Fresh"
      SetTreatment <- "Fresh"
      SetSubs <- "Fresh" #Subway... Eat Fresh...
      
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
      
        if (is.na(UTreatments) != TRUE) {
          for (i in 1:length(UTreatments)) {
            if (substr(UTreatments[i], 1, 1) == " ") {
              UTreatments[i] <- substr(UTreatments[i], 2, nchar(UTreatments[i])) 
            } 
            if (substr(UTreatments[i], nchar(UTreatments[i]), nchar(UTreatments[i])) == " ") {
              UTreatments[i] <- substr(UTreatments[i], 1, (nchar(UTreatments[i]) - 1)) 
            } 
          }
        }
        
        if(sum(is.na(SetTreatment)) != length(SetTreatment)) {
          if (length(SetTreatment)!= 0) {
            if (SetTreatment != "Fresh") {
              for (i in 1:length(UTreatments)) {
                NTreatments <- UTreatments[UTreatments %in% SetTreatment]
              }
            } else {
              NTreatments <- UTreatments
            }
          } else {
            NTreatments <- NULL
          }
        } else {
          NTreatments <- NULL
        }
  
        
        
        

        UMedia <- c()
        NMedia <- c()
        for (i in 1:length(VALUES$CFSeqStrain)) {
          UMedia[i] <- VALUES$CFSeqStrain[[i]][[3]]$Medium
        }
        if(length(na.omit(UMedia)) != 0) {
          UMedia <- unlist(strsplit(na.omit(UMedia), ","))
        }
        
        if (is.na(UMedia) != TRUE) {
          for (i in 1:length(UMedia)) {
            if (substr(UMedia[i], 1, 1) == " ") {
              UMedia[i] <- substr(UMedia[i], 2, nchar(UMedia[i])) 
            } 
            if (substr(UMedia[i], nchar(UMedia[i]), nchar(UMedia[i])) == " ") {
              UMedia[i] <- substr(UMedia[i], 1, (nchar(UMedia[i]) - 1)) 
            } 
          }
        }
        
        if(sum(is.na(SetMedium)) != length(SetMedium)) {
          if (length(SetMedium)!= 0) {
            if (SetMedium != "Fresh") {
              for (i in 1:length(UMedia)) {
                NMedia <- UMedia[UMedia %in% SetMedium]
              }
            } else {
              NMedia <- UMedia
            }
          } else {
            NMedia <- NULL
          }
        } else {
          NMedia <- NULL
        }
        
        

        UGenotype <- c()
        NGenotype <- c()
        for (i in 1:length(VALUES$CFSeqStrain)) {
          UGenotype[i] <- VALUES$CFSeqStrain[[i]][[3]]$Genotype
        }
        if(length(na.omit(UGenotype)) != 0) {
          UGenotype <- unlist(strsplit(na.omit(UGenotype), ","))
        }
      
        if (is.na(UGenotype) != TRUE) {
          for (i in 1:length(UGenotype)) {
            if (substr(UGenotype[i], 1, 1) == " ") {
              UGenotype[i] <- substr(UGenotype[i], 2, nchar(UGenotype[i])) 
            } 
            if (substr(UGenotype[i], nchar(UGenotype[i]), nchar(UGenotype[i])) == " ") {
              UGenotype[i] <- substr(UGenotype[i], 1, (nchar(UGenotype[i]) - 1)) 
            } 
          }
        }
        
        if(sum(is.na(SetGene)) != length(SetGene)) {
          if (length(SetGene)!= 0) {
            if (SetGene != "Fresh") {
              for (i in 1:length(UGenotype)) {
                NGenotype <- UGenotype[UGenotype %in% SetGene]
              }
            } else {
              NGenotype <- UGenotype
            }
          } else {
            NGenotype <- NULL
          }
        } else {
          NGenotype <- NULL
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
          
          if (is.na(USubSpecies) != TRUE) {
          for (i in 1:length(USubSpecies)) {
            if (substr(USubSpecies[i], 1, 1) == " ") {
              USubSpecies[i] <- substr(USubSpecies[i], 2, nchar(USubSpecies[i])) 
            } 
            if (substr(USubSpecies[i], nchar(USubSpecies[i]), nchar(USubSpecies[i])) == " ") {
              USubSpecies[i] <- substr(USubSpecies[i], 1, (nchar(USubSpecies[i]) - 1)) 
            } 
          }
          }
          
          if(sum(is.na(SetSubs)) != length(SetSubs)) {
            if (length(SetSubs)!= 0) {
              if (SetSubs != "Fresh") {
                for (i in 1:length(USubSpecies)) {
                  NSubSpecies <- USubSpecies[USubSpecies %in% SetSubs]
                }
              } else {
                NSubSpecies <- USubSpecies
              }
            } else {
              NSubSpecies <- NULL
            }
          } else {
            NSubSpecies <- NULL
          }
            
          SetSubs <<- NSubSpecies
          
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
          
        
        if (input$Select_Strain != "[Select A Strain]") {
          CompleteS <<- TRUE
        }
      
        
        StrainSelect <<- input$Select_Strain
        
        SetTreatment <<- NTreatments
        SetMedium <<- NMedia
        SetGene <<- NGenotype
        
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
              
              if (is.na(UStrains) != TRUE) {
              for (i in 1:length(UStrains)) {
                if (substr(UStrains[i], 1, 1) == " ") {
                  UStrains[i] <- substr(UStrains[i], 2, nchar(UStrains[i])) 
                } 
                if (substr(UStrains[i], nchar(UStrains[i]), nchar(UStrains[i])) == " ") {
                  UStrains[i] <- substr(UStrains[i], 1, (nchar(UStrains[i]) - 1)) 
                } 
              }
              }
              
              if(sum(is.na(SetStrain)) != length(SetStrain)) {
                if (length(SetStrain)!= 0) {
                  if (SetStrain != "Fresh") {
                    for (i in 1:length(UStrains)) {
                      NStrains <- UStrains[UStrains %in% SetStrain]
                    }
                  } else {
                    NStrains <- UStrains
                  }
                } else {
                  NStrains <- NULL
                }
              } else {
                NStrains <- NULL
              }
              
              

              UMedia <- c()
              NMedia <- c()
              for (i in 1:length(VALUES$CFSeqTreatment)) {
                UMedia[i] <- VALUES$CFSeqTreatment[[i]][[3]]$Medium
              }
              if(length(na.omit(UMedia)) != 0) {
                UMedia <- unlist(strsplit(na.omit(UMedia), ","))
              }
              
              if (is.na(UMedia) != TRUE) {
                for (i in 1:length(UMedia)) {
                  if (substr(UMedia[i], 1, 1) == " ") {
                    UMedia[i] <- substr(UMedia[i], 2, nchar(UMedia[i])) 
                  } 
                  if (substr(UMedia[i], nchar(UMedia[i]), nchar(UMedia[i])) == " ") {
                    UMedia[i] <- substr(UMedia[i], 1, (nchar(UMedia[i]) - 1)) 
                  } 
                }
              }
              
              if(sum(is.na(SetMedium)) != length(SetMedium)) {
                if (length(SetMedium)!= 0) {
                  if (SetMedium != "Fresh") {
                    for (i in 1:length(UMedia)) {
                      NMedia <- UMedia[UMedia %in% SetMedium]
                    }
                  } else {
                    NMedia <- UMedia
                  }
                } else {
                  NMedia <- NULL
                }
              } else {
                NMedia <- NULL
              }
              
              

              UGenotype <- c()
              NGenotype <- c()
              for (i in 1:length(VALUES$CFSeqTreatment)) {
                UGenotype[i] <- VALUES$CFSeqTreatment[[i]][[3]]$Genotype
              }
              if(length(na.omit(UGenotype)) != 0) {
                UGenotype <- unlist(strsplit(na.omit(UGenotype), ","))
              }
              
              if (is.na(UGenotype) != TRUE) {
              for (i in 1:length(UGenotype)) {
                if (substr(UGenotype[i], 1, 1) == " ") {
                  UGenotype[i] <- substr(UGenotype[i], 2, nchar(UGenotype[i])) 
                } 
                if (substr(UGenotype[i], nchar(UGenotype[i]), nchar(UGenotype[i])) == " ") {
                  UGenotype[i] <- substr(UGenotype[i], 1, (nchar(UGenotype[i]) - 1)) 
                } 
              }
              }
              
              if(sum(is.na(SetGene)) != length(SetGene)) {
                if (length(SetGene)!= 0) {
                  if (SetGene != "Fresh") {
                    for (i in 1:length(UGenotype)) {
                      NGenotype <- UGenotype[UGenotype %in% SetGene]
                    }
                  } else {
                    NGenotype <- UGenotype
                  }
                } else {
                  NGenotype <- NULL
                }
              } else {
                NGenotype <- NULL
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
              
              if (is.na(USubSpecies) != TRUE) {
                for (i in 1:length(USubSpecies)) {
                  if (substr(USubSpecies[i], 1, 1) == " ") {
                    USubSpecies[i] <- substr(USubSpecies[i], 2, nchar(USubSpecies[i])) 
                  } 
                  if (substr(USubSpecies[i], nchar(USubSpecies[i]), nchar(USubSpecies[i])) == " ") {
                    USubSpecies[i] <- substr(USubSpecies[i], 1, (nchar(USubSpecies[i]) - 1)) 
                  } 
                }
              }
              
              if(sum(is.na(SetSubs)) != length(SetSubs)) {
                if (length(SetSubs)!= 0) {
                  if (SetSubs != "Fresh") {
                    for (i in 1:length(USubSpecies)) {
                      NSubSpecies <- USubSpecies[USubSpecies %in% SetSubs]
                    }
                  } else {
                    NSubSpecies <- USubSpecies
                  }
                } else {
                  NSubSpecies <- NULL
                }
              } else {
                NSubSpecies <- NULL
              }
              
              
              SetSubs <<- NSubSpecies
              
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
          
          SetStrain <<- NStrains
          SetMedium <<- NMedia
          SetGene <<- NGenotype
          
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
            
            if (is.na(UStrains) != TRUE) {
              for (i in 1:length(UStrains)) {
                if (substr(UStrains[i], 1, 1) == " ") {
                  UStrains[i] <- substr(UStrains[i], 2, nchar(UStrains[i])) 
                } 
                if (substr(UStrains[i], nchar(UStrains[i]), nchar(UStrains[i])) == " ") {
                  UStrains[i] <- substr(UStrains[i], 1, (nchar(UStrains[i]) - 1)) 
                } 
              }
            }
            

            if(sum(is.na(SetStrain)) != length(SetStrain)) {
              if (length(SetStrain)!= 0) {
                if (SetStrain != "Fresh") {
                  for (i in 1:length(UStrains)) {
                    NStrains <- UStrains[UStrains %in% SetStrain]
                  }
                } else {
                  NStrains <- UStrains
                }
              } else {
                NStrains <- NULL
              }
            } else {
              NStrains <- NULL
            }
            
            
            UTreatments <- c()
            NTreatments <- c()
            for (i in 1:length(VALUES$CFSeqMedium)) {
              UTreatments[i] <- VALUES$CFSeqMedium[[i]][[3]]$Treatment
            }
            if(length(na.omit(UTreatments)) != 0) {
              UTreatments <- unlist(strsplit(na.omit(UTreatments), ","))
            }
            
            if (is.na(UTreatments) != TRUE) {
              for (i in 1:length(UTreatments)) {
                if (substr(UTreatments[i], 1, 1) == " ") {
                  UTreatments[i] <- substr(UTreatments[i], 2, nchar(UTreatments[i])) 
                } 
                if (substr(UTreatments[i], nchar(UTreatments[i]), nchar(UTreatments[i])) == " ") {
                  UTreatments[i] <- substr(UTreatments[i], 1, (nchar(UTreatments[i]) - 1)) 
                } 
              }
            }
            
            if(sum(is.na(SetTreatment)) != length(SetTreatment)) {
              if (length(SetTreatment)!= 0) {
                if (SetTreatment != "Fresh") {
                  for (i in 1:length(UTreatments)) {
                    NTreatments <- UTreatments[UTreatments %in% SetTreatment]
                  }
                } else {
                  NTreatments <- UTreatments
                }
              } else {
                NTreatments <- NULL
              }
            } else {
              NTreatments <- NULL
            }
            
            
            
            UGenotype <- c()
            NGenotype <- c()
            for (i in 1:length(VALUES$CFSeqMedium)) {
              UGenotype[i] <- VALUES$CFSeqMedium[[i]][[3]]$Genotype
            }
            if(length(na.omit(UGenotype)) != 0) {
              UGenotype <- unlist(strsplit(na.omit(UGenotype), ","))
            }
            
            if (is.na(UGenotype) != TRUE) {
            for (i in 1:length(UGenotype)) {
              if (substr(UGenotype[i], 1, 1) == " ") {
                UGenotype[i] <- substr(UGenotype[i], 2, nchar(UGenotype[i])) 
              } 
              if (substr(UGenotype[i], nchar(UGenotype[i]), nchar(UGenotype[i])) == " ") {
                UGenotype[i] <- substr(UGenotype[i], 1, (nchar(UGenotype[i]) - 1)) 
              } 
            }
            }
            
            if(sum(is.na(SetGene)) != length(SetGene)) {
              if (length(SetGene)!= 0) {
                if (SetGene != "Fresh") {
                  for (i in 1:length(UGenotype)) {
                    NGenotype <- UGenotype[UGenotype %in% SetGene]
                  }
                } else {
                  NGenotype <- UGenotype
                }
              } else {
                NGenotype <- NULL
              }
            } else {
              NGenotype <- NULL
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
            
            if (is.na(USubSpecies) != TRUE) {
              for (i in 1:length(USubSpecies)) {
                if (substr(USubSpecies[i], 1, 1) == " ") {
                  USubSpecies[i] <- substr(USubSpecies[i], 2, nchar(USubSpecies[i])) 
                } 
                if (substr(USubSpecies[i], nchar(USubSpecies[i]), nchar(USubSpecies[i])) == " ") {
                  USubSpecies[i] <- substr(USubSpecies[i], 1, (nchar(USubSpecies[i]) - 1)) 
                } 
              }
            }
            
            if(sum(is.na(SetSubs)) != length(SetSubs)) {
              if (length(SetSubs)!= 0) {
                if (SetSubs != "Fresh") {
                  for (i in 1:length(USubSpecies)) {
                    NSubSpecies <- USubSpecies[USubSpecies %in% SetSubs]
                  }
                } else {
                  NSubSpecies <- USubSpecies
                }
              } else {
                NSubSpecies <- NULL
              }
            } else {
              NSubSpecies <- NULL
            }
            
            
            SetSubs <<- NSubSpecies
            
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
          
          SetTreatment <<- NTreatments
          SetStrain <<- NStrains
          SetGene <<- NGenotype
          
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
          
          if (is.na(UStrains) != TRUE) {
            for (i in 1:length(UStrains)) {
              if (substr(UStrains[i], 1, 1) == " ") {
                UStrains[i] <- substr(UStrains[i], 2, nchar(UStrains[i])) 
              } 
              if (substr(UStrains[i], nchar(UStrains[i]), nchar(UStrains[i])) == " ") {
                UStrains[i] <- substr(UStrains[i], 1, (nchar(UStrains[i]) - 1)) 
              } 
            }
          }

          if(sum(is.na(SetStrain)) != length(SetStrain)) {
            if (length(SetStrain)!= 0) {
              if (SetStrain != "Fresh") {
                for (i in 1:length(UStrains)) {
                  NStrains <- UStrains[UStrains %in% SetStrain]
                }
              } else {
                NStrains <- UStrains
              }
            } else {
              NStrains <- NULL
            }
          } else {
            NStrains <- NULL
          }
          
          
          
          UTreatments <- c()
          NTreatments <- c()
          for (i in 1:length(VALUES$CFSeqGenotype)) {
            UTreatments[i] <- VALUES$CFSeqGenotype[[i]][[3]]$Treatment
          }
          if(length(na.omit(UTreatments)) != 0) {
            UTreatments <- unlist(strsplit(na.omit(UTreatments), ","))
          }
          
          if (is.na(UTreatments) != TRUE) {
            for (i in 1:length(UTreatments)) {
              if (substr(UTreatments[i], 1, 1) == " ") {
                UTreatments[i] <- substr(UTreatments[i], 2, nchar(UTreatments[i])) 
              } 
              if (substr(UTreatments[i], nchar(UTreatments[i]), nchar(UTreatments[i])) == " ") {
                UTreatments[i] <- substr(UTreatments[i], 1, (nchar(UTreatments[i]) - 1)) 
              }
            }
          }
          
          if(sum(is.na(SetTreatment)) != length(SetTreatment)) {
            if (length(SetTreatment)!= 0) {
              if (SetTreatment != "Fresh") {
                for (i in 1:length(UTreatments)) {
                  NTreatments <- UTreatments[UTreatments %in% SetTreatment]
                }
              } else {
                NTreatments <- UTreatments
              }     
            } else {
              NTreatments <- NULL
            }
          } else {
            NTreatments <- NULL
          }
          
          
          
          UMedia <- c()
          NMedia <- c()
          for (i in 1:length(VALUES$CFSeqGenotype)) {
            UMedia[i] <- VALUES$CFSeqGenotype[[i]][[3]]$Medium
          }
          if(length(na.omit(UMedia)) != 0) {
            UMedia <- unlist(strsplit(na.omit(UMedia), ","))
          }
          
          if (is.na(UMedia) != TRUE) {
            for (i in 1:length(UMedia)) {
              if (substr(UMedia[i], 1, 1) == " ") {
                UMedia[i] <- substr(UMedia[i], 2, nchar(UMedia[i])) 
              } 
              if (substr(UMedia[i], nchar(UMedia[i]), nchar(UMedia[i])) == " ") {
                UMedia[i] <- substr(UMedia[i], 1, (nchar(UMedia[i]) - 1)) 
              }
            }
          }
          
          if(sum(is.na(SetMedium)) != length(SetMedium)) {
            if (length(SetMedium)!= 0) {
              if (SetMedium != "Fresh") {
                for (i in 1:length(UMedia)) {
                  NMedia <- UMedia[UMedia %in% SetMedium]
                }
              } else {
                NMedia <- UMedia
              }
            } else {
              NMedia <- NULL
            }
          } else {
            NMedia <- NULL
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
          
          if (is.na(USubSpecies) != TRUE) {
            for (i in 1:length(USubSpecies)) {
              if (substr(USubSpecies[i], 1, 1) == " ") {
                USubSpecies[i] <- substr(USubSpecies[i], 2, nchar(USubSpecies[i])) 
              } 
              if (substr(USubSpecies[i], nchar(USubSpecies[i]), nchar(USubSpecies[i])) == " ") {
                USubSpecies[i] <- substr(USubSpecies[i], 1, (nchar(USubSpecies[i]) - 1)) 
              }
            }
          }
          
          if(sum(is.na(SetSubs)) != length(SetSubs)) {
            if (length(SetSubs)!= 0) {
              if (SetSubs != "Fresh") {
                for (i in 1:length(USubSpecies)) {
                  NSubSpecies <- USubSpecies[USubSpecies %in% SetSubs]
                }
              } else {
                NSubSpecies <- USubSpecies
              }
            } else {
              NSubSpecies <- NULL
            }
          } else {
            NSubSpecies <- NULL
          }
            
          
          
          SetSubs <<- NSubSpecies
          
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
        
        SetTreatment <<- NTreatments
        SetStrain <<- NStrains
        SetMedium <<- NMedia

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
              
              if (is.na(UStrains) != TRUE) {
                for (i in 1:length(UStrains)) {
                  if (substr(UStrains[i], 1, 1) == " ") {
                    UStrains[i] <- substr(UStrains[i], 2, nchar(UStrains[i])) 
                  } 
                  if (substr(UStrains[i], nchar(UStrains[i]), nchar(UStrains[i])) == " ") {
                    UStrains[i] <- substr(UStrains[i], 1, (nchar(UStrains[i]) - 1)) 
                  }
                }
              }
              
              if(sum(is.na(SetStrain)) != length(SetStrain)) {
                if (length(SetStrain)!= 0) {
                  if (SetStrain != "Fresh") {
                    for (i in 1:length(UStrains)) {
                      NStrains <- UStrains[UStrains %in% SetStrain]
                    }
                  } else {
                    NStrains <- UStrains
                  }
                } else {
                  NStrains <- NULL
                }    
              } else {
                NStrains <- NULL
              }
              
              
              
              UTreatments <- c()
              NTreatments <- c()
              
              for (i in 1:length(VALUES$CFSeqSubSpecies)) {
                UTreatments[i] <- VALUES$CFSeqSubSpecies[[i]][[3]]$Treatment
              }
              if(length(na.omit(UTreatments)) != 0) {
                UTreatments <- unlist(strsplit(na.omit(UTreatments), ","))
              }
              
              if (is.na(UTreatments) != TRUE) {
                for (i in 1:length(UTreatments)) {
                  if (substr(UTreatments[i], 1, 1) == " ") {
                    UTreatments[i] <- substr(UTreatments[i], 2, nchar(UTreatments[i])) 
                  } 
                  if (substr(UTreatments[i], nchar(UTreatments[i]), nchar(UTreatments[i])) == " ") {
                    UTreatments[i] <- substr(UTreatments[i], 1, (nchar(UTreatments[i]) - 1)) 
                  }
                }
              }
              
              if(sum(is.na(SetTreatment)) != length(SetTreatment)) {
                if (length(SetTreatment)!= 0) {
                  if (SetTreatment != "Fresh") {
                    for (i in 1:length(UTreatments)) {
                      NTreatments <- UTreatments[UTreatments %in% SetTreatment]
                    }
                  } else {
                    NTreatments <- UTreatments
                  }
                } else {
                  NTreatments <- NULL
                }      
              } else {
                NTreatments <- NULL
              }
              
              
              UMedia <- c()
              NMedia <- c()
              for (i in 1:length(VALUES$CFSeqSubSpecies)) {
                UMedia[i] <- VALUES$CFSeqSubSpecies[[i]][[3]]$Medium
              }
              if(length(na.omit(UMedia)) != 0) {
                UMedia <- unlist(strsplit(na.omit(UMedia), ","))
              }
              
              if (is.na(UMedia) != TRUE) {
                for (i in 1:length(UMedia)) {
                  if (substr(UMedia[i], 1, 1) == " ") {
                    UMedia[i] <- substr(UMedia[i], 2, nchar(UMedia[i])) 
                  } 
                  if (substr(UMedia[i], nchar(UMedia[i]), nchar(UMedia[i])) == " ") {
                    UMedia[i] <- substr(UMedia[i], 1, (nchar(UMedia[i]) - 1)) 
                  }
                }
              }
              
              if(sum(is.na(SetMedium)) != length(SetMedium)) {
                if (length(SetMedium)!= 0) {
                  if (SetMedium != "Fresh") {
                    for (i in 1:length(UMedia)) {
                      NMedia <- UMedia[UMedia %in% SetMedium]
                    }
                  } else {
                    NMedia <- UMedia
                  }
                } else {
                  NMedia <- NULL
                }
              } else {
                NMedia <- NULL
              }
                
              
              UGenotype <- c()
              NGenotype <- c()
              for (i in 1:length(VALUES$CFSeqSubSpecies)) {
                UGenotype[i] <- VALUES$CFSeqSubSpecies[[i]][[3]]$Genotype
              }
              if(length(na.omit(UGenotype)) != 0) {
                UGenotype <- unlist(strsplit(na.omit(UGenotype), ","))
              }
              
              if (is.na(UGenotype) != TRUE) {
              for (i in 1:length(UGenotype)) {
                if (substr(UGenotype[i], 1, 1) == " ") {
                  UGenotype[i] <- substr(UGenotype[i], 2, nchar(UGenotype[i])) 
                } 
                if (substr(UGenotype[i], nchar(UGenotype[i]), nchar(UGenotype[i])) == " ") {
                  UGenotype[i] <- substr(UGenotype[i], 1, (nchar(UGenotype[i]) - 1)) 
                }
              }
              }
              
              if(sum(is.na(SetGene)) != length(SetGene)) {
                if (length(SetGene)!= 0) {
                  if (SetGene != "Fresh") {
                    for (i in 1:length(UGenotype)) {
                      NGenotype <- UGenotype[UGenotype %in% SetGene]
                    }
                  } else {
                    NGenotype <- UGenotype
                  }
                } else {
                  NGenotype <- NULL
                }
              } else {
                NGenotype <- NULL
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
            
            SetTreatment <<- NTreatments
            SetStrain <<- NStrains
            SetMedium <<- NMedia
            SetGene <<- NGenotype
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
          SetStrain <<- "Fresh"
          SetTreatment <<- "Fresh"
          SetGene <<- "Fresh"
          SetMedium <<- "Fresh"
          SetSubs <<- "Fresh"
          
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
    
  #3. Action Button: More Details Dialogue + Run Analysis Window
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
                    paste("Strain(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Strain)
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
                    paste("Treatment(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Treatment)
                  })
                ),
                br(),
                br(),
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Gene Perturbation(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Genotype)
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
                tags$a(href = FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Link, "Study Link - Go to Gene Expression Omnibus (GEO) Record"),
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
    
  #4. Press Show Analysis Button --> Show Analysis Screen with default comparison 
  {
      shiny::observeEvent(input$ShowAnalysis, {
        
        output$VolcanoPlot <- NULL
        output$MAPlot <- NULL
        
        #Update Comparison list and pathways for specific study 
        {
          #Comparisons (Again)
          Comparison <- names(FilterInput[[as.integer(substr(input$current_id,11,15))]][4:length(FilterInput[[as.integer(substr(input$current_id,11,15))]])])
          Comparison <- Comparison[grep("Full Results", Comparison)]
          Comparison <- gsub(" Full Results", "", Comparison)

          #Pathways
          PathwaySelect <- "Empty"
          TestGenes <- FilterInput[[as.integer(substr(input$current_id,11,15))]][[1]][,1][1:50]

          for (j in 1:length(TestGenes)) {
            for (i in 1:length(GenePathway_Data)) {
              if (length(grep(TestGenes[j], GenePathway_Data[[i]]$Gene)) != 0) {
                PathwaySelect <- unique(GenePathway_Data[[i]]$Pathway)
                PathwayFile <- GenePathway_Data[[i]]
              }
            }
          }

          observe({
            updateSelectInput(session, inputId = "Select_Comparison", "Choose Your Comparison of Interest*", choices = c("[Select A Comparison]", Comparison),
                              selected = "[Select A Comparison]")
            updateSelectInput(session, inputId = "Find_Pathway", "Highlight a Pathway", choices = c("[Select a Pathway]", PathwaySelect),
                              selected = "[Select a Pathway]")
          })
          
          
        }
        
        #Show / Hide outputs based on study type
        {
          if (length(Comparison) == 0) {
            
            shinyjs::hide(id = "Select_Comparison")
            shinyjs::hide(id = "ComparisonTable")
            shinyjs::show(id = "No_Analysis")
            shinyjs::hide(id = "FA2")
            shinyjs::hide(id = "FA3")
            shinyjs::show(id = "View_Metadata2")
            shinyjs::show(id = "downloadCountDT2")
            
          } else {
            
            shinyjs::hide(id = "No_Analysis")
            shinyjs::show(id = "ComparisonTable")
            shinyjs::show(id = "Select_Comparison")
            shinyjs::show(id = "FA2")
            shinyjs::show(id = "FA3")
            shinyjs::show(id = "Find_Gene")
            shinyjs::show(id = "Find_Pathway")
            shinyjs::show(id = "downloadDataDT")
            shinyjs::show(id = "downloadDataVP")
            shinyjs::show(id = "downloadDataMA")
            shinyjs::show(id = "Analysis_Plots")
            shinyjs::show(id = "CLEAR_DATA1")
            shinyjs::hide(id = "CLEAR_DATA2")
            shinyjs::hide(id = "View_Metadata2")
            shinyjs::hide(id = "downloadCountDT2")
            
          }
          
          if (PathwaySelect == "Empty") {
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
        
        #Show study metadata
        {
          observeEvent(input$View_Metadata, {
            
            showModal(modalDialog(
              title = "Study Details",
              div(
                class = "text-center",
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Strain(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Strain)
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
                    paste("Treatment(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Treatment)
                  })
                ),
                br(),
                br(),
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Gene Perturbation(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Genotype)
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
                tags$a(href = FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Link, "Study Link - Go to Gene Expression Omnibus (GEO) Record"),
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
                  inputId = "dismiss_modal",
                  label = "Close",
                  class = "btn-danger"
                )
              )
            ))
            
          })
          
          observeEvent(input$View_Metadata2, {
            
            showModal(modalDialog(
              title = "Study Details",
              div(
                class = "text-center",
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Strain(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Strain)
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
                    paste("Treatment(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Treatment)
                  })
                ),
                br(),
                br(),
                div(
                  style = "display: inline-block;",
                  shiny::renderText({
                    paste("Gene Perturbation(s): ", FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Genotype)
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
                tags$a(href = FilterInput[[as.integer(substr(input$current_id,11,15))]][[3]]$Link, "Study Link - Go to Gene Expression Omnibus (GEO) Record"),
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
                  inputId = "dismiss_modal",
                  label = "Close",
                  class = "btn-danger"
                )
              )
            ))
            
          })
          
        }

        if (length(Comparison) != 0) {
          
          #Generate Comparison Table to view
          {
            
            DFC <- data.frame(stringsAsFactors=FALSE)
            
            if (ncol(FilterInput[[as.integer(substr(input$current_id,11,15))]][[2]]) > 2) {
              DFC <- FilterInput[[as.integer(substr(input$current_id,11,15))]][[2]][,c(2:ncol(FilterInput[[as.integer(substr(input$current_id,11,15))]][[2]]))]
            } else {
              DFC <- as.data.frame(FilterInput[[as.integer(substr(input$current_id,11,15))]][[2]][,2])
              colnames(DFC) <- names(FilterInput[[as.integer(substr(input$current_id,11,15))]][[2]])[2]
            }
            
            for (i in 1:ncol(DFC)) {
              for (j in 1:length(DFC[,i])) {
                
                if (substr(DFC[,i][j], 1, 1) == " ") {
                  DFC[,i][j] <- substr(DFC[,i][j], 2, nchar(DFC[,i][j])) 
                } 
                if (substr(DFC[,i][j], nchar(DFC[,i][j]), nchar(DFC[,i][j])) == " ") {
                  DFC[,i][j] <- substr(DFC[,i][j], 1, (nchar(DFC[,i][j]) - 1)) 
                } 
                
              }
            }
            
            DFC <- distinct(DFC)
            
            valuesC <- reactiveValues(
              df = DFC,
              dt_row = NULL,
              keep_track_id = nrow(df) + 1
              
            )
            
            output$ComparisonTable <- DT::renderDataTable(
              {
                shiny::isolate(valuesC$df)
              },
              escape = F,
              rownames = FALSE,
              options = list(dom = 't')
              #or options = list(Searching = FALSE)
            )
            
          }
          
          #User Selects a Comparison
          {
            observeEvent(input$Select_Comparison, {
              
              if (input$Select_Comparison != "[Select A Comparison]") {
                
                #Create data table, volcano plot, and MA plot
                {  
                  DF2 = data.frame(log2FC=double(),
                                   log2CPM = double(),
                                   FValue=double(),
                                   PValue=double(),
                                   stringsAsFactors=FALSE)
                  
                  DF2 <- FilterInput[[as.integer(substr(input$current_id,11,15))]][[which(names(FilterInput[[as.integer(substr(input$current_id,11,15))]]) == paste0(input$Select_Comparison, " Full Results"))]]
                  colnames(DF2) <- c("log2FC", "log2CPM", "FValue", "PValue")
                  
                  values3 <<- reactiveValues(
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
                      
                      plot_ly(values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                              text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                            -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                        layout(xaxis = list(title="log2FC"),
                               yaxis = list(title="-log10(PValue)"))
                      

                      
                    })
                  }
                  
                  #Create + Output MA Plot in App
                  {
                    
                    output$MAPlot <- renderPlotly({
                      
                      plot_ly(values3$df, x = values3$df$log2CPM, y = values3$df$log2FC, 
                              text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:', 
                                            values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                        layout(xaxis = list(title="log2CPM"),
                               yaxis = list(title="log2FC"))
                  
                      
                    })
                  }
                  
                  #Download Plots
                  {
                    VP <- plot_ly(values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                         text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                       -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                      layout(xaxis = list(title="log2FC"),
                             yaxis = list(title="-log10(PValue)"))
                    
                    MA <- plot_ly(values3$df, x = values3$df$log2CPM, y = values3$df$log2FC, 
                                    text = ~paste("log2CPM: ", values3$df$logCPM, '$<br>log2FC:', 
                                                  values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                      layout(xaxis = list(title="log2CPM"),
                             yaxis = list(title="log2FC"))
                    
                    output$fullAnalysisOutput <- downloadHandler(
                      filename = function() {
                        paste0(input$Select_Comparison, ".zip")
                      },
                      content <- function(file) {
                        UserDir = tempdir()
                        write.csv(DF2, paste0(UserDir, "/Analysis_Table.csv"))
                        MyWidget = combineWidgets(VP, MA, ncol = 2)
                        saveWidget(MyWidget, paste0(UserDir, "/Analysis_Plots.html"))
                        zip(file, paste0(UserDir, "/Analysis_Table.csv"), paste0(UserDir, "/Analysis_Plots.html"), flags = "-r9Xj")
                      },
                      contentType = "application/zip")
                  }
                  
                  
                }
                
                #Update the gene list and pathways
                {
                  GeneList <<- rownames(values3$df)
                  
                  observe({
                    updateSelectizeInput(session, inputId = "Find_Gene", "Find a Gene", choices = c(Choose = "", GeneList),
                                         options = list(maxOptions = 15000))
                    updateSelectizeInput(session, inputId = "Find_Pathway", "Find a Pathway", choices = c(Choose = "", PathwaySelect),
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
                          
                          VP1 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                         text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                       -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(showlegend = FALSE, xaxis = list(title="log2FC"),
                                   yaxis = list(title="-log10(PValue)"))
                          
                          VP1 <- add_trace(VP1, x = Found_Gene$log2FC, y = -log10(Found_Gene$PValue), type = "scatter", mode = "markers", size = 50, color = I("red"), inherit = FALSE, text = ~paste("log2FC: ", Found_Gene$log2FC, '$<br>-log10(PValue):',
                                                                                                                                                                                                      -log10(Found_Gene$PValue), '$<br>Gene Name:', row.names(Found_Gene)))
                          
                          VP1
                          
                          
                        })
                      }
                      
                      #MAPlot
                      {
                        
                        output$MAPlot <- renderPlotly({
                          
                          MA1 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2CPM, y = values3$df$log2FC,
                                         text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:',
                                                       values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="log2CPM"),
                                   yaxis = list(title="log2FC"))
                          
                          MA1 <- add_trace(MA1, x = Found_Gene$log2CPM, y = Found_Gene$log2FC, type = "scatter", mode = "markers", size = 50, color = I("red"), inherit = FALSE, text = ~paste("log2FC: ", Found_Gene$log2CPM, '$<br>-log10(PValue):',
                                                                                                                                                                                               Found_Gene$log2FC, '$<br>Gene Name:', row.names(Found_Gene)))
                          MA1 
                          
                         
                          
                        })
                      }
                      
                      #Download Plots
                      {
                        VP1 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                       text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                     -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                          layout(showlegend = FALSE, xaxis = list(title="log2FC"),
                                 yaxis = list(title="-log10(PValue)"))
                        
                        VP1 <- add_trace(VP1, x = Found_Gene$log2FC, y = -log10(Found_Gene$PValue), type = "scatter", mode = "markers", size = 50, color = I("red"), inherit = FALSE, text = ~paste("log2FC: ", Found_Gene$log2FC, '$<br>-log10(PValue):',
                                                                                                                                                                                                    -log10(Found_Gene$PValue), '$<br>Gene Name:', row.names(Found_Gene)))
                          
                        MA1 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2CPM, y = values3$df$log2FC,
                                       text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:',
                                                     values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                          layout(xaxis = list(title="log2CPM"),
                                 yaxis = list(title="log2FC"))
                        
                        MA1 <- add_trace(MA1, x = Found_Gene$log2CPM, y = Found_Gene$log2FC, type = "scatter", mode = "markers", size = 50, color = I("red"), inherit = FALSE, text = ~paste("log2FC: ", Found_Gene$log2CPM, '$<br>-log10(PValue):',
                                                                                                                                                                                             Found_Gene$log2FC, '$<br>Gene Name:', row.names(Found_Gene)))
                        
                        output$fullAnalysisOutput <- downloadHandler(
                          filename = function() {
                            paste0(input$Select_Comparison, ".zip")
                          },
                          content <- function(file) {
                            UserDir = tempdir()
                            write.csv(DF2, paste0(UserDir, "/Analysis_Table.csv"))
                            MyWidget = combineWidgets(VP1, MA1, ncol = 2)
                            saveWidget(MyWidget, paste0(UserDir, "/Analysis_Plots.html"))
                            zip(file, paste0(UserDir, "/Analysis_Table.csv"), paste0(UserDir, "/Analysis_Plots.html"), flags = "-r9Xj")
                          },
                          contentType = "application/zip")
                      }
                      
                      #Adjust Table
                      {
                        if (input$Find_Gene != "[Select a Gene]") {
                          
                          DFAdjG <- values3$df[which(rownames(values3$df) == input$Find_Gene),]
                          
                          valuesAdjG <- reactiveValues(
                            df = DFAdjG,
                            dt_row = NULL,
                            keep_track_id = nrow(df) + 1
                          )
                          
                          
                          output$ResultTable <- DT::renderDataTable(
                            {
                              shiny::isolate(valuesAdjG$df)
                            },
                            escape = F,
                            rownames = TRUE,
                            options = list(processing = FALSE)
                          )
                        } else {
                          
                          output$ResultTable <- DT::renderDataTable(
                            {
                              shiny::isolate(values3$df)
                            },
                            escape = F,
                            rownames = TRUE,
                            options = list(processing = FALSE)
                          )
                        }
                      }
                      
                    }
                    
                  })
                }
                
                #Highlight a pathway on the Volcano or MA Plot
                {
                  observeEvent(input$Find_Pathway, {
                    
                    if (input$Find_Pathway != "" & input$Find_Pathway != "[Select a Pathway]") {
                      
                      length(grep(TestGenes[j], GenePathway_Data[[i]]$Gene)) != 0
                      
                      Selected_Pathway <- PathwayFile[which(PathwayFile$Pathway == input$Find_Pathway),]
                      Selected_Pathway <- Selected_Pathway %>% separate(Gene, c("Gene1", "Gene2"), sep = " ")

                      Selected_Genes <- subset(values3$df, rownames(values3$df) %in% Selected_Pathway$Gene1)
                      Selected_Genes <- rbind(Selected_Genes, subset(values3$df, rownames(values3$df) %in% Selected_Pathway$Gene2))

                      #Volcano Plot
                      {
                        
                        output$VolcanoPlot <- renderPlotly({
                          
                          VP2 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                         text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                       -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="log2FC"),
                                   yaxis = list(title="-log10(PValue)"))
                          
                          VP2 <- add_trace(VP2, x = Selected_Genes$log2FC, y = -log10(Selected_Genes$PValue), type = "scatter", mode = "markers", size = 50, color = I("green"), inherit = FALSE, text = ~paste("log2FC: ", Selected_Genes$log2FC, '$<br>-log10(PValue):',
                                                                                                                                                                                                                --log10(Selected_Genes$PValue), '$<br>Gene Name:', row.names(Selected_Genes)))
                          VP2
                          
                          
                        })
                      }
                      
                      #MA Plot
                      {
                        
                        output$MAPlot <- renderPlotly({
                          
                          MA2 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2CPM, y = values3$df$log2FC,
                                         text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>logFC:',
                                                       values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="log2CPM"),
                                   yaxis = list(title="log2FC"))
                          
                          MA2 <- add_trace(MA2, x = Selected_Genes$log2CPM, y = Selected_Genes$log2FC, type = "scatter", mode = "markers", size = 50, color = I("green"), inherit = FALSE, text = ~paste("log2CPM: ", Selected_Genes$log2CPM, '$<br>logFC:',
                                                                                                                                                                                                         Selected_Genes$log2FC, '$<br>Gene Name:', row.names(Selected_Genes)))
                          MA2
                          
                          
                        })
                      }
                      
                      #Download Plots
                      {
                        VP2 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                       text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                     -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                          layout(xaxis = list(title="log2FC"),
                                 yaxis = list(title="-log10(PValue)"))
                        
                        VP2 <- add_trace(VP2, x = Selected_Genes$log2FC, y = -log10(Selected_Genes$PValue), type = "scatter", mode = "markers", size = 50, color = I("green"), inherit = FALSE, text = ~paste("log2FC: ", Selected_Genes$log2FC, '$<br>-log10(PValue):',
                                                                                                                                                                                                              --log10(Selected_Genes$PValue), '$<br>Gene Name:', row.names(Selected_Genes)))
                        
                        MA2 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2CPM, y = values3$df$log2FC,
                                       text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:',
                                                     values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                          layout(xaxis = list(title="log2CPM"),
                                 yaxis = list(title="log2FC"))
                        
                        MA2 <- add_trace(MA2, x = Selected_Genes$log2CPM, y = Selected_Genes$log2FC, type = "scatter", mode = "markers", size = 50, color = I("green"), inherit = FALSE, text = ~paste("log2CPM: ", Selected_Genes$log2CPM, '$<br>logFC:',
                                                                                                                                                                                                       Selected_Genes$log2FC, '$<br>Gene Name:', row.names(Selected_Genes)))
                        
                        output$fullAnalysisOutput <- downloadHandler(
                          filename = function() {
                            paste0(input$Select_Comparison, ".zip")
                          },
                          content <- function(file) {
                            UserDir = tempdir()
                            write.csv(Selected_Genes, paste0(UserDir, "/Analysis_Table.csv"))
                            MyWidget = combineWidgets(VP2,MA2, ncol = 2)
                            saveWidget(MyWidget, paste0(UserDir, "/Analysis_Plots.html"))
                            zip(file, paste0(UserDir, "/Analysis_Table.csv"), paste0(UserDir, "/Analysis_Plots.html"), flags = "-r9Xj")
                          },
                          contentType = "application/zip")
                      }
                      
                      #Adjust Table
                      {
                        if (input$Find_Gene != "[Select a Gene]") {
                          
                          DFAdjPath <- Selected_Genes
                          
                          valuesAdjG <- reactiveValues(
                            df = DFAdjPath,
                            dt_row = NULL,
                            keep_track_id = nrow(df) + 1
                          )
                          
                          
                          output$ResultTable <- DT::renderDataTable(
                            {
                              shiny::isolate(valuesAdjG$df)
                            },
                            escape = F,
                            rownames = TRUE,
                            options = list(processing = FALSE)
                          )
                        } else {
                          
                          output$ResultTable <- DT::renderDataTable(
                            {
                              shiny::isolate(values3$df)
                            },
                            escape = F,
                            rownames = TRUE,
                            options = list(processing = FALSE)
                          )
                        }
                      }
                      
                    }
                  })
                }
                
                #Highlight genes under P-Value cutoff + Adjust Table
                {
                  observeEvent(input$Choose_P, {
                    
                    if (input$Choose_P != "[Select a Cutoff]") {
                      
                      PCutoff <- values3$df[which(values3$df$PValue < as.double(input$Choose_P)),]

                      #Volcano Plot
                      {
                        
                        output$VolcanoPlot <- renderPlotly({
                          
                          VP3 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                         text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                       -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="log2FC"),
                                   yaxis = list(title="-log10(PValue)"))
                          
                          VP3 <- add_trace(VP3, x = PCutoff$log2FC, y = -log10(PCutoff$PValue), type = "scatter", mode = "markers", size = 50, color = I("orange"), inherit = FALSE, text = ~paste("log2FC: ", PCutoff$log2FC, '$<br>-log10(PValue):',
                                                                                                                                                                                                   -log10(PCutoff$PValue), '$<br>Gene Name:', row.names(PCutoff)))
                          VP3

                          
                        })
                      }
                      
                      #MA Plot
                      {
                        
                        output$MAPlot <- renderPlotly({
                          
                          MA3 <- plot_ly(type = "scatter", values3$df, x = values3$df$logCPM, y = values3$df$log2FC,
                                         text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:',
                                                       values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="log2CPM"),
                                   yaxis = list(title="log2FC"))
                          
                          MA3 <- add_trace(MA3, x = PCutoff$log2CPM, y = PCutoff$log2FC, type = "scatter", mode = "markers", size = 50, color = I("orange"), inherit = FALSE, text = ~paste("log2CPM: ", PCutoff$log2CPM, '$<br>log2FC:',
                                                                                                                                                                                            PCutoff$log2FC, '$<br>Gene Name:', row.names(PCutoff)))
                          MA3
                          
                          
                        })
                      }
                
                      
                      #Adjust Table
                      {
                        if (input$Choose_P != "[Select a Cutoff]") {
                          
                          DFAdj <- values3$df[which(values3$df$PValue < as.double(input$Choose_P)),]
                          
                          valuesAdj <- reactiveValues(
                            df = DFAdj,
                            dt_row = NULL,
                            keep_track_id = nrow(df) + 1
                          )
                          
                          
                          output$ResultTable <- DT::renderDataTable(
                            {
                              shiny::isolate(valuesAdj$df)
                            },
                            escape = F,
                            rownames = TRUE,
                            options = list(processing = FALSE)
                          )
                        } else {
                          
                          output$ResultTable <- DT::renderDataTable(
                            {
                              shiny::isolate(values3$df)
                            },
                            escape = F,
                            rownames = TRUE,
                            options = list(processing = FALSE)
                          )
                        }
                      }
                      
                      #Download Plots
                      {
                        VP3 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                       text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                     -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                          layout(xaxis = list(title="log2FC"),
                                 yaxis = list(title="-log10(PValue)")) 
                        VP3 <- add_trace(VP3, x = PCutoff$log2FC, y = -log10(PCutoff$PValue), type = "scatter", mode = "markers", size = 50, color = I("orange"), inherit = FALSE, text = ~paste("log2FC: ", PCutoff$log2FC, '$<br>-log10(PValue):',
                                                                                                                                                                                                 -log10(PCutoff$PValue), '$<br>Gene Name:', row.names(PCutoff)))
                        
                        MA3 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2CPM, y = values3$df$log2FC,
                                       text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:',
                                                     values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                          layout(xaxis = list(title="log2CPM"),
                                 yaxis = list(title="log2FC"))
                        MA3 <- add_trace(MA3, x = PCutoff$log2CPM, y = PCutoff$log2FC, type = "scatter", mode = "markers", size = 50, color = I("orange"), inherit = FALSE, text = ~paste("log2CPM: ", PCutoff$log2CPM, '$<br>log2FC:',
                                                                                                                                                                                          PCutoff$log2FC, '$<br>Gene Name:', row.names(PCutoff)))
                        
                        output$fullAnalysisOutput <- downloadHandler(
                          filename = function() {
                            paste0(input$Select_Comparison, ".zip")
                          },
                          content <- function(file) {
                            UserDir = tempdir()
                            write.csv(DFAdj, paste0(UserDir, "/Analysis_Table.csv"))
                            MyWidget = combineWidgets(VP3, MA3, ncol = 2)
                            saveWidget(MyWidget, paste0(UserDir, "/Analysis_Plots.html"))
                            zip(file, paste0(UserDir, "/Analysis_Table.csv"), paste0(UserDir, "/Analysis_Plots.html"), flags = "-r9Xj")
                          },
                          contentType = "application/zip")
                      }
                      
                      }
                      
                      
                  })
                }
                
                #Highlight genes over Fold-Change cutoff + Adjust Table
                {
                  observeEvent(input$Choose_FC, {
                    
                    if (input$Choose_FC != "[Select a Cutoff]") {
                      
                      FCCutoff <- values3$df[which(abs(values3$df$log2FC) >= as.double(input$Choose_FC)),]
          
                      
                      #Create Volcano Plot
                      {
                        
                        output$VolcanoPlot <- renderPlotly({
                          
                          VP4 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                         text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                       -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="log2FC"),
                                   yaxis = list(title="-log10(PValue)"))
                          
                          VP4 <- add_trace(VP4, x = FCCutoff$log2FC, y = -log10(FCCutoff$PValue), type = "scatter", mode = "markers", size = 50, color = I("orange"), inherit = FALSE, text = ~paste("log2FC: ", FCCutoff$log2FC, '$<br>-log10(PValue):',
                                                                                                                                                                                                     -log10(FCCutoff$PValue), '$<br>Gene Name:', row.names(FCCutoff)))
                          VP4
                          
                          
                        })
                      }
                      
                      #Create MA Plot
                      {
                        
                        output$MAPlot <- renderPlotly({
                        
                          MA4 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2CPM, y = values3$df$log2FC,
                                         text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:',
                                                       values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                            layout(xaxis = list(title="log2CPM"),
                                   yaxis = list(title="log2FC"))
                          
                          MA4 <- add_trace(MA4, x = FCCutoff$log2CPM, y = FCCutoff$log2FC, type = "scatter", mode = "markers", size = 50, color = I("orange"), inherit = FALSE, text = ~paste("log2CPM: ", FCCutoff$log2CPM, '$<br>log2FC:',
                                                                                                                                                                                              FCCutoff$log2FC, '$<br>Gene Name:', row.names(FCCutoff)))
                          MA4
                          
                          
                        })
                      }
                      
                      #Adjust Table
                      {
                        if (input$Choose_FC != "[Select a Cutoff]") {
                          
                          DFAdj2 <- values3$df[which(abs(values3$df$log2FC) >= as.double(input$Choose_FC)),]
                          
                          valuesAdj2 <- reactiveValues(
                            df = DFAdj2,
                            dt_row = NULL,
                            keep_track_id = nrow(df) + 1
                          )
                          
                          
                          output$ResultTable <- DT::renderDataTable(
                            {
                              shiny::isolate(valuesAdj2$df)
                            },
                            escape = F,
                            rownames = TRUE,
                            options = list(processing = FALSE)
                          )
                        } else {
                          
                          output$ResultTable <- DT::renderDataTable(
                            {
                              shiny::isolate(values3$df)
                            },
                            escape = F,
                            rownames = TRUE,
                            options = list(processing = FALSE)
                          )
                        }
                      }
                      
                      #Download Plots
                      {
                        VP4 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                                       text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                                     -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                          layout(xaxis = list(title="log2FC"),
                                 yaxis = list(title="-log10(PValue)"))
                        
                        VP4 <- add_trace(VP4, x = FCCutoff$log2FC, y = -log10(FCCutoff$PValue), type = "scatter", mode = "markers", size = 50, color = I("orange"), inherit = FALSE, text = ~paste("log2FC: ", FCCutoff$log2FC, '$<br>-log10(PValue):',
                                                                                                                                                                                                   -log10(FCCutoff$PValue), '$<br>Gene Name:', row.names(FCCutoff)))
                        
                        MA4 <- plot_ly(type = "scatter", values3$df, x = values3$df$log2CPM, y = values3$df$log2FC,
                                       text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:',
                                                     values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                          layout(xaxis = list(title="log2CPM"),
                                 yaxis = list(title="log2FC"))
                        
                        MA4 <- add_trace(MA4, x = FCCutoff$log2CPM, y = FCCutoff$log2FC, type = "scatter", mode = "markers", size = 50, color = I("orange"), inherit = FALSE, text = ~paste("log2CPM: ", FCCutoff$log2CPM, '$<br>log2FC:',
                                                                                                                                                                                            FCCutoff$log2FC, '$<br>Gene Name:', row.names(FCCutoff)))
                        
                        output$fullAnalysisOutput <- downloadHandler(
                          filename = function() {
                            paste0(input$Select_Comparison, ".zip")
                          },
                          content <- function(file) {
                            UserDir = tempdir()
                            write.csv(DFAdj2, paste0(UserDir, "/Analysis_Table.csv"))
                            MyWidget = combineWidgets(VP4, MA4, ncol = 2)
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
          
          #Reset Analysis Outputs
          {
            observeEvent(input$RESET_ANALYSIS_OUTPUTS, {
              
              observe({
                updateSelectizeInput(session, inputId = "Choose_P", choices = c("[Select a Cutoff]", 0.001, 0.01, 0.05, 0.10, 0.15, 0.20),
                                     selected = "[Select a Cutoff]")
              })
              observe({
                updateSelectizeInput(session, inputId = "Choose_FC", choices = c("[Select a Cutoff]", 2, 3, 4, 5, 10, 50),
                                     selected = "[Select a Cutoff]")
              })
              observe({
                updateSelectizeInput(session, inputId = "Find_Gene", choices = c(Choose = "", GeneList),
                                     options = list(maxOptions = 15000))
              })
              observe({
                updateSelectizeInput(session, inputId = "Find_Pathway", choices = c(Choose = "", PathwaySelect),
                                     options = list(maxOptions = 15000))
              })
              
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
                  
                  plot_ly(values3$df, x = values3$df$log2FC, y = -log10(values3$df$PValue),
                          text = ~paste("log2FC: ", values3$df$log2FC, '$<br>-log10(PValue):',
                                        -log10(values3$df$PValue), '$<br>Gene Name:', row.names(values3$df))) %>%
                    layout(xaxis = list(title="log2FC"),
                           yaxis = list(title="-log10(PValue)"))
                  
                  
                  
                })
              }
              
              #Create + Output MA Plot in App
              {
                
                output$MAPlot <- renderPlotly({
                  
                  plot_ly(values3$df, x = values3$df$log2CPM, y = values3$df$log2FC, 
                          text = ~paste("log2CPM: ", values3$df$log2CPM, '$<br>log2FC:', 
                                        values3$df$log2FC, '$<br>Gene Name:', row.names(values3$df))) %>%
                    layout(xaxis = list(title="log2CPM"),
                           yaxis = list(title="log2FC"))
                  
                  
                })
              }
              
            })
            
          }
          
        } else {
          
          shiny::observeEvent(input$ShowAnalysis, {
            
            if (input$Select_Comparison != "[Select A Comparison]") {
              
            output$ResultTable <- DT::renderDataTable(
              {
                shiny::isolate(VALUES$CFSeq_Data[[SpeciesID]][[as.integer(substr(input$current_id,11,15))]]$Count_Table)
              },
              escape = F,
              rownames = TRUE,
              options = list(processing = FALSE)
            )
            
            }
            
          })
          
        }
        
      })
  }

  #5. Reset Analysis Window
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
          observe({
            updateSelectizeInput(session, inputId = "Choose_P", choices = c("[Select a Cutoff]", 0.001, 0.01, 0.05, 0.10, 0.15, 0.20),
                                 selected = "[Select a Cutoff]")
          })
          observe({
            updateSelectizeInput(session, inputId = "Choose_FC", choices = c("[Select a Cutoff]", 2, 3, 4, 5, 10, 50),
                                 selected = "[Select a Cutoff]")
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
          output$Comparison <-  NULL
          DF2 <- NULL
          Values3 <- NULL
          GeneList <- NULL
        }
        
      })
      
    }
    
  #6. Dismiss Modal Dialog
  {
      shiny::observeEvent(input$dismiss_modal, {
        shiny::removeModal()
      })
    }
    
  #7. Download Data
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
    
      output$downloadCountDT2 <- downloadHandler(
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
  
  #8. Reset Species
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


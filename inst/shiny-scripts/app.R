library(shiny)
library(shinyalert)
library(plotly)

# Define UI
ui <- fluidPage(
  
  # Title
  titlePanel("Antibody Loop Comparison"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Description and instructions
      tags$p("This is a Shiny App that is part of the ReadAb package in R. 
              Its purpose is to illustrate the similarity of binding loops 
              between antibodies in the ReadAb dataset."),
      
      tags$b("Description:"),
      tags$p("ReadAb is an R package for the analysis of antibody binding loops. 
              This Shiny App allows users to select antibodies from the ReadAb 
              dataset, specify binding loops, and compare the sequences."),
      
      # Break elements for visual clarity
      br(),
      br(),
      
      
      tags$p("Instructions: Upload antibody PDB files, then choose 
              renumbering schemes and specify weights for each loop."),
      
      br(),
      
      # Allow the user to upload up to 4 Antibody PDBs
      
      # Upload the each antibody PDB, specifying the PDB ID, renumbering scheme, and
      # the chain ids for the heavy and light chains
      
      #Antibody 1
      
      # Antibody PDB path
      fileInput(
        inputId = "file1",
        label = "Upload antibody 1 PDB file",
        accept = c(".pdb")
      ),
      # Renumbering scheme
      selectInput(
        inputId = "scheme1",
        label = "Select Renumbering Scheme for Antibody 1",
        choices = c("Chothia", "AHo/Honneger", "IMGT", "Kabat")
      ),
      # PDB ID
      textInput(inputId = "pdb1", label = "PDB ID for Antibody 1"),
      # Chain IDs
      tags$p("Specify ID of the heavy and light chains in the antibody 1 PDB"),
      textInput(inputId = "heavy1", label = "Heavy chain IDs separated by a comma (ex. A,B,C)"),
      textInput(inputId = "light1", label = "Light chain ID separated by a comma (ex. A,B,C)"),
      br(),
      
      # Antibody 2
      
      # Antibody PDB path
      fileInput(
        inputId = "file2",
        label = "Upload the PDB file for antibody 2",
        accept = c(".pdb")
      ),
      # Renumbering scheme
      selectInput(
        inputId = "scheme2",
        label = "Select Renumbering Scheme for Antibody 2",
        choices = c("Chothia", "AHo/Honneger", "IMGT", "Kabat")
      ),
      # PDB ID
      textInput(inputId = "pdb2", label = "PDB ID for Antibody 2"),
      # Chain IDs
      tags$p("Specify ID of the heavy and light chains in the antibody 2 PDB"),
      textInput(inputId = "heavy2", label = "Heavy chain IDs separated by a comma ex. (A,B,C)"),
      textInput(inputId = "light2", label = "Light chain ID separated by a comma ex. (A,B,C)"),
      br(),
     
      # Antibody 3
      
      # Antibody PDB path
      fileInput(
        inputId = "file3",
        label = "Upload the PDB file for antibody 3 (optional)",
        accept = c(".pdb")
      ),
      # Renumbering scheme
      selectInput(
        inputId = "scheme3",
        label = "Select Renumbering Scheme for Antibody 3",
        choices = c("Chothia", "AHo/Honneger", "IMGT", "Kabat")
      ),
      # PDB ID
      textInput(inputId = "pdb3", label = "PDB ID for Antibody 3"),
      # Chain IDs
      tags$p("Specify ID of the heavy and light chains in the antibody 3 PDB"),
      textInput(inputId = "heavy3", label = "Heavy chain IDs separated by a comma ex. (A,B,C)"),
      textInput(inputId = "light3", label = "Light chain ID separated by a comma ex. (A,B,C)"),
      br(),
      
      # Antibody 4
      
      # Antibody PDB path
      fileInput(
        inputId = "file4",
        label = "Upload the PDB file for antibody 4 (optional)",
        accept = c(".pdb")
      ),
      # Renumbering scheme
      selectInput(
        inputId = "scheme4",
        label = "Select Renumbering Scheme for Antibody 4",
        choices = c("Chothia", "AHo/Honneger", "IMGT", "Kabat")
      ),
      # PDB ID
      textInput(inputId = "pdb4", label = "PDB ID for Antibody 4"),
      # Chain IDs
      tags$p("Specify ID of the heavy and light chains in the antibody 4 PDB"),
      textInput(inputId = "heavy4", label = "Heavy chain IDs separated by a comma ex. (A,B,C)"),
      textInput(inputId = "light4", label = "Light chain ID separated by a comma ex. (A,B,C)"),
      br(),
      
      # Loop weights
      tags$p("Specify weights for each loop (values between 0 and 1, summing to 1)"),
      textInput(inputId = "wH1", label = "Weight for the H1 loop", value = "0.167"),
      textInput(inputId = "wH2", label = "Weight for the H2 loop", value = "0.167"),
      textInput(inputId = "wH3", label = "Weight for the H3 loop", value = "0.167"),
      textInput(inputId = "wL1", label = "Weight for the L1 loop", value = "0.167"),
      textInput(inputId = "wL2", label = "Weight for the L2 loop", value = "0.167"),
      textInput(inputId = "wL3", label = "Weight for the L3 loop", value = "0.167"),
      
      # Run button
      actionButton(inputId = "runButton", label = "Run")
    ),
    
    # Main panel for outputs
    mainPanel(
      tabsetPanel(
        tabPanel("H1 Comparison", plotlyOutput("H1Plot"), tableOutput("H1Table")),
        tabPanel("H2 Comparison", plotlyOutput("H2Plot"), tableOutput("H2Table")),
        tabPanel("H3 Comparison", plotlyOutput("H3Plot"), tableOutput("H3Table")),
        tabPanel("L1 Comparison", plotlyOutput("L1Plot"), tableOutput("L1Table")),
        tabPanel("L2 Comparison", plotlyOutput("L2Plot"), tableOutput("L2Table")),
        tabPanel("L3 Comparison", plotlyOutput("L3Plot"), tableOutput("L3Table")),
        tabPanel("Overall Loop Comparison", plotlyOutput("AllPlot"), tableOutput("AllTable"))
      )
    )
  )
)

server <- function(input, output) {
  
  observeEvent(eventExpr = input$runButton, {
    
    # Initialize an empty list to store antibodies
    antibodies <- list()
    
    # Loop over the antibody inputs (1 to 4)
    for(i in 1:4) {
      file <- paste0("file", i)
      
      # Process only if a file is uploaded
      if (!is.null(input[[file]])) {
        pdb <- paste0("pdb", i)
        scheme <- paste0("scheme", i)
        heavy <- paste0("heavy", i)
        light <- paste0("light", i)
        
        # Parse heavy and light chain IDs
        heavyIDs <- strsplit(input[[heavy]], ",")[[1]]
        lightIDs <- strsplit(input[[light]], ",")[[1]]
        
        # Read the antibody data using the ReadAb package
        antibody <- ReadAb::ReadAntibody(pdbPath = input[[file]]$datapath, 
                                         numbering = input[[scheme]],
                                         heavy = heavyIDs,
                                         light = lightIDs)
        
        # Store antibody in the list by PDB ID
        antibodies[[input[[pdb]]]] = antibody
      }
    }
    
    # Assess similarity for each loop type
    h1Similarity <- ReadAb::AssessLoopSimilarity(antibodies, loop = 'H1')
    h2Similarity <- ReadAb::AssessLoopSimilarity(antibodies, loop = 'H2')
    h3Similarity <- ReadAb::AssessLoopSimilarity(antibodies, loop = 'H3')
    l1Similarity <- ReadAb::AssessLoopSimilarity(antibodies, loop = 'L1')
    l2Similarity <- ReadAb::AssessLoopSimilarity(antibodies, loop = 'L2')
    l3Similarity <- ReadAb::AssessLoopSimilarity(antibodies, loop = 'L3')
    
    # Assess overall similarity with weights
    overallSimilarity <- ReadAb::AssessOverallLoopSimilarity(antibodies,
                                                             wH1 = as.numeric(input$wH1),
                                                             wH2 = as.numeric(input$wH2),
                                                             wH3 = as.numeric(input$wH3),
                                                             wL1 = as.numeric(input$wL1),
                                                             wL2 = as.numeric(input$wL2),
                                                             wL3 = as.numeric(input$wL3))
    
    # Render plots and tables for each loop
    
    # H1 Loop
    output$H1Plot <- renderPlotly({
      ReadAb::DisplaySimilarityPlot(similarityMatrix = h1Similarity, loop = "H1")
    })
    output$H1Table <- renderTable({
      h1Similarity
    })
    
    # H2 Loop
    output$H2Plot <- renderPlotly({
      ReadAb::DisplaySimilarityPlot(h2Similarity, loop = "H2")
    })
    output$H2Table <- renderTable({
      h2Similarity
    })
    
    # H3 Loop
    output$H3Plot <- renderPlotly({
      ReadAb::DisplaySimilarityPlot(h3Similarity, loop = "H3")
    })
    output$H3Table <- renderTable({
      h3Similarity
    })
    
    # L1 Loop
    output$L1Plot <- renderPlotly({
      ReadAb::DisplaySimilarityPlot(l1Similarity, loop = "L1")
    })
    output$L1Table <- renderTable({
      l1Similarity
    })
    
    # L2 Loop
    output$L2Plot <- renderPlotly({
      ReadAb::DisplaySimilarityPlot(l2Similarity, loop = "L2")
    })
    output$L2Table <- renderTable({
      l2Similarity
    })
    
    # L3 Loop
    output$L3Plot <- renderPlotly({
      ReadAb::DisplaySimilarityPlot(l3Similarity, loop = "L3")
    })
    output$L3Table <- renderTable({
      l3Similarity
    })
    
    # All Loops
    output$AllPlot <- renderPlotly({
      ReadAb::DisplaySimilarityPlot(overallSimilarity, loop = 'all')
    })
    output$AllTable <- renderTable({
      overallSimilarity
    })
  })
}


# Create Shiny app ----
shiny::shinyApp(ui, server)

# [END]

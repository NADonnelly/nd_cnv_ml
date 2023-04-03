#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load Up ========

#Load packages
library(shiny)
library(shinydashboard)
library(tibble)
library(ranger)
library(glmnet)
library(dplyr)
library(magrittr)
library(workflows)
library(recipes)

conflicted::conflict_prefer("box", "shinydashboard")
conflicted::conflict_prefer("observe", "shiny")

#Load the model
m_30 <- readRDS("rf_model_30.rds")
m_5  <- readRDS("en_model_5.rds")

# Define UI ========

ui <- dashboardPage(
  
  
  ## Header =====
  
  #Set the header for our dashboard
  dashboardHeader(title = "ND-GC Prediction Model"),
  
  
  ## Sidebar =====
  
  #This is where you define the contents of the sidebar - we use this for
  #the model names
  dashboardSidebar(
    
    sidebarMenu(
      menuItem("30-Item RF Model", tabName = "model30" , icon = icon("code")),
      menuItem("5-Item EN Model" , tabName = "model5"  , icon = icon("code"))
    )
    
  ),
  
  
  ## Body =====
  
  #This is where we define the contents of the body of the dashboard
  dashboardBody(
    
    tabItems(
      
      ### 30 item model UI =====
      
      #First tab: 30 item model
      tabItem(
        
        #Identify the tab
        tabName = "model30",
        
        #Set tab name
        h2("30 Item Random Forest Model"),
        
        
        #Now we need to define our user interface
        fluidRow(
          
          #We need to put the output predicted probability somewhere, so lets do that here
          valueBox(subtitle = "Predicted Probability ND-GC",
                   value = textOutput("text30"),
                   icon = shiny::icon("info"),
                   color = "light-blue",
                   width = 6),
          
          
          box(
            title = "Reset to defaults",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            actionButton("reset_30","Reset to defaults"),
            width = 6)),
        
        fluidRow(
          
          #We need to gather our inputs for each of the 30 variables (yawn):
          
          #CAPA 
          # 1  pcb2i01 Agoraphobia intensity = AGO      
          box(
            title = "AGO",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_ago",
                         "Does agoraphobia (fear of open or public places) lead to a restricted lifestyle for the young person?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)),
          
          # 2  pcc0i01 Situational anxious affect intensity = SIT
          box(
            title = "SIT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_sit",
                         "Does anxiety in particular situations lead to a restricted lifestyle for the young person?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)),
          
          # 3  pbf4i01          Avoidance of being alone intensity = ALO
          box(
            title = "ALO",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_alo",
                         "Does the young person try to avoid being on their own?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)),
          
          # 4  pbf5i01          Anticipatory distress intensity = ANT
          box(
            title = "ANT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_ant",
                         "Is the young person distressed when they think you might be going to leave them? Or when they have to leave you?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)),
          
          # 5  pfb7i02          Sleep problems -  Initial insomnia intensity = INI *** 
          box(
            title = "InI",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_ini",
                         "Does the young person experience initial insomnia (Does it take more than an hour to get to sleep at night?)",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)),  
          
          # 6  prb8i01          Often blurts out answers to questions = BLT ***
          box(
            title = "BLT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_blt",
                         "Does the young person tend to blurt out the answers before the person's finished asking the question?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)), 
          
          # 7  pge2i01          Vandalism intensity = VAN
          box(
            title = "VAN",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_van",
                         "Has the young person ever written on walls? Have they ever damaged or broken or smashed up anything?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)), 
          
          
          # Health & Development
          
          # 8  P_Preg_23        How much did your child weigh at birth? = WGT
          box(
            title = "WGT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            sliderInput("m30_wgt",
                         "How much did your child weigh at birth (in kg)?",
                         min = 0, max = 10, value = 3.3,step = 0.1)),
          
          # 9  P_Health_dev_6   Is your child clumsy? = CLM
          box(
            title = "CLM",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_clm",
                         "Is the young person clumsy?",
                         c("Not at all"    = 1,
                           "Just a little" = 2,
                           "Pretty Much"   = 3,
                           "Very Much"     = 4),
                         selected = 1)),
          
          # 10 P_Health_dev_9   Is your child behind in reading = REA ***
          box(
            title = "REA",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_rea",
                         "Is the young person behind in reading?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)), 
          
          # 11 P_Health_dev_11  Is your child educationally statemented = EST ***
          box(
            title = "EST",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_est",
                         "Does the young person have an educational health care plan?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)), 
          
          # 12 P_Health_dev_12  Was your child talking by the age of 2 = SP2
          box(
            title = "SP2",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_sp2",
                         "Was the young person talking by the age of two?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)), 
          
          # 13 P_Health_dev_15  Has your child had speech therapy = SLT ***
          box(
            title = "SLT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_slt",
                         "Has the young person had speech therapy?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)), 
          
          # 14 P_Health_dev_20  Frequent infections of the chest/airways = RTI
          box(
            title = "RTI",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_rti",
                         "Does the young person get frequent infections of the chest/airways?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)), 
          
          # SDQ
          # 15 P_SDQ_1          Considerate of other's feelings = CNS
          box(
            title = "CNS",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_cns",
                         "Is the young person consider of other people's feelings?",
                         c("Not true"       = 0,
                           "Somewhat true"  = 1,
                           "Certainly true" = 2),
                         selected = 0)),
          
          
          # 16 P_SDQ_6          Rather solitary, tends to play alone = SOL
          box(
            title = "SOL",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_sol",
                         "Is the young person rather solitary and tends to play alone?",
                         c("Not true"       = 0,
                           "Somewhat true"  = 1,
                           "Certainly true" = 2),
                         selected = 0)),
          
          # 17 P_SDQ_9          Helpful if someone is hurt = HRT
          box(
            title = "HRT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_hrt",
                         "Is the young person helpful if someone is hurt?",
                         c("Not true"       = 0,
                           "Somewhat true"  = 1,
                           "Certainly true" = 2),
                         selected = 0)),
          
          # 18 P_SDQ_16         Easily distracted, concentration wanders = DIS
          box(
            title = "DIS",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_dis",
                         "Is the young person easily distracted, does their concentration wander?",
                         c("Not true"       = 0,
                           "Somewhat true"  = 1,
                           "Certainly true" = 2),
                         selected = 0)),
          
          # 19 P_SDQ_19         Often tells lies = LIE
          box(
            title = "LIE",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_lie",
                         "Does the young person often tell lies?",
                         c("Not true"       = 0,
                           "Somewhat true"  = 1,
                           "Certainly true" = 2),
                         selected = 0)),
          
          # 20 P_SDQ_20         Often cheats = CHT
          box(
            title = "CHT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_cht",
                         "Does the young person often cheat?",
                         c("Not true"       = 0,
                           "Somewhat true"  = 1,
                           "Certainly true" = 2),
                         selected = 0)),
          
          # 21 P_SDQ_29         Inconsiderate of others = INC
          box(
            title = "INC",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_inc",
                         "Is the young person often inconsiderate of others?",
                         c("Not true"       = 0,
                           "Somewhat true"  = 1,
                           "Certainly true" = 2),
                         selected = 0)),
          
          #SCQ
          # 22 P_ASQ_14         Special interests unusual in their intensity = SII
          box(
            title = "SII",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_sii",
                         "Does the young person have special interests that are unusual in their intensity?",
                         c("No"  = 0,
                           "Yes" = 1),
                         selected = 0)),
          
          # 23 P_ASQ_39         Imaginative play with another child = PIM
          box(
            title = "PIM",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_pim",
                         "When they were 4 to 5 did they ever play imaginative games with another child in such a way that you could tell they understood what each other was pretending?",
                         c("No"  = 0,
                           "Yes" = 1),
                         selected = 0)),
          
          # 24 P_ASQ_40         Play cooperatively = PCO
          box(
            title = "PCO",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_pco",
                         "When they were 4 to 5 did they play co-operatively in games that need some form of joining in with a group of other children, such as hide and seek or ball games?",
                         c("No"  = 0,
                           "Yes" = 1),
                         selected = 0)),
          
          
          #DCDQ
          
          # 25 CMD_2            Catches a small ball thrown from 6-8ft = CBA ***
          box(
            title = "CBA",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_cba",
                         "Can your child catch a small ball (e.g. tennis ball size) thrown from a distance of 6-8 feet (1.8-2.4 metres) ",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1)),
          
          # 26 CMD_5            Runs as fast and easily as other children = CRF ***
          box(
            title = "CRF",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 400,
            radioButtons("m30_crf",
                         "Does your child runs as fast and in a similar way to other children of the same age and gender?",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1)),
          
          # 27 CMD_6            Can organise her body to do a planned motor activity = COB ***
          box(
            title = "COB",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 400,
            radioButtons("m30_cob",
                         "If your child has a plan to do a motor activity, they can organise their body to follow the plan and effectively complete the task (e.g., building a cardboard or cushion ‘fort’, moving on playground equipment, building a house or a structure with blocks, or using craft materials)?",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1)),
          
          # 28 CMD_11           Likes participating in games requiring good motor skills = CGM
          box(
            title = "CGM",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 400,
            radioButtons("m30_cgm",
                         "Is the young person interested in and likes participating in sports or active games requiring good motor skills?",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1)),
          
          # 29 CMD_14           Child would never be described as a bull in a china shop (clumsy and breaks things) = CBL
          box(
            title = "CBL",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 400,
            radioButtons("m30_cbl",
                         "The young person would never be described as a bull in a china shop (clumsy and breaks things)",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1)),
          
          # 30 CMD_15           Child does not fatigue easily or appear to slouch and "fall out" of the chair       = CSL
          box(
            title = "CSL",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 400,
            radioButtons("m30_csl",
                         "The young person does not fatigue easily or appear to slouch and 'fall out' of the chair",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1)) 

        
      )),
      
      
      ### 5 item model UI ======
      
      #Second Tab: 5 item model
      tabItem(
        
        #Tab name
        tabName = "model5",
        
        # Tab title
        h2("5 Item Elastic Net Model"),
        
        #Now we need to define our interface
        fluidRow(
          
          #We need to put the output predicted probability somewhere, so lets do that here
          valueBox(subtitle = "Predicted Probability ND-GC",
                  value = textOutput("text5"),
                  icon = shiny::icon("info"),
                  color = "light-blue",
                  width = 6),
          
          box(
            title = "Reset to defaults",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            actionButton("reset_5","Reset to defaults"),
            width = 6
          )),
        
        
        fluidRow(
          
          #We need to gather our inputs for each of the 5 variables:
          box(
            title = "AGO",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m5_ago",
                         "Does agoraphobia (fear of open or public places) lead to a restricted lifestyle for the young person?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)),
          
          box(
            title = "ANT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m5_ant",
                         "Is the young person distressed when they think you might be going to leave them? Or when they have to leave you?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)),
          
          
          box(
            title = "BLT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m5_blt",
                         "Does the young person tend to blurt out the answers before the person's finished asking the question?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)), 
          
          
          box(
            title = "SP2",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m5_sp2",
                         "Was the young person talking by the age of two?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 1)), 
          
          
          box(
            title = "CGM",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 400,
            radioButtons("m5_cgm",
                         "Is the young person interested in and likes participating in sports or active games requiring good motor skills?",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1))

        )
        
      )

  )
)
)

# Server backend =====

#This does the backend stuff
server <- function(input, output) {
  
  #This is the output for the 30 item model
  
  output$text30 <- 
    
    
    ## 30 item output =====
    renderText({ 
      
      d_input_30 = tibble(
        #CAPA 
        "pcb2i01" = input$m30_ago |> as.integer(),
        "pcc0i01" = input$m30_sit |> as.integer(),
        "pbf4i01" = input$m30_alo |> as.integer(),
        "pbf5i01" = input$m30_ant |> as.integer(),
        "pfb7i02" = input$m30_ini |> as.integer(),
        "prb8i01" = input$m30_blt |> as.integer(),
        "pge2i01" = input$m30_van |> as.integer(),
        
        # Health & Development
        "P_Preg_23"       = input$m30_wgt |> as.double(),
        "P_Health_dev_6"  = input$m30_clm |> as.integer(),
        "P_Health_dev_9"  = input$m30_rea |> as.integer(), 
        "P_Health_dev_11" = input$m30_est |> as.integer(), 
        "P_Health_dev_12" = input$m30_sp2 |> as.integer(),
        "P_Health_dev_15" = input$m30_slt |> as.integer(), 
        "P_Health_dev_20" = input$m30_rti |> as.integer(),
        
        # SDQ
        "P_SDQ_1"  = input$m30_cns  |> as.integer(),
        "P_SDQ_6"  = input$m30_sol  |> as.integer(),
        "P_SDQ_9"  = input$m30_hrt  |> as.integer(),
        "P_SDQ_16" = input$m30_dis  |> as.integer(),
        "P_SDQ_19" = input$m30_lie  |> as.integer(),
        "P_SDQ_20" = input$m30_cht  |> as.integer(),
        "P_SDQ_29" = input$m30_inc  |> as.integer(),
        
        #SCQ
        "P_ASQ_14"  = input$m30_sii  |> as.integer(),
        "P_ASQ_39"  = input$m30_pim  |> as.integer(),
        "P_ASQ_40"  = input$m30_pco  |> as.integer(),
        
        #DCDQ
        "CMD_2"  = input$m30_cba |> as.integer(),
        "CMD_5"  = input$m30_crf |> as.integer(),
        "CMD_6"  = input$m30_cob |> as.integer(), 
        "CMD_11" = input$m30_cgm |> as.integer(),
        "CMD_14" = input$m30_cbl |> as.integer(),
        "CMD_15" = input$m30_csl |> as.integer())
      
      m_30 |>
        predict(d_input_30,type = "prob") |>
        pull(`.pred_ND-GC`) |>
        round(digits = 3) |>
        paste()
      
    })
  
  #Reset button on 30 item model
  observe({
    
    x <- input$reset_30
    
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_ago",selected = 0)
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_sit",selected = 0)
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_alo",selected = 0)
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_ant",selected = 0)
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_ini",selected = 0)
  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_blt",selected = 0)
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_van",selected = 0)
  updateSliderInput( session = getDefaultReactiveDomain(),"m30_wgt",value    = 3.3)
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_clm",selected = 1)
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_rea",selected = 0)

  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_est",selected = 0)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_sp2",selected = 1)   
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_slt",selected = 0)   
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_rti",selected = 0)   
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cns",selected = 0)   

  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_sol",selected = 0)   
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_hrt",selected = 0)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_dis",selected = 0)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_lie",selected = 0) 
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cht",selected = 0)  
  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_inc",selected = 0)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_sii",selected = 0) 
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_pim",selected = 0)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_pco",selected = 0)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cba",selected = 1)  
  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_crf",selected = 1)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cob",selected = 1)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cgm",selected = 1)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cbl",selected = 1)  
  updateRadioButtons(session = getDefaultReactiveDomain(),"m30_csl",selected = 1)  

  })
  
  
  
  ## 5 item model =====
  
  #This is the output for the 5 item model
  output$text5 <- 
    renderText({ 
      
      d_input_5 = tibble(
        "pcb2i01"         = input$m5_ago |> as.integer(),
        "pbf5i01"         = input$m5_ant |> as.integer(),
        "prb8i01"         = input$m5_blt |> as.integer(),
        "P_Health_dev_12" = input$m5_sp2 |> as.integer(),
        "CMD_11"          = input$m5_cgm |> as.integer())
      

      m_5 |>
        predict(d_input_5,type = "prob") |>
        pull(`.pred_ND-GC`) |>
        round(digits = 4) |>
        paste()
      
    })
  
  
  #Reset button on 5 item model
  observe({
    
    x <- input$reset_5

    updateRadioButtons(session = getDefaultReactiveDomain(),"m5_ago",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m5_ant",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m5_blt",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m5_sp2",selected = 1)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m5_cgm",selected = 1)
    
  })
 
  
}


# Run the application 
shinyApp(ui = ui, server = server)




#CAPA 
# pcb2i01 = AGO
# pcc0i01 = SIT
# pbf4i01 = ALO
# pbf5i01 = ANT
# pfb7i02 = INI 
# prb8i01 = BLT 
# pge2i01 = VAN

# Health & Development
# P_Preg_23 = WGT
# P_Health_dev_6 = CLM
# P_Health_dev_9 = REA
# P_Health_dev_11 = EST 
# P_Health_dev_12 = SP2
# P_Health_dev_15 = SLT 
# P_Health_dev_20 = RTI

# SDQ
#P_SDQ_1 = CNS
#P_SDQ_6 = SOL
#P_SDQ_9 = HRT
#P_SDQ_16 = DIS
#P_SDQ_19 = LIE
#P_SDQ_20 = CHT
#P_SDQ_29 = INC

#SCQ
#P_ASQ_14 = SII
#P_ASQ_39 = PIM
#P_ASQ_40 = PCO

#DCDQ
#CMD_2 = CBA
#CMD_5 = CRF 
#CMD_6 = COB 
#CMD_11 = CGM
#CMD_14 = CBL
#CMD_15 = CSL


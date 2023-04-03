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
library(dplyr)
library(magrittr)
library(workflows)
library(recipes)

conflicted::conflict_prefer("box", "shinydashboard")
conflicted::conflict_prefer("observe", "shiny")

#Load the model
m_30 <- readRDS("rf_model_30.rds")
m_5  <- readRDS("en_model_5.rds")
m_lr <- readRDS("lr_model_30.rds")



# Define UI ========

ui <- dashboardPage(
  
  
  ## Header =====
  
  #Set the header for our dashboard
  dashboardHeader(title = "CNV Prediction Model"),
  
  
  ## Sidebar =====
  
  #This is where you define the contents of the sidebar - we might use this for
  #the 30 item model versus the 4 item model
  dashboardSidebar(
    
    sidebarMenu(
      menuItem("30-Item RF Model", tabName = "model30" , icon = icon("code")),
      menuItem("5-Item EN Model" , tabName = "model5"  , icon = icon("code")),
      menuItem("30-Item LR Model", tabName = "model_lr", icon = icon("code"))
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
        h2("30 Item Prediction Model"),
        
        
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
            width = 6)
          ),
        
        fluidRow(
          
          #We need to gather our inputs for each of the 30 variables (yawn):
          
          # CAPA Items
          
          # 1  pcb1i01         Anxious affect -  Fear of activities in public                                              = FPB                 
          box(
            title = "FPB",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_fpb",
                         "Does fear of activities in public lead to a restricted lifestyle for the child?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),
          
          # 2  pcb3i01         Anxious affect -  Agoraphobia                                                               = AGO     
          box(
            title = "AGO",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_ago",
                         "Does agoraphobia (fear of open or public places) lead to a restricted lifestyle for the child?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),
          
          # 3  pce0i01         Anxious affect -  Fear of blood/injections                                                  = FBI   
          box(
            title = "FBI",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_fbi",
                         "Does the child experiencing fear of blood/injections intrude into at least two activities, and is uncontrollable at least some of the time?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),
          
          # 4  pcd2i01         Rumination, Obsessions and Compulsions - ?Rumination Intensity                              = RUM     
          box(
            title = "RUM",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_rum",
                         "Does rumination (unproductive dwelling on particular themes) occur at least 1 hour daily and intrude in at least two activities and be at least sometimes uncontrollable by the child",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),          
          
          # 5  pda0i02         Depression  -   Episode of Depressed mood intensity                                         = DeI   
          box(
            title = "DeI",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_dei",
                         "Was there a week when the child felt miserable most days?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),    
          
          
          # 6  pda0i03         Depression  -   Period of 2 continuous months without depressed mood in last year           = D2M   
          box(
            title = "D2M",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_d2m",
                         "Has there been a period of two months in the last year when the participant did not feel depressed in mood?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),   
          
          # 7  pda1i01         Depression  -   Distinct quality of depressed mood                                          = DeQ  
          box(
            title = "DeQ",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_deq",
                         "Does depressed mood have a subjectively different quality from sadness, contrasted with an experience that caused sadness, such as loss of a pet or watching a sad film",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 8  pbe1i01         Physical symptoms on separation from caregiver (separation anxiety symptom)                 = SAP  
          box(
            title = "SAP",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_sap",
                         "Does the child get aches and pains, feel sick, get headaches etc. on school days, or at other times when separated from a parent/carer?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 9  pbe7i01         Separation Anxiety -  Intensity of separation worries/anxiety (across multiple activities)  = SAI   
          box(
            title = "SAI",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_sai",
                         "Does the child experience excessive worries or fear concerning separation from the persons to whom the child is attached",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 10 pbf8i01         Separation Anxiety -  Frequency of separation anxiety                                       = SAF  
          box(
            title = "SAF",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_saf",
                         "Does the child sleep with a family member because of persistent refusal to sleep through the night without being near a major attachment figure?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 11 pbf2i01         Separation Anxiety -  Avoidance of sleeping away from family                                = SAS  
          box(
            title = "SAS",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_sas",
                         "Does the child avoid, or attempt to avoid, sleeping away from family, as a result of worrying or anxiety about separation from home or family?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 12 pfb7i01         Sleep probs -  Total Insomnia intensity                                                     = InT 
          box(
            title = "InT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_int",
                         "Does the child experience overall insomnia greater than 1 hour per night?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 13 pfb7i02         Sleep probs -  Initial insomnia intensity                                                   = InI  
          box(
            title = "InI",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_ini",
                         "Does the child experience initial insomnia (Does it take more than an hour to get to sleep at night?)",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 14 prc8i01         hyperactivity -  Forgetful in daily activities intensity                                    = FGT 
          box(
            title = "FGT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_fgt",
                         "Is the child forgetful enough that forgetfulness intrudes into at least two activities and is at least sometimes uncontrollable?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 15 prb8i01         Often blurts out answers to questions (ADHD Impulsivity/Hyperactivity symptom)              = BLT 
          box(
            title = "BLT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_blt",
                         "Does the child often blurt out answers to questions?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 16 pgc3i01         oppositional / conduct disorder  - Lying intensity                                          = LIE 
          box(
            title = "LIE",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_lie",
                         "Does the child lie or distort the truth with intent to deceive others? To be present - lies are told for gain, or to escape punishment, in at least two activities that do not result in others getting in trouble.",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 17 pgc5i01         oppositional / conduct disorder  - Cheating intensity                                       = CHT  
          box(
            title = "CHT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_cht",
                         "Does the child ever cheat in attempts to gain increased marks at school or increased success in other settings by unfair means. To be considered present - takes place in at least two activities, and at least sometimes not responsive to admonition if caught.",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          #Your child's health and development questionairre items
          
          # 18 P_Health_dev_9  Health and Development - Is your child behind in reading                                    = REA 
          box(
            title = "REA",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_rea",
                         "Is your child behind in reading?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 19 P_Health_dev_11 Health and Development - educationally statemented                                          = EST  
          box(
            title = "EST",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_est",
                         "Is your child educationally statemented?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 20 P_Health_dev_13 Health and Development - did your child walk by 18 months                                   = W18 
          box(
            title = "W18",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_w18",
                         "Did your child walk by 18 months of age?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 21 P_Health_dev_15 Health and Development - has your child had speech therapy                                  = SLT    
          box(
            title = "SLT",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_slt",
                         "Has your child had speech therapy?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 22 P_Health_dev_21 Health and Development - other problems with airways/lungs                                  = LUN 
          box(
            title = "LUN",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_lun",
                         "Has your child had problems with their airways or lungs?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 23 P_Health_dev_24 Health and Development - heart problems                                                     = CAR    
          box(
            title = "CAR",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_car",
                         "Has your child had heart problems?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 24 P_Health_dev_27 Health and Development - skeletal or muscular problems                                      = MSK 
          box(
            title = "MSK",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_msk",
                         "Has your child had skeletal or muscular problems?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          #SCQ items
          
          # 25 P_ASQ_7         ASQ Behav and social comm - invented words, odd indirect, metaphorical ways                 = AWM   
          box(
            title = "AWM",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_awm",
                         "Have your child ever used words that they seem to have invented or made up or ever put things in indirect or metaphorical ways?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          # 26 P_ASQ_8         ASQ Behav and social comm - say the same thing over and over                                = ARP
          box(
            title = "ARP",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 300,
            radioButtons("m30_arp",
                         "Has your child ever said the same thing over and over again in exactly the same way, or insist that you say things over and over in the exact same way?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0)
          ),  
          
          
          #DCDQ Items (note these are five choice ordinal, not binary, and that we reverse coded so that higher values mean more impaired)
          
          # 27 CMD_2           Coordination and motor development - catches a small ball thrown from 6-8ft                 = CBA  
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
                         selected = 1)
          ),
          
          # 28 CMD_5           Coordination and motor development - runs as fast and easily as other children              = CRF 
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
                         selected = 1)
          ),
          
          
          # 29 CMD_6           Coordination and motor development - can organise her body to do a planned motor activity   = COB
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
                         selected = 1)
          ),
          
          # 30 CMD_10          Coordination and motor development - cuts pictures and shapes accurately                    = CCP  
          box(
            title = "CCP",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 4,
            height = 400,
            radioButtons("m30_ccp",
                         "Does your child cut pictures and shapes accurately?",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1)
          )
        )
        
      ),
      
      
      ### 5 item model UI ======
      
      #Second Tab: 5 item model
      tabItem(
        
        #Tab name
        tabName = "model5",
        
        # Tab title
        h2("4 Item Prediction Model"),
        
        #Now we need to define our interface
        fluidRow(
          
          #We need to put the output predicted probability somewhere, so lets do that here
          valueBox(subtitle = "Predicted Probability ND-GC",
                  value = textOutput("text4"),
                  icon = shiny::icon("info"),
                  color = "light-blue",
                  width = 6),
          
          box(
            title = "Reset to defaults",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            actionButton("reset_4","Reset to defaults"),
            width = 6
          )),
        
        fluidRow(
          
          #We need to gather our inputs for each of the 4 variables:
          
          # pbe1i01 = SAP = Separation Anxiety: do they get aches and pains, feel sick, get headaches etc. on school days, 
          #                 or at other times when separated from a parent/carer (Binary Variable)
          
          # P_Health_dev_15 = SLT = Health and Development: Has your child had speech therapy? (Binary Variable)
          
          # pfb7i02 = InI = Sleep probs - Insomnia: Initial insomnia present (Binary Variable)
          
          # pda0i02  = DeI = Depression - Was there a week when the participant felt miserable most days (Binary Variable)                                    
          
          box(
            title = "SAP",
            solidHeader = FALSE,
            collapsible = TRUE,
            radioButtons("m4_sap",
                         "Does your child get aches and pains, feel sick, get headaches etc. on school days, or at other times when separated from a parent/carer?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            width = 6,
            height = 200
          ),
          
          box(
            title = "SLT",
            solidHeader = FALSE,
            collapsible = TRUE,
            radioButtons("m4_slt",
                         "Has your child had speech therapy?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            width = 6,
            height = 200
          ),
          
          box(
            title = "InI",
            solidHeader = FALSE,
            collapsible = TRUE,
            radioButtons("m4_ini",
                         "Does your child experience initial insomnia (Does it take more than an hour to get to sleep at night?)",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            width = 6,
            height = 200
          ),
          
          box(
            title = "DeI",
            solidHeader = FALSE,
            collapsible = TRUE,
            radioButtons("m4_dei",
                         "Was there a week when your child felt miserable most days?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            width = 6,
            height = 200)          
          
          
        )
        
      ),
      
      
      
      ### 30 Item LR model UI ======
      
      #Second Tab: lr model
      tabItem(
        
        #Tab name
        tabName = "model_lr",
        
        # Tab title
        h2("Logistic Regression Prediction Model"),
        
        #Now we need to define our interface
        
        fluidRow(
          
          #We need to put the output predicted probability somewhere, so lets do that here
          valueBox(subtitle = "Predicted Probability ND-GC",
                   value = textOutput("textLR"),
                   icon = shiny::icon("info"),
                   color = "light-blue",
                   width = 6),
          
          box(
            title = "Reset to defaults",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            actionButton("reset_lr","Reset to defaults"),
            width = 6)
          ),
        
        fluidRow(
          
          # CAPA Items
          box(
            title = "Please answer all questions",
            solidHeader = FALSE,
            collapsible = TRUE,
            width = 12,

            # 1  pcb1i01         Anxious affect -  Fear of activities in public                                              = FPB   
            radioButtons("mLR_fpb",
                         "Does fear of activities in public lead to a restricted lifestyle for the child?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 2  pcb3i01         Anxious affect -  Agoraphobia                                                               = AGO                 
            radioButtons("mLR_ago",
                         "Does agoraphobia (fear of open or public places) lead to a restricted lifestyle for the child?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 3  pce0i01         Anxious affect -  Fear of blood/injections                                                  = FBI        
            radioButtons("mLR_fbi",
                         "Does the child experiencing fear of blood/injections intrude into at least two activities, and is uncontrollable at least some of the time?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 4  pcd2i01         Rumination, Obsessions and Compulsions - ?Rumination Intensity                              = RUM  
            radioButtons("mLR_rum",
                         "Does rumination (unproductive dwelling on particular themes) occur at least 1 hour daily and intrude in at least two activities and be at least sometimes uncontrollable by the child",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 5  pda0i02         Depression  -   Episode of Depressed mood intensity                                         = DeI   
            radioButtons("mLR_dei",
                         "Was there a week when the child felt miserable most days?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 6  pda0i03         Depression  -   Period of 2 continuous months without depressed mood in last year           = D2M   
            radioButtons("mLR_d2m",
                         "Has there been a period of two months in the last year when the participant did not feel depressed in mood?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 7  pda1i01         Depression  -   Distinct quality of depressed mood                                          = DeQ              
            radioButtons("mLR_deq",
                         "Does depressed mood have a subjectively different quality from sadness, contrasted with an experience that caused sadness, such as loss of a pet or watching a sad film",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),    
            
            # 8  pbe1i01         Physical symptoms on separation from caregiver (separation anxiety symptom)                 = SAP              
            radioButtons("mLR_sap",
                         "Does the child get aches and pains, feel sick, get headaches etc. on school days, or at other times when separated from a parent/carer?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 9  pbe7i01         Separation Anxiety -  Intensity of separation worries/anxiety (across multiple activities)  = SAI   
            radioButtons("mLR_sai",
                         "Does the child experience excessive worries or fear concerning separation from the persons to whom the child is attached",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 10 pbf8i01         Separation Anxiety -  Frequency of separation anxiety                                       = SAF  
            radioButtons("mLR_saf",
                         "Does the child sleep with a family member because of persistent refusal to sleep through the night without being near a major attachment figure?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 11 pbf2i01         Separation Anxiety -  Avoidance of sleeping away from family                                = SAS              
            radioButtons("mLR_sas",
                         "Does the child avoid, or attempt to avoid, sleeping away from family, as a result of worrying or anxiety about separation from home or family?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 12 pfb7i01         Sleep probs -  Total Insomnia intensity                                                     = InT 
            radioButtons("mLR_int",
                         "Does the child experience overall insomnia greater than 1 hour per night?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 13 pfb7i02         Sleep probs -  Initial insomnia intensity                                                   = InI  
            radioButtons("mLR_ini",
                         "Does the child experience initial insomnia (Does it take more than an hour to get to sleep at night?)",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 14 prc8i01         hyperactivity -  Forgetful in daily activities intensity                                    = FGT 
            radioButtons("mLR_fgt",
                         "Is the child forgetful enough that forgetfulness intrudes into at least two activities and is at least sometimes uncontrollable?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 15 prb8i01         Often blurts out answers to questions (ADHD Impulsivity/Hyperactivity symptom)              = BLT 
            radioButtons("mLR_blt",
                         "Does the child often blurt out answers to questions?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 16 pgc3i01         oppositional / conduct disorder  - Lying intensity                                          = LIE 
            radioButtons("mLR_lie",
                         "Does the child lie or distort the truth with intent to deceive others? To be present - lies are told for gain, or to escape punishment, in at least two activities that do not result in others getting in trouble.",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 17 pgc5i01         oppositional / conduct disorder  - Cheating intensity                                       = CHT  
            radioButtons("mLR_cht",
                         "Does the child ever cheat in attempts to gain increased marks at school or increased success in other settings by unfair means. To be considered present - takes place in at least two activities, and at least sometimes not responsive to admonition if caught.",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            #Your child's health and development questionnaire items
            
            # 18 P_Health_dev_9  Health and Development - Is your child behind in reading                                    = REA 
            radioButtons("mLR_rea",
                         "Is your child behind in reading?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 19 P_Health_dev_11 Health and Development - educationally statemented                                          = EST  
            radioButtons("mLR_est",
                         "Is your child educationally statemented?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 20 P_Health_dev_13 Health and Development - did your child walk by 18 months                                   = W18 
            radioButtons("mLR_w18",
                         "Did your child walk by 18 months of age?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 21 P_Health_dev_15 Health and Development - has your child had speech therapy                                  = SLT 
            radioButtons("mLR_slt",
                         "Has your child had speech therapy?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 22 P_Health_dev_21 Health and Development - other problems with airways/lungs                                  = LUN 
            radioButtons("mLR_lun",
                         "Has your child had problems with their airways or lungs?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 23 P_Health_dev_24 Health and Development - heart problems                                                     = CAR   
            radioButtons("mLR_car",
                         "Has your child had heart problems?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 24 P_Health_dev_27 Health and Development - skeletal or muscular problems                                      = MSK 
            radioButtons("mLR_msk",
                         "Has your child had skeletal or muscular problems?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            #SCQ items
            
            # 25 P_ASQ_7         ASQ Behav and social comm - invented words, odd indirect, metaphorical ways                 = AWM  
            radioButtons("mLR_awm",
                         "Have your child ever used words that they seem to have invented or made up or ever put things in indirect or metaphorical ways?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            # 26 P_ASQ_8         ASQ Behav and social comm - say the same thing over and over                                = ARP
            radioButtons("mLR_arp",
                         "Has your child ever said the same thing over and over again in exactly the same way, or insist that you say things over and over in the exact same way?",
                         c("Yes" = 1,
                           "No"  = 0),
                         selected = 0),
            
            #DCDQ Items (note these are five choice ordinal, not binary, and that we reverse coded so that higher values mean more impaired)
            
            # 27 CMD_2           Coordination and motor development - catches a small ball thrown from 6-8ft                 = CBA
            radioButtons("mLR_cba",
                         "Can your child catch a small ball (e.g. tennis ball size) thrown from a distance of 6-8 feet (1.8-2.4 metres) ",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1),
            
            # 28 CMD_5           Coordination and motor development - runs as fast and easily as other children              = CRF 
            radioButtons("mLR_crf",
                         "Does your child runs as fast and in a similar way to other children of the same age and gender?",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1),
            
            # 29 CMD_6           Coordination and motor development - can organise her body to do a planned motor activity   = COB
            radioButtons("mLR_cob",
                         "If your child has a plan to do a motor activity, they can organise their body to follow the plan and effectively complete the task (e.g., building a cardboard or cushion ‘fort’, moving on playground equipment, building a house or a structure with blocks, or using craft materials)?",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1),
            
            # 30 CMD_10          Coordination and motor development - cuts pictures and shapes accurately                    = CCP  
            radioButtons("mLR_ccp",
                         "Does your child cut pictures and shapes accurately?",
                         c("Not at all like your child"  = 5,
                           "A bit like your child"       = 4,
                           "Moderately like your child"  = 3,
                           "Quite a bit like your child" = 2,
                           "Extremely like your child"   = 1),
                         selected = 1)
          )
          
          
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
                          #CAPA variables
                          "pcb1i01" = input$m30_fpb  |> as.integer(),
                          "pcb3i01" = input$m30_ago  |> as.integer(),
                          "pce0i01" = input$m30_fbi  |> as.integer(),
                          "pcd2i01" = input$m30_rum  |> as.integer(),
                          "pda0i02" = input$m30_dei  |> as.integer(),
                          "pda0i03" = input$m30_d2m  |> as.integer(),
                          "pda1i01" = input$m30_deq  |> as.integer(),
                          "pbe1i01" = input$m30_sap  |> as.integer(),
                          "pbe7i01" = input$m30_sai  |> as.integer(),
                          "pbf8i01" = input$m30_saf  |> as.integer(),
                          "pbf2i01" = input$m30_sas  |> as.integer(),
                          "pfb7i01" = input$m30_ini  |> as.integer(),
                          "pfb7i02" = input$m30_int  |> as.integer(),
                          "prc8i01" = input$m30_fgt  |> as.integer(),
                          "prb8i01" = input$m30_blt  |> as.integer(),
                          "pgc3i01" = input$m30_lie  |> as.integer(),
                          "pgc5i01" = input$m30_cht  |> as.integer(),
                          
                          #Health and Development Qs
                          "P_Health_dev_9"  = input$m30_rea  |> as.integer(),
                          "P_Health_dev_11" = input$m30_est  |> as.integer(),
                          "P_Health_dev_13" = input$m30_w18  |> as.integer(),
                          "P_Health_dev_15" = input$m30_slt  |> as.integer(),
                          "P_Health_dev_21" = input$m30_lun  |> as.integer(),
                          "P_Health_dev_24" = input$m30_car  |> as.integer(),
                          "P_Health_dev_27" = input$m30_msk  |> as.integer(),
                          
                          #SCQ
                          "P_ASQ_7" = input$m30_awm  |> as.integer(),
                          "P_ASQ_8" = input$m30_arp  |> as.integer(),
                          
                          #DCDQ
                          "CMD_2"  = input$m30_cba  |> as.integer(),
                          "CMD_5"  = input$m30_crf  |> as.integer(),
                          "CMD_6"  = input$m30_cob  |> as.integer(),
                          "CMD_10" = input$m30_ccp  |> as.integer()
      )
      

      
      m_30 |>
        predict(d_input_30,type = "prob") |>
        pull(`.pred_ND-CNV`) |>
        round(digits = 3) |>
        paste()
      
    })
  
  #Reset button on 30 item model
  observe({
    
    x <- input$reset_30
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_fpb",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_ago",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_fbi",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_rum",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_dei",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_d2m",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_deq",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_sap",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_sai",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_saf",selected = 0)

    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_sas",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_ini",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_int",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_fgt",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_blt",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_lie",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cht",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_rea",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_est",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_w18",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_slt",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_lun",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_car",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_msk",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_awm",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_arp",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cba",selected = 1)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_crf",selected = 1)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_cob",selected = 1)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m30_ccp",selected = 1)
  })
  
  
  
  ## 5 item model =====
  
  #This is the output for the 5 item model
  output$text4 <- 
    renderText({ 
      
      d_input_5 = tibble("pbe1i01"         = input$m4_sap |> as.integer(),
                         "P_Health_dev_15" = input$m4_slt |> as.integer(),
                         "pfb7i02"         = input$m4_ini |> as.integer(),
                         "pda0i02"         = input$m4_dei |> as.integer()
      )
      
      # m_5 |>
      #   predict(d_input_4,type = "prob") |>
      #   pull(`.pred_ND-CNV`) |>
      #   round(digits = 3) %>%
      #   paste("predicted probability ND-GC = ",.,sep = "")
      
      m_5 |>
        predict(d_input_5,type = "prob") |>
        pull(`.pred_ND-CNV`) |>
        round(digits = 3) |>
        paste()
      
    })
  
  
  #Reset button on 5 item model
  observe({
    
    x <- input$reset_5

    updateRadioButtons(session = getDefaultReactiveDomain(),"m4_sap",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m4_slt",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m4_ini",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"m4_dei",selected = 0)
    
  })
  
  
  #LR model =====

  output$textLR <- 
    renderText({ 
    
    d_input_LR = tibble(
      #CAPA variables
      "pcb1i01" = input$mLR_fpb  |> as.integer(),
      "pcb3i01" = input$mLR_ago  |> as.integer(),
      "pce0i01" = input$mLR_fbi  |> as.integer(),
      "pcd2i01" = input$mLR_rum  |> as.integer(),
      "pda0i02" = input$mLR_dei  |> as.integer(),
      "pda0i03" = input$mLR_d2m  |> as.integer(),
      "pda1i01" = input$mLR_deq  |> as.integer(),
      "pbe1i01" = input$mLR_sap  |> as.integer(),
      "pbe7i01" = input$mLR_sai  |> as.integer(),
      "pbf8i01" = input$mLR_saf  |> as.integer(),
      "pbf2i01" = input$mLR_sas  |> as.integer(),
      "pfb7i01" = input$mLR_ini  |> as.integer(),
      "pfb7i02" = input$mLR_int  |> as.integer(),
      "prc8i01" = input$mLR_fgt  |> as.integer(),
      "prb8i01" = input$mLR_blt  |> as.integer(),
      "pgc3i01" = input$mLR_lie  |> as.integer(),
      "pgc5i01" = input$mLR_cht  |> as.integer(),
      
      #Health and Development Qs
      "P_Health_dev_9"  = input$mLR_rea  |> as.integer(),
      "P_Health_dev_11" = input$mLR_est  |> as.integer(),
      "P_Health_dev_13" = input$mLR_w18  |> as.integer(),
      "P_Health_dev_15" = input$mLR_slt  |> as.integer(),
      "P_Health_dev_21" = input$mLR_lun  |> as.integer(),
      "P_Health_dev_24" = input$mLR_car  |> as.integer(),
      "P_Health_dev_27" = input$mLR_msk  |> as.integer(),
      
      #SCQ
      "P_ASQ_7" = input$mLR_awm  |> as.integer(),
      "P_ASQ_8" = input$mLR_arp  |> as.integer(),
      
      #DCDQ
      "CMD_2"  = input$mLR_cba  |> as.integer(),
      "CMD_5"  = input$mLR_crf  |> as.integer(),
      "CMD_6"  = input$mLR_cob  |> as.integer(),
      "CMD_10" = input$mLR_ccp  |> as.integer()
    )
    
    
    
    m_lr |>
      predict(d_input_LR,type = "prob") |>
      pull(`.pred_ND-CNV`) |>
      round(digits = 3) |>
      paste()
    
  })
  
  #Reset button on 30 item model
  observe({
    
    x <- input$reset_lr
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_fpb",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_ago",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_fbi",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_rum",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_dei",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_d2m",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_deq",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_sap",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_sai",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_saf",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_sas",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_ini",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_int",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_fgt",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_blt",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_lie",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_cht",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_rea",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_est",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_w18",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_slt",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_lun",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_car",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_msk",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_awm",selected = 0)
    
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_arp",selected = 0)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_cba",selected = 1)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_crf",selected = 1)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_cob",selected = 1)
    updateRadioButtons(session = getDefaultReactiveDomain(),"mLR_ccp",selected = 1)
  })  
  
}


# Run the application 
shinyApp(ui = ui, server = server)

# Introduction ======

#In this script we load the raw data from the main database and clean it up

# Load packages and set up ======

#Required packages
pacman::p_load(tidyverse,readxl,patchwork)

#Set our threshold for proportion of similar responses
prop_similar_thres = 0.95

#Set our threshold of proportion missing for removing variables and
#participants 
prop_missing_thres = 0.25

# Load Data ======

#This makes a big list of the sheets used within the master database, each sheet being a different 
#psychiatric/health assessment

D <- lapply(c("IDs", 
              "CAPAADHD", 
              "CAPAWorries", 
              "CAPAPhobiasAnxAff",
              "CAPASepAnx", 
              "CAPAOCD", 
              "CAPAODDCD1",
              "CAPAODDCD2", 
              "CAPADepression",
              "CAPAHypomania",
              "CAPA TobAlcDrugs",
              "CAPASleep",
              "CAPASuicide",
              "CAPAOCD",
              "CAPAPsychosis",
              "CAPATrichTic",
              "EPQ01-6 FamEnvHealth QA",
              "EPQ12 DCDQ",
              "EPQ10 Social Comunication QA (…",
              "EPQ09 Strengths & Difficulties…"), 
            read_excel, 
            path = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/MASTERDATABASE_BE_09_11_18.xlsx", na=c("#NULL!","888","777","666","555","999"))


#Merge list items together to create dataframe and wrangle variable names by replacing "#" with "." to match names properly
D0 = 
  Reduce(function(x, y) merge(x, y, all = TRUE), D) |>
  as_tibble() |>
  mutate(X1 = -1) |>
  rename_with(~gsub("#",".",.x), .cols = everything()) 


#Read in Master participant list for confirmed genotypes
MPL <- read_excel("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/MASTER PARTICIPANT LIST.xlsx", 
                  sheet = "MASTER PARTICIPANT LIST")

#Select IMAGINE_ID individuals 
MPL = 
  MPL |>
  select( `Primary ID`, `Genotype code`, Confirmed,`IMAGINE (yes(1)/no(0)/F2Fconsent only(2))` ) |> 
  filter(`IMAGINE (yes(1)/no(0)/F2Fconsent only(2))` == 1 ) |>
  select(-c(`IMAGINE (yes(1)/no(0)/F2Fconsent only(2))`,Confirmed)) |>
  rename("IDs"="Primary ID", 
         "GenotypeCode"="Genotype code")

# table(MPL$GenotypeCode)


#Merge confirmed individuals with data to create our dataset
DF0 = inner_join(D0,MPL,by = "IDs")

# table(DF$GenotypeCode)

#Lightly wrangle by removing some variables
DF0 = 
  DF0 |>
  select(-c(X1,P_Health_dev_10,opp_incap,rum_incap)) |>
  mutate(P_Health_dev_4 = map(P_Health_dev_4,~ifelse(.x == "60-72",NA,.x))) |>
  mutate(P_Health_dev_4 = P_Health_dev_4 |> as.numeric())


#Reorder variables so they look clearer when printed in the command line
DF0 = 
  DF0 |>
  relocate(IDs,Gender,GenotypeCode)


rm(list = c("D","D0", "MPL"))



# Removing non-numeric variables ======

# Remove various missing codes 
DF = 
  DF0 |>
  mutate(across(.cols = everything(), ~ifelse(.x == 888,0,.x))) |>
  mutate(across(where(is.character) | where(is.numeric), 
                ~ifelse(.x %in% c("555","666", "777","999"),NA,.x))) 


#Drop CAPA variables other than intensity
toDrop_CAPA  = c("f0","v0","x0","d0","o0","O0","00")

#Drop other variables that need to be searched for
toDrop_Other = c("Imp","_symp","Fam","Life","FES","incap","Incap","Edu","diagnosis","severity","stata")

#Specific named variables to remove
toDrop_Named = c("IMAGINE_ID","Name","DOB",
                 "CAPAinterviewdateW1","CAPAinterviewdateW2","CAPAinterviewdateW3",
                 "IDs_L","IDs_T",
                 "mextrainfos",
                 "Other_CNV","OtherCNVcoded",
                 "Subject","Parent_DOB",
                 "D_Comp",
                 "enter1c","enter2c",
                 "Notes","DOB_Comment","Carer_DOB",
                 "P_Health_dev_8a","P_Health_dev_14","P_Health_dev_21a","P_Health_dev_22a", "P_Health_dev_23a","P_Health_dev_24a",
                 "P_Health_dev_25a","P_Health_dev_26a","P_Health_dev_27a","P_Health_dev_28a","P_Health_dev_29a","P_Health_dev_34fa")

DF = 
  DF |>
  
  #Drop CAPA variables other than intensity
  select(-contains(toDrop_CAPA)) |>
  
  #Drop SDQ impairment questions  
  select(-contains(toDrop_Other)) |>
  
  #In the handover document, they only take 2 pregnancy related variables, so lets do that directly
  #And they don't actually recode those variables either
  select(-(contains("Preg") & !contains(c("19","23")))) |>
  
  #Drop extra variables
  select(-all_of(toDrop_Named))


# table(DF$GenotypeCode)


# Recode variables =====

#At this stage, my understanding is that all our variables should be numeric, but some are
#showing up as characters; lets deal with that here
DF = 
  DF |>
  mutate(across(where(is.character) & !IDs,as.numeric))


#Take a little look at our variables
# skimr::skim(DF)


## ASQ variables ====

DF = 
  DF |>
  mutate(across(contains("ASQ"),~recode(.x, `2` = 0))) |>
  mutate(across(contains("ASQ"),~na_if(.x,99))) |>
  mutate(across(c(P_ASQ_2,P_ASQ_3,P_ASQ_10,P_ASQ_20:P_ASQ_40), 
                ~ifelse(.x ==1,0,1)))


## DCDQ Variables =====

#Reverse code DCDQ
DF = 
  DF |>
  mutate(across(contains("CMD"),
                ~recode(.x, `1`= 5, `2` = 4, `3` = 3, `4`= 2,`5`= 1))) 


## SDQ Variables ======

## SDQ variables recoding 
DF =
  DF %>% 
  mutate(across(contains("SDQ") & !contains("Imp"), ~recode(.x, `1`=0, `2`=1, `3`=2)))


## Recode Health and Development Variables =====

## Health and development variables recoding to 1/0 - apart from the Health-dev variables
#that are continuous or ordinal

DF =
  DF %>% 
  mutate(across(c(P_Health_dev_7:P_Health_dev_29,P_Health_dev_34f),~recode(.x, `1`=1, `2`=0, `22`=0))) |>
  mutate(P_Health_dev_6 = ifelse(P_Health_dev_6==6, NA, P_Health_dev_6))


## CAPA variables ====

# CAPA intensity variables recoding to 1/0
DF = 
  DF |>
  mutate(across(contains("i0"),~recode(.x, `3` = 1, `2` = 1, `0` = 0, .default = 999 ))) |>
  mutate(across(contains("i0"),~na_if(.x,999)))


# Recode another variable that has a strange coding where its only either 0 or 2, to 0 /1
DF = 
  DF |>
  mutate(pbd7e01 = recode(pbd7e01, `2` = 1))


# Filter Items =====


## Filter Variables based on responses ------

## Get the values of responses to all items - tidy way
Variables = 
  
  left_join(
    
    #This gets the names of every possible response
    DF |>
      summarise(across(everything(), ~.x |> unique() |> list() )) |>
      pivot_longer(everything(),names_to = "Variable",values_to = "Responses") |>
      mutate(n = map_dbl(Responses,length)),
    
    #This works out the proportion of each response out of all responses
    DF |>
      summarise(across(everything(), ~.x |> table() |> prop.table() |> list() )) |>
      pivot_longer(everything(),names_to = "Variable",values_to = "props") |>
      mutate(max_prop = map_dbl(props,max)), 
    
    by = c("Variable")) |>
  
  #Arrange the results nicely
  relocate(Variable,Responses,n,props,max_prop)


#We now have a tibble, "Variables", which contains all the variables in DF, the number of different responses 
#given to that item and the most common response to that item

#Look at the distribution of these
Variables |>
  ggplot(aes(max_prop)) +
  geom_histogram(bins = 20) +
  theme_bw()

#Now we do some selection of these items

# select only variables from DF based on those where most common response rate is  <=.95
DF1 = 
  DF |>
  select(Variables |> filter(max_prop <= prop_similar_thres) |> pull(Variable)) 

#At this stage we have 589 participants and 271 variables


## Filter items and participants based on proportional missingness ========

### Variable Missingness ====

## Calculate missing value rate for each variable we retain

Variables2 = 
  Variables |>
  filter(max_prop <= prop_similar_thres) |>
  left_join(
    DF1 |>
      summarise(across(everything(), ~.x |> is.na() |> sum()/nrow(DF1) )) |>
      pivot_longer(everything(),names_to = "Variable",values_to = "prop_miss") ,
    by = "Variable") 



#Visualise the proportions missing
p1 = 
  Variables2 |> 
  select(Variable,prop_miss) |>
  mutate(Variable = factor(Variable)) |>
  mutate(Variable = fct_reorder(Variable,prop_miss,min)) |>
  ggplot(aes(x = Variable, y = prop_miss)) +
  geom_col() +
  coord_flip()

p2 =
  Variables2 |> 
  select(Variable,prop_miss) |>
  ggplot(aes(prop_miss)) +
  geom_histogram(breaks = seq(0,1,0.01) )

p1 + p2

#So actually, we don't have many variables we lose, and also we can threshold at 0.25 

#Select only those variables which are < threshold missing
DF2 = 
  DF1 |>
  select(Variables2 |> 
           filter(prop_miss < prop_missing_thres) |> 
           pull(Variable))


### Participant Missingness =====

#Now we filter based on the proportion of missing data within 
#each individual participant
DF3 = 
  DF2 |>
  
  #Apply a function to each row
  rowwise() |>
  
  #Make sure all non-ID variables are numeric
  mutate(across(where(is.character) & !IDs,as.numeric)) |>
  
  #Work out proportion of na values
  mutate(prop_miss = c_across(where(is.numeric)) |> is.na() |> sum()/(ncol(DF2)-1)) |>
  relocate(IDs,Gender,GenotypeCode,prop_miss) |>
  ungroup()

#Plot the participant missingess
DF3 |>
  select(IDs,prop_miss) |>
  ggplot(aes(prop_miss)) +
  geom_histogram(breaks = seq(0,1,0.01) )


#Finally remove the participants with more than the threshold amount of missing data
DF3  = 
  DF3 |>
  filter(prop_miss < prop_missing_thres)

# table(DF3$GenotypeCode)

#This leaves us with 489 individuals, 236 variables (233 informative) and 113 controls




# Save Up =====


## write our filtered data for use in further analysis
write_csv(DF3 |> select(-prop_miss), "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/CleanedData.csv")


## Save a variable of only IDs and genotypes 
ID_Gens = 
  DF3 |> 
  select(IDs, GenotypeCode) |>
  add_count(GenotypeCode) |>
  select(-n)

#Tabulate the genotypes of all the included participants
#Controls are code 16
ID_Gens |> 
  count(GenotypeCode) |> 
  print(n = 23)

## Write out separate file of ids and genotypes
write_csv(ID_Gens, "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/ParticipantGenotypes.csv")
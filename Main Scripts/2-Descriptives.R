# Introduction =========

#This is the script that makes tables of the final selected participants and their demographics 


## load packages and data ======

#Get packages
pacman::p_load(tidyverse,tidymodels,readxl,gtsummary,lubridate)

#Set preferences
`%nin%` = Negate(`%in%`)


#Load the full dataset
DF = read_csv("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/CleanedData.csv")

#Make the class labels
D = 
  DF |>
  mutate(group = map_chr(GenotypeCode,~ifelse(.x == 16,"Control","ND-CNV"))) |>
  mutate(group = factor(group, levels = c("Control","ND-CNV"))) |>
  relocate(group)


# Participant Genotype Details ======


#Load the meaning of each genotype code
Gens <- read_excel("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/MASTER PARTICIPANT LIST.xlsx", sheet = 3)

# Load the genotype for each participant
MPL <- read_excel("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/MASTER PARTICIPANT LIST.xlsx", sheet = 1)


#Join the data and the genotype data together
D2 <-
  D |>
  left_join(Gens |> rename(GenotypeCode = Code),
            by = "GenotypeCode") |>
  relocate(IDs,group,GenotypeCode,CNV)
  

table_genotype_detail = 
  D2 |> 
  select(IDs, CNV) |>
  left_join(  MPL |>
                select(`Primary ID`, `Genotype as far as we know`,`Recruitment CNV`) |>
                rename(IDs = `Primary ID`),
              by = "IDs") |>
  mutate(CNV_detail = case_when(
    CNV == "More than 1 priority CNV" ~ `Genotype as far as we know`,
    CNV == "Other (non-priority)"     ~ `Genotype as far as we know`,
    TRUE ~ CNV
  )) |>
  mutate(CNV_detail = case_when(
    is.na(CNV_detail) ~ `Recruitment CNV`,
    TRUE ~ CNV_detail
  )) |>
  count(CNV_detail)|>
  arrange(-n) 


#Tabulate in the viewer
table_genotype_detail |>
  knitr::kable(col.names = c("ND-CNV","N"),format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


#Save for posterity
write_csv(table_genotype_detail, "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/FinalIncludedGenotypes.csv")


#Lets use the new spreadsheet from Sam Chawner

mpl  <- read_excel("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/ParticipantGenotypes.xlsx", sheet = 1)
gens <- read_excel("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/ParticipantGenotypes.xlsx", sheet = 2)

#First of all, lets confirm that all the controls are the same genotype in both spreadsheets - 
DF |>
  select(IDs,GenotypeCode)  |>
  janitor::clean_names() |>
  rename(id = i_ds) |>
  left_join(mpl |> 
              select(`Primary ID`,`Genotype as far as we know`, `Genotype code`) |>
              janitor::clean_names() |>
              rename(id = primary_id),
            by = "id")|>
  relocate(id,genotype_as_far_as_we_know, genotype_code.x,genotype_code.y) |>
  filter(str_detect(genotype_as_far_as_we_know,"Control")) |>
  print(n = 500)

#Yes they are!

#Now lets inspect our counts of genotypes
gens |> print(n = 30)

d_ids = DF |>
  select(IDs)  |>
  rename(id = IDs) |>
  pull(id)

d_genotypes = 
  mpl |> 
  select(`Primary ID`,`Genotype as far as we know`, `Genotype code`) |>
  janitor::clean_names() |>
  rename(id = primary_id) |>
  filter(id %in% d_ids) |>
  left_join(gens |> 
              janitor::clean_names() |> 
              rename(genotype_code = code) |>
              mutate(genotype_code = as.double(genotype_code))|>
              drop_na(genotype_code),
            by = "genotype_code") 


#OK so we now need to do a few extra things. We need to convert any groups is n < 5 into "other"
#We need to make a list of all the genotypes collapsed into the "other" group
#We need to decide what to do with the "more than 1 cnv" group (code 21)

#First lets look into group 21

d_genotypes |>
  filter(genotype_code == 21) |>
  group_by(genotype_as_far_as_we_know) |>
  count() |>
  arrange(-n)

#So all the participants flagged as multiple CNVs are n < 5 so they can go into the other group

d_genotypes |>
  
  #Lets add a count of all the participants in each group
  add_count(genotype_code) |>
  
  #We make a new variable and record when n < 5
  mutate(new_code = case_when(n < 5 ~ "Other (non-priority)",
                              genotype_code == 21 ~ "Other (non-priority)",
                              TRUE ~ cnv)) |>
  group_by(new_code) |>
  count() |>
  arrange(-n) |>
  knitr::kable(col.names = c("Genetic Condition","N"),format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)
  

#And lets make a big list of all the genotypes subsumed into that "other"
#group

d_g_other = 
  d_genotypes |>
  
  #Lets add a count of all the participants in each group
  add_count(genotype_code) |>
  
  #We make a new variable and record when n < 5
  mutate(new_code = case_when(n < 5 ~ "Other (non-priority)",
                              genotype_code == 21 ~ "Other (non-priority)",
                              TRUE ~ cnv)) |>
  group_by(new_code) |>
  filter(new_code == "Other (non-priority)") |>
  ungroup() |>
  distinct(genotype_as_far_as_we_know)

#How many deletions/duplications/mixtures/others?

d_g_other |> 
  mutate(type = case_when(
    str_detect(genotype_as_far_as_we_know,"deletion") & str_detect(genotype_as_far_as_we_know,"duplication") ~ "mixed",
    str_detect(genotype_as_far_as_we_know,"deletion") ~ "del",
    str_detect(genotype_as_far_as_we_know,"duplication") ~ "dup",

    TRUE ~ "other")
  ) |>
  count(type)
         
#Print all the free text
d_g_other|>
  pull(genotype_as_far_as_we_know) |>
  paste(sep = "; ",collapse = "; ")
  



# Participant Demographics =====


## read info from our database xlsx files
DBIDs <-
  read_xlsx("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/MASTERDATABASE_BE_09_11_18.xlsx", sheet = "IDs")

Fam <- 
  read_xlsx("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/MASTERDATABASE_BE_09_11_18.xlsx", sheet = "EPQ01-6 FamEnvHealth QA")



## Age and Gender =====


DBIDs <-
  DBIDs |> 
  select(IDs, DOB, Gender, CAPAinterviewdateW1) |>
  mutate(DOB = ymd(DOB),
         CAPAinterviewdateW1 = ymd(CAPAinterviewdateW1),
         Age = interval(DOB,CAPAinterviewdateW1) |> 
           as.period() |>
           as.numeric("years"),
         Gender = as.numeric(Gender)) 


d_demographic = 
  D2 |> 
  select(IDs,group,GenotypeCode,CNV,Gender) |>
  left_join(DBIDs, by = c("IDs","Gender")) |>
  mutate(Gender = case_when(
    Gender == 1 ~ "Male",
    Gender == 2 ~ "Female"
  ))

# table(d_demographic$group,d_demographic$Gender)


## Educational level ====


#Now we work out educational levels, income and ethnic background

d_family =
  d_demographic |> 
  select(-c(DOB,CAPAinterviewdateW1)) |>
  left_join(Fam |> select(IDs,
                          P_Edu_1a, P_Edu_1b, P_Edu_1c, P_Edu_1d, P_Edu_1e,
                          P_Edu_2,
                          P_Fam_Back_1), 
            by = "IDs") |>
  mutate(Family = map_chr(IDs,~str_split(.x,pattern = "-",n = 4,simplify = TRUE) %>% .[1])) |>
  relocate(Family,.before = group)


#This is all wrangling stuff
Education <- 
  d_family |>
  select(IDs,Family,group,P_Edu_1a, P_Edu_1b, P_Edu_1c, P_Edu_1d, P_Edu_1e) |>
  
  #We want to know if all the variables are NA

  #Replace values with logicals based on the rule that only "1" values are true
  mutate(across(contains("P_Edu"),~case_when(
    . == "1" ~ T,
    . == "2" ~ F,
    is.na(.) ~ NA,
    TRUE     ~ NA
  ))) |>
  
  #Apply some conditions to determine the highest educational level
  mutate(Highest = case_when( P_Edu_1a & !P_Edu_1b & !P_Edu_1c & !P_Edu_1d & !P_Edu_1e ~ "Low",
                              P_Edu_1a & is.na(P_Edu_1b) & is.na(P_Edu_1c) & is.na(P_Edu_1d) & is.na(P_Edu_1e) ~ "Low",
                              P_Edu_1b & !P_Edu_1c & !P_Edu_1d & !is.na(P_Edu_1e)  ~ "Middle",
                              P_Edu_1a & P_Edu_1e ~ "Middle",
                             (P_Edu_1c | P_Edu_1d) ~ "High",
                              P_Edu_1e ~ "Middle",
                             !P_Edu_1a & !P_Edu_1b & !P_Edu_1c & !P_Edu_1d & !P_Edu_1e ~ "No School Leaving Exams",
                             TRUE ~ "Unknown"))




#So we address controls often have missing info by duplicating across from the affected siblings
Education = 
  d_family |>
  select(IDs,Family,group) |>
  left_join(Education |>
              select(Family, Highest) |>
              mutate(DupFam = ifelse(duplicated(Family), 1, 0)) |>
              filter(DupFam==0),
              by = c("Family")) |>
  select(-DupFam)

table(Education$group, Education$Highest)

## Income =====

Income <- 
  d_family |>
  select(IDs,Family,group,P_Edu_2) |>
  rename(Income = P_Edu_2) |>
  mutate(Income = case_when(
    Income == "1"| Income == "2" | Income == "8" ~ "<=£19,999",
    Income == "3"| Income == "4" ~ "£20,000 - £39,999",
    Income == "5"| Income == "6" ~ "£40,000 - £59,999",
    Income == "7" ~ "£60,000 + ",
    TRUE ~ "Unknown"
  ))

#Values not covered by the case_when conditions are replaced with NAs

Income = 
  d_family |>
  select(IDs,Family,group) |>
  left_join(Income |>
              select(Family, Income) |>
              mutate(DupFam = ifelse(duplicated(Family), 1, 0)) |>
              filter(DupFam==0),
            by = c("Family")) |>
  select(-DupFam)

table(Income$group, Income$Income)


## Ethnicity =====

Background<-
  d_family |>
  select(IDs,Family,group,P_Fam_Back_1) |>
  rename(Ethnicity = P_Fam_Back_1) |>
  mutate(Ethnicity = case_when(
    (Ethnicity == "1" | Ethnicity == "2")  ~ "European",
    (Ethnicity == "4" | Ethnicity == "8" |Ethnicity == "10" )  ~ "Other",
    TRUE ~ "Unknown"
  )) |>
  
  #There are some manual replacements
  mutate(Ethnicity = case_when(
    Family=="E143"  ~ "European",
    Family=="E154"  ~ "Other",
    Family=="E268"  ~ "European",
    Family=="E269"  ~ "European",
    Family=="E272"  ~ "European",
    Family=="E298"  ~ "European",
    Family=="E332"  ~ "European",
    Family=="E344"  ~ "European",
    Family=="E272"  ~ "European",
    Family=="I101"  ~ "European",
    Family=="I129"  ~ "European",
    Family=="IM172" ~ "European",
    Family=="IM185" ~ "European",
    Family=="IM188" ~ "European",
    Family=="IM228" ~ "European",
    Family=="IM243" ~ "European",
    TRUE ~ Ethnicity
  )) |>
  
  #Set factor levels for our possible ethnicity categories
  mutate(Ethnicity = factor(Ethnicity, levels = c("European", "Other", "Unknown")))



Background = 
  d_family |>
  select(IDs,Family,group) |>
  left_join(Background |>
              select(Family, Ethnicity) |>
              mutate(DupFam = ifelse(duplicated(Family), 1, 0)) |>
              filter(DupFam==0),
            by = c("Family")) |>
  select(-DupFam)

table(Background$group, Background$Ethnicity)



# Combine datasets and tabulate =====



#Mash together our datasets
d_table = 
  d_demographic |> 
  
  #Grab our basic demograhpic details
  select(IDs,group,Age,Gender) |>
  
  #Glue on educational detail
  left_join(
    Education |>
      select(IDs,  Highest),
    by = "IDs") |>
  
  #Glue on family income
  left_join(Income |>
              select(IDs, Income),
            by = "IDs") |>
  
  # Glue on ethnicity data
  left_join(Background |>
              select(IDs, Ethnicity),
            by = "IDs")



# Make a nice summary table ====

d_table |>
  mutate(Highest = factor(Highest, levels = c("No School Leaving Exams","Low","Middle","High","Unknown") ),
         Income  = factor(Income , levels = c("<=£19,999","£20,000 - £39,999","£40,000 - £59,999","£60,000 + ","Unknown"))) |>
  rename(`Highest Educational Level` = Highest,Group = group) |>
  select(-IDs) |>
  tbl_summary(by = c(Group)) %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Group**") %>%
  bold_labels() %>%
  add_p()



# Carrier/control family pairs =====

#We need to work out how many participants with a CNV also have a sibling included in the study versus who does not
Education |>
  select(-Highest) |>
  group_by(Family) |>
  tally() |>
  arrange(-n) |>
  ggplot(aes(n)) +
  geom_histogram(binwidth = 1)

#What proportion of CNV carriers have a sibling included?
Education |>
  select(-Highest) |>
  left_join(Education |>
              select(-Highest) |>
              group_by(Family) |>
              tally(),
            by = "Family") |>
  filter(group == "ND-CNV") %>%
  group_by(n) %>%
  summarise(f_count = n(),
            f_prop = (n()/nrow(.)) |> round(digits = 2))

#What about the other way around?
Education |>
  select(-Highest) |>
  left_join(Education |>
              select(-Highest) |>
              group_by(Family) |>
              tally(),
            by = "Family") |>
  filter(group == "Control") %>%
  group_by(n) %>%
  summarise(f_count = n(),
            f_prop = (n()/nrow(.)) |> round(digits = 2))


# Study dates =====

#We can get the dates that the questionairres were completed

q_dates <- read_excel("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/MASTERDATABASE_BE_09_11_18.xlsx", 
                  sheet = "IDs")

q_dates = 
  D2 |> 
  select(IDs,group,CNV,Gender) |>
  left_join(q_dates |>
              select(IDs,CAPAinterviewdateW1), 
            by = "IDs")

#Get the earliest and latest dates-
q_dates |> 
  summarise(min_date = min(CAPAinterviewdateW1),
            max_date = max(CAPAinterviewdateW1))

# Figure 1 In R ======


# Make figure 1 as a R plot?

pacman::p_load(DiagrammeR,DiagrammeRsvg,xml2)

#Try our own
my_graph = 
grViz("
digraph cnv_classication_flow {

  # a 'graph' statement
  graph [layout = dot, compound = true, 
        rankdir = TB,
         nodesep = 1, ranksep = .25,
         color = crimson]

  # several 'node' statements
  node [fontname = Arial, fontcolor = darkslategray,
        shape = rectangle, fixedsize = true, width = 3,
        color = darkslategray]


  # several 'edge' statements
  edge [color = grey]
  
  # subgraph for Training Data Information
  subgraph cluster0 {
    node [fixedsize = true, width = 3, height = 1]
    '@@3' -> '@@4' -> '@@5' -> '@@6' -> '@@7'
    
  }

  # subgraph for Testing Data 
  subgraph cluster1 {
    node [fixedsize = true, width = 3, height = 1]
    '@@1' -> '@@2'
  }
  
  node [fixedsize = true, width = 3, height = 1]
  a[label = 'Starting Dataset: 589 participants, \\n 1450 variables']   
  b[label = 'Remove administrative, free text, \\n data and time variables, \\n leaving 478 variables']   
  c[label = 'Remove variables where > 95% \\n of responses are identifcal, \\n leaving 271 variables']   
  d[label = 'Remove variables and participants \\n with > 25% missing data \\n leaving 489 participants and \\n 233 variables']

  a -> b -> c -> d


  d -> '@@3'       [lhead = cluster0]
  d -> '@@1'       [lhead = cluster1]
  
  '@@7' -> '@@2'
  
  node [fixedsize = true, width = 3]
  '@@7' -> '@@8' 
  '@@8' -> '@@9'  
  

}

[1]: 'Test Data: 99 participants'
[2]: 'Evaluation of best model  \\n performance on test data'
[3]: 'Training Data: 390 participants'
[4]: 'Nested CV with all variables: \\n  imputation, ML model fitting \\n and performance evaluation'
[5]: 'Select top 30 variables for each \\n model plus variables \\n selected by multiple models'
[6]: 'Nested CV with variable subsets'
[7]: 'Evaluation of performance with \\n best variable subsets'
[8]: 'Variable Importance for \\n best variable set'
[9]: 'Exploratory Graph Analysis \\n of best variable set'

") 


#And save
my_graph %>%
  export_svg() %>%
  read_xml() %>%
  write_xml("figure_1.svg")


#And then we wrangle this lightly in Inkscape to get it perfect


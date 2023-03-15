# Introduction =========

#This is the script that makes tables of the final selected participants and their demographics 


## load packages and data ======

#Get packages
pacman::p_load(tidyverse,tidymodels,readxl,gtsummary,lubridate)

#Set preferences
`%nin%` = Negate(`%in%`)


#Set our data directory
if(str_detect(Sys.info()[['nodename']],'IT088825')){
  data_dir = 'C://Users/pyxnd/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data'
  
}else if(str_detect(Sys.info()[['nodename']],'AVOCADO')){
  
  data_dir = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/"
}

setwd(data_dir)



#Load the full dataset
DF = read_csv("CleanedData.csv")


# Participant Genotype Details ======

#Load the latest genotype information from Sam Chawner

mpl  <- read_excel("ParticipantGenotypes.xlsx", sheet = 1)
gens <- read_excel("ParticipantGenotypes.xlsx", sheet = 2)

d_ids = 
  DF |>
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


## Table 2: Genotypes ======

#Now we make the Table - currently Table 2
d_genotypes |>
  
  #Lets add a count of all the participants in each group
  add_count(genotype_code) |>

  #We make a new variable and record when n < 5
  mutate(new_code = case_when(n < 5 ~ "Other (non-priority)",
                              genotype_code == 21 ~ "Other (non-priority)",
                              TRUE ~ cnv)) |>
  group_by(new_code) |>
  count() |>
  mutate(new_code = case_when(new_code == "Other (non-priority)" ~ "Other",
                              TRUE ~ new_code)) |>
  arrange(-n) |>
  knitr::kable(col.names = c("Genomic Condition","N"),format = "html", booktabs = TRUE) |>
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

#How many deletions/duplication/mixtures/others?
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
  
#We then put a sorted form of this into the Table 2 caption


#Including the others and the non-controls, how many genotypes do we have?
d_g_other |> 
  mutate(type = case_when(
    str_detect(genotype_as_far_as_we_know,"deletion") & str_detect(genotype_as_far_as_we_know,"duplication") ~ "mixed",
    str_detect(genotype_as_far_as_we_know,"deletion") ~ "del",
    str_detect(genotype_as_far_as_we_know,"duplication") ~ "dup",
    
    TRUE ~ "other")
  ) |>
  count(type) |>
  summarise(sumn = sum(n) + 14)

#We add 14 as there are 14 non contorl, non-other genotypes in our table above


# Participant Demographics =====


## read info from our database xlsx files
DBIDs <-
  read_xlsx("MASTERDATABASE_BE_09_11_18.xlsx", 
            sheet = "IDs")

Fam <- 
  read_xlsx("MASTERDATABASE_BE_09_11_18.xlsx", 
            sheet = "EPQ01-6 FamEnvHealth QA")



## Age and Gender =====


DBIDs <-
  DBIDs |> 
  select(IDs, DOB, Gender, CAPAinterviewdateW1) |>
  mutate(DOB = ymd(DOB),
         CAPAinterviewdateW1 = ymd(CAPAinterviewdateW1),
         Age = interval(DOB,CAPAinterviewdateW1) |> 
           as.period() |>
           as.numeric("years"),
         Gender = as.numeric(Gender)) |>
  janitor::clean_names() |>
  rename(id = i_ds,
         capa_interviewdate_w1 = cap_ainterviewdate_w1) |>
  select(-dob)


d_demographic = 
  d_genotypes  |> 
  select(id,genotype_code) |>
  mutate(group = if_else(genotype_code == 16,"control","ND-GC")) |>

  left_join(DBIDs, by = c("id")) |>
  mutate(gender = case_when(
    gender == 1 ~ "Male",
    gender == 2 ~ "Female"
  )) |>
  select(-capa_interviewdate_w1)

# table(d_demographic$group,d_demographic$gender)



## Educational level ====

#Now we work out educational levels, income and ethnic background

d_family =
  d_demographic |> 
  left_join(Fam |> select(IDs,
                          P_Edu_1a, P_Edu_1b, P_Edu_1c, P_Edu_1d, P_Edu_1e,
                          P_Edu_2,
                          P_Fam_Back_1) |>
              rename(id = IDs), 
            by = "id") |>
  mutate(Family = map_chr(id,~str_split(.x,pattern = "-",n = 4,simplify = TRUE) %>% .[1])) |>
  relocate(Family,.before = group)


#This is all wrangling stuff
Education <- 
  d_family |>
  select(id,Family,group,P_Edu_1a, P_Edu_1b, P_Edu_1c, P_Edu_1d, P_Edu_1e) |>
  
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
  Education |>
  select(id,Family,group,Highest) |>
  group_by(Family) |>
  nest() |>
  mutate(highest_new = map_chr(data,~if_else(dim(.x)[1] > 1,
                                             .x$Highest[!str_detect(.x$Highest,'Unknown')][1],
                                             .x$Highest[[1]]))) |>
  unnest(data) |>
  select(-Highest) |>
  rename(Highest = highest_new) |> 
  ungroup() |>
  mutate(Highest = case_when(is.na(Highest) ~ "Unknown",
                              TRUE ~ Highest )) |>
  mutate(Highest = factor(Highest,levels = c("No School Leaving Exams","Low","Middle","High","Unknown")))

table(Education$group, Education$Highest)



## Income =====

Income <- 
  d_family |>
  select(id,Family,group,P_Edu_2) |>
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
  Income |>
  arrange(Family) |>
  select(id,Family,group,Income) |>
  group_by(Family) |>
  nest() |>
  mutate(highest_new = map_chr(data,~if_else(dim(.x)[1] > 1,
                                             .x$Income[!str_detect(.x$Income,'Unknown')][1],
                                             .x$Income[[1]]))) |>
  unnest(data) |>
  select(-Income) |>
  rename(Income = highest_new) |>
  
  ungroup() |>
  mutate(Income = case_when(is.na(Income) ~ "Unknown",
                             TRUE ~ Income )) |>
  
  #Set factor levels for our possible ethnicity categories
  mutate(Income = factor(Income, levels = c("<=£19,999", "£20,000 - £39,999","£40,000 - £59,999","£60,000 + ", "Unknown")))


table(Income$group, Income$Income)


## Ethnicity =====

Background<-
  d_family |>
  select(id,Family,group,P_Fam_Back_1) |>
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
  )) 



Background = 
  Background |>
  arrange(Family) |>
  select(id,Family,group,Ethnicity) |>
  group_by(Family) |>
  nest() |>
  mutate(highest_new = map_chr(data,~if_else(dim(.x)[1] > 1,
                                             .x$Ethnicity[!str_detect(.x$Ethnicity,'Unknown')][1],
                                             .x$Ethnicity[[1]]))) |>
  unnest(data) |>
  select(-Ethnicity) |>
  rename(Ethnicity = highest_new) |>
  ungroup() |>
  mutate(Ethnicity = case_when(is.na(Ethnicity) ~ "Unknown",
                             TRUE ~ Ethnicity )) |>
  
  #Set factor levels for our possible ethnicity categories
  mutate(Ethnicity = factor(Ethnicity, levels = c("European", "Other", "Unknown")))


table(Background$group, Background$Ethnicity)


# Combine datasets and tabulate =====


#Mash together our datasets
d_table = 
  
  d_demographic |> 
  
  #Grab our basic demograhpic details
  select(id,group,age,gender) |>
  
  #Glue on educational detail
  left_join(
    Education |>
      select(id,  Highest) ,
    by = "id") |>
  
  #Glue on family income
  left_join(Income|>
              select(id, Income),
            by = "id") |>
  
  # Glue on ethnicity data
  left_join(Background|>
              select(id, Ethnicity),
            by = "id") |>
  arrange(id)



# Make a nice summary table ====

d_table |>
  rename(`Highest Educational Level` = Highest,
         Group = group,Gender = gender,Age = age) |>
  select(-id) |>
  tbl_summary(by = c(Group)) %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Group**") %>%
  bold_labels()

#We are not going to do stats tests as there is literature suggesting this isnt a good idea
#%>%  add_p()


# Save ======

write_csv(d_table, "DemographicTable.csv")



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
  filter(group == "ND-GC") %>%
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
  filter(group == "control") %>%
  group_by(n) %>%
  summarise(f_count = n(),
            f_prop = (n()/nrow(.)) |> round(digits = 2))


# Study dates =====

#We can get the dates that the questionnaires were completed

q_dates <- read_excel("MASTERDATABASE_BE_09_11_18.xlsx", 
                  sheet = "IDs")

q_dates = 
  d_table |> 
  select(id,group,gender) |>
  left_join(q_dates |>
              select(IDs,CAPAinterviewdateW1) |>
              rename(id = IDs), 
            by = "id")

#Get the earliest and latest dates-
q_dates |> 
  summarise(min_date = min(CAPAinterviewdateW1),
            max_date = max(CAPAinterviewdateW1))






# Figure 1 In R ======


# Make figure 1 as a R plot
pacman::p_load(DiagrammeR,DiagrammeRsvg,xml2)

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
  a[label = 'Starting Dataset: 589 participants, \\n 1451 variables']   
  b[label = 'Remove administrative, free text, \\n data and time variables, \\n leaving 474 variables']   
  c[label = 'Remove variables where > 90% \\n of responses are identifcal, \\n leaving 211 variables']   
  d[label = 'Remove variables and participants \\n with > 25% missing data \\n leaving 493 participants and \\n 192 variables']
  e[label = 'Remove variables with > 0.80 \\n correlation with other variables \\n leaving 493 participants and \\n 176 variables']

  a -> b -> c -> d -> e


  e -> '@@3'       [lhead = cluster0]
  e -> '@@1'       [lhead = cluster1]
  
  '@@7' -> '@@2'
  
  node [fixedsize = true, width = 3]
  '@@7' -> '@@8' 
  '@@8' -> '@@9'  
  

}

[1]: 'Test Data: 100 participants'
[2]: 'Evaluation of best model  \\n performance on test data'
[3]: 'Training Data: 393 participants'
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


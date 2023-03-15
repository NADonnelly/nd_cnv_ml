# Introduction ======

#Lets load up
pacman::p_load(tidyverse,tidymodels,
               gtsummary,patchwork,
               readxl,lubridate, crayon)

conflicted::conflicts_prefer(crayon::`%+%`)

#Set our data directory
if(str_detect(Sys.info()[['nodename']],'IT088825')){
  data_dir = 'C://Users/pyxnd/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data'
  
}else if(str_detect(Sys.info()[['nodename']],'AVOCADO')){
  
  data_dir = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/"
}

setwd(data_dir)

#Load the full dataset
DF = read_csv("CleanedData.csv")


#Load the demographic information for use in the stratification
d_demo = read_csv("DemographicTable.csv")


#Set some preferences
tidymodels_prefer()
`%nin%` = Negate(`%in%`)


# Initial data splits ========

#Make the class labels
DF0 = 
  DF |>
  mutate(group = map_chr(GenotypeCode,~ifelse(.x == 16,"Control","ND-GC"))) |>
  mutate(group = factor(group, levels = c("Control","ND-GC"))) |>
  relocate(group,.after = GenotypeCode)

#Have as quick look at all variables
# skimr::skim(DF0)

# DF0 |>
#   summarytools::dfSummary() |>
#   summarytools::view()


# Force ordinal variables to be integers (necessary for correct imputation)

## work out if a variable is all integers
is_int_var = 
  DF0 |> 
  pivot_longer(-c(IDs,GenotypeCode,group)) |> 
  group_by(name) |> 
  nest() |>
  mutate(all_int = map_lgl(data, ~ .x |>
                             drop_na() |> 
                             mutate(frac = map_lgl(value,~.x%%1 == 0)) |> 
                             summarise(is_int = all(frac)) |>
                             pull())) |>
  select(name,all_int) 

#Now change to integers wherever appropriate
DF1 = 
  DF0 |>
  mutate(across(is_int_var |> filter(all_int == T) |> pull(name), as.integer)) 



#Set a random number seed for reproducibility
set.seed(09062022)


#Set our data splitting ratios and number of folds
split_prop = 4/5


#Prepare a variable that combines multiple demographic features for stratification

d_demo = 
  d_demo |>
  mutate(age_strata = rsample::make_strata(age,breaks = 4)) |>
  mutate(multi_demo = interaction(group,age_strata,gender)) 


#Prepare our initial data splits                                  
d_split = 
  DF1 |>
  left_join(d_demo |>
              select(id,multi_demo) |>
              rename(IDs = id),
            by = "IDs") |>
  relocate(multi_demo,.after = group) |>
  initial_split(strata = multi_demo,prop = split_prop)


#Save the split data for use elsewhere and consistency
write_rds(d_split ,"nested_cv_split_data.rds")
#d_split = read_rds("nested_cv_split_data.rds")

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)


#Get the IDs of the participants in the testing split
d_test_id = d_split |> complement()


# Split data table =======

#Check the characteristics of the test vs training data for balance
bind_rows(d_train|>
            select(IDs) |>
            left_join(d_demo |>
                        select(-c(age_strata,multi_demo))|>
                        rename(IDs = id),
                      by = "IDs") |>
            mutate(set = "train"),
          d_test|>
            select(IDs) |>
            left_join(d_demo |>
                        select(-c(age_strata,multi_demo))|>
                        rename(IDs = id),
                      by = "IDs") |>
            mutate(set = "test")) |>
  select(-IDs) |>
  tbl_summary(by = c(set)) %>%
  add_overall() %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Dataset**") %>%
  bold_labels()



# Impute full training dataset =-======

#Set up a recipe with all the variables all at once, with imputation
rec_impute =
  d_train |>
  select(-c(IDs,Gender,GenotypeCode,multi_demo)) |>
  recipe(group ~ .) |>
  step_zv(all_predictors()) |>
  step_impute_bag(all_predictors())   


#Prep a model for imputing missing data
rec_impute_prep = 
  rec_impute |>
  prep()

#Get our imputed data back
d_train_impute = 
  rec_impute_prep |>
  juice()


#And make a test set using the same imputation model, but the test (assessment) 
#data
d_test_impute = 
  rec_impute_prep |>
  bake(new_data = d_test)


#And combine - we rebuild a split object with the imputed training and test data
d_last_fit = 
  make_splits(d_train_impute,d_test_impute)

#Lets save this so we can use it again
write_rds(d_last_fit ,"nested_cv_imputed_data.rds")
# d_last_fit = read_rds("nested_cv_imputed_data.rds")

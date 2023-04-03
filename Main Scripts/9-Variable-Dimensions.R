# Introduction ========

# Lets explore the final selected set of variables


#Lets load up
pacman::p_load(tidyverse,tidymodels,
               naniar,rgl,polycor,Matrix,
               parallel,doSNOW,tcltk,
               ggcorrplot,
               qgraph, EGAnet, EFATools ,
               patchwork)

library(EFAtools)

tidymodels_prefer()



#Set our data directory
if(str_detect(Sys.info()[['nodename']],'IT088825')){
  data_dir = 'C://Users/pyxnd/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data'
  
}else if(str_detect(Sys.info()[['nodename']],'AVOCADO')){
  
  data_dir = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/"
}

setwd(data_dir)



#Load the imputed dataset
d_last_fit = read_rds("nested_cv_imputed_data.rds")



## Prepare variable definitions =====


#Load our final variables
d_var = read_rds("nested_cv_selected_var_definitions.rds")


#We can try and get variable definitions from the spreadsheet found on the 
#imagine-id website?

#We can look these items up in the data dictionary....
vd <- readxl::read_excel("DATA DICTIONARY.xlsx", 
                         sheet = "Sheet1")

vl <- readxl::read_excel("DATA ENTRY CODEBOOK .xlsx", 
                         sheet = "VARIABLE LIST")

#Wrangle lightly
vl = 
  vl |>
  janitor::clean_names() |>
  select(-x4) |>
  filter(variable %in% d_var$variable)

vd = 
  vd |> 
  janitor::clean_names() |>
  select(question_number,question_description,assessment,section,question_description_me) |>
  drop_na(question_description) |>
  mutate(question_number = map_chr(question_number,str_remove,"_L")) |>
  filter(question_number %in% d_var$variable) |>
  rename(variable = question_number) |>
  mutate(variable_definition_dict  = case_when(is.na(question_description_me) ~ question_description,
                                           TRUE ~ question_description_me)) |>
  select(variable,assessment,section,variable_definition_dict) 


variable_table = 
  left_join(vl,vd,by = "variable") |>
  select(variable,data_type,assessment,section,variable_definition,variable_definition_dict) 



#We might benefit from make short names for each variable.
variable_table |>
  print(n = 30)


#What do we got? We have 30 variables in total. Note that these are somewhat different to the variables first identified
#using the SVM model (which is not unexpected given they are rather different models)

#CAPA 
# 1  pcb2i01          Agoraphobia intensity                                                               = AGO
# 2  pcc0i01          Situational anxious affect intensity                                                = SIT
# 3  pbf4i01          Avoidance of being alone intensity                                                  = ALO
# 4  pbf5i01          Anticipatory distress intensity                                                     = ANT
# 5  pfb7i02          Sleep problems -  Initial insomnia intensity                                        = INI ***  
# 6  prb8i01          Often blurts out answers to questions                                               = BLT ***
# 7  pge2i01          Vandalism intensity                                                                 = VAN

# Health & Development
# 8  P_Preg_23        How much did your child weigh at birth?                                             = WGT
# 9  P_Health_dev_6   Is your child clumsy?                                                               = CLM
# 10 P_Health_dev_9   Is your child behind in reading                                                     = REA ***
# 11 P_Health_dev_11  Is your child educationally statemented                                             = EST ***
# 12 P_Health_dev_12  Was your child talking by the age of 2                                              = SP2
# 13 P_Health_dev_15  Has your child had speech therapy                                                   = SLT ***
# 14 P_Health_dev_20  Frequent infections of the chest/airways                                            = RTI

# SDQ
# 15 P_SDQ_1          Considerate of other's feelings                                                     = CNS
# 16 P_SDQ_6          Rather solitary, tends to play alone                                                = SOL
# 17 P_SDQ_9          Helpful if someone is hurt                                                          = HRT
# 18 P_SDQ_16         Easily distracted, concentration wanders                                            = DIS
# 19 P_SDQ_19         Often tells lies                                                                    = LIE
# 20 P_SDQ_20         Often cheats                                                                        = CHT
# 21 P_SDQ_29         Inconsiderate of others                                                             = INC

#SCQ
# 22 P_ASQ_14         Special interests unusual in their intensity                                        = SII
# 23 P_ASQ_39         Imaginative play with another child                                                 = PIM
# 24 P_ASQ_40         Play cooperatively                                                                  = PCO

#DCDQ
# 25 CMD_2            Catches a small ball thrown from 6-8ft                                              = CBA ***
# 26 CMD_5            Runs as fast and easily as other children                                           = CRF ***
# 27 CMD_6            Can organise her body to do a planned motor activity                                = COB ***
# 28 CMD_11           Likes participating in games requiring good motor skills                            = CGM
# 29 CMD_14           Child would never be described as a bull in a china shop (clumsy and breaks things) = CBL
# 30 CMD_15           Child does not fatigue easily or appear to slouch and "fall out" of the chair       = CSL

# The *** indicates the variables that were selected in the first iteration of the manuscript

#Lets smush these together
variable_table = 
  variable_table|>
  
  #This is my rather arbitrary variable naming
  mutate(short_name = case_when(variable == "pcb2i01" ~ "AGO",
                                variable == "pcc0i01" ~ "SIT",
                                variable == "pbf4i01" ~ "ALO",
                                variable == "pbf5i01" ~ "ANT",
                                variable == "pfb7i02" ~ "INI",
                                variable == "prb8i01" ~ "BLT",
                                variable == "pge2i01" ~ "VAN",
                                variable == "P_Preg_23" ~ "WGT",
                                variable == "P_Health_dev_6" ~ "CLM",
                                variable == "P_Health_dev_9" ~ "REA",
                                variable == "P_Health_dev_11" ~ "EST",                              
                                variable == "P_Health_dev_12" ~ "SP2",
                                variable == "P_Health_dev_15" ~ "SLT",
                                variable == "P_Health_dev_20" ~ "RTI",
                                variable == "P_SDQ_1" ~ "CNS",
                                variable == "P_SDQ_6" ~ "SOL",
                                variable == "P_SDQ_9" ~ "HRT",
                                variable == "P_SDQ_16" ~ "DIS",
                                variable == "P_SDQ_19" ~ "LIE",
                                variable == "P_SDQ_20" ~ "CHT",
                                variable == "P_SDQ_29" ~ "INC",
                                variable == "P_ASQ_14" ~ "SII",
                                variable == "P_ASQ_39" ~ "PIM",
                                variable == "P_ASQ_40" ~ "PCO",
                                variable == "CMD_2" ~ "CBA",
                                variable == "CMD_5" ~ "CRF",
                                variable == "CMD_6" ~ "COB",
                                variable == "CMD_11" ~ "CGM",
                                variable == "CMD_14" ~ "CBL",
                                variable == "CMD_15" ~ "CSL")) |>
  
  mutate(variable_definition_dict = case_when(is.na(variable_definition_dict) ~ variable_definition,
                                              TRUE ~ variable_definition_dict)) |>
  select(-variable_definition)


#Lets make a column that tells us if a variable is binary, ordinal or continuous (we include this information
#in the table in the paper)
variable_table = 
  d_last_fit |> 
  training() |> 
  select(all_of(variable_table$variable)) |> 
  pivot_longer(everything()) |> 
  nest_by(name) |> 
  mutate(nd = dim(unique(data))[1]) |>
  mutate(var_type = case_when(nd == 2 ~ "binary",
                              nd >2 & nd <=5 ~ "ordinal",
                              nd > 5 ~ "continuous")) |>
  rename(variable = name) |>
  select(variable,var_type) |>
  left_join(variable_table,by = "variable") |>
  select(variable,short_name,data_type,var_type,assessment,section,variable_definition_dict) |>
  arrange(variable)


#Now we save this
write_csv(variable_table,"./VariableDefinitionsExpanded.csv")


#The idea is that we then get the full variable definitions from Jess and prepare a nice sheet of these definitions and the data  
#types
vd2 <- readxl::read_excel("VariableDefinitions.xlsx", 
                          sheet = "Sheet1")




# Correlation between variables ========

# Lets look at correlation between variables

#Do we use the training set only or both the training and test sets?
DF = bind_rows(d_last_fit |> analysis(),
               d_last_fit |> assessment())

DF = DF |> relocate(group)


## Do Polychor Correlations ======

#Do the polychoric crrelation between all pairs of variables using the polychor package

## Set up dataframe for polychoric correlation calculation
D3 = 
  DF |>
  select(-group)

#Or we can use the hetcor package??
r3.df = 
  D3 |>
  mutate(across(.cols = everything(),~as.numeric(.x))) |>
  as.matrix() |>
  polycor::hetcor(use = "pairwise.complete.obs",
                  parallel = TRUE,
                  ncores = 12)

ggcorrplot::ggcorrplot(r2.df)
ggcorrplot::ggcorrplot(r3.df$correlations)

#Definitely similar 

#Save this
write_csv(r3.df, "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_corr_mat.csv")
# r3.df = read_csv( "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_corr_mat.csv")


## Extract our variables of interest =====

#We did the correlation on all variables
corr = 
  r3.df |>
  as_tibble() |>
  mutate(variable = colnames(D3)) |>
  relocate(variable) |>
  select(variable,all_of(d_var$variable)) |>
  filter(variable %in% d_var$variable) |>
  mutate(variable = factor(variable,levels = d_var$variable)) |>
  arrange(variable)

## Correlation plot
corr.mat = 
  corr  |> 
  select(-variable) |> 
  as.matrix()

row.names(corr.mat) = variable_table$short_name
colnames(corr.mat)  = variable_table$short_name

ggcorrplot(corr.mat, hc.order = TRUE, outline.color = "white")



# Exploratory Graph Analysis ======


#As in the original draft of the paper, we can do EGA

#Can we switch out our variable names for our short names
name_swap        = variable_table$variable
names(name_swap) = variable_table$short_name

d_ega = 
  bind_rows(d_last_fit |> analysis(),
            d_last_fit |> assessment()) |> relocate(group) |>
  select(group,all_of(d_var$variable)) |>
  mutate(id = 1:n()) |>
  pivot_longer(-c(id,group)) |>
  mutate(name = fct_recode(name, !!!name_swap)) |>
  pivot_wider(names_from = name, values_from = value) |>
  select(-id)


#Do a quick correlation plot where the variables are clumped by correlations, with no grouping 
qgraph(cor_auto(d_ega), layout="spring")

## Calculate EGA =====

# exploratory graph analysis
EGA_0 <- 
  EGA(d_ega |> select(-group), 
      plot.EGA = TRUE,
      model = "glasso",
      algorithm = "walktrap",
      corr = "spearman")

#So this produces 6 dimensions

# 1 - Health/Development
# 2 - Co-ordination
# 3 - Pro-social behaviour
# 4 - Conduct problems
# 5 - Situation anxiety + insomnia
# 6 - anxiety related to being left alone

# Plus weight at birth which doesn't fit in any

#Every time we seem to run this we get somewhat different results, which suggests a rather unstable structure

#So we run the bootstrapped EGA method to test this out
set.seed(08092022)

parallel:::setDefaultClusterOptions(setup_strategy = "parallel")

EGA_boot_0 <- 
  bootEGA(d_ega |> select(-group), 
          iter = 10000,
          type = "resampling",
          seed = 15082022,
          uni.method = "LE",
          corr = "spearman",
          model = "glasso",
          algorithm = "walktrap",
          typicalStructure = TRUE,
          plot.typicalStructure = TRUE,
          plot.type = "GGally",
          ncores = 16)

#This produces a 6 dimension solution - but the graph has 5?
EGA_boot_0$summary.table

EGA_boot_0$boot.ndim |> 
  as_tibble() |> 
  ggplot(aes(N.Dim)) + 
  geom_histogram(binwidth = 1)

#This is the bit of the object that contains the useful data
# EGA_boot_0$typicalGraph


## Plot Graph =======

# Plot the empirical EGA and bootstrap alongside each other
EGA_0 |> plot() |EGA_boot_0 |> plot()

#So some of the names are transposed and the conduct thing gets merged

# #Plot our bootstrapped graph
p_EGA_boot_0 =
  EGA_boot_0 |>
  plot(plot.args = list(vsize = 5,
                        label.size = 4,
                        legend.names = c("1: Conduct",
                                         "2: Health/Development",
                                         "3: Co-ordination",
                                         "4: Situation Anxiety/Sleep",
                                         "5: Separation Anxiety") )) +
  labs(title = "EGA Plot")



## Item Stability =====

#This is confusing, so we will look at dimensional stability 
# (we are drawing on the package tutorial here https://www.mdpi.com/2624-8611/3/3/32/htm)

EGA_boot_0_stab = dimensionStability(EGA_boot_0)

EGA_boot_0_stab$dimension.stability$structural.consistency

#We can see that our small dimensions are very stable, but our dimensions with
#more variables are very unstable (explains why the bootstrapping keeps giving 
#different results)

#Proportion of times each item replicated within the EGA defined dimension
EGA_boot_0_stab$item.stability$item.stability$empirical.dimensions |> 
  as_tibble(rownames = "variable") |>
  arrange(-value) |>
  print(n = 30)

#So you can see the situational anxiety and conduct items are very consistently assigned to
#their dimension over bootstraps, but 10 variables have a proportion < 0.75, suggesting
#they are rather unstable in their dimensional affiliation


#Item stability across all variables
EGA_boot_0_stab$item.stability$item.stability$all.dimensions |>
  as_tibble(rownames = "Variable") |>
  mutate(across(where(is.double),~if_else(.x < 0.15,0,.x))) |>
  print(n = 30)

#So some variables are popping in between dimensions much more than others

#Mean network loadings
EGA_boot_0_stab$item.stability$mean.loadings |>
  as_tibble(rownames = "Variable") |>
  
  #We are going to zero out anything < 0.15 to clarify when printing
  mutate(across(where(is.double),~if_else(.x < 0.15,0,.x))) |>
  print(n = 30)

#So here we can see that no variables are really cross-loading
#strongly on a realistic number of dimensions, but some variables
#are not really loading on any dimensions strongly

#In the EGA tutorial the authors suggest we could improve the consistency of 
#our dimension fit by removing the less stable variables

#Here we repeat the bootstrap EGA, removing the variables with an item stability < 0.75


# Fit final EGA model ======

parallel:::setDefaultClusterOptions(setup_strategy = "parallel")

EGA_boot_1 <- 
  bootEGA(d_ega |>
            select(all_of(EGA_boot_0_stab$item.stability$item.stability$empirical.dimensions |> 
                            as_tibble(rownames = "variable") |>
                            filter(value > 0.75) |>
                            pull(variable))), 
          iter=9999,
          type = "resampling",
          seed = 15082022,
          uni.method = "LE",
          corr = "spearman",
          model = "glasso",
          algorithm = "walktrap",
          typicalStructure = TRUE,
          plot.typicalStructure = TRUE,
          plot.type = "GGally",
          ncores = 16)


#It takes ages to fit this model, so lets save it
write_rds(EGA_boot_1,'final_ega_model.rds')
# EGA_boot_1 = read_rds('final_ega_model.rds')



#Now lets look at the stability again
EGA_boot_1_stab = dimensionStability(EGA_boot_1)

EGA_boot_1_stab$dimension.stability$structural.consistency

#Now we have 5 very consistent dimensions

# 1 - Communication/Play
# 2 - Conduct
# 3 - Co-ordination
# 4 - Situational anxiety and sleep
# 5 - Separation Anxiety

#Item stability across all variables
EGA_boot_1_stab$item.stability$item.stability$all.dimensions |>
  as_tibble(rownames = "Variable") |>
  mutate(across(where(is.double),~if_else(.x < 0.15,0,.x))) |>
  print(n = 30)

#Item loadings
EGA_boot_1_stab$item.stability$mean.loadings |>
  as_tibble(rownames = "Variable") |>
  
  #We are going to zero out anything < 0.15 to clarify when printing
  mutate(across(where(is.double),~if_else(.x < 0.15,0,.x))) |>
  print(n = 30)



## Figure 3 ========


#Make the plot that will become figure 3
p_EGA_boot_1 =
  EGA_boot_1 |> 
  plot(plot.args = list(vsize = 5,
                        label.size = 4,
                        node.color = "black",
                        edge.color = c("grey20","grey20"),
                        legend.names = c("1: Conduct",
                                         "2: Separation Anxiety",
                                         "3: Situational Anxiety/Sleep",
                                         "4: Communication/Play",
                                         "5: Movement/Co-ordination") ))

p_EGA_boot_1
  

ggsave("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Figures/Figure Parts/figure_3_ega.pdf",
       plot = p_EGA_boot_1,width = 6,height = 3)




## Extract dimensions =====

## extract typical domains from bootstrap analysis and write to excel file to create table
# Typical <- EGA_boot_0$typicalGraph$typical.dim.variables

Typical <- EGA_boot_1$typicalGraph$typical.dim.variables


#We need to think of some names for the dimensions

# Dim 1 is Conduct
# Dim 2 is Separation Anxiety
# Dim 3 is Situational Anxiety and sleep
# Dim 4 is Communication and play
# Dim 5 is Co-ordination

#Glue that together with the variable information
var_dims = 
  Typical |>
  as_tibble() |>
  rename(short_name = items) |>
  left_join(variable_table,by = "short_name")

var_dims = 
  var_dims |> 
  mutate(dim_name = case_when(
    dimension == 1 ~ "1: Conduct",
    dimension == 2 ~ "2: Separation Anxiety",
    dimension == 3 ~ "3: Situational Anxiety/Sleep",
    dimension == 4 ~ "4: Communication/Play",
    dimension == 5 ~ "5: Movement/Co-ordination"))
  
var_dims = 
  var_dims |> 
  select(variable,short_name,dim_name,variable_definition_dict) |>
  left_join(d_var |> select(variable,dropout_loss,.lower,.upper),
            by = "variable")  

var_dims|>
  print(n = 30)

#Save as a spreadsheet
write_csv(var_dims, "nested_cv_variable_dimensions.csv")



#The package also contains a function that makes a method section from the object (amazing!)
methods.section(EGA_boot_1)



## Supplementary Table 7 ======

#Tibbulate the full table of variable definitions
var_dims     = read_csv("nested_cv_variable_dimensions.csv")
var_expanded = read_csv("VariableDefinitionsExpanded.csv")
d_var_full   = readxl::read_excel("VariableDefinitions.xlsx", sheet = "Sheet1")

#Make a table for the manuscript
left_join(var_dims |>
            select(variable,short_name,dim_name),
          d_var_full |>
            select(Code, `Paper Description`) |>
            janitor::clean_names() |>
            rename(variable = code),
          by = "variable") |>
  left_join(var_expanded |> 
              select(short_name,var_type),
            by = "short_name") |>
  relocate(variable,short_name,paper_description,var_type,dim_name)|>
  rename(Variable = short_name,
         `Variable Definition` = paper_description,
         `Variable Type` = var_type,
         `Dimension Name` = dim_name) |>
  select(-variable) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)



# Make a variable importance plot, but coloured by dimension
var_dims |>
  ggplot(aes(x = dropout_loss, y = forcats::fct_reorder(short_name,dropout_loss,max),
             xmin = .lower, xmax = .upper, colour = dim_name)) +
  geom_point(size = 4) +
  geom_linerange() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Dropout Loss (1 - AUC)", y = "Variable", 
       title = "Best RF Model") +
  facet_wrap(~dim_name, ncol = 1)

#So there is no real pattern in terms of different dimensions being more predictive -
#maybe the development/movement dimension is more predictive



## CFA on the graph structure =====

# Do CFA on the domains identified by EGA, using an EGA built in function
EGA_CFA <- 
  CFA(EGA_boot_1$EGA, d_ega |> select(-group), estimator = "WLSMV")

#Get a summary of the model
summary(EGA_CFA)

#Interpret the indices of model fit
effectsize::interpret_cfi(EGA_CFA$fit.measures[4] )
effectsize::interpret_rmsea(EGA_CFA$fit.measures[5])


#Make a plot of the factors
EGA_CFA |> 
  plot()



# Factor analysis ======


#This is an alternative to the EGA but we do not include it in the manuscript


## N Factors ======

# We can also do exploratory factor analysis


#We can also use the EFAtools toolbox to ask how many factors we should extract
EFAtools::N_FACTORS(d_ega |> select(-group))

#This produces a range of results between 2 and 11,  One eigenvalue is particularly big, only 2 are > 1

#This is the parallel analysis result showing 5 factors with significant eigenvalues after PCA
EFAtools::PARALLEL(d_ega |> select(-group), eigen_type = "PCA")

## Fit EFA Models =====

#We can do the FA in EFA tools - lets do 2 and 5 and then compare

#Do a 2 factor model - using an oblique rotation (given this is the default of psych and also having correlated factors makes sense)
efa_2 = EFAtools::EFA(d_ega |> select(-group),n_factors = 2,rotation = "oblimin", method = "ML")
efa_2

#Do a 5 factor model, using an oblique rotation 
efa_5 = EFAtools::EFA(d_ega |> select(-group),n_factors = 5,rotation = "oblimin", method = "ML")
efa_5


#Compare fitting methods with a 5 factor solution
EFAtools::EFA_AVERAGE(d_ega |> select(-group),n_factors = 5, N = 500,
                      method = c("PAF","ML","ULS"), rotation = c("oblique"),
                      show_progress = TRUE)


#The plot shows that we are identifying four factors (although some are clearer than others)

# 1 loads heavily on the co-ordination variables
# 2 loads on speech + reading - so communication
# 3 loads on conduct-related variables
# 4 loads on lying/cheating
# 5 loads on separation anxiety


## Tabulate EFA model fit =====

#We can compare the 2, and 5 factor model fits like this:
tab_efa =
  bind_rows(
  efa_2$fit_indices |> as_tibble() |> mutate(Model = "2 Factor"),
  efa_5$fit_indices |> as_tibble() |> mutate(Model = "5 Factor")
) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  mutate(`RMSEA (90% CI)` = paste(round(RMSEA,digits = 3) ," (",round(RMSEA_LB,digits = 3),", ", round(RMSEA_UB,digits = 3), ")",sep = "")) |>
  select(-c(RMSEA_LB:RMSEA_UB)) |>
  relocate(Model, CAF,CFI,`RMSEA (90% CI)`,AIC,BIC,Fm,chi,chi_null,df,df_null,p_null,p_chi) 

tab_efa |>
  pull(CFI) |>
  effectsize::interpret_cfi(rules = "default")

tab_efa |>
  pull(RMSEA) |>
  effectsize::interpret_rmsea(rules = "default")


# So 5 is better than two, based on the CFI and RMSEA


## Tabulate factor loadings =====

#Lets tabulate our 5 factor model loadings
efa_5$rot_loadings[,] |> 
  as_tibble(rownames = "Measure") |>
  rename(motor         = V1,          
         communication = V2,
         considerate   = V3,
         conduct       = V4,
         anxiety       = V5) |>
  left_join(efa_5$h2 |> as_tibble(rownames = "Measure"),by = "Measure") |>
  rename("h2" = value) |>
  mutate(across(where(is.double),round, digits = 2),
         Measure = str_to_title(Measure)) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


## Plot factor loadings ====

p_factor_loading = 
  efa_5$rot_loadings[,] |> 
  as_tibble(rownames = "Measure") |>
  rename(motor         = V1,          
         communication = V2,
         considerate   = V3,
         conduct       = V4,
         anxiety       = V5) |>
  pivot_longer(-Measure,names_to = "factor",values_to = "loadings") |>
  mutate(factor = factor(factor, levels = c("motor","communication","considerate","conduct","anxiety"))) |>
  ggplot(aes(x = loadings,y = Measure)) +
  geom_col(aes(fill = loadings)) +
  geom_vline(xintercept = 0) +
  facet_wrap(~factor,nrow = 1) +
  scale_fill_gradient2(name = "Loading",
                       high = "red", mid = "white", low = "blue",
                       midpoint=0, guide=F)  +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.title = element_text(size  = 8,family = "sans"),
        axis.text  = element_text(size  = 8,family = "sans"),
        strip.text = element_text(size  = 9,family = "sans")) +
  labs(x = "Factor Loading", x = "Measure", title = "5 Factor Loadings") +
  coord_cartesian(xlim = c(-0.9,0.9))

p_factor_loading


p_corr = 
  ggcorrplot(corr.mat, hc.order = FALSE, outline.color = "white") +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.title = element_text(size  = 6,family = "sans"),
        axis.text  = element_text(size  = 6,family = "sans"),
        axis.text.x.bottom = element_text(size  = 6,family = "sans"),
        axis.text.y.left =   element_text(size  = 6,family = "sans",angle = 45),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6)) 

p_loading_stack = 
  efa_5$rot_loadings[,] |> 
  as_tibble(rownames = "Measure") |>
  rename(motor         = V1,          
         communication = V2,
         considerate   = V3,
         conduct       = V4,
         anxiety       = V5) |>
  pivot_longer(-Measure,names_to = "factor",values_to = "loadings") |>
  mutate(factor = factor(factor, levels = c("motor","communication","considerate","conduct","anxiety"))) |>
  mutate(Measure = factor(Measure,levels = corr.mat |>
                            as_tibble(rownames = "variable") |> 
                            pull(variable) |> rev())) |>
  ggplot(aes(x = loadings,y = fct_rev(Measure))) +
  geom_col(aes(fill = factor),position = "stack") +
  scale_fill_brewer(palette = "Dark2") +
  geom_vline(xintercept = 0,size = 0.25, lty = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        axis.title = element_text(size  = 6,family = "sans"),
        axis.text  = element_text(size  = 6,family = "sans"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6)) +
  annotate(
    geom = 'segment',
    y = 0,
    yend = 0,
    x = -Inf,
    xend = Inf)


p_corr_fa = (p_corr | p_loading_stack) + plot_layout(ncol = 2,widths = c(4,1))

# ggsave("./Figures/figure_4_corr_fa.pdf",plot = p_corr_fa,width = 10,height = 6)




# Prepare Minimal Variable Set =====

#What if we took the variable from each factor with the highest importance and did a ML model with just those variables?
#Then we would have a 5 item screener, which would be very handy and have some semi theoretical basis

#Lets do this = 
vars_EGA = 
  var_dims |>
  group_by(dim_name) |>
  slice_max(dropout_loss, n = 1) 

#Now we can do away and do ML on these - but we can keep it simple and just use the linear SVM Model
form_list = 
  vars_EGA|> 
  mutate(type = "EGA")|>
  group_by(type) |>
  nest() |>
  
  #Make our formulas
  mutate(formulas = map_chr(data,~paste("group ~",.x |> 
                                          pull(variable) |> 
                                          paste(collapse = " + ")))) |>
  
  #Make recipes which we will append our formulas into
  mutate(recipes = map(formulas, ~recipe(formula = as.formula(.x),
                                         data    = d_last_fit |> analysis()) |>
                         step_zv(all_predictors())))

#Extract the recipes and convert them into a named list
final_recipes = 
  form_list |>
  pull(recipes) 

names(final_recipes) <-  
  form_list |> pull(type)


#Save these variable sets
write_rds(final_recipes,"final_selected_vars.rds")


## Tabulate the variable descriptions ========

final_recipes = read_rds("final_selected_vars.rds")

#We can look these items up in the data dictionary
vars_final  = 
  final_recipes$EGA$var_info |> 
  filter(variable != "group") |>
  mutate("var_set" = "EGA")

d_var = read_csv("VariableDefinitionsExpanded.csv")

d_var |> 
  filter(variable %in% vars_final$variable) |>
  pull(paper_description)

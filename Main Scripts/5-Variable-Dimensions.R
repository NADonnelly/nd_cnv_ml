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


#Load the imputed dataset
d_last_fit = read_rds("nested_cv_imputed_data.rds")


## Prepare variable defintions =====

#Load our final variables
d_var = read_rds("nested_cv_selected_var_definitions.rds")


#We might benefit from make short names for each variable.

#What do we got? We have 30 variables in total:

#CAPA 
# 1  pcb1i01         Anxious affect - ?Fear of activities in public                                              = FPB                                                                      
# 2  pcb3i01         Anxious affect - ?Agoraphobia                                                               = AGO           
# 3  pce0i01         Anxious affect - ?Fear of blood/injections                                                  = FBI                          
# 4  pcd2i01         Rumination, Obsessions and Compulsions - ?Rumination Intensity                              = RUM                       
# 5  pda0i02         Depression  - ? Episode of Depressed mood intensity                                         = DeI                                       
# 6  pda0i03         Depression  - ? Period of 2 continuous months without depressed mood in last year           = D2M                                                                         
# 7  pda1i01         Depression  - ? Distinct quality of depressed mood                                          = DeQ                                       
# 8  pbe1i01         Physical symptoms on separation from caregiver (separation anxiety symptom)                 = SAP                                                                           
# 9  pbe7i01         Separation Anxiety - ?Intensity of separation worries/anxiety (across multiple activities?) = SAI                                                                       
# 10 pbf8i01         Separation Anxiety - ?Frequency of separation anxiety                                       = SAF                                
# 11 pbf2i01         Separation Anxiety - ?Avoidance of sleeping away from family                                = SAS                                       
# 12 pfb7i01         Sleep probs - ?Total Insomnia intensity                                                     = InT                         
# 13 pfb7i02         Sleep probs - ?Initial insomnia intensity                                                   = InI                            
# 14 prc8i01         hyperactivity - ?Forgetful in daily activities intensity                                    = FGT                                               
# 15 prb8i01         Often blurts out answers to questions (ADHD Impulsivity/Hyperactivity symptom)              = BLT                                                                             
# 16 pgc3i01         oppositional / conduct disorder  - Lying intensity                                          = LIE              
# 17 pgc5i01         oppositional / conduct disorder  - Cheating intensity                                       = CHT                   

#Your child's health and development 
# 18 P_Health_dev_9  Health and Development - Is your child behind in reading                                    = REA                  
# 19 P_Health_dev_11 Health and Development - educationally statemented                                          = EST                                  
# 20 P_Health_dev_13 Health and Development - did your child walk by 18 months                                   = W18                            
# 21 P_Health_dev_15 Health and Development - has your chld had speech therapy                                   = SLT                        
# 22 P_Health_dev_21 Health and Development - other problems with airways/lungs                                  = LUN               
# 23 P_Health_dev_24 Health and Develpoment - heart problems                                                     = CAR                                     
# 24 P_Health_dev_27 Health and Development - skeletal or muscular problems                                      = MSK                      

#SCQ
# 25 P_ASQ_7         ASQ Behav and social comm - invented words, odd indirect, metaphorical ways                 = AWM      
# 26 P_ASQ_8         ASQ Behav and social comm - say the same thing over and over                                = ARP

#DCDQ
# 27 CMD_2           Coordination and motor development - catches a small ball thrown from 6-8ft                 = CBA      
# 28 CMD_5           Coordination and motor development - runs as fast and easily as other children              = CRF 
# 29 CMD_6           Coordination and motor development - can organise her body to do a planned motor activity   = COB
# 30 CMD_10          Coordination and motor development - cuts pictures and shapes accurately                    = CCP   
 

#We can try and get variable definitions from the spreadsheet found on the 
#imagine-id website?

#We can look these items up in the data dictionary....
vd <- readxl::read_excel("Copy-of-F2f-data-dictionnary-assessment-questions-removed.xlsx", 
                         sheet = "Sheet1")

#We have made a nice sheet of these definitions
vd2 <- readxl::read_excel("Variable_definitions.xlsx", 
                          sheet = "Sheet1")


vd = 
  vd |> 
  select(`Question number`,`Question description`,Assessment,Section,`Question description (ME)`) |>
  drop_na(`Question description`) |>
  mutate(`Question number` = map_chr(`Question number`,str_remove,"_L")) |>
  filter(`Question number` %in% d_var$VARIABLE) |>
  rename(VARIABLE = `Question number`) |>
  right_join(d_var,by = c("VARIABLE")) |>
  mutate(var_def = if_else(is.na(`Question description (ME)`),`VARIABLE DEFINITION`,`Question description (ME)`)) |>
  select(VARIABLE,Assessment,Section,var_def) 

#Lets smush these together
d_var = 
  d_var |>
  rename(variable  = VARIABLE,
         data_type = `DATA TYPE`,
         var_def   = `VARIABLE DEFINITION`) |>
  left_join(vd |> 
              select(VARIABLE, var_def,Assessment,Section) |> 
              rename(variable = VARIABLE,var_def_long = var_def), 
            by = c("variable")) |>
  left_join(vd2 |> rename(variable = Code), by = c("variable")) |>
  
  #This is my rather arbitrary variable naming
  mutate(short_name = case_when(variable == "pcb1i01" ~ "FPB",
                                variable == "pcb3i01" ~ "AGO",
                                variable == "pce0i01" ~ "FBI",
                                variable == "pcd2i01" ~ "RUM",
                                variable == "pda0i02" ~ "DeI",
                                variable == "pda0i03" ~ "D2M",
                                variable == "pda1i01" ~ "DeQ",
                                variable == "pbe1i01" ~ "SAP",
                                variable == "pbe7i01" ~ "SAI",
                                variable == "pbf8i01" ~ "SAF",
                                variable == "pbf2i01" ~ "SAS",                              
                                variable == "pfb7i01" ~ "InT",
                                variable == "pfb7i02" ~ "InI",
                                variable == "prc8i01" ~ "FGT",
                                variable == "prb8i01" ~ "BLT",
                                variable == "pgc3i01" ~ "LIE",
                                variable == "pgc5i01" ~ "CHT",
                                variable == "P_Health_dev_9" ~ "REA",
                                variable == "P_Health_dev_11" ~ "EST",
                                variable == "P_Health_dev_13" ~ "W18",
                                variable == "P_Health_dev_15" ~ "SLT",
                                variable == "P_Health_dev_21" ~ "LUN",
                                variable == "P_Health_dev_24" ~ "CAR",
                                variable == "P_Health_dev_27" ~ "MSK",
                                variable == "P_ASQ_7" ~ "AWM",
                                variable == "P_ASQ_8" ~ "ARP",
                                variable == "CMD_2" ~ "CBA",
                                variable == "CMD_5" ~ "CRF",
                                variable == "CMD_6" ~ "COB",
                                variable == "CMD_10" ~ "CCP")) 






#Tabulate our variable definitions
d_var  |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)




d_var|>
  select(variable, Assessment, Section, short_name,var_def_long,`Paper Description`) |>
  rename(short_def = var_def_long,
         long_def   = `Paper Description`) |>
  janitor::clean_names() |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


#Save this
write_rds(d_var |>
            janitor::clean_names(), 
          "nested_cv_selected_var_definitions_expanded.rds")




# Correlation between variables ========

# Lets look at correlation between variables

#Do we use the training set only or both the training and test sets?
DF = bind_rows(d_last_fit |> analysis(),
               d_last_fit |> assessment())

DF = DF |> relocate(group)

DF2 = 
  DF |>
  select(all_of(d_var$variable)) 

# - we can use the polychoric correlation method as in the original version of the paper


## Do Polychor Correlations ======

#The code I inherited does a lot of fancy stuff to make a polychoric correlation matrix
#with our variables, and then does PCA with that
#
# But we can do polychoric correlation directly within the pca function from psych, so why
#wouldn't we do that?
#
#Well it turns out that doesn't work, annoyingly. I think from experimentation that this is
#directly because there are too many variables?

#Or is there a problem with one variable? pbd7e01 is coded as either 0 or 2, but its dichotomous,
# not continuous/polytomous. I think this is causing problems. But there are other issues

#e.g. 
# 
# mixedCor(DF[,53:77] |>
#             mutate(pbd7e01 = recode(pbd7e01, `2` = 1))|> 
#            as.matrix(),
#          global = F,smooth = T)
#
#Which I think is due to zero value cells when you crosstable variables.
#
#So instead we are going to use a loop via map to do each correlation via


## Set up dataframe for polychoric correlation calculation
D3 = 
  DF |>
  select(-group)


## Set up a matrix of combinations to iterate through
Is   <- expand_grid(x = 1:ncol(D3),y = 1:ncol(D3))
cors <- matrix(nrow = ncol(D3), ncol = ncol(D3))


## Set up parallel processing cluster

numCores <- 16
cl       <- makeSOCKcluster(numCores)

registerDoSNOW(cl)

parallel::clusterExport(cl, "Is")
parallel::clusterExport(cl, "D3")

parallel::clusterEvalQ(cl, library(polycor))
parallel::clusterEvalQ(cl, library(dplyr))

## set up progress bar
pb       <- txtProgressBar(max = nrow(Is), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts     <- list(progress = progress)

## calculate polychoric correlation between each variable
r <-
  
  foreach(i = 1:nrow(Is), .options.snow = opts, .combine = "c") %dopar% {
    
    vs <- (Is[i,])
    v1 <- pull(vs[1])
    v2 <- pull(vs[2])
    
    cors[v1,v2] <- polychor(D3[,v1] |> pull(),
                            D3[,v2] |> pull(), 
                            ML = FALSE)
    
  }

conflicted::conflict_prefer("stopCluster", "parallel")

stopCluster(cl)

#Recover our correlation matrix from the parallel thing
r2 <- 
  matrix(r, 
         nrow     = sqrt(length(r)), 
         ncol     = sqrt(length(r)), 
         dimnames = list(names(D3), names(D3)))


#Not sure what this thing about positive definite matrixes  - does it turn the correlations in covariances?
r2 <- nearPD(r2, corr = T)
r2 <- as.matrix(r2$mat)

#Save up
r2.df <- as.data.frame(r2)

write_csv(r2.df, "nested_cv_corr_mat.csv")

r2.df = read_csv( "nested_cv_corr_mat.csv")


## Extract our variables of interest =====

#We did the correlation on all variables
corr = 
  r2.df |>
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

row.names(corr.mat) = d_var$short_name
colnames(corr.mat) = d_var$short_name

ggcorrplot(corr.mat, hc.order = TRUE, outline.color = "white")



# Exploratory Graph Analysis ======


#As in the original draft of the paper, we can do EGA

#Can we switch out our variable names for our short names
name_swap        = d_var$variable
names(name_swap) = d_var$short_name

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
EGA_0 <- EGA(d_ega |> select(-group), 
             plot.EGA = TRUE,
             model = "glasso",
             algorithm = "walktrap",
             corr = "spearman")

#So this produces 4 dimensions

# 1 - Anxiety/Conduct
# 2 - Health/Development
# 3 - Insomnia
# 4 - Depression

#Every time we seem to run this we get somewhat different results, which suggests a rather unstable structure
set.seed(08092022)

#Now we run the bootstrapped EGA method
parallel:::setDefaultClusterOptions(setup_strategy = "parallel")

EGA_boot_0 <- bootEGA(d_ega |> select(-group), 
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

#This actually produces a 4 dimension solution - but the modal network size is 5. This is all very odd - why is it doing that
EGA_boot_0$summary.table
EGA_boot_0$boot.ndim |> as_tibble() |> ggplot(aes(N.Dim)) + geom_histogram(binwidth = 1)

#This is the bit of the object that contains the useful data
# EGA_boot_0$typicalGraph


## Plot Graph =======

# Plot the empirical EGA and bootstrap alongside each other
EGA_0 |> plot() |EGA_boot_0 |> plot()

#So dimensions 1 and 2 are transposed

# #Plot our bootstrapped graph
p_EGA_boot_0 =
  EGA_boot_0 |>
  plot(plot.args = list(vsize = 5,
                        label.size = 4,
                        legend.names = c("1: Development/Health",
                                         "2: Anxiety/Hyperactivity",
                                         "3: Sleep",
                                         "4: Depression") )) +
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

#So you can see the insomnia and depression items are very consistently assigned to
#their dimension over bootstraps, but 1- variables have a proportion < 0.75, suggesting
#they are rather unstable in their dimensional affiliation


#Item stability across all variables
EGA_boot_0_stab$item.stability$item.stability$all.dimensions |>
  as_tibble(rownames = "Variable") |>
  mutate(across(where(is.double),~if_else(.x < 0.15,0,.x))) |>
  print(n = 30)

#So some variables are popping in between dimensions 1 and 5 or 2 and 5, or 1 and 2 and 5

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

parallel:::setDefaultClusterOptions(setup_strategy = "parallel")

EGA_boot_1 <- bootEGA(d_ega |> select(all_of(EGA_boot_0_stab$item.stability$item.stability$empirical.dimensions |> 
                                               as_tibble(rownames = "variable") |>
                                               filter(value > 0.75) |>
                                               pull(variable))), 
                      iter=10000,
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

#Now lets look at the stability again
EGA_boot_1_stab = dimensionStability(EGA_boot_1)

EGA_boot_1_stab$dimension.stability$structural.consistency

#Now we have 4 very consistent dimensions

# 1 - Anxiety (separation anxiety, agoraphobia)
# 2 - Co-ordination and development
# 3 - Insomnia
# 4 - Depression

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



## Figure 4 ========


#Make the plot that will become figure 4
p_EGA_boot_1 =
  EGA_boot_1 |>
  plot(plot.args = list(vsize = 5,
                        label.size = 4,
                        node.color = "black",
                        edge.color = c("grey20","grey20"),
                        legend.names = c("1: Anxiety",
                                         "2: Movement/Development",
                                         "3: Insomnia",
                                         "4: Depression") ))

p_EGA_boot_1
  

ggsave("./Figures/figure_4_ega.pdf",plot = p_EGA_boot_1,width = 6,height = 3)



#The package also contains a function that makes a method section from the object (amazing!)
methods.section(EGA_boot_1)


## Extract dimensions =====

## extract typical domains from bootstrap analysis and write to excel file to create table
# Typical <- EGA_boot_0$typicalGraph$typical.dim.variables

Typical <- EGA_boot_1$typicalGraph$typical.dim.variables


#We need to think of some names for the dimensions

# Dim 1 is Anxiety
# Dim 2 is Development and Co-ordination
# Dim 3 is Insomnia
# Dim 4 is Depression

#Glue that together with the variable information
var_dims = 
  Typical |>
  as_tibble() |>
  rename(short_name = items) |>
  left_join(d_var,by = "short_name")

var_dims = 
  var_dims |> 
  mutate(dim_name = case_when(
    dimension == 1 ~ "Anxiety",
    dimension == 2 ~ "Development/Co-ordination",
    dimension == 3 ~ "Insomnia",
    dimension == 4 ~ "Depression"))
  
var_dims |> 
  select(variable,short_name,var_def,dim_name,dropout_loss,.lower,.upper) |>
  print(n = 30)

#Save as a spreadsheet

write_csv(var_dims |>
            select(variable,short_name,var_def,dim_name,dropout_loss,.lower,.upper), 
          "nested_cv_variable_dimensions.csv")


## Table 6 ======

var_dims = read_csv("nested_cv_variable_dimensions.csv")
d_var    = read_rds("nested_cv_selected_var_definitions_expanded.rds")

#Make a table as in the manuscript

left_join(var_dims |>
            select(variable,dim_name),
          d_var|>
            select(variable, short_name,var_def_long,paper_description) |>
            rename(short_def = var_def_long,
                   long_def  = paper_description) ,
          by = "variable") |>
  relocate(variable,short_name,short_def,long_def,dim_name)|>
  select(-short_def) |>
  rename(Variable = variable,
         `Variable Name` = short_name,
         `Variable Definition` = long_def,
         `Dimension Name` = dim_name) |>
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
       title = "Best SVM Model") +
  facet_wrap(~dim_name, ncol = 1)

#So there is no real pattern in terms of different dimensions being more predictive -
#maybe the development/movement dimension is more predictive




## CFA on the graph structure =====


# Do CFA on the domains identified by EGA, using an EGA built in function

#This is odd because it is getting different numbers of dimensions

EGA_CFA <- CFA(EGA_boot_1$EGA, d_ega |> select(-group), estimator = "ML")

effectsize::interpret(EGA_CFA$fit)

summary(EGA_CFA)

p_cfa = 
  EGA_CFA |> plot()




# UVA =====

# The EGA package also has a method that tries to identify redundant variables

uva_0 = UVA(d_ega |> select(-group))

#Apparently we have no redundant variables, which is nice


# Factor analysis ======

## N Factors ======

# We can also do exploratory factor analysis

#We can use the psych toolbox for this
nf_psych <-psych::fa.parallel(d_ega |> select(-group), cor="poly", fm="ml")

#Psych suggests 6 components / 9 factors

#We can also use the EFAtools toolbox to ask how many factors we should extract
EFAtools::N_FACTORS(d_ega |> select(-group))

#This produces a range of results between 2 and 12,  One eigenvalue is particularly big, only 2 are > 1

#This is the parallel analysis result showing 6 factors with significant eigenvalues after PCA
# - although two are very close to the line

EFAtools::PARALLEL(d_ega |> select(-group), eigen_type = "PCA")
EFAtools::PARALLEL(d_ega |> select(-group), eigen_type = "SMC")


#We can do the FA in EFA tools - lets do 2, 4 and 6 and then compare

#Do a 2 factor model - using an oblique rotation (given this is the default of psych and also having correlated factors makes sense)
efa_2 = EFAtools::EFA(d_ega |> select(-group),n_factors = 2,rotation = "oblimin", method = "ML")
efa_2

#Do a 4 factor model, using an oblique rotation 
efa_4 = EFAtools::EFA(d_ega |> select(-group),n_factors = 4,rotation = "oblimin", method = "ML")
efa_4


#Do a 6 factor model, using an oblique rotation 
efa_6 = EFAtools::EFA(d_ega |> select(-group),n_factors = 6,rotation = "oblimin", method = "ML")
efa_6

#Use PAF - basically similar results
efa_6_paf = EFAtools::EFA(d_ega |> select(-group),n_factors = 6,rotation = "oblimin", method = "PAF")
efa_6_paf

#Do a 
efa_6_bi = EFAtools::EFA(d_ega |> select(-group),n_factors = 6,rotation = "bifactorQ", method = "ML")

#Compare rotations
EFAtools::COMPARE(
  EFAtools::EFA(d_ega |> select(-group),n_factors = 6,rotation = "oblimin", method = "ML")$rot_loadings,
  # EFAtools::EFA(d_ega |> select(-group),n_factors = 6,rotation = "oblimin", method = "PAF")$rot_loadings,
  EFAtools::EFA(d_ega |> select(-group),n_factors = 6,rotation = "bifactorQ", method = "ML")$rot_loadings,
  x_labels = c("ML and oblimin","ML and bifactorQ")
)

#Differences are only small


#Compare fitting methods with a 6 factor solution
EFAtools::EFA_AVERAGE(d_ega |> select(-group),n_factors = 6, N = 500,
                      method = c("PAF","ML","ULS"), rotation = c("oblique"),
                      show_progress = TRUE)

#The plot shows that we are identifying four factors (although some are clearer than others)

# 1 loads heavily on the co-ordination/development variables
# 2 loads on a set of CAPA variables for depression
# 3 loads on CAPA sleep
# 4 loads on separation anxiety
# 5 loads on conduction/hyperactivity
# 6 loads on ASQ/speech

### Tabulate EFA model fit =====

#We can compare the 2, 4 factor and 6 factor model fits like this:
tab_efa =
  bind_rows(
  efa_2$fit_indices |> as_tibble() |> mutate(Model = "2 Factor"),
  efa_4$fit_indices |> as_tibble() |> mutate(Model = "4 Factor"),
  efa_6$fit_indices |> as_tibble() |> mutate(Model = "6 Factor")
) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  mutate(`RMSEA (90% CI)` = paste(round(RMSEA,digits = 3) ," (",round(RMSEA_LB,digits = 3),", ", round(RMSEA_UB,digits = 3), ")",sep = "")) |>
  select(-c(RMSEA:RMSEA_UB)) |>
  relocate(Model, CAF,CFI,`RMSEA (90% CI)`,AIC,BIC,Fm,chi,chi_null,df,df_null,p_null,p_chi) 

tab_efa |>
  pull(CFI) |>
  effectsize::interpret_cfi(rules = "default")

tab_efa |>
  pull(RMSEA) |>
  effectsize::interpret_rmsea(rules = "default")


# So 6 is probably best, based on the CFI and RMSEA


### Tabulate factor loadings =====

#Lets tabulate our 6 factor model loadings
efa_6$rot_loadings[,] |> 
  as_tibble(rownames = "Measure") |>
  rename(motor      = V1,          
         depression = V2,
         insomnia   = V3,
         anxiety    = V4,
         language   = V5,
         conduct    = V6) |>
  left_join(efa_6$h2 |> as_tibble(rownames = "Measure"),by = "Measure") |>
  rename("h2" = value) |>
  mutate(across(where(is.double),round, digits = 2),
         Measure = str_to_title(Measure)) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


### Plot factor loadings ====

p_factor_loading = 
  efa_6$rot_loadings[,] |> 
  as_tibble(rownames = "Measure") |>
  rename(motor      = V1,          
         depression = V2,
         insomnia   = V3,
         anxiety    = V4,
         language   = V5,
         conduct    = V6) |>
  pivot_longer(-Measure,names_to = "factor",values_to = "loadings") |>
  mutate(factor = factor(factor, levels = c("motor","depression","insomnia","anxiety","language","conduct"))) |>
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
  labs(x = "Factor Loading", x = "Measure", title = "6 Factor Loadings") +
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
  efa_6$rot_loadings[,] |> 
  as_tibble(rownames = "Measure") |>
  rename(motor      = V1,          
         depression = V2,
         insomnia   = V3,
         anxiety    = V4,
         language   = V5,
         conduct    = V6) |>
  pivot_longer(-Measure,names_to = "factor",values_to = "loadings") |>
  mutate(factor = factor(factor, levels = c("motor","depression","insomnia","anxiety","language","conduct"))) |>
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

ggsave("./Figures/figure_4_corr_fa.pdf",plot = p_corr_fa,width = 10,height = 6)




# What do we do with this? =====

#What if we took the variable from each factor with the highest importance and did a ML model with just those variables?
#Then we would have a 6 item screener, which would be very handy and have some semi theoretical basis

#Lets do this = 

#Or ever the four items from the EGA

vars_EGA = 
  var_dims |>
  group_by(dim_name) |>
  slice_max(dropout_loss, n = 1) 


#To do this with the factor loads we need to select the factor each variable loads most on I guess
vars_FA = 
  efa_6$rot_loadings[,] |> 
  as_tibble(rownames = "Measure") |>
  rename(motor      = V1,          
         depression = V2,
         sleep      = V3,
         sep        = V4,
         conduct    = V6,
         speech     = V5)  |>
  pivot_longer(-Measure,names_to = "factor",values_to = "loadings") |>
   mutate(loadings = abs(loadings)) |>
  group_by(Measure) |>
  slice_max(loadings,n = 1) |>
  rename(VARIABLE = Measure) |>
  left_join(d_var,by = "VARIABLE") |>
  ungroup() |>
  group_by(factor) |>
  slice_max(dropout_loss, n = 1) 

#Now we can do away and do ML on these - but we can keep it simple and just use the linear SVM Model
form_list = 
  bind_rows(vars_FA |> 
              mutate(type = "FA") ,
            vars_EGA|> 
              mutate(type = "EGA")) |>
  group_by(type) |>
  nest() |>
  
  #Make our formulas
  mutate(formulas = map_chr(data,~paste("group ~",.x |> 
                                          pull(VARIABLE) |> 
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
write_rds(final_recipes,"ega_fa_selected_vars.rds")


#Tabulate the variable descriptions
final_recipes = read_rds("ega_fa_selected_vars.rds")


vars_final  = 
  bind_rows(
    final_recipes$FA$var_info |> 
      filter(variable != "group") |>
      mutate("var_set" = "FA"),
    final_recipes$EGA$var_info |> 
      filter(variable != "group") |>
      mutate("var_set" = "EGA"))


#We can look these items up in the data dictionary....
#Read in Master participant list for confirmed genotypes
VL <- readxl::read_excel("DATA ENTRY CODEBOOK .xlsx", 
                         sheet = "VARIABLE LIST")

#Wrangle lightly
VL = 
  VL |>
  select(VARIABLE:`VARIABLE DEFINITION`) |>
  filter(VARIABLE %in% best_vars)

vars_final = 
  vars_final |>
  mutate(var_def = map_chr(variable,~ VL |>
                             select(VARIABLE:`VARIABLE DEFINITION`) |>
                             filter(VARIABLE %in% .x) |>
                             pull(`VARIABLE DEFINITION`)))





# An extra something: EGA on all variables (yes all) =====



#Lets load up
pacman::p_load(tidyverse,tidymodels,
               naniar,rgl,polycor,Matrix,
               parallel,doSNOW,tcltk,
               ggcorrplot,
               qgraph, EGAnet, EFATools  )

library(EFAtools)
tidymodels_prefer()


#Load the imputed dataset
d_last_fit = read_rds("nested_cv_imputed_data.rds")



d_ega = 
  bind_rows(d_last_fit |> analysis(),
            d_last_fit |> assessment()) |> relocate(group) |>
  mutate(id = 1:n()) |>
  select(-id)


#Do a quick correlation plot where the variables are clumped by correlations, with no grouping 

qgraph(cor_auto(d_ega), layout="spring")




# exploratory graph analysis
EGA_all <- EGA(d_ega |> select(-group), plot.EGA = TRUE)

#It just looks like total chaos


#Now we run the bootstrapped EGA method
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

EGA_boot_all <- bootEGA(d_ega |> select(-group), iter=20000)


## Plot Graph =======

#Plot our bootstrapped graph
EGA_boot_all |> 
  plot(plot.args = list(legend.names = c("1: social/communication",
                                         "2: hyperactivity/conduct",
                                         "3: speech/education",
                                         "4: motor/coordination") )) +
  labs(title = "EGA Plot")

#This is the bit of the object that contains the useful data
# EGA_boot_all$typicalGraph

#The package also contains a function that makes a method section from the object (amazing!)
methods.section(EGA_boot_all)


## Extract dimensions =====

## extract typical domains from bootstrap analysis and write to excel file to create table
Typical <- EGA_boot_all$typicalGraph$typical.dim.variables


#We need to think of some names for the dimensions

# Dim 1 is social/communication
# Dim 2 is hyperactivity/conduct
# Dim 3 is speech/education
# Dim 4 is motor/coordination

#Glue that together with the variable information
var_dims = 
  Typical |>
  as_tibble() |>
  rename(`SHORT NAME` = items) |>
  left_join(d_var,by = "SHORT NAME")


var_dims = 
  var_dims |> 
  mutate(`DIM NAME` = case_when(
    dimension == 1 ~ "social/communication",
    dimension == 2 ~ "hyperactivity/conduct",
    dimension == 3 ~ "speech/education",
    dimension == 4 ~ "motor/coordination"))

var_dims |> 
  print(n = 29)

#Save as a spreadsheet

write_csv(var_dims |>
            select(VARIABLE,`SHORT NAME`,`VARIABLE DEFINITION`,`DATA TYPE`,dimension,`DIM NAME`,
                   dropout_loss,.lower,.upper), 
          "nested_cv_variable_dimensions.csv")


#Make a table as in the manuscript

var_dims |>
  select(`SHORT NAME`,`VARIABLE DEFINITION`,dimension,`DIM NAME`) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

# Make a variable importance plot, but coloured by dimension

var_dims |>
  ggplot(aes(x = dropout_loss, y = forcats::fct_reorder(`SHORT NAME`,dropout_loss,max),
             xmin = .lower, xmax = .upper, colour = `DIM NAME`)) +
  geom_point(size = 4) +
  geom_linerange() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Dropout Loss (1 - AUC)", y = "Variable", 
       title = "Best SVM Model") +
  facet_wrap(~`DIM NAME`, ncol = 1)

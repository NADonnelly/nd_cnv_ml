

# Introduction =====

#Sample size calculation

#We would like to know whethere we have power to generate an accurate prediction model

#Riley et al https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/sim.7992 produced an R tool 
#for doing power calculation in this situation, although it appears ot be mostly around lgoistc regression models

#The model is based on a number of different parameters:

# - The number of participants
# - The number of outcome events (in our case, the prevalence)
# - The number of predictors

# Prevalence in our sample is a bit peculiar because our cohort is deliberate over-representative
#of ND-CNV carriers compared to the general population (and I would expect also of the CAMHS/
#community paeds population). Therefore our model may not be weell calibrated for use in the general 
#population

# - Prevalence in the general population is ~ 0.5% https://pubmed.ncbi.nlm.nih.gov/32778765/ 
# - Prevalence in a review is 15~20% https://www.sciencedirect.com/science/article/pii/S0002929710002089
# - Prevalence in a large sample with DD from the US/Canada is ~ 14.2% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3171215/
# - Prevalence in Chinese patients with neurodevelopmental disorders is ~ 21% https://pubmed.ncbi.nlm.nih.gov/33402738/

# Regarding the number of predictors, we have a balance to strike - we have lots of potential
#predictor variables, and we would like to reduce this number to those most usefully predictive
#of our outcome, but we also would like to retain enough to gain some insight into underlying
#dimensionality without our big set of predictors

#Regarding an estimate of the C statistic (AUC) to use in the model, I'm not sure what to use. We want
#the ability of behavioural data to predict CNV status. There is data on using genotype to predict e.g. psychiatric
#disorders https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7610853/ - which is very variable e.g AUC range  .48-0.95 

#From some work looking in pubmed it does not appear there is a substantial literature on predicting who is most
#likely to have a cnv using behavioural/clinical psychiatric data. We we will have to make a guess

library(pmsampsize)



#Lets explore this  - how does the minimum sample size change with the number of parameters and c-statistic
pow_data = crossing(par_range = seq(from = 2, to = 50, by = 2),
                    cstat = seq(from = 0.55, to = 0.95, by = 0.05))

pow_data = 
  pow_data |>
  mutate(pow_data = map2(par_range,cstat,~pmsampsize(type = "b",shrinkage = 0.9,
                                                     #We are going to assume the prevalence is the clinical 
                                                     #population prevalence level
                                                     prevalence = 0.15,
                                                     parameters = .x,cstatistic = .y)))

#Lets save this given it takes so long to make
write_rds(pow_data,"./power_calculation_simulations.rds")

pow_data = read_rds("./power_calculation_simulations.rds")


#We can make plots of the number of participants needed to have a well-powered study, based
#on our ranges of potential parameters and model c-statistics
p_pow = 
  pow_data |>
  mutate(sample_size = map_dbl(pow_data, ~.x$sample_size),
         EPP         = map_dbl(pow_data, ~.x$EPP)) |>
  filter(par_range < 30) |>
  ggplot(aes(x = cstat, y = sample_size, colour = factor(par_range))) +
  geom_line(aes(lty  = factor(par_range)),size = 0.5) +
  theme_bw() +
  theme(panel.background = element_blank(),
        # panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  annotate(
    "text", label = "study n = 489",
    x = 0.55, y = 530, size = 4, colour = "black"
  ) +
  labs(x = "Model AUC", y = "Number of Samples Needed", title = "Samples by AUC and n parameters",
       colour = "Number of Parameters",lty = "Number of Parameters") +
  coord_cartesian(ylim = c (100,1000),xlim = c(0.5,0.95)) +
  scale_y_continuous(breaks = seq(100, 1000, by = 100)) +
  geom_hline(yintercept = 489, size = 0.5, colour = "grey10")


#From this plot we can see that even with an AUC of > 0.9, we can only generate a model with 16-18 
#predictors given our sample size

#We can save this as a plot, perhaps for the supplement?
ggsave(plot = p_pow,filename = "./Figures/Figure_Supplement_Power_Calculation.pdf",width = 8,height = 4)

#Our EPP close to this point is ~ 5
pow_data |>
  mutate(sample_size = map_dbl(pow_data, ~.x$sample_size),
         EPP         = map_dbl(pow_data, ~.x$EPP)) |>
  filter(par_range == 16)


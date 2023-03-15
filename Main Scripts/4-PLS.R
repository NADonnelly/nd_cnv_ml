
# Introduction =======

#We are going to do partial least squares on the full imputed dataset

#Load packages
pacman::p_load(tidyverse,tidymodels,patchwork, mixOmics)

#Set package preferences
tidymodels_prefer()


#Set our data directory
if(str_detect(Sys.info()[['nodename']],'IT088825')){
  data_dir = 'C://Users/pyxnd/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data'
  
}else if(str_detect(Sys.info()[['nodename']],'AVOCADO')){
  
  data_dir = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/"
}

setwd(data_dir)


#Load the imputed data - note that this dataset is created in script 3-Data-Split.R
d_last_fit = read_rds("nested_cv_imputed_data.rds")


#Select only the training dataset
d = d_last_fit |>  analysis()

#Clean up
rm(list = c("d_last_fit"))


# Set up =====

#Set a seed for the random number generator for reproducibility
set.seed(05082022)

#The PLS code likes the exposure and outcome variables in separate arrays, so we extract
#from our main data tibble
X <-  d |> select(-group) |> as.matrix()
Y <-  d |> select(group)  |> pull()

# check the dimensions of the X dataframe (participants X variables)
# dim(X)  
# length(Y)
# summary(Y)

# Initial PCA =====

# PCA on the independent variables - asking for 10 components
pca.cnv = pca(X, ncomp = 10, center = TRUE, scale = TRUE) 

# base R barplot of the eigenvalues (explained variance per component)
plot(pca.cnv) 

#So it looks a lot like there is one big component and maybe another 2 which explain
#most of the variance in the data

#Plot scatter plots of the PCA results for the first 2 components, coloured by genotype,
#using the mixedomics package
plotIndiv(pca.cnv, group = Y, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on CNV, comp 1 - 2') # onto the PCA subspace


#So we see a tighter cluster of controls and a wider distribution of CNV carriers


# Basic Sparse PLSDA =====

#Now we do (sparse) PLS - DA
cnv.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later


# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(cnv.splsda , comp = 1:2, 
          group = Y, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')


# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(cnv.splsda, 
                                comp.predicted=2, 
                                dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(cnv.splsda, comp = 1:2,
          group = Y, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")


#So this is doing an imperfect job of classification, we need to tune our PLSDA model


# Tune our (s)PLS-DA =====

## Select number of components =====

#Select the right number of components by a cross-validation approach

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.cnv <- perf(cnv.splsda, 
                        validation = "Mfold", 
                        folds = 5, 
                        nrepeat = 50,        # use repeated cross-validation
                        progressBar = TRUE, 
                        auc = TRUE,          # include AUC values
                        cpus = 16) 

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.cnv, 
     col = color.mixo(5:7),
     sd = TRUE,
     legend.position = "horizontal")


# print the optimal value of components according to perf()
perf.splsda.cnv$choice.ncomp

#Either 2 or ??8 depending on the measure used


## Select the number of variables per component =====

#Next we can determine the number of variables used to construct each latent component
#This is iterative and you need to use a cross-validation approach

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 200, 10))


# undergo the tuning process to determine the optimal number of variables
tune.splsda.cnv <- tune.splsda(X, Y, 
                               ncomp = 4,         # calculate for first 4 components given the result of perf
                               validation = 'Mfold',
                               folds = 5, 
                               nrepeat = 50,      # use repeated cross-validation
                               dist = 'max.dist', # use max.dist measure
                               measure = "BER",   # use balanced error rate of dist measure
                               test.keepX = list.keepX,
                               cpus = 16)         # allow for paralleliation to decrease runtime

# plot output of variable number tuning
plot(tune.splsda.cnv, col = color.jet(4)) 


# what is the optimal number of components according to tune.splsda()?
tune.splsda.cnv$choice.ncomp$ncomp 

#One, apparently? We need two for plotting at the very least


# what are the optimal number of variables per component according to tune.splsda()
tune.splsda.cnv$choice.keepX

#50 for comp 1, 40 for comp 2

#Store these
# optimal.ncomp <- tune.splsda.cnv$choice.ncomp$ncomp
optimal.ncomp <- 2
optimal.keepX <- tune.splsda.cnv$choice.keepX[1:optimal.ncomp]


# Fit final (s)PLS-DA model ======

# Now we can fit a final (s)PLS-DA model

# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)


#Variance explained
final.splsda$prop_expl_var

## Plot Final model =====

#This is based on the example code from the mixOmics package and are not included in the final manuscript


#Plot our results
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Y, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on ND-GC Dataset, comp 1 & 2')



# Variable stability plots - these show how important a given variable is for a component


# form new perf() object which utilises the final model
perf.splsda.cnv <- perf(final.splsda, 
                        folds = 5, nrepeat = 10, # use repeated cross-validation
                        validation = "Mfold", dist = "max.dist",  # use max.dist measure
                        progressBar = TRUE)


# plot the stability of each feature for the first two components, 'h' type refers to histogram
par(mfrow=c(1,2))
plot(perf.splsda.cnv$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.cnv$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

#Component 2 has rather less stability than 1


# Correlation circle plot
plotVar(final.splsda, comp = c(1,2), cex = 3) # generate correlation circle plot


# Correlation circle plot representing the variables selected by sPLS-DA performed on the 
#CNV data. Only the variables selected by sPLS-DA are shown in components 1 and 2.

#We can see some clustering of variables with similar loadings within questionnaires - the CMD
#variables are in similar locations, as are the SCQ and pra variables for example




# Make Supplementary Figure 1 ======

# Now we make supplementary figure 1 - we will plot the PCA results to illustrate that one component appears to explain
#a particularly large proportion of variation in our dataset, and we also plot the first 3 PLS components
#from our final tuned PLSDA model

#Plot the PCA components
p_pca =
  pca.cnv$prop_expl_var$X |>
  as_tibble(rownames = "Component") |>
  mutate(value = round(value * 100,digits = 2)) |>
  rename(`Variance Explained (%)` = value) |>
  mutate(Component = fct_reorder(Component,`Variance Explained (%)`,max)) |>
  mutate(Component = fct_rev(Component)) |>
  ggplot(aes(x = Component, y = `Variance Explained (%)`)) +
  geom_col(fill = "#7D3232") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = -45)) +
  labs(x = "PCA Component")

#Note the %explaiend variance for the components
pca.cnv$prop_expl_var$X

#Plot the sPLS-DA


#Plot component 1 vs 2
p_pls1 = 
  plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
            group = Y, ind.names = FALSE, # colour by class label
            ellipse = TRUE, legend = FALSE, # include 95% confidence ellipse
            col = c("#E1AF64","#32324B"),
            pch = c(20,20),
            alpha = 0.3,
            size.xlabel = 8,size.ylabel = 8, size.axis = 6,
            point.lwd = 0.2,
            title = ' (a) sPLS-DA on CNV, comp 1 & 2')



#Modify the plots using ggplot
p_pls1 = 
  p_pls1$graph + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  labs(x = "PLSDA Component 1", y = "PLSDA Component 2")
        

#Form the figure using patchwork
sf_1 = p_pca + p_pls1 +  plot_annotation(tag_levels = "A")


#Save the figure
ggsave(filename = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Figures/Figure Parts/sf_1.pdf",
       height = 3,
       width  = 6,
       plot = sf_1)


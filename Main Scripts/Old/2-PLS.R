
# Introduction =======

#We are going to do partial least squares on the full imputed dataset

#Lets load up
pacman::p_load(tidyverse,tidymodels,patchwork, mixOmics)

tidymodels_prefer()


d_last_fit = read_rds("nested_cv_imputed_data.rds")

d = d_last_fit |>  analysis()



# PLSDA casestudy ======


set.seed(5249) # for reproducibility, remove for normal use

data(srbct) # extract the small round bull cell tumour data

X <- srbct$gene # use the gene expression data as the X matrix
Y <- srbct$class # use the class data as the Y matrix

dim(X) # check the dimensions of the X dataframe
length(Y)
summary(Y)

# run pca method on data
pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE) 
plot(pca.srbct)  # barplot of the eigenvalues (explained variance per component)


plotIndiv(pca.srbct, group = srbct$class, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') # onto the PCA subspace



#Now we do (sparse) PLS - DA
srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later


# plot the samples projected onto the first two components of the PLS-DA subspace



plotIndiv(srbct.splsda , comp = 1:2, 
          group = srbct$class, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = srbct$class, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")



# We can do tuning on our sPLS-DA

#Select the right number of components by a cross-validation approach

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = TRUE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")


perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()


# Selecting the number of variables

#Next we can determine the number of variables used to construct each latent component
#This is iterative and you need to use a cross-validation approach

# I notice this approach is kind of similar to the PCA based approach used in the original
#version of the manuscript

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 6, # calculate for first 6 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 10, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 10) # allow for paralleliation to decrease runtime

plot(tune.splsda.srbct, col = color.jet(6)) # plot output of variable number tuning


# what is the optimal value of components according to tune.splsda()?
tune.splsda.srbct$choice.ncomp$ncomp 

# what are the optimal values of variables according to tune.splsda()
tune.splsda.srbct$choice.keepX

#Store these
optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]


# Now we can fit a final (s)PLS-DA model
# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)


#Plot our results
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = srbct$class, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = srbct$class, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')


# set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Tumour Type", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), 
           legend = legend)

cim

#Clustered Image Map of the genes selected by sPLS-DA on the SRBCT gene expression data 
#across all 3 components. A hierarchical clustering based on the gene expression levels 
#of the selected genes, with samples in rows coloured according to their tumour subtype 
#(using Euclidean distance with Complete agglomeration method).


# Variable stability plots - these show how important a given variable is for a component


# form new perf() object which utilises the final model
perf.splsda.srbct <- perf(final.splsda, 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.srbct$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.srbct$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.srbct$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)


# Correlation circle plot
var.name.short <- substr(srbct$gene.name[, 2], 1, 10) # form simplified gene names

plotVar(final.splsda, comp = c(1,2), var.names = list(var.name.short), cex = 3) # generate correlation circle plot

# Correlation circle plot representing the genes selected by sPLS-DA performed on the 
#SRBCT gene expression data. Gene names are truncated to the first 10 characters. 
#Only the genes selected by sPLS-DA are shown in components 1 and 2.


## You can use the PLSDA model directly for prediction =====

#You need to split into test and training sets - 

train <- sample(1:nrow(X), 50) # randomly select 50 samples in training
test <- setdiff(1:nrow(X), train) # rest is part of the test set

# store matrices into training and test set:
X.train <- X[train, ]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]

# train the model
train.splsda.srbct <- splsda(X.train, Y.train, ncomp = optimal.ncomp, keepX = optimal.keepX)


# use the model on the Xtest set
predict.splsda.srbct <- predict(train.splsda.srbct, X.test, 
                                dist = "mahalanobis.dist")


# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda.srbct$class$mahalanobis.dist[,2]
table(factor(predict.comp2, levels = levels(Y)), Y.test)

# evaluate the prediction accuracy for the first three components
predict.comp3 <- predict.splsda.srbct$class$mahalanobis.dist[,3]
table(factor(predict.comp3, levels = levels(Y)), Y.test)


auc.splsda = auroc(final.splsda, roc.comp = 1, print = FALSE) # AUROC for the first component

auc.splsda = auroc(final.splsda, roc.comp = 3, print = FALSE) # AUROC for all three components


# Lets do the same thing with our data =====

## Set up =====



set.seed(05082022) # for reproducibility


X <-  d |> select(-group) |> as.matrix()
Y <-  d |> select(group) |> pull()

dim(X) # check the dimensions of the X dataframe
length(Y)
summary(Y)

## Initial PCA =====

# run pca method on data
pca.cnv = pca(X, ncomp = 10, center = TRUE, scale = TRUE) 
plot(pca.cnv)  # barplot of the eigenvalues (explained variance per component)

#So it looks a lot like there is one big component and maybe another 2 which explain
#most of the variance in the data


plotIndiv(pca.cnv, group = Y, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on CNV, comp 1 - 2') # onto the PCA subspace


#So we see a cluster of controls and a wider distribution of CNV carriers

## Basic Sparse PLSDA =====

#Now we do (sparse) PLS - DA
cnv.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later


# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(cnv.splsda , comp = 1:2, 
          group = Y, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(cnv.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(cnv.splsda, comp = 1:2,
          group = Y, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")


#So this is not doing a perfect job of classification, we need to tune our PLSDA model


## Tune our (s)PLS-DA =====

# We can do tuning on our sPLS-DA

### Select number of components =====

#Select the right number of components by a cross-validation approach

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.cnv <- perf(cnv.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 50, # use repeated cross-validation
                          progressBar = TRUE, auc = TRUE, cpus = 4) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.cnv, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")


perf.splsda.cnv$choice.ncomp # what is the optimal value of components according to perf()

#Either 2 or 3 depending on the measure used


### Select the number of variables per component =====

#Next we can determine the number of variables used to construct each latent component
#This is iterative and you need to use a cross-validation approach

# I notice this approach is kind of similar to the PCA based approach used in the original
#version of the manuscript

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.cnv <- tune.splsda(X, Y, ncomp = 3, # calculate for first 3 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 10) # allow for paralleliation to decrease runtime

plot(tune.splsda.cnv, col = color.jet(3)) # plot output of variable number tuning


# what is the optimal value of components according to tune.splsda()?
tune.splsda.cnv$choice.ncomp$ncomp 

# what are the optimal values of variables according to tune.splsda()
tune.splsda.cnv$choice.keepX

#Store these
optimal.ncomp <- tune.splsda.cnv$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.cnv$choice.keepX[1:optimal.ncomp]


## Fit final (s)PLS-DA model ======

# Now we can fit a final (s)PLS-DA model

# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)



## Plot Final model =====

#Plot our results
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Y, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = Y, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')

plotIndiv(final.splsda, comp = c(2,3), # plot samples from final model
          group = Y, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 2 & 3')



# Variable stability plots - these show how important a given variable is for a component


# form new perf() object which utilises the final model
perf.splsda.cnv <- perf(final.splsda, 
                        folds = 5, nrepeat = 10, # use repeated cross-validation
                        validation = "Mfold", dist = "max.dist",  # use max.dist measure
                        progressBar = TRUE)


# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.cnv$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.cnv$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.cnv$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)

#This doesn't really work

# Correlation circle plot
plotVar(final.splsda, comp = c(1,2), cex = 3) # generate correlation circle plot

# Correlation circle plot representing the variables selected by sPLS-DA performed on the 
#CNV data. Only the variables selected by sPLS-DA are shown in components 1 and 2.



## You can use the PLSDA model directly for prediction =====

#You need to split into test and training sets - 
X.train <- X
Y.train <- Y

X.test  <- d_last_fit |> assessment() |> select(-group) |> as.matrix()
Y.test  <- d_last_fit |> assessment() |> select(group)  |> pull()

# train the model
train.splsda.cnv <- splsda(X.train, Y.train, ncomp = optimal.ncomp, keepX = optimal.keepX)


# use the model on the Xtest set
predict.splsda.cnv <- predict(train.splsda.cnv, X.test, 
                                dist = "mahalanobis.dist")


# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda.cnv$class$mahalanobis.dist[,2]
table(factor(predict.comp2, levels = levels(Y)), Y.test)

# evaluate the prediction accuracy for the first three components
predict.comp3 <- predict.splsda.cnv$class$mahalanobis.dist[,3]
table(factor(predict.comp3, levels = levels(Y)), Y.test)


auc.splsda = auroc(final.splsda, roc.comp = 1, print = FALSE) # AUROC for the first component

auc.splsda = auroc(final.splsda, roc.comp = 3, print = FALSE) # AUROC for all three components


# Make Figure 2 ======

# Now lets reflect on what we want from a plot. The original version of the manuscript had the PCA 
# eigenvalues/explained variance - we could have that, and also have the PLS component plots
# after tuning

#Plot the PCA componenets
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


#Plot the sPLS-DA


#Plot our results
p_pls1 = plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
                   group = Y, ind.names = FALSE, # colour by class label
                   ellipse = TRUE, legend = FALSE, # include 95% confidence ellipse
                   col = c("#E1AF64","#32324B"),
                   pch = c(20,20),alpha = 0.3,
                   size.xlabel = 8,size.ylabel = 8, size.axis = 6,
                   point.lwd = 0.2,
                   title = ' (a) sPLS-DA on CNV, comp 1 & 2')

p_pls2 = plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
                   group = Y, ind.names = FALSE, # colour by class label
                   ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
                   col = c("#E1AF64","#32324B"),
                   pch = c(20,20),alpha = 0.3,
                   size.xlabel = 8,size.ylabel = 8, size.axis = 6,
                   point.lwd = 0.2,
                   title = '(b) sPLS-DA on CNV, comp 1 & 3')


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
        
p_pls2 = 
  p_pls2$graph + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  labs(x = "PLSDA Component 1", y = "PLSDA Component 3")

#Form the figure
figure_2 = p_pca + p_pls1 + p_pls2 + plot_annotation(tag_levels = "A")

ggsave(filename = "./Figures/figure_2.pdf",
       height = 3,width = 8,
       plot = figure_2)


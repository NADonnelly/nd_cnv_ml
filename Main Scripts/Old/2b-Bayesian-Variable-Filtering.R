# This time, with Bayes ======

#Can we acheive regularisation with the use of the horseshoe prior?
library(tidyverse)
library(tidybayes)
library(rstanarm)
library(tidymodels)

#Use our somewhat naughty imputed dataset
DF <- read_csv("ImputedData2.csv")

# options(mc.cores = 12)

set.seed(1234)

DF = 
  DF |>
  relocate(group) |>
  mutate(group = factor(group, levels = c("Control","ND-CNV")))

#Make our horseshoe prior
n <- nrow(DF)
D <- ncol(DF)-1

#How many variables do you think will be useful predictors? Lets say 30 as a 
#starting value

p0 <- 30

tau0 <- p0/(D - p0) * 1/sqrt(n)
prior_coeff <- hs(global_scale = tau0, slab_scale = 1)


#Right well, lets fit it
fit_1 = 
  stan_glm(group ~ ., 
           data = DF,
           family = binomial(link = "logit"),
           prior           = prior_coeff, 
           prior_intercept = normal(0, 1.5),
           cores = 2, seed = 12345,
           chains = 4, iter = 4000, warmup = 1000,
           refresh = 1,
           adapt_delta = 0.999)



#Plot our top 30 predictors

left_join(
  fit_1 |>
    coef() |>
    as_tibble(rownames = "variable") |>
    rename(estimate = value),
  fit_1 |>
    posterior_interval(prob = 0.90) |>
    as_tibble(rownames = "variable"),
  by = "variable") |>
  filter(variable != "(Intercept)") |>
  slice_max(abs(estimate), n = 30) |>
  mutate(
    variable = fct_reorder(variable, abs(estimate))
  ) |>
  ggplot(aes(x = estimate, y = variable, xmin = `5%`,xmax = `95%`)) +
  geom_errorbar() +
  geom_point() +
  labs(y = NULL) +
  theme_bw() +
  geom_vline(xintercept = 0,lty = 3)


#We can save the variables we extracted thus
left_join(
  fit_1 |>
    coef() |>
    as_tibble(rownames = "variable") |>
    rename(estimate = value),
  fit_1 |>
    posterior_interval(prob = 0.90) |>
    as_tibble(rownames = "variable"),
  by = "variable") |>
  filter(variable != "(Intercept)") |>
  slice_max(abs(estimate), n = 30) |>
  mutate(
    variable = fct_reorder(variable, abs(estimate))
  )

bayes_horseshoe_vars = 
  left_join(
    fit_1 |>
      coef() |>
      as_tibble(rownames = "variable") |>
      rename(estimate = value),
    fit_1 |>
      posterior_interval(prob = 0.90) |>
      as_tibble(rownames = "variable"),
    by = "variable") |>
  filter(variable != "(Intercept)") |>
  slice_max(abs(estimate), n = 30) |>
  mutate(
    variable = fct_reorder(variable, abs(estimate))
  ) |>
  mutate(type = "Horseshoe")

## Save selected variables ====

# we could save this
write_csv(bayes_horseshoe_vars, "selected_variables_horseshoe.csv")


# Now lets look at this projpred thing ======

## Vignette example ======

# library(projpred)
# library(rstanarm)
# library(posterior)
# library(bayesplot)
# bayesplot_theme_set(ggplot2::theme_bw())
# 
# data("df_gaussian", package = "projpred")
# 
# dat_gauss <- 
#   data.frame(y = df_gaussian$y, df_gaussian$x)
# 
# # Number of regression coefficients:
# ( D <- sum(grepl("^X", names(dat_gauss))) )
# 
# 
# # Prior guess for the number of relevant (i.e., non-zero) regression
# # coefficients:
# p0 <- 5
# # Number of observations:
# N <- nrow(dat_gauss)
# # Hyperprior scale for tau, the global shrinkage parameter (note that for the
# # Gaussian family, 'rstanarm' will automatically scale this by the residual
# # standard deviation):
# tau0 <- p0 / (D - p0) * 1 / sqrt(N)
# 
# 
# ncores <- 8
# 
# ### Only for technical reasons in this vignette (you can omit this when running
# ### the code yourself):
# 
# ###
# options(mc.cores = ncores)
# 
# refm_fit <- stan_glm(
#   y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 +
#     X15 + X16 + X17 + X18 + X19 + X20,
#   family = gaussian(),
#   data = dat_gauss,
#   prior = hs(global_scale = tau0),
#   chains = 4, iter = 2000,
#   ###
#   seed = 2052109, QR = TRUE, refresh = 0
# )
# 
# 
# refm_fit |> summary()
# 
# #Do the cross validation variable selection thing
# cvvs <- 
#   refm_fit |>
#   cv_varsel(
#     validate_search = TRUE,
#     nclusters_pred = 20,
#     ###
#     nterms_max = 9,
#     seed = 411183
# )
# 
# plot(cvvs, stats = c("elpd", "rmse"), deltas = TRUE)
# 
# ( soltrms <- solution_terms(cvvs) )
# 
# modsize_decided <- 6
# 
# 
# ( soltrms_final <- head(soltrms, modsize_decided) )
# 
# 
# #Now we can project the refernece model onto the final submodel
# prj <- project(
#   refm_fit,
#   solution_terms = soltrms_final,
#   seed = 15705533
# )
# 
# 
# # Get the projected posterior draws
# prj_mat <- as.matrix(prj)
# 
# 
# #Then process them
# prj_drws <- as_draws_matrix(prj_mat)
# 
# # In the following call, as.data.frame() is used only because pkgdown
# # versions > 1.6.1 don't print the tibble correctly.
# as.data.frame(summarize_draws(
#   prj_drws,
#   "median", "mad", function(x) quantile(x, probs = c(0.025, 0.975))
# ))
# 
# #Plot marginal posteriors
# mcmc_intervals(prj_mat) +
#   ggplot2::coord_cartesian(xlim = c(-1.5, 1.6))
# 
# # Do an equivalent of pp_check
# prj_predict <- proj_predict(prj, .seed = 762805)
# # Using the 'bayesplot' package:
# ppc_dens_overlay(y = dat_gauss$y, yrep = prj_predict, alpha = 0.9, bw = "SJ")


## Apply our own dataset ======

# Load packages
library(projpred)
library(rstanarm)
library(posterior)
library(bayesplot)
library(tidymodels)
library(tidyverse)

bayesplot_theme_set(ggplot2::theme_bw())

# Same as above

#Use our somewhat naughty imputed dataset
DF <- read_csv("ImputedData2.csv")

DF = 
  DF |>
  relocate(group) |>
  mutate(group = factor(group, levels = c("Control","ND-CNV")))

options(mc.cores = 2)

set.seed(1234)


#Make our horseshoe prior
n <- nrow(DF)
D <- ncol(DF)-1

#How many variables do you think will be useful predictors? Lets aim for 30 max
p0 <- 30

tau0 <- p0/(D - p0) * 1/sqrt(n)
prior_coeff <- hs(global_scale = tau0, slab_scale = 1)


#Right well, lets fit it
fit_1 = 
  stan_glm(group ~ ., 
           data = DF1,
           family = binomial(link = "logit"),
           prior           = prior_coeff, 
           prior_intercept = normal(0, 1.5),
           cores = 8, seed = 12345,
           chains = 4, iter = 4000, warmup = 1000,
           refresh = 0,
           adapt_delta = 0.999)



fit_1 |> summary()

#Do the cross validation variable selection thing - this is very slow but doesn't use
#a lot of CPU power, which is very frustrating
cvvs <- 
  fit_1 |>
  cv_varsel(
    cv_method = "LOO",
    validate_search = TRUE,
    nclusters_pred = 30,
    nterms_max = 30,
    seed = 411183
  )

#Plot our results
plot(cvvs, stats = c("auc","acc"), deltas = TRUE)

#How many terms should we take?
suggest_size(cvvs)


#The package itself suggests we take 11, I am going to take 25 as well as that 
#is the point at which accuracy hits the same as the full model
modsize_decided <- 25

( soltrms <- solution_terms(cvvs) )
( soltrms_final <- head(soltrms, modsize_decided) )


#Now we can project the reference model onto the final submodel
prj <- 
  project(
  fit_1,
  solution_terms = soltrms_final,
  seed = 15705533
)


# Get the projected posterior draws
prj_mat <- as.matrix(prj)


#Then process them
prj_drws <- posterior::as_draws_matrix(prj_mat)

# In the following call, as.data.frame() is used only because pkgdown
# versions > 1.6.1 don't print the tibble correctly.

prj_drws |>
  posterior::summarize_draws("median", "mad", 
                  function(x) quantile(x, probs = c(0.025, 0.975))) |>
  filter(variable != "(Intercept)") 



#Plot marginal posteriors
mcmc_intervals(prj_mat) +
  ggplot2::coord_cartesian(xlim = c(-1.5, 1.6))

# Do an equivalent of pp_check
prj_predict <- proj_predict(prj, .seed = 762805)


# We can save the variables we identified here

bayes_projpred_vars = 
  prj_drws |>
  posterior::summarize_draws("median", "mad", 
                             function(x) quantile(x, probs = c(0.025, 0.975))) |>
  filter(variable != "(Intercept)") |>
  mutate(type = paste("projpred_",modsize_decided,sep = ""))

write_csv(bayes_projpred_vars, paste("selected_variables_projpred_",modsize_decided, ".csv",sep = ""))

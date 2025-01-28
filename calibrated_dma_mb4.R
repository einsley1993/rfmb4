

library(tidyverse) 
library(knitr)
library(here)
library(tableone)
library(corrr)
library(robumeta)
library(MetaUtility)
library(fastDummies)
library(weightr)
library(PublicationBias)
library(xtable)
library(boot)
library(testthat)
library(data.table)
library(ggplot2)
library(Replicate)
library(clubSandwich) 
library(devtools)


getwd()
setwd("/Users/yili/Desktop/Meta Analysis and Replication Forecasting/MB4")

# read in data
dma = fread("mb4_dma.csv")
dr = fread("mb4_dr.csv")

# remove one study with NA value for n_excluded_1
dr <- dr[!is.na(dr$n_excluded_1),]


dma_yi_sei <- dma %>% 
  select(study_ID, yi, vi) %>% 
  mutate(sei = sqrt(vi))


###########################################################################
# 1. For uncorrected MA: use Wang & Lee's method with calibrated estimates
###########################################################################
### First, reproduce tau2 hat and mu hat using Wang & Lee, Section 2
# Use metafor package:
library(metafor)
meta.re = rma.uni( yi = dma$yi,
                   vi = dma$vi,
                   method = "REML") # default: REML. Previously used "DL"
meta.re
#  estimate: 0.3296

t2hat.naive = meta.re$tau2
t2hat.naive
# 0.02343389

t2hat.naive_dma <- t2hat.naive


# Use robumeta packgae:
dmaMod = robu( yi ~ 1,
               data = dma, 
               studynum = as.factor(study_ID),
               var.eff.size = vi,
               modelweights = "HIER",
               small = TRUE)
dmaMod
# Estimate StdErr t-value  dfs        P(|t|>) 95% CI.L 95% CI.U Sig
# 1 X.Intercept.    0.333 0.0391    8.53 32.9 0.000000000759    0.254    0.413 ***
# Tau.sq = 0.02389638 
# tau2 different, because here we specify modelweights

dmaMod$b.r
# 0.3331773

k <-  nrow(dma)
k

dma_cop1 <- dma %>% select(yi, vi) %>% 
  mutate(
         Qi = (yi - sum(yi/vi)/sum(1/vi))^2/vi
         ) 
sum_inv_vi <-  sum(1/dma_cop1$vi)
sum_inv_vi2 <-  sum(1/dma_cop1$vi^2)
Q <-  sum(dma_cop1$Qi)

sum_inv_vi
sum_inv_vi2
Q

## tau2 hat:
tau2 <-  (Q - (k-1)) / (sum_inv_vi - sum_inv_vi2/sum_inv_vi)
tau2
# 0.02343389

## mu hat: 
dma_cop2 <- dma_cop1 %>% 
  mutate(
    mui_num = yi/(vi + tau2),
    mui_denom = 1/(vi + tau2)
  ) 
mu <- sum(dma_cop2$mui_num) / sum(dma_cop2$mui_denom)
mu
# 0.3329477




### Wang & Lee, Section 2.2 Conventional method:
# Point estimate: mu
# Var(mu):
var_mu <- 1 / sum(dma_cop2$mui_denom)
var_mu
# 0.001164185
sqrt(var_mu)
# 0.03412016


### Wang & Lee, Section 2.3 Calibrated estimate Eq(1):
dma_yi_cali_wang <- dma %>% 
  select(yi, vi) %>% 
  mutate(yi_cali_wang = mu + sqrt(tau2/(tau2+vi))*(yi-mu))
var(dma_yi_cali_wang$yi_cali_wang)
# 0.02348293

# Original var(yi): larger than tau2, because of overdispersion
var(dma$yi)
# 0.2200958

# write.csv(dma_yi_cali_wang, "dma_yi_cali_wang.csv")



### Wang & Lee, Section 2.4 Prediction interval
# Q(i/{k+1}) = hat theta_(i)*
# Sort yi_cali_wang, assign p = i/(k+1)
dma_order <- dma_yi_cali_wang %>% 
  select(yi, vi, yi_cali_wang) %>% 
  arrange(yi_cali_wang) %>% 
  mutate(i = 1:k,
         p = i/(k+1),
         yi_cali_wang_lead = case_when(i < k ~ lead(yi_cali_wang),
                                  TRUE ~ 999))
dma_order <- dma_order %>% 
  mutate(Q_p = case_when(p < 1/(k+1) ~ mu + (dma_order$yi_cali_wang[1] - mu)/qnorm(1/(k+1)) * qnorm(p),
                         p < (i+1)/(k+1) ~ yi_cali_wang + (yi_cali_wang_lead - yi_cali_wang) / (qnorm((i+1)/(k+1)) - qnorm(i/(k+1))) * (qnorm(p) - qnorm(i/(k+1))),
                         p > k/(k+1) ~ mu + (yi_cali_wang - mu)/qnorm(k/(k+1)) * qnorm(p))
         )

hist(dma_order$Q_p)
mean(dma_order$Q_p)
quantile(dma_order$Q_p, c(0.025, 0.975))
# 2.5%      97.5% 
# 0.06675066 0.58583667 



### Wang & Lee, Section 2.5 Prediction distribution
range(dma$yi)
# -0.8873296  1.7746593

mean(dma_order$yi_cali_wang)
# 0.3480867
range(dma_order$yi_cali_wang)
k

# bandwidth: 0.1 - 0.5
# x: [-1, 5], step size = 0.01
h <- 0.2

dma_pred_dist <- tibble(x = seq(-2, 2, 0.005), fh = NA)
for(i in 1:nrow(dma_pred_dist)){
  dma_pred_dist$fh[i] <- sum(dnorm((dma_pred_dist$x[i] - dma_order$yi_cali_wang)/h))/(k*h)
}

ggplot(dma_pred_dist, aes(x = x, y = fh)) + 
  geom_point(size = 0.5, color = "blue") +
  # plot calibrated estimates as well as points:
  geom_point(data = dma_yi_cali_wang,
             aes(x = yi_cali_wang, y = 0),
             size = 3,
             shape = 18,
             color =  "orange",
             alpha = 1 ) + 
  
  # Plot the density curve if using conventional method (Wang & Lee Section 2.2)
  stat_function(fun = dnorm, n = 601, color = "darkgray", size = 1,
                args = list(mean = mu, sd = sqrt(tau2 + var_mu))) +
  scale_x_continuous(breaks = seq(-2, 4, 1)) +
  scale_y_continuous(breaks = seq(0, 0.9, 0.1)) +
  ggtitle("Predictive distributions (Wang & Lee 2018 method). Dataset = dma.")

# ggsave("Wang_Lee_calibrated_pred_dists_dma.pdf", height = 6, width = 6)





######################################################################################
# 2. For uncorrected and corrected MR, prediction intervals: Mathur & VanderWeele 2020
######################################################################################
### Section 2.2 Extension to MR
# Calculate calibrated theta_i(Z = z).
# Z is set to the Z_i levels in dr, for each replication study.
# In dr:
summary(dr$mean_age_1)
# mean = 253.1
hist(dr$mean_age_1)

# Take from R01 code: Meta-Regression (Has Mods), no pub bias
# In dma:
summary(dma$mean_age_1)
# mean = 401.3 
table(dma$method)
# help reach 
# 20   118 

# common variables in dma and dr:
intersect(names(dma), names(dr))

intersect_vars <- intersect(names(dma), names(dr))[3:16]
intersect_vars
# [1] "response_mode"      "exposure_phase"     "method"             "dependent_measure" 
# [5] "participant_design" "infant_type"        "mean_age_1"         "scenarios"         
# [9] "stimuli"            "stimuli_type"       "choice_object"      "Hamlin_Lab"        
# [13] "intent_val"         "outcome_val"       

# variables that have more than one level in dma:
# method, scenarios, stimuli, stimuli_type, 
# choice_object, Hamlin_Lab, intent_val, outcome_val
table(dma$method)
table(dr$method) # reach

table(dma$scenarios) 
table(dr$scenarios) # help_hinder

table(dma$stimuli)
table(dr$stimuli) # movies

table(dma$stimuli_type)
table(dr$stimuli_type) # real

table(dma$choice_object)
table(dr$choice_object) # shapes

table(dma$Hamlin_Lab)
table(dr$Hamlin_Lab)

table(dma$intent_val)
table(dr$intent_val) # opp

table(dma$outcome_val)
table(dr$outcome_val) # opp


### MR conditional on Z, uncorrected for pub bias
mrMod2 = robu( yi ~ 
                 mean_age_1 + n_excluded_1 +
                 method + scenarios + stimuli + stimuli_type + 
               choice_object + Hamlin_Lab + intent_val + outcome_val,
               data = dma, 
               studynum = as.factor(study_ID),
               var.eff.size = vi,
               modelweights = "HIER", 
               small = TRUE )

mrMod2

tau2_mrMod2 <- mrMod2$mod_info$tau.sq
tau2_mrMod2
# 0.00644501

# In Eq(2.3), hat theta_i is y_i here. 
# The z level to be conditioned on are specified as below: for now, condition on mean and modal
# Mode in dma:
table(dma$method) 
# help reach 
# 20   118 
table(dma$scenarios) 
# fair_unfair   give_take help_hinder       other 
# 29          14          84          11 
table(dma$stimuli)

# live_show    movies 
# 85        53 
table(dma$stimuli_type)
# cartoon    real 
# 24     114 
table(dma$choice_object)
# experimenters        people       puppets        shapes 
# 29             2            69            38 
table(dma$Hamlin_Lab)
# FALSE  TRUE 
# 91    47 
table(dma$intent_val)
# neg opp pos 
# 9 112  17 
table(dma$outcome_val)
# neut  opp  rev same 
# 15  104    1   18 


# Below is for plotting PDF for the 1st REP study:
Z_mod = as.matrix( data.frame( X.Intercept. = 1,
                               mean_age_1 = dr$mean_age_1[1], 
                               n_excluded_1 = dr$n_excluded_1[1],
                               methodreach = 0,
                               scenariosgive_take = 0,
                               scenarioshelp_hinder = 1,
                               scenariosother = 0,
                               stimulimovies = 0,
                               stimuli_typereal = 1,
                               choice_objectpeople = 0,
                               choice_objectpuppets = 1,
                               choice_objectshapes = 0,
                               Hamlin_LabTRUE = as.numeric(dr$Hamlin_Lab[1] == TRUE),
                               intent_valopp = 1,
                               intent_valpos = 0,
                               outcome_valopp = 1,
                               outcome_valrev = 0,
                               outcome_valsame = 0
                              ) 
                        )

preds_Z_mod <- predict(mrMod2, Z_mod)[1]
preds_Z_mod
# This is hat E[theta | Z = z] = \beta_0 + \beta_1 Z.

# Now compute \beta_0 + \beta_1 Z_i. Still use the predict function.
Zi_mat = as.matrix( data.frame( X.Intercept. = 1,
                                mean_age_1 = mean(dma$mean_age_1),
                                n_excluded_1 = mean(dma$n_excluded_1),
                                methodreach = dma$method == "reach",
                                scenariosgive_take = dma$scenarios == "give_take",
                                scenarioshelp_hinder = dma$scenarios == "help_hinder",
                                scenariosother = dma$scenarios == "other",
                                stimulimovies = dma$stimuli == "movies",
                                stimuli_typereal = dma$stimuli_type == "real",
                                choice_objectpeople = dma$choice_object == "people",
                                choice_objectpuppets = dma$choice_object == "puppets",
                                choice_objectshapes = dma$choice_object == "shapes",
                                Hamlin_LabTRUE = dma$Hamlin_Lab == "TRUE",
                                intent_valopp = dma$intent_val == "opp",
                                intent_valpos = dma$intent_val == "pos",
                                outcome_valopp = dma$outcome_val == "opp",
                                outcome_valrev = dma$outcome_val == "rev",
                                outcome_valsame = dma$outcome_val == "same"
                                ) 
                    )

preds_Zi = apply( Zi_mat,
                  MARGIN = 1,
                  FUN = function(x) predict( object = mrMod2, pred.vector = x) )

preds_Zi = as.data.frame( t(preds_Zi) )


dma_cop3 <- dma %>% 
  select(V1, study_ID,
         mean_age_1, n_excluded_1,
         method, scenarios, stimuli, stimuli_type,
         choice_object, Hamlin_Lab, intent_val, outcome_val,
         yi, vi) %>% 
  mutate(pred_Zi = preds_Zi$prediction) %>% 
  mutate(tilde_yi_Z = preds_Z_mod + 
           sqrt(t2hat.naive_dma/(t2hat.naive_dma + vi)) * (yi - pred_Zi))



var(dma_cop3$tilde_yi_Z)
# 0.01566326

# Different from this:
t2hat.naive_dma
# 0.01721693


### Calculate joint density
# f_Y|X (t):
sigma_y_given_x <- sqrt(dr$vi[1]) # use vi in dr study, one at a time

# Generate a range of values to be plotted (on X-axis later)
w_vec <- seq(-5, 5, by = 0.01)

# Calculate the conditional PDF f_Y|X(t|w). 
# Here, Y: theta hat, X: true theta (from theta tilde's distribution)
pdf_y_given_x <- function(t, w) {
  dnorm(t, mean = w, sd = sigma_y_given_x)
}


# Calculate f_Y(t) = int f_X(w) * f_Y|X(t) dw
t_vec <- seq(-5, 5, by = 0.01)
density_matrix <- matrix(NA, nrow = length(w_vec), 
                         ncol = length(t_vec))

dma_cop3_order <- dma_cop3 %>% 
  select(V1, yi, vi, tilde_yi_Z) %>% 
  arrange(tilde_yi_Z)

# Calculate the joint density f_X(w) * f_Y|X(t) over w
for(j in 1:length(t_vec)){
  t <- t_vec[j]
  
  for(i in 1:length(w_vec)){
    w <- w_vec[i]
    
    # f(theta): true theta's density
    pdf_x <- sum(dnorm((w - dma_cop3_order$tilde_yi_Z)/h))/(nrow(dr)*h)
    
    # f(theta hat):  
    joint_density <- pdf_x * pdf_y_given_x(t, w) #f_X(w) * f_Y|X(t|w)
    
    density_matrix[i, j] <- joint_density
  }
  
}

dim(density_matrix)

# integrated density (integrate over w)
int_dens <- colSums(density_matrix)



# Normalize above to make it a proper probability density
pdf_y <- int_dens / sum(int_dens)

# Dataset for plotting:
plot_dat <- data.frame(t = t_vec,
                       f_y = pdf_y)
sum(plot_dat$f_y)


# Plot the densities
ggplot(plot_dat, aes(x = t, y = f_y)) +
  geom_line() +
  ylab("density f_Y(t)") +
  theme_minimal()

# ggsave("Density of theta hat, uncorrected (1st replication study in dr).pdf", 
       # heigh = 4, width = 4)

t_vec[which.max(pdf_y)]

ecdf_dat <- plot_dat %>% mutate(ecdf = cumsum(f_y))



# Find 95% prediction interval
# 2.5%:
q1 <- 0.025
for (i in 1:nrow(ecdf_dat)) {
  if (ecdf_dat$ecdf[i] > q1) {
    break
  }
}

# Row number: i - 1. The last row less than 2.5%
ecdf_dat$t[i - 1]


# 97.5%:
q2 <- 0.975

for (j in 1:nrow(ecdf_dat)) {
  if (ecdf_dat$ecdf[j] >= q2) {
    break
  }
}

# Row number: j + 1. The first row greater than 97.5%
ecdf_dat$t[j + 1]




### MR conditional on Z, CORRECTED for pub bias
mods = c( "mean_age_1","n_excluded_1",
          "method", "scenarios", "stimuli", "stimuli_type", 
          "choice_object", "Hamlin_Lab", "intent_val", "outcome_val")

gotError = TRUE  # initialize so the while-loop is entered
while ( gotError == TRUE ) {
  
  linPred = paste( "yi ~", paste(mods, collapse = " + ") )
  
  tryCatch({
    mrMod2 = robu( eval( parse(text = linPred) ), 
                   data = dma, 
                   studynum = as.factor(study_ID),
                   var.eff.size = vi,
                   modelweights = "HIER",
                   small = TRUE)
    
    gotError = FALSE
    
  }, error = function(err){
    gotError <<- TRUE
    
    # remove one moderator from the end of the list
    cat( paste( "\n Removing ", 
                mods[ length(mods) ],
                " from moderators and trying again" ) )
    mods <<- mods[ 1 : ( length(mods) - 1 ) ]
    
  })
  
}

mods




eta = 4.70
pvals = dma$p_value 
dma$affirm = dma$p_value < 0.05
# weight for model
weights = rep( 1, length(pvals) )
weights[ dma$affirm == FALSE ] = eta # if non-affirmative, set weights to selection ratio
weights

tau2_mrMod2

mrMod3 = robu( eval( parse( text = linPred ) ),
               studynum = as.factor(study_ID),
               data = dma,
               userweights = weights / (dma$vi + tau2_mrMod2),
               var.eff.size = vi,
               small = TRUE )
mrMod3

Xmat = as.matrix( data.frame( X.Intercept. = 1,
                              mean_age_1 = dr$mean_age_1,
                              n_excluded_1 = dr$n_excluded_1,
                              methodreach = dr$method == "reach",
                              scenariosgive_take = dr$scenarios == "give_take",
                              scenarioshelp_hinder = dr$scenarios == "help_hinder",
                              scenariosother = dr$scenarios == "other",
                              stimulimovies = dr$stimuli == "movies",
                              stimuli_typereal = dr$stimuli_type == "real",
                              choice_objectpeople = dr$choice_object == "people",
                              choice_objectpuppets = dr$choice_object == "puppets",
                              choice_objectshapes = dr$choice_object == "shapes",
                              Hamlin_LabTRUE = dr$Hamlin_Lab == "TRUE",
                              intent_valopp = dr$intent_val == "opp",
                              intent_valpos = dr$intent_val == "pos",
                              outcome_valopp = dr$outcome_val == "opp",
                              outcome_valrev = dr$outcome_val == "rev",
                              outcome_valsame = dr$outcome_val == "same"
                              ) 
                  )
  

# sanity check
expect_equal( attr(Xmat, "dimnames")[[2]], mrMod3$labels )



### Upload unredacted dr dataset:
dr_full = fread("REP_full.csv")

# remove one study with NA value for n_excluded_1
dr_full <- dr_full[!is.na(dr_full$n_excluded_1),]

# yi in dr:
dr_full$yi <- dr_full$d_calc
dr_full$vi <- dr_full$d_var_calc


# Next, will write a function that outputs densities. 
# Input is the model, i.e. mrMod2 or mrMod3. 

######################################
### Now do this for all REPs in dr ### ### This needs yi in dr
######################################

# Generate a range of values for density curves
w_vec <- seq(-5, 5, by = 0.01)
t_vec <- seq(-5, 5, by = 0.01)

# Function of conditional PDF f_Y|X(t|w)
pdf_y_given_x <- function(t, w) {
  dnorm(t, mean = w, sd = sigma_y_given_x)
}


q1 <- 0.025
q2 <- 0.975


# SAPB Function:
SAPB_interval_func <- function(fitted.model){
  
  dr_pred <- dr_full %>% # need REP point estimates 
    select(lab_id, 
           yi, vi) %>% 
    mutate(pred_yi_peak = NA, # predicted yi with highest density
           pred_yi_lo = NA,
           pred_yi_hi = NA)
  
  # fitted.model = mrMod2
  # Get preds_Zi:
  preds_Zi = apply( Zi_mat,
                    MARGIN = 1,
                    FUN = function(x) predict( object = fitted.model, pred.vector = x) )
  
  preds_Zi = as.data.frame( t(preds_Zi) )
  
  dma_reduc <- dma %>%
    select(study_ID, yi, vi, mean_age_1, n_excluded_1, method, scenarios, 
           stimuli, stimuli_type, choice_object, Hamlin_Lab, intent_val, outcome_val) %>% 
    mutate(pred_Zi = preds_Zi$prediction)
  
  # Get prediction interval for each row
  for(k in 1:nrow(dr_pred)){
    # Pre-process Z_i variables to dummy variables
    Z_mod = as.matrix( data.frame( X.Intercept. = 1,
                                        mean_age_1 = dr_full$mean_age_1[k],
                                        n_excluded_1 = dr_full$n_excluded_1[k],
                                        methodreach = 0,
                                        scenariosgive_take = 0,
                                        scenarioshelp_hinder = 1,
                                        scenariosother = 0,
                                        stimulimovies = 0,
                                        stimuli_typereal = 1,
                                        choice_objectpeople = 0,
                                        choice_objectpuppets = 1,
                                        choice_objectshapes = 0,
                                        Hamlin_LabTRUE = as.numeric(dr_full$Hamlin_Lab == TRUE)[k],
                                        intent_valopp = 1,
                                        intent_valpos = 0,
                                        outcome_valopp = 1,
                                        outcome_valrev = 0,
                                        outcome_valsame = 0 ) 
    )
    
    preds_Z_mean_mod <- predict(fitted.model, Z_mod)[1]
    # This is hat E[theta | Z = z] = \beta_0 + \beta_1 Z.
    
    dma_reduc <- dma_reduc %>% 
      mutate(tilde_yi_Z = preds_Z_mod + 
               sqrt(tau2_mrMod2/(tau2_mrMod2 + vi)) * (yi - pred_Zi))
    
    
    # f_Y|X (t):
    sigma_y_given_x <- sqrt(dr$vi[k]) # use vi in dr study, one at a time
    
    
    density_matrix <- matrix(NA, nrow = length(w_vec), ncol = length(t_vec))
    
    # Calculate the joint density f_X(w) * f_Y|X(t) over w
    for(j in 1:length(t_vec)){
      t <- t_vec[j]
      
      for(i in 1:length(w_vec)){
        w <- w_vec[i]
        
        # f(theta): true theta's density
        pdf_x <- sum(dnorm((w - dma_reduc$tilde_yi_Z)/h))/(nrow(dr)*h)
        
        # f(theta hat):  
        joint_density <- pdf_x * pdf_y_given_x(t, w) #f_X(w) * f_Y|X(t|w)
        
        density_matrix[i, j] <- joint_density
      }
      
    }
    
    # integrated density (integrate over w)
    int_dens <- colSums(density_matrix)
    
    # Normalize above to make it a proper probability density
    pdf_y <- int_dens / sum(int_dens)
    
    # Dataset for plotting:
    plot_dat <- data.frame(t = t_vec,
                           f_y = pdf_y)
    
    # Predicted yi with highest density:
    dr_pred$pred_yi_peak[k] <- t_vec[which.max(pdf_y)]
    
    
    # Find 95% prediction interval:
    ecdf_dat <- plot_dat %>% mutate(ecdf = cumsum(f_y))
    
    
    # 2.5%:
    for (i in 1:nrow(ecdf_dat)) {
      if (ecdf_dat$ecdf[i] > q1) {
        break
      }
    }
    dr_pred$pred_yi_lo[k] <- ecdf_dat$t[i - 1]
    
    
    # 97.5%:
    for (j in 1:nrow(ecdf_dat)) {
      if (ecdf_dat$ecdf[j] >= q2) {
        break
      }
    }
    dr_pred$pred_yi_hi[k] <- ecdf_dat$t[j + 1]
  }
  
  # % of prediction intervals that cover the point estimate yi: 
  dr_pred$cover_yi <- with(dr_pred, 
                           ifelse((yi - pred_yi_lo)*(yi - pred_yi_hi) <= 0, 1, 0))
  
  SAPB_interval_out = list(mean(dr_pred$cover_yi), dr_pred)
  
  return(SAPB_interval_out)
}

# Need point estimates yi in dr:
dr_pred_mrMod2_out <- SAPB_interval_func(fitted.model = mrMod2)
# To Do:
dr_pred_mrMod3_out <- SAPB_interval_func(fitted.model = mrMod3)




# Merge SAPB results for mrMod2 and mrMod3
SAPB_pred_mrMod23 <- dr_pred_mrMod2_out[[2]] %>% 
  rename(mrMod2_pred_yi_peak = pred_yi_peak,
         mrMod2_pred_yi_lo = pred_yi_lo,
         mrMod2_pred_yi_hi = pred_yi_hi,
         mrMod2_cover_yi = cover_yi) %>% 
  left_join(dr_full, by = join_by(lab_id, yi, vi))









############################################
##### Make predictions under 4 models ######
############################################

##### 1. From dmaMod (No Pub Bias) #####
# same estimate for everyone
dmaMod$b.r

crit = qnorm(.975)
yio = dmaMod$b.r
vio = dmaMod$reg_table$SE^2
t2 = dmaMod$mod_info$tau.sq
# one for each study:
totalVar = t2 + vio + dr$vi

dr_full$pred.dmaMod = dmaMod$b.r
dr_full$dmaMod_pred_yi_lo = c(yio) - c(crit) * sqrt(totalVar)
dr_full$dmaMod_pred_yi_hi = c(yio) + c(crit) * sqrt(totalVar)





##### 2. From metaCorr (Corrected for Pub Bias) #####
# same estimate for everyone
eta = 4.70 # selection ratio
metaCorr = pubbias_meta( yi = dma$yi,
                           vi = dma$vi, 
                           selection_ratio = eta, # selection ratio
                           model_type = "robust",
                           favor_positive = TRUE )
metaCorr
metaCorr$stats$estimate
# 0.2672472

# Prediction interval
crit = qnorm(.975)
yio = metaCorr$stats$estimate
vio = metaCorr$stats$se^2
t2 = dmaMod$mod_info$tau.sq # use the naive tau2
# one for each study:
totalVar = t2 + vio + dr$vi

dr_full$pred.metaCorr = metaCorr$stats$estimate
dr_full$metaCorr_pred_yi_lo = c(yio) - c(crit) * sqrt(totalVar)
dr_full$metaCorr_pred_yi_hi = c(yio) + c(crit) * sqrt(totalVar)



##### 3. From mrMod2 (No Pub Bias) #####
# design matrix with dummy variables for categorical vars columns need to match mrMod's covariate order
# sanity check: are the columns in correct order?
expect_equal( attr(Xmat, "dimnames")[[2]], mrMod2$labels )

# make predictions
# x = Xmat[1,]
# use robumeta's fn to easily get inference
# predict( object = mrMod2, pred.vector = x)
mrMod2.preds = apply( Xmat,
                      MARGIN = 1,
                      FUN = function(x) predict( object = mrMod2, pred.vector = x) )

mrMod2.preds = as.data.frame( t(mrMod2.preds) )

# sanity check: OK
# make own predictions
Bhat = as.matrix( mrMod2$b.r, ncol = 1 )
myPreds = Xmat %*% Bhat
expect_equal( as.numeric(myPreds), mrMod2.preds$prediction )

pred.mrMod2 <- mrMod2.preds$prediction

### *** Main Results 1: Point estimate prediction, and 95% prediction interval *** ###
dr_full$pred.mrMod2 = mrMod2.preds$prediction
dr_full$predLo.mrMod2 = mrMod2.preds$lower
dr_full$predHi.mrMod2 = mrMod2.preds$upper

summary(dr_full$pred.mrMod2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.03138  0.23401  0.29395  0.26401  0.31710  0.35666 
summary(dr_full$predLo.mrMod2)
summary(dr_full$predHi.mrMod2)



##### 4. From mrMod3 (Corrected for Pub Bias) #####
# OK!
expect_equal( attr(Xmat, "dimnames")[[2]], mrMod3$labels )
# Already ran before
mods
pvals = dma$p_value
# weight for model
weights = rep( 1, length(pvals) )
weights[ dma$affirm == FALSE ] = eta # if non-affirmative, set weights to selection ratio
weights

userweights_dma = weights / (dma$vi + mrMod2$mod_info$tau.sq)
length(userweights_dma)
# 138
dim(dma)
# 138 37

### STOP HERE:
mrMod3 = robu( eval( parse( text = linPred ) ),
               studynum = as.factor(study_ID),
               data = dma,
               userweights = userweights_dma, # use tau2 from mrMod2
               var.eff.size = vi,
               small = TRUE )
mrMod3


# Unweighted:
mrMod3.preds = apply( Xmat,
                      MARGIN = 1,
                      FUN = function(x) predict( object = mrMod3,
                                                 pred.vector = x,
                                                 userweights = userweights_dr) )
# mrMod3.preds = as.data.frame( t(mrMod3.preds) )
# pred.mrMod3 <- mrMod3.preds$prediction
# summary(pred.mrMod3)

Bhat = as.matrix( mrMod3$b.r, ncol = 1:18 )
dim(Bhat)
dim(Xmat)
# Calculate by hand:
pred.mrMod3 = Xmat %*% Bhat
# expect_equal( as.numeric(myPreds), mrMod3.preds$prediction )

### *** Main Results 2: Point estimate prediction, and 95% CI *** ###
dr_full$pred.mrMod3 = pred.mrMod3
# dr$predLo.mrMod3 = mrMod3.preds$lower
# dr$predHi.mrMod3 = mrMod3.preds$upper





































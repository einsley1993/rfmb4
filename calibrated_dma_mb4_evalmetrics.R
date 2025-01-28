
####################################################
####### Put all predictions into the dataset #######
####################################################
# All point estimates here:
dr_pred_est <- dr_full %>% 
  select(lab_id, yi, vi,
         pred.dmaMod, 
         pred.metaCorr,
         pred.mrMod2,
         pred.mrMod3
         ) 


# All prediction intervals here:
dr_pred_interval <- dr_full %>% 
  select(lab_id, yi, vi, 
         dmaMod_pred_yi_lo, dmaMod_pred_yi_hi,
         metaCorr_pred_yi_lo, metaCorr_pred_yi_hi,
         predLo.mrMod2, predHi.mrMod2) %>% 
  mutate(dmaMod_pred_yi_lo = dmaMod_pred_yi_lo,
         dmaMod_pred_yi_hi = dmaMod_pred_yi_hi,
         metaCorr_pred_yi_lo = metaCorr_pred_yi_lo,
         metaCorr_pred_yi_hi = metaCorr_pred_yi_hi,
         mrMod2_pred_yi_lo = predLo.mrMod2,
         mrMod2_pred_yi_hi = predHi.mrMod2
         # mrMod3:
         ) 




#### Outcome metric #1: Raw difference between MA and REP point estimates (just subtraction)
## (1) dmaMod: not condition on moderators
# dr
dmaMod_dr = robu( yi ~ 1,
                  data = dr_full,
                  studynum = as.factor(lab_id),
                  var.eff.size = vi,
                  modelweights = "CORR",
                  small = TRUE)
dmaMod_dr$reg_table$b.r
# 0.03106345

# dma
dmaMod$b.r
# 0.3331773

dmaMod_dr$reg_table$b.r - dmaMod$b.r
# -0.3021138



## (2) mdMod2: condition on moderators
# Can only include variabels with >= 2 levels:
mrMod2_dr = robu( yi ~ mean_age_1 + n_excluded_1 + Hamlin_Lab,
                  data = dr_full, 
                  studynum = as.factor(lab_id),
                  var.eff.size = vi,
                  modelweights = "CORR",
                  small = TRUE)
mrMod2_dr


# dma
mrMod2_dma = robu( yi ~ mean_age_1 + n_excluded_1 +
                     method + scenarios + stimuli + stimuli_type + 
                     choice_object + Hamlin_Lab + intent_val + outcome_val,
                   data = dma, 
                   studynum = as.factor(study_ID),
                   var.eff.size = vi,
                   modelweights = "HIER",
                   small = TRUE)
mrMod2_dma


# Condition on moderators' mean levels
mrMod2_dr_mean = predict(mrMod2_dr, 
                         as.matrix(data.frame(X.Intercept. = 1,
                                              mean_age_1 = mean(dr_full$mean_age_1),
                                              n_excluded_1 = mean(dr_full$n_excluded_1, na.rm = TRUE),
                                              Hamlin_LabTRUE = 0
                                              )) )
mrMod2_dr_mean

mrMod2_dma_mean = predict(mrMod2_dma, 
                         as.matrix(data.frame(X.Intercept. = 1,
                                              mean_age_1 = mean(dma$mean_age_1),
                                              n_excluded_1 = mean(dma$n_excluded_1,na.rm = TRUE),
                                              methodreach = 0,
                                              scenariosgive_take = 0,
                                              scenarioshelp_hinder = 1,
                                              scenariosother = 0,
                                              stimulimovies = 0,
                                              stimuli_typereal = 1,
                                              choice_objectpeople = 0,
                                              choice_objectpuppets = 1,
                                              choice_objectshapes = 0,
                                              Hamlin_LabTRUE = 0,
                                              intent_valopp = 1,
                                              intent_valpos = 0,
                                              outcome_valopp = 1,
                                              outcome_valrev = 0,
                                              outcome_valsame = 0)) )
mrMod2_dma_mean

mrMod2_dr_mean[1] - mrMod2_dma_mean[1]
# -0.1297366 



## (3) metaCorr: not condition on moderators
# dr: Won't correct for publication bias. eta = 1.
metaCorr_dr = corrected_meta( yi = dr_full$yi,
                              vi = dr_full$vi, 
                              eta = 1, 
                              cluster = dr_full$lab_id, 
                              model = "robust",
                              favor.positive = TRUE)
metaCorr_dr$stats$estimate
# 0.03108822

# dma
metaCorr_dma = corrected_meta( yi = dma$yi,
                               vi = dma$vi, 
                               eta = 4.70, # selection ratio
                               cluster = dma$study_ID, 
                               model = "robust",
                               favor.positive = TRUE)
metaCorr_dma$stats$estimate
# 0.2672472

metaCorr_dr$stats$estimate - metaCorr_dma$stats$estimate
# -0.236159




## (4) mrMod3: condition on moderators
meta.re_dr = rma.uni( yi = dr_full$yi, vi = dr_full$vi)
t2hat.naive_dr = meta.re_dr$tau2
t2hat.naive_dr

mrMod3_dr = robu( yi ~ mean_age_1 + n_excluded_1 + Hamlin_Lab,
                  data = dr_full, 
                  studynum = as.factor(lab_id),
                  # userweights = weights, 
                  var.eff.size = vi,
                  modelweights = "CORR",
                  small = TRUE)
mrMod3_dr



# dma
meta.re_dma = rma.uni( yi = dma$yi, vi = dma$vi)
meta.re_dma
t2hat.naive_dma = meta.re_dma$tau2
t2hat.naive_dma

# weight for model
weights = rep( 1, length(dma$p_value) )
weights[ dma$affirm == FALSE ] = 4.70 # if non-affirmative, set weights to selection ratio
weights

mrMod3_dma = robu( yi ~ mean_age_1 + n_excluded_1 +
                     method + scenarios + stimuli + stimuli_type + 
                     choice_object + Hamlin_Lab + intent_val + outcome_val,
                   data = dma, 
                   studynum = as.factor(study_ID),
                   userweights = weights / (vi + t2hat.naive_dma), # Or use tau2 from mrMod2?
                   var.eff.size = vi,
                   modelweights = "HIER",
                   small = TRUE)
mrMod3_dma


# Condition on moderators' mean levels
mrMod3_dr_mean = predict(mrMod3_dr, 
                         as.matrix(data.frame(X.Intercept. = 1,
                                              mean_age_1 = mean(dr_full$mean_age_1),
                                              n_excluded_1 = mean(dr_full$n_excluded_1),
                                              Hamlin_LabTRUE = 0
                                              )) )
mrMod3_dr_mean

mrMod3_dma$reg_table$labels

# Gives error because of userweights: ask Zach Fisher
mrMod3_dma_mean = predict(mrMod3_dma, 
                          as.matrix(data.frame(X.Intercept. = 1,
                                               mean_age_1 = mean(dma$mean_age_1),
                                               n_excluded_1 = mean(dma$n_excluded_1, na.rm = TRUE),
                                               methodreach = 0,
                                               scenariosgive_take = 0,
                                               scenarioshelp_hinder = 1,
                                               scenariosother = 0,
                                               stimulimovies = 0,
                                               stimuli_typereal = 1,
                                               choice_objectpeople = 0,
                                               choice_objectpuppets = 1,
                                               choice_objectshapes = 0,
                                               Hamlin_LabTRUE = 0,
                                               intent_valopp = 1,
                                               intent_valpos = 0,
                                               outcome_valopp = 1,
                                               outcome_valrev = 0,
                                               outcome_valsame = 0)) )
mrMod3_dma_mean

mrMod3_dr_mean[1] - mrMod3_dma_mean[1]
#







#### Outcome metric #2: p_orig from Maya’s JRSSA paper (statistical inconsistency between MA and REP)
#### P_orig: the probability that the effect estimate reported in the REP study would be as extreme or more extreme
#### than it actually was if, in fact, the REP study and the MA study were statistically consistent in the sense of being
#### drawn from the same distribution
## (1) dmaMod:
dma_reduc_ls <- list()
p_orig_dmaMod <- rep(NA, nrow(dr_pred))
for(k in 1:nrow(dr_pred)){
  preds_Z_mod <- dmaMod$b.r 
  tau2 <- dmaMod$mod_info$tau.sq
  
  dma_reduc_ls[[k]] <- dma_reduc %>% 
    mutate(tilde_yi_Z = preds_Z_mod + 
             sqrt(tau2/(tau2 + vi)) * (yi - pred_Zi))
  
  mu <- mean(dma_reduc_ls[[k]]$tilde_yi_Z)
  var_mu <- var(dma_reduc_ls[[k]]$tilde_yi_Z)
  
  
  p_orig_dmaMod[k] <- 2*(1 - pnorm(abs(dr_full$yi[k] - mu)/sqrt(tau2 + var_mu + dr_full$vi[k])))
  
}

hist(p_orig_dmaMod)
summary(p_orig_dmaMod)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000383 0.1338332 0.2891191 0.3703269 0.5563114 0.9618430 


## (2) mrMod2:
dma_reduc_ls <- list()
p_orig_mrMod2 <- rep(NA, nrow(dr_pred))
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
                                      Hamlin_LabTRUE = as.numeric(dr_full$Hamlin == TRUE)[k],
                                      intent_valopp = 1,
                                      intent_valpos = 0,
                                      outcome_valopp = 1,
                                      outcome_valrev = 0,
                                      outcome_valsame = 0) 
  )
  
  preds_Z_mean_mod <- predict(mrMod2, Z_mod)[1]
  tau2 <- mrMod2$mod_info$tau.sq
  
  dma_reduc_ls[[k]] <- dma_reduc %>% 
    mutate(tilde_yi_Z = preds_Z_mod + 
             sqrt(tau2/(tau2 + vi)) * (yi - pred_Zi))
  
  # Maya's JRSSA paper, Eq (4.1)
  # theta_orig = REP's point estimate (yi), SE_orig = sqrt(vi)
  
  # mu = mean(calibrated estimates from MA), SE_mu = SD(calibrated estimates from MA)
  mu <- mean(dma_reduc_ls[[k]]$tilde_yi_Z)
  var_mu <- var(dma_reduc_ls[[k]]$tilde_yi_Z)
  
  
  p_orig_mrMod2[k] <- 2*(1 - pnorm(abs(dr_full$yi[k] - mu)/sqrt(tau2 + var_mu + dr_full$vi[k])))
  
  summary_output <- p_orig_mrMod2
}

hist(p_orig_mrMod2)
summary(p_orig_mrMod2)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0000181 0.1271733 0.2647787 0.3552030 0.5453965 0.9696490 


## (3) metaCorr:
dma_reduc_ls <- list()
p_orig_metaCorr <- rep(NA, nrow(dr_pred))
for(k in 1:nrow(dr_pred)){
  preds_Z_mod <- metaCorr$stats$estimate
  tau2 <- t2hat.naive
  
  dma_reduc_ls[[k]] <- dma_reduc %>% 
    mutate(tilde_yi_Z = preds_Z_mean_mod + 
             sqrt(tau2/(tau2 + vi)) * (yi - pred_Zi))
  
  mu <- mean(dma_reduc_ls[[k]]$tilde_yi_Z)
  var_mu <- var(dma_reduc_ls[[k]]$tilde_yi_Z)
  
  
  p_orig_metaCorr[k] <- 2*(1 - pnorm(abs(dr_full$yi[k] - mu)/sqrt(tau2 + var_mu + dr_full$vi[k])))
  
  
}

hist(p_orig_metaCorr)
summary(p_orig_metaCorr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0001526 0.0948651 0.4090636 0.3798962 0.5266054 0.8571484 



## (4) mrMod3: error, ask Zach Fisher
dma_reduc_ls <- list()
p_orig_mrMod3 <- rep(NA, nrow(dr_pred))
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
                                 Hamlin_LabTRUE = 0,
                                 intent_valopp = 1,
                                 intent_valpos = 0,
                                 outcome_valopp = 1,
                                 outcome_valrev = 0,
                                 outcome_valsame = 0)  ) 
  
  preds_Z_mean_mod <- predict(mrMod3, Z_mod)[1]
  
  # Use the tau2 from mrMod2 ?
  tau2 <- mrMod2$mod_info$tau.sq
  
  dma_reduc_ls[[k]] <- dma_reduc %>% 
    mutate(tilde_yi_Z = preds_Z_mod + 
             sqrt(tau2/(tau2 + vi)) * (yi - pred_Zi))

  mu <- mean(dma_reduc_ls[[k]]$tilde_yi_Z)
  var_mu <- var(dma_reduc_ls[[k]]$tilde_yi_Z)
  
  
  p_orig_mrMod3[k] <- 2*(1 - pnorm(abs(dr_full$yi[k] - mu)/sqrt(tau2 + var_mu + dr_full$vi[k])))
  
  
}

hist(p_orig_mrMod3)
summary(p_orig_mrMod3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 




dr_pred_est_porig <- dr_pred_est %>% mutate(p_orig_dmaMod = p_orig_dmaMod,
                                            p_orig_mrMod2 = p_orig_mrMod2,
                                            p_orig_metaCorr = p_orig_metaCorr
                                            # Will update:
                                            # p_orig_mrMod3 = p_orig_mrMod3
                                            )


write.csv(dr_pred_est_porig, "dr_pred_est_porig.csv")

apply(dr_pred_est_porig[,8:10], 2, mean)



#### Outcome metric #3: MAE
dr_pred_est %>% summarize(MAE_pred.dmaMod = mean(abs(yi - pred.dmaMod)), 
                          MAE_pred.mrMod2 = mean(abs(yi - pred.mrMod2)),
                          MAE_pred.metaCorr = mean(abs(yi - pred.metaCorr)), 
                          MAE_pred.mrMod3 = mean(abs(yi - pred.mrMod3)))

#### Outcome metric #4: RMSE
dr_pred_est %>% summarize(RMSE_pred.dmaMod = sqrt(mean((yi - pred.dmaMod)^2)), 
                          RMSE_pred.mrMod2 = sqrt(mean((yi - pred.mrMod2)^2)),
                          RMSE_pred.metaCorr = sqrt(mean((yi - pred.metaCorr)^2)), 
                          RMSE_pred.mrMod3 = sqrt(mean((yi - pred.mrMod3)^2)))

#### Outcome metric #5: % Overestimate
dr_pred_est %>% summarize(percOverest_pred.dmaMod = mean(yi < pred.dmaMod), 
                          percOverest_pred.mrMod2 = mean(yi < pred.mrMod2),
                          percOverest_pred.metaCorr = mean(yi < pred.metaCorr), 
                          percOverest_pred.mrMod3 = mean(yi < pred.mrMod3))




#### Outcome metric #6: Percentage of REP studies whose estimates were inside their corresponding prediction intervals
dr_pred_interval %>% summarize(percCover_yi_pred.dmaMod = mean(ifelse((yi - dmaMod_pred_yi_hi)*(yi - dmaMod_pred_yi_lo) <= 0, 1, 0)), 
                               percCover_yi_pred.mrMod2 = mean(ifelse((yi - mrMod2_pred_yi_hi)*(yi - mrMod2_pred_yi_lo) <= 0, 1, 0)),
                               percCover_yi_pred.metaCorr = mean(ifelse((yi - metaCorr_pred_yi_hi)*(yi - metaCorr_pred_yi_lo) <= 0, 1, 0))
                               # Will update:
                               # percCover_yi_pred.mrMod3 = mean(ifelse((yi - mrMod3_pred_yi_hi)*(yi - mrMod3_pred_yi_lo) <= 0, 1, 0))
                               )

# percCover_yi_pred.dmaMod percCover_yi_pred.mrMod2 percCover_yi_pred.metaCorr
#               0.8484848                0.6363636                  0.8484848


# Not reported: width of prediction interval
dr_pred_interval %>% summarize(width_yi_pred.dmaMod = mean(dmaMod_pred_yi_hi - dmaMod_pred_yi_lo), 
                               width_yi_pred.mrMod2 = mean(mrMod2_pred_yi_hi - mrMod2_pred_yi_lo),
                               width_yi_pred.metaCorr = mean(metaCorr_pred_yi_hi - metaCorr_pred_yi_lo)
                               # Will update:
                               # width_yi_pred.mrMod3 = mean(mrMod3_pred_yi_hi - mrMod3_pred_yi_lo)
                               )
# width_yi_pred.dmaMod width_yi_pred.mrMod2 width_yi_pred.metaCorr
# 1             2.684151             1.432742               2.682219





########### Plot the intervals ###########
dp_pred_est <- dr_pred_est %>% arrange( desc(yi) )
dp_pred_interval <- dr_pred_interval %>% arrange( desc(yi) )


ggplot( ) + 
  # raw yi in dr:
  geom_point( data = dp_pred_est, 
              aes( x = yi,
                   y = 1:nrow(dp_pred_est) ),
              color = "blue", size = 1 ) +
  
  # predicted yi with the highest density:
  # geom_point( data = dp_pred_interval, 
  #             aes( x = mrMod2_pred_yi_peak,
  #                  y = 1:nrow(dp_pred_est) ),
  #             color = "pink", size = 2.5, shape = 18 ) +
  
  # predicted point estimate from mrMod2:
  geom_point( data = dp_pred_est, 
              aes( x = pred.mrMod2,
                   y = 1:nrow(dp_pred_est) ),
              color = "purple", shape = 17, size = 1.5 ) +
  
  # Replication mean:
  # geom_vline( xintercept = drMod$b.r,
  #             color = "black" ) +
  
  # Naive meta-analysis estimate:
  geom_vline( xintercept = dmaMod$b.r,
              color = "gray" ) +
  
  # Corrected meta-analysis estimate:
  geom_vline( xintercept = metaCorr$stats$estimate,
              color = "dark green" ) +
  
  # Prediction intervals (conditional on Z_i):
  geom_errorbarh( data = dp_pred_interval,
                  aes( xmax = mrMod2_pred_yi_hi, 
                       xmin = mrMod2_pred_yi_lo,
                       y = 1:nrow(dp_pred_est)),
                  color = "purple")  +
  
  # Prediction intervals (from dmaMod):
  geom_errorbarh( data = dp_pred_interval,
                  aes( xmax = dmaMod_pred_yi_hi, 
                       xmin = dmaMod_pred_yi_lo,
                       y = 1:nrow(dp_pred_est)),
                  color = "gray")  +
  
  theme_bw() + 
  theme(axis.title = element_blank(),
        plot.subtitle = element_text(vjust = -1 )) + 
  labs(title = "Predicted point estimates and intervals",
       subtitle = expression(
         paste(
           "Blue dots: raw yi.",
           "\n Purple dots: predicted point estimates from mrMod2.",
           # "\n Black line: drMod (unconditional) estimate.",
           "\n Gray line: dmaMod estimate.",
           "\n Green line: metaCorr estimate.",
           "\n Purple bar: intervals from mrMod2.",
           "\n Gray bar: intervals from dmaMod.",
           sep = " "
         )
       ))
  
getwd()
ggsave("Prediction intervals_MR_dr.pdf", height = 13, width = 16)  

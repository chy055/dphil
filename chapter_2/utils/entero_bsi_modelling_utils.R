
###############################################################################

# functions used in baseline_models_enterobacterales.Rmd

## operation functions
winsorise_covar <- function(x, lowerlim = 0.025, upperlim = 0.975){
  quantiles <- quantile( x, c(lowerlim, upperlim) , na.rm = TRUE)
  x[ x < quantiles[1] ] <- quantiles[1]
  x[ x > quantiles[2] ] <- quantiles[2]
  x
}


# factorise categorical variables:
factorise_categorical <- function(dataframe, factorise_all = TRUE, list_of_categorical_variables, factorise_multilevel = TRUE){
  
  df <- dataframe
  
  if (factorise_all == TRUE){
    
    df[list_of_categorical_variables] <- lapply(df[list_of_categorical_variables], factor)
    
  }
  
  if (factorise_multilevel == TRUE){
    
    # set multi-level categorical reference levels
    df$AdmissionDayOfWeek <- factor(df$AdmissionDayOfWeek, levels = c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'))
    df$AdmissionTimeOfDay <- factor(df$AdmissionTimeOfDay, levels = c('morning', 'afternoon', 'evening', 'overnight'))
    df$AVPU <- factor(df$AVPU, levels = c('Alert', 'Verbal', 'Pain / Unresponsive'))
    df$Ethnicity <- factor(df$Ethnicity, levels = c('white', 'other', 'unrecorded'))
    df$Source <- factor(df$Source, levels = c('urinary', 'chest', 'abdominal', 'other', 'unknown source', 'multiple potential sources'))
    df$Specialty <- factor(df$Specialty, levels = c('Acute and general medicine', 'Medical subspecialty', 'Acute and general surgery', 'Other'))
    df$Active_Partial_Inactive <- factor(df$Active_Partial_Inactive, levels = c('Active baseline', 'Partial baseline activity', 'Inactive baseline'))
    
    return(df)
    
  }else{
    
    return(df)
    
  }
}


###############################################################################

# make df for survival models
make_survival_df <- function(pre_analysis_df, winsorise_continuous = FALSE, sym_wins_vars, sym_limits = c(0.025, 0.975), asym_upper_vars, asym_upper_limits = c(0, 0.95), asym_lower_vars, asym_lower_limits = c(0.05, 1)){
  
  edited_df <- pre_analysis_df %>%
    select(ClusterID, NewCulturePeriodNumber,
           
           # categorical variables:
           LinkedSex, Ethnicity, Specialty, AdmissionTimeOfDay, AdmissionDayOfWeek, OutOfHoursAdmission, PriorHospitalisation, AVPU, SupplementaryOxygen, diabetes, dialysis, immunosuppression, palliative, CommunityOnset, PriorBacteraemiaEpisode, Polymicrobial,
           ESBL, Ecoli, OtherESCAPPM, NonEnterobacterales, Other, Enterobacter, Klebsiella,
           PriorAmoxicillin, PriorCoamoxiclav, PriorAnyAntibiotic, PriorBetaLactam, GivenMultiAntibioticBaseline, GivenBroadSpectrum, ActiveCoamoxiclav, ActiveNonCoamox, ActiveAminoglycoside, ActiveOther, ActivePreCultureAntibiotics, ActivePreCultureAminoglycoside, GivenTrimethoOrNitro, 
           InactiveCoamoxiclav, EUCASTbreakpoints, GivenBaselineAminoglycoside,
           InactiveBaseline, PartialBaselineActivity,
           Active_Partial_Inactive,
           # Source
           urinary_source, chest_source, abdominal_source, other_source, surgical_source, Multi_source,
           UnadjustedPostIndexPeriod,
           
           
           # continuous variables:
           Weight, Height, BMI, SBP, DBP, HeartRate, RespiratoryRate, OxygenSaturation, Temperature, albumin, alkaline_phosphatase, alt, aptt, ast, bilirubin, calcium_adjusted, creatinine, egfr, eosinophils, haematocrit, haemoglobin, immature_granulocytes, 
           lymphocytes, magnesium, mch, mcv, monocytes, neutrophils, phosphate, platelets, potassium, pt, sodium, urea, white_cells, crp,
           AgeAtAdmission, PriorAdmissionsTotalDuration, IMDPercentile, Charlson_sum, Elixhauser_sum,
           UnadjustedBacteraemiaPeriodDays,
           
           # survival data
           OverallSurvival_days_adjusted, OverallSurvival_days_study_end, CensorStatus_day30, CensorStatus_day365, CensorStatus_study_end, LengthOfStayFromCulture) %>%
    
    mutate(OverallSurvival_days_capped = case_when(OverallSurvival_days_adjusted >= 30 ~ 30,
                                                   TRUE ~ OverallSurvival_days_adjusted))
  
  if(winsorise_continuous == TRUE){
    
    if(!missing(sym_wins_vars) & !missing(asym_upper_vars)){
      
      if(!missing(asym_lower_vars)){
        
        winsorised_df <- edited_df %>%
          mutate_at(sym_wins_vars, ~winsorise_covar(., lowerlim = sym_limits[1]), upperlim = sym_limits[2]) %>%
          mutate_at(asym_upper_vars, ~winsorise_covar(., lowerlim = asym_upper_limits[1]), upperlim = asym_upper_limits[2]) %>%
          mutate_at(asym_lower_vars, ~winsorise_covar(., lowerlim = asym_lower_limits[1]), upperlim = asym_lower_limits[2])
        
        return(winsorised_df)
        
      }else{
        
        winsorised_df <- edited_df %>%
          mutate_at(sym_wins_vars, ~winsorise_covar(., lowerlim = sym_limits[1]), upperlim = sym_limits[2]) %>%
          mutate_at(asym_upper_vars, ~winsorise_covar(., lowerlim = asym_upper_limits[1]), upperlim = asym_upper_limits[2])
        
        return(winsorised_df)
        
      }
      
      
      
    }else{
      
      cat('\nError: missing asymmetrical and/or symmetrical lists, even if empty.\n')
      
    }
    
  }else{
    
    return(edited_df)
    
  }
}



# forward selection using multiple fractional polynomials (requires mfp package)
add_MFP_term <- function(mfp_fit_object, data_df, categorical_variables_list, continuous_covariates_list, not_to_scale_variables_list){ 
  
  base_formula <- mfp_fit_object$fit$formula
  
  pre_backelim_covariates_list <- c(continuous_covariates_list, categorical_variables_list)
  
  backselection_covar <- as_tibble(mfp_fit_object$coefficients, rownames = "id") %>%
    mutate(new_id = str_extract(id, paste0(pre_backelim_covariates_list, collapse = '|'))) %>%
    .$new_id %>% unique(.)
  
  dropped_covariates_list <- setdiff(pre_backelim_covariates_list, backselection_covar)
  dropped_covariates_list <- dropped_covariates_list[dropped_covariates_list %notin% c('Group', 'AVPU', 'AdmissionDayOfWeek')]
  
  print(dropped_covariates_list)
  
  list_of_df_by_covar <- list()
  add_results <- list()
  for(i in seq_along(dropped_covariates_list)){
    
    if(dropped_covariates_list[[i]] %in% continuous_covariates_list){
      
      if(dropped_covariates_list[[i]] %in% not_to_scale_variables_list){
        
        covar_to_add <- paste0('fp(', dropped_covariates_list[[i]],', df = 4, scale=FALSE, select = 1,  alpha = 0.05)')
        
      }else{
        
        covar_to_add <- paste0('fp(', dropped_covariates_list[[i]],', df = 4, scale=TRUE, select = 1,  alpha = 0.05)')  
        
      }
      
    }else{
      
      covar_to_add <- dropped_covariates_list[[i]]
      
    }
    
    pre_surv_df <- data_df[, c(c(backselection_covar, dropped_covariates_list[[i]]), 'OverallSurvival_days_capped', 'CensorStatus_day30')]
    # sanity check that subset worked for right number of covariates:
    surv_df <- pre_surv_df[complete.cases(pre_surv_df), ]
    list_of_df_by_covar[[dropped_covariates_list[[i]]]] <- surv_df
    cat(paste("\nafter adding", dropped_covariates_list[[i]], "no. of complete cases:", nrow(surv_df)))
    
    pre_add.formula <- as.formula(paste0('Surv(OverallSurvival_days_capped, CensorStatus_day30) ~ (', as.character(base_formula[3]), ')'))
    add.formula <- as.formula(paste0('Surv(OverallSurvival_days_capped, CensorStatus_day30) ~ (', as.character(base_formula[3]), " + ", covar_to_add,  ')'))
    
    add_mfp_cox <- mfp(add.formula, family = cox, data = surv_df,
                       select = 1, method = 'efron',
                       keep = c("Group"),       #maintaining the exposure group given its significance for this study
                       verbose = FALSE, maxits = 20)
    add_summary <- summary(add_mfp_cox)
    
    add_tibble <- as_tibble(add_summary$coefficients, rownames = 'transformed_covariate')
    
    add_tibble <- add_tibble[grep(dropped_covariates_list[[i]], add_tibble$transformed_covariate), ]
    
    # calculate log likelihood test
    add_loglik <- signif(anova(coxph(pre_add.formula, data=surv_df), coxph(add_mfp_cox$formula, data = surv_df))$"Pr(>|Chi|)"[2],3)
    add_tibble$p_loglik <- add_loglik
    
    # calculate AIC
    add_aic <- extractAIC(add_mfp_cox)[[2]]
    add_tibble$AIC <- add_aic
    
    add_results[[i]] <- add_tibble
    
  }
  
  add_results_df <- do.call('rbind', add_results)%>%
    arrange(AIC, 'Pr(>|z|)')
  
  base_formula_AIC <- extractAIC(coxph(base_formula, data=surv_df))[[2]]
  
  add_results_df <- rbind(add_results_df,c('<<none>>', 0, 0, 0, 0, 1, 1, base_formula_AIC)) %>%
    mutate(AIC = as.numeric(AIC)) %>%
    arrange(AIC)
  
  
  
  add_results_df <- add_results_df %>%
    rename(HR = 'exp(coef)',
           se = 'se(coef)',
           p_wald = 'Pr(>|z|)') %>%
    mutate(HR = as.numeric(HR),
           se = as.numeric(se),
           p_wald = as.numeric(p_wald)) %>%
    mutate(upper_95 = exp(log(HR) + 1.96*se),
           lower_95 = exp(log(HR) - 1.95*se))
  
  is.num <- sapply(add_results_df, is.numeric)
  add_results_df[is.num] <- lapply(add_results_df[is.num], round, 3)
  
  add_results_df <- add_results_df %>%
    mutate("95% CI" = paste0("(", lower_95, "-", upper_95, ")")) %>%
    mutate(covariate = str_extract(transformed_covariate, paste0(pre_backelim_covariates_list, collapse = '|'))) %>%
    select(covariate, transformed_covariate, HR, "95% CI", p_wald, p_loglik, AIC) %>%
    arrange(AIC, p_wald)
  
  
  
  return(list(add_results_df,list_of_df_by_covar))
  
}


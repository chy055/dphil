###############################################################################
# utilities and functions for cluster_sepsis.Rmd

# make cluster df with only Seymour et al's 29 variables
make_survival_df_sepsis <- function(pre_analysis_df, winsorise_continuous = FALSE, sym_wins_vars, sym_limits = c(0.025, 0.975), asym_upper_vars, asym_upper_limits = c(0, 0.95), asym_lower_vars, asym_lower_limits = c(0.05, 1)){
  
  edited_df <- pre_analysis_df %>%
    select(ClusterID, SpellID,
           
           # categorical variables:
           LinkedSex, CulturePositive,
           
           # continuous variables:
           AgeAtAdmission, albumin, alt, ast, bicarbonate, bilirubin, urea, chloride, crp, creatinine,
           Elixhauser, esr, GCS, glucose, HeartRate, haemoglobin, immature_granulocytes, inr, lactate, OxygenSaturation, po2,
           platelets, RespiratoryRate, sodium, SBP, Temperature, troponin, white_cells,
           
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



# Make cluster df with Seymour and all extra IORD variables
make_survival_df_sepsis_plus <- function(pre_analysis_df, winsorise_continuous = FALSE, sym_wins_vars, sym_limits = c(0.025, 0.975), asym_upper_vars, asym_upper_limits = c(0, 0.95), asym_lower_vars, asym_lower_limits = c(0.05, 1)){
  
  edited_df <- pre_analysis_df %>%
    select(ClusterID, SpellID,
           
           # categorical variables:
           LinkedSex, Ethnicity, SpecialtyGroup, AdmissionTimeOfDay, AdmissionDayOfWeek, OutOfHoursAdmission, PriorHospitalisation, SupplementaryOxygen, AVPU, diabetes, dialysis, immunosuppression, palliative, 
           # Source
           urinary_source, respiratory_source, abdominal_source, other_source, cns_source, skin_soft_tissue_orthopedic_source, unspecific_source, unknown_source,
           
           # continuous variables:
           AgeAtAdmission, Weight, Height, BMI, 
           albumin, alkaline_phosphatase, alt, aptt, ast, bicarbonate, bilirubin, urea, calcium_adjusted, Charlson, chloride, crp, creatinine,
           DBP, Elixhauser, eosinophils, esr, GCS, glucose, haematocrit, HeartRate, haemoglobin, IMDPercentile, immature_granulocytes, inr, lactate, lymphocytes, magnesium, mch, mcv, monocytes, neutrophils, 
           OxygenSaturation, phosphate, po2, platelets, potassium, pt, RespiratoryRate, sodium, SBP, Temperature, troponin, white_cells,
           preAdmissionDays,
           
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



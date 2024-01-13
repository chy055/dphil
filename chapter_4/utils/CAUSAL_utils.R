###############################################################################


fisher.test.simulate.p.values <- function(data, variable, by, ...) {
  result <- list()
  test_results <- stats::fisher.test(data[[variable]], data[[by]], simulate.p.value = TRUE)
  result$p <- test_results$p.value
  result$test <- test_results$method
  result
}

my_negbin_tidy <- function(x, exponentiate = FALSE, ...) {
  
  # custom tidying function for gtsummary::tbl_regression(model, exponentiate=T, tidy_fun my_negbin_tidy)
  # enables log(OR) and CI to be exponentiated without throwing an error
  
  df_tidy <- broom::tidy(x, ...)
  
  # exponentiate coef and CI if requested
  if (exponentiate) {
    df_tidy <-
      df_tidy %>%
      mutate_at(vars(any_of(c("estimate", "conf.low", "conf.high"))), exp)
  }
  
  df_tidy
}



make_univ_formula <- function(selected_variable, linear_or_polynomial, max_permitted_degree = 2, scale_variable = FALSE){
  
  variable_str <- as.character(selected_variable)
  df_str <- as.character(max_permitted_degree)
  
  if(linear_or_polynomial == 'linear') {
    
    univ_age_cox <- as.formula(paste('Surv(OS_14, CensorStatus_day14)~', variable_str))
    
  }else if(linear_or_polynomial == 'polynomial' & scale_variable == FALSE){
    
    univ_age_cox <- as.formula(paste0('Surv(OS_14, CensorStatus_day14) ~ fp(', variable_str,', df = ', df_str, ', scale=FALSE)'))
    
  }else if(linear_or_polynomial == 'polynomial' & scale_variable == TRUE){
    
    univ_age_cox <- as.formula(paste0('Surv(OS_14, CensorStatus_day14) ~ fp(', variable_str,', df = ', df_str, ', scale=TRUE)'))
    
  }
  
  return(univ_age_cox)
  
}


categorical_univariate_fit <- function(selected_variable, survival_dataframe){
  variable_str <- as.character(selected_variable)
  
  pre_surv_df <<- survival_dataframe
  
  pre_surv_df <- pre_surv_df[, c(variable_str, 'OS_14', 'CensorStatus_day14')]
  # sanity check that subset worked for right number of covariates:
  surv_df <- pre_surv_df[complete.cases(pre_surv_df), ]
  
  cat(paste(nrow(surv_df),'\n'))
  
  linear_formula <- make_univ_formula(selected_variable = selected_variable,
                                      linear_or_polynomial = 'linear')
  
  fit <- mfp(formula = linear_formula,
             family=cox, method = 'efron', alpha = 0.05,
             data = surv_df)
  
  return(fit)
  
}



compute_iptw <- function(df, num_model, denom_model){
  
  # takes DF of the treatment model
  # uses treat_up as the treatment variable
  # output: df containing original columns & iptw-related ones
  
  output <- df
  output$propensity_num = predict(num_model, df, type='response')
  output$propensity_denom = predict(denom_model, df, type='response')
  
  output = output %>%
    
    # Probability of observed outcome
    mutate(propensity_num_outcome = ifelse(treat_up == 1, propensity_num, 1 - propensity_num),
           propensity_denom_outcome = ifelse(treat_up == 1, propensity_denom, 1 - propensity_denom)) %>% 
    # Numerator / denominator
    mutate(weights_no_time = propensity_num_outcome / propensity_denom_outcome) %>% 
    group_by(id) %>% 
    mutate(ipw_manual = cumprod(weights_no_time)) %>% 
    mutate(ipw_manual = round_half_up(ipw_manual, digits = 6)) %>%
    ungroup() 
  
  return(output)
  
}


fit_splines <- function(df, continuous_vars, knots = c(0.1, 0.5, 0.9)){
  
  var_quants = list()
  for (i in seq_along(continuous_vars)){
    var = continuous_vars[[i]]
    var_str = as.character(var)
    
    qs = unname(quantile(df[[var]], probs = knots))
    
    var_quants[[var_str]] <- qs
    
  }
  
  return(var_quants)
  
}



compute_iptw_means <- function(iptw_df){
  
  # takes DF including IPTWs (e.g., from compute_iptw function above)
  # output: df containing summary iptw_means by time interval
  
  output_means = iptw_df %>% group_by(time) %>%
    mutate(mean = mean(ipw_manual, na.rm=TRUE)) %>%
    ungroup() %>%
    distinct(time, .keep_all=TRUE)
  
  return(output_means)
  
}



plot_iptws <- function(iptw_df, iptw_means_df){
  
  # takes DF incl. IPTWs, and IPTW means computed using compute_iptw_means
  # output: plot of IPTWs by time
  
  p = iptw_df %>%
    ggplot(aes(x=as.factor(time), y=ipw_manual)) + 
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(color='black', size=0.1, alpha=0.9) +
    geom_point(data = iptw_means_df, mapping = aes(x = as.factor(time), y=mean), shape = 4, color = 'red') +
    geom_hline(yintercept = 1, color = 'red', alpha = 0.25)
  
  return(p)
  
}


robust_se <- function(model_obj){
  
  # takes a model obj
  # output: robust standard errors using Stata's default HC1 - heteroskedasticity-consistent
  output = model_obj %>%
    coeftest(., vcov = vcovCL, type = 'HC1', cluster = ~id) %>%
    tidy(., conf.int = TRUE) %>%
    mutate(OR = exp(estimate),
           ll = exp(conf.low),
           ul = exp(conf.high)) %>%
    dplyr::select(term, OR, ll, ul, everything()) %>%
    filter(!grepl("time",term)) %>%
    filter(!grepl("Intercept",term)) %>%
    dplyr::select(term, OR, ll, ul, p.value) %>%
    mutate(across(.cols = c(OR, ll, ul, p.value), ~ round_half_up(.x, digits=3)))
  
  return(output)
  
}



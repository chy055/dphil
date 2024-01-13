###############################################################################
# utilities and functions for cluster_enterobacterales.Rmd

# rounding numbers to specified number of digits
round.off <- function (x, digits=0) 
{
  posneg = sign(x)
  z = trunc(abs(x) * 10 ^ (digits + 1)) / 10
  z = floor(z * posneg + 0.5) / 10 ^ digits
  return(z)
}


# perform yeo-johnson transformation on selected continuous variables of df
yeo_transform <- function(data, continuous_vars, standardise = TRUE){
  
  new_df <- data
  for (i in seq_along(continuous_vars)){
    
    variable_str <- as.character(continuous_vars[[i]])
    
    if(standardise == TRUE){
      
      trans_results <- yeojohnson(data[[variable_str]], standardize = TRUE)
      
      new_vector <- trans_results$x.t
      new_df[[variable_str]] = new_vector
      
    }else{
      
      trans_results <- yeojohnson(data[[variable_str]], standardize = FALSE)
      
      new_vector <- trans_results$x.t
      new_df[[variable_str]] = new_vector
      
    }
    
    
  }
  
  return(new_df)
  
}



# compare selected clusters (all vars)
compare_clusters <- function(df, cluster_1, cluster_2){
  
  one_v_two = df %>%
    filter(cluster %in% c(cluster_1, cluster_2))
  
  overall_means = one_v_two %>%
    summarise(across(everything(), mean, na.rm=TRUE)) %>%
    pivot_longer(cols = -cluster) %>%
    select(-cluster) %>%
    rename('overall_mean' = 'value')
  
  cluster_means = one_v_two %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean, na.rm=TRUE)) %>%
    pivot_longer(cols = -cluster)
  
  cluster_plot_df = cluster_means %>%
    left_join(overall_means, by='name') %>%
    mutate(diff_mean = value - overall_mean) %>%
    arrange(cluster, desc(diff_mean))
  
  variate_order = unique(cluster_plot_df$name)
  
  group.colors <- c("1" = "#F8766D", "2" = "#7DAE01", "3" ="#00BFC4", "4" = "#C77CFF")
  
  p <- cluster_plot_df %>%
    ggplot(mapping = aes(x=diff_mean, y=name, group=cluster, color=cluster)) + 
    geom_line() +
    geom_hline(yintercept = 0, linetype='dashed') +
    aes(y=factor(name, levels=variate_order)) + 
    labs(title = paste('Hard k-means:',cluster_1,'vs', cluster_2), 
         y = 'variate',
         x = 'difference of cluster mean and mean of both clusters') + 
    scale_x_continuous(limits = c(-1.1, 1.1), breaks = seq(-1, 1, by = 0.5)) +
    scale_color_manual(values=group.colors) +
    theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), 
          axis.text.x = element_text(angle = 0))
  
  return(p)
  
}



# compare clusters in terms of continuous variables
compare_continuous_clusters <- function(standardised_df, continuous_var, reduced_vars, type_clustering, n_clusters, variate_order){
  
  if(!missing(reduced_vars)){
    
    cluster_medians <- standardised_df %>%
      select(c(cluster, continuous_var[continuous_var %in% reduced_vars])) %>%
      group_by(cluster) %>%
      summarise(across(everything(), median, na.rm=TRUE)) %>%
      pivot_longer(cols = -cluster)
    
    all_medians <- standardised_df %>%
      select(c(cluster, continuous_var[continuous_var %in% reduced_vars])) %>%
      summarise(across(-cluster, median, na.rm=TRUE)) %>%
      pivot_longer(cols = everything()) %>%
      rename(total_median = 'value')
    
  }else{
    
    cluster_medians <- standardised_df %>%
      select(c(cluster, all_of(continuous_var))) %>%
      group_by(cluster) %>%
      summarise(across(everything(), median, na.rm=TRUE)) %>%
      pivot_longer(cols = -cluster)
    
    all_medians <- standardised_df %>%
      select(c(cluster, all_of(continuous_var))) %>%
      summarise(across(-cluster, median, na.rm=TRUE)) %>%
      pivot_longer(cols = everything()) %>%
      rename(total_median = 'value')
    
  }

  
  compare_var_plot <- cluster_medians %>%
    mutate(cluster = as.character(cluster)) %>%
    left_join(all_medians, by='name') %>%
    arrange(name, cluster) %>%
    
    mutate(ratio = value/total_median,
           abs_ratio = abs(value/total_median),
           diff = value - total_median,
           abs_diff = abs(value - total_median)) %>%
    
    group_by(name) %>%
    mutate(max_diff = max(abs_diff),
           max_ratio = max(abs_ratio)) %>%
    ungroup() %>%
    arrange(desc(max_ratio), name) %>%
    
    mutate(variate = factor(name))

  
  if(missing(variate_order)){

  p1 <- ggplot(data = compare_var_plot, mapping = aes(x=variate, y=ratio, group=cluster, color=cluster)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype='dashed')+
    aes(x=fct_inorder(variate)) +
    labs(title = paste0(type_clustering,': ratio of medians'),
         x = 'variate',
         y = 'ratio of cluster median to overall median') +
    theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 30, vjust=1, hjust=1))

    p2 <- ggplot(data = compare_var_plot, mapping = aes(x=variate, y=abs_ratio, group=cluster, color=cluster)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype='dashed')+
      aes(x=fct_inorder(variate)) +
      labs(title = paste0(type_clustering,': ratio of absolute medians'),
           x = 'variate',
           y = 'absolute ratio of cluster median to overall median') +
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 30, vjust=1, hjust=1))

    compare_var_plot <- compare_var_plot %>%
      group_by(variate) %>%
      mutate(max_val_cluster = case_when(cluster==n_clusters ~ diff)) %>%
      fill(max_val_cluster, .direction = 'downup') %>%
      ungroup() %>%
      arrange(max_val_cluster, desc(cluster), name)

    p3 <- ggplot(data = compare_var_plot, mapping = aes(x=variate, y=diff, group=cluster, color=cluster)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype='dashed') +
      aes(x=fct_inorder(variate)) +
      labs(title = paste0(type_clustering,': differences of median'),
           x = 'variate',
           y = 'difference of cluster median and overall median') +
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 60, vjust=1, hjust=1))

    p4 <- compare_var_plot %>%
      group_by(variate) %>%
      mutate(max_abs_cluster = case_when(cluster==3 ~ abs_diff)) %>%
      fill(max_abs_cluster, .direction = 'downup') %>%
      ungroup() %>%
      arrange(desc(max_abs_cluster), desc(cluster), name) %>%
      ggplot(mapping = aes(x=variate, y=abs_diff, group=cluster, color=cluster)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype='dashed') +
      aes(x=fct_inorder(variate)) +
      labs(title = paste0(type_clustering,': absolute differences of median'),
           x = 'variate',
           y = 'absolute difference of cluster median and overall median') +
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 60, vjust=1, hjust=1))

  }else{

    p1 <- ggplot(data = compare_var_plot, mapping = aes(x=variate, y=ratio, group=cluster, color=cluster)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype='dashed')+
      aes(x=factor(variate, levels=variate_order)) +
      labs(title = paste0(type_clustering,': ratio of medians'),
           x = 'variate',
           y = 'ratio of cluster median to overall median') +
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 30, vjust=1, hjust=1))

    p2 <- ggplot(data = compare_var_plot, mapping = aes(x=variate, y=abs_ratio, group=cluster, color=cluster)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype='dashed')+
      aes(x=factor(variate, levels=variate_order)) +
      labs(title = paste0(type_clustering,': ratio of absolute medians'),
           x = 'variate',
           y = 'absolute ratio of cluster median to overall median') +
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 30, vjust=1, hjust=1))

    compare_var_plot <- compare_var_plot %>%
      group_by(variate) %>%
      mutate(max_val_cluster = case_when(cluster==n_clusters ~ diff)) %>%
      fill(max_val_cluster, .direction = 'downup') %>%
      ungroup() %>%
      arrange(max_val_cluster, desc(cluster), name)

    p3 <- ggplot(data = compare_var_plot, mapping = aes(x=variate, y=diff, group=cluster, color=cluster)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype='dashed') +
      aes(x=factor(variate, levels=variate_order)) +
      labs(title = paste0(type_clustering,': differences of median'),
           x = 'Variate',
           y = 'Difference between standardised cluster median and overall median') +
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

    p4 <- compare_var_plot %>%
      group_by(variate) %>%
      mutate(max_abs_cluster = case_when(cluster==3 ~ abs_diff)) %>%
      fill(max_abs_cluster, .direction = 'downup') %>%
      ungroup() %>%
      arrange(desc(max_abs_cluster), desc(cluster), name) %>%
      ggplot(mapping = aes(x=variate, y=abs_diff, group=cluster, color=cluster)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype='dashed') +
      aes(x=factor(variate, levels=variate_order)) +
      labs(title = paste0(type_clustering,': absolute differences of median'),
           x = 'variate',
           y = 'absolute difference of cluster median and overall median') +
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

  }

  return(list(p1,p2,p3,p4, compare_var_plot))
}



# compare clusters by categorical variables
compare_categorical_clusters <- function(standardised_df, continuous_var, type_clustering, n_clusters, variate_order){
  
  column_names <- colnames(standardised_df)
  categorical_cols <- column_names[column_names %notin% continuous_var]
  
  cluster_categorical_means <- standardised_df %>%
    select(all_of(categorical_cols)) %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean, na.rm=TRUE)) %>%
    pivot_longer(cols = -cluster)
  
  all_categorical_means <- standardised_df %>%
    select(all_of(categorical_cols)) %>%
    summarise(across(-cluster, mean, na.rm=TRUE)) %>%
    pivot_longer(cols = everything()) %>%
    rename(total_mean = 'value')
  
  if(missing(variate_order)){
  
    # plot proportions
    df <- cluster_categorical_means %>%
      mutate(cluster = as.character(cluster)) %>%
      left_join(all_categorical_means, by='name') %>%
      arrange(name, cluster) %>%
      
      mutate(ratio = value/total_mean,
             diff = value - total_mean) %>%
      
      group_by(name) %>%
      mutate(max_diff = max(diff),
             max_ratio = max(ratio)) %>%
      ungroup() %>%
      
      mutate(variate = factor(name)) %>%
      
      # order x-axis according to last group
      group_by(variate) %>%
      mutate(max_cluster = case_when(cluster==n_clusters ~ diff)) %>%
      fill(max_cluster, .direction = 'downup') %>%
      ungroup() %>%
      arrange(max_cluster, desc(cluster), name)
      
    p1 <- df %>% ggplot(mapping = aes(x=variate, y=diff, group=cluster, color=cluster)) + 
      geom_line() +
      aes(x=fct_inorder(variate)) +
      labs(title = paste0(type_clustering,': differences of cluster proportions'), 
           x = 'variate',
           y = 'difference of cluster proportions to overall proportions') + 
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 60, vjust=1, hjust=1))
  }else{
    
    df <- cluster_categorical_means %>%
      mutate(cluster = as.character(cluster)) %>%
      left_join(all_categorical_means, by='name') %>%
      arrange(name, cluster) %>%
      
      mutate(ratio = value/total_mean,
             diff = value - total_mean) %>%
      
      group_by(name) %>%
      mutate(max_diff = max(diff),
             max_ratio = max(ratio)) %>%
      ungroup() %>%
      
      mutate(variate = factor(name)) %>%
      
      # order x-axis according to last group
      group_by(variate) %>%
      mutate(max_cluster = case_when(cluster==n_clusters ~ diff)) %>%
      fill(max_cluster, .direction = 'downup') %>%
      ungroup() %>%
      arrange(max_cluster, desc(cluster), name) 
      
    p1 <- df %>%  ggplot(mapping = aes(x=variate, y=diff, group=cluster, color=cluster)) + 
      geom_line() +
      aes(x=factor(variate, levels = variate_order)) +
      labs(title = paste0(type_clustering,': differences of cluster proportions'), 
           x = 'variate',
           y = 'difference of cluster proportions to overall proportions') + 
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "grey95"), axis.text.x = element_text(angle = 60, vjust=1, hjust=1))
    
  }
  
  return(list(p1, df))
  
}


# compute accuracy of imputations by RMSE
rmse_results <- function(imputed_df, iterations, selected_var){
  var = as.character(selected_var)
  
  results_list = list()
  for(i in 1:iterations){
    
    x <- complete(imputed_df, i)
    x = x[var]
    colnames(x) <- 'predicted'
    
    x$true = pre_impute_df[[var]]
    
    rmse = sqrt(mean((x$true - x$predicted)^2))
    
    iter_results <- as_tibble(rmse) %>% rename('rmse' = 'value')
    iter_results$iter = i
    
    results_list[[i]] = iter_results
    
  }
  
  results <- do.call('rbind', results_list)
  
  return(results)  
}



# plot simple correlation plots between desired (fixed_var) and other continuous variables 
plot_corr <- function(df, fixed_var, cont_vars, which_cluster, colors){
  
  remaining_vars = cont_vars[cont_vars %notin% c(fixed_var)]
  
  corr_plots = list()
  
  if(!missing(which_cluster) && !missing(colors)){
    
    for(i in seq_along(remaining_vars)){
      var = remaining_vars[[i]]
      
      corr_plots[[i]] = df %>%
        filter(cluster == which_cluster) %>%
        ggplot(aes(x = .data[[var]], y = .data[[fixed_var]], group=cluster, color=cluster)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "loess", se = TRUE, formula = y ~ x) +
        labs(title = paste(var,'vs.', fixed_var), x=var, y=fixed_var) +
        scale_color_manual(values = colors) +
        theme_minimal()
    }
    
  }else{
    
    for(i in seq_along(remaining_vars)){
      var = remaining_vars[[i]]
      
      corr_plots[[i]] = df %>%
        ggplot(aes(x = .data[[var]], y = .data[[fixed_var]])) +
        geom_point(alpha = 0.2) +
        geom_smooth(method = "loess", se = TRUE, formula = y ~ x) +
        labs(title = paste(var,'vs.', fixed_var), x=var, y=fixed_var) +
        theme_minimal()
    }
  }
  
  return(corr_plots)
  
}




# prepare categorical cluster differences for word cloud
prep_categorical_wordcloud <- function(cluster_df, cont_vars, n_clusters = 3){
  
  column_names <- colnames(cluster_df)
  categorical_cols <- column_names[column_names %notin% cont_vars]
  
  cluster_categorical_means <- cluster_df %>%
    select(all_of(categorical_cols)) %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean, na.rm=TRUE)) %>%
    pivot_longer(cols = -cluster)
  
  all_categorical_means <- cluster_df %>%
    select(all_of(categorical_cols)) %>%
    summarise(across(-cluster, mean, na.rm=TRUE)) %>%
    pivot_longer(cols = everything()) %>%
    rename(total_mean = 'value')
  
  output <- cluster_categorical_means %>%
    mutate(cluster = as.character(cluster)) %>%
    left_join(all_categorical_means, by='name') %>%
    arrange(name, cluster) %>%
    
    mutate(ratio = value/total_mean,
           diff = value - total_mean) %>%
    
    group_by(name) %>%
    mutate(max_diff = max(diff),
           max_ratio = max(ratio)) %>%
    ungroup() %>%
    
    mutate(variate = factor(name)) %>%
    
    # order x-axis according to last group
    group_by(variate) %>%
    mutate(max_cluster = case_when(cluster==n_clusters ~ diff)) %>%
    fill(max_cluster, .direction = 'downup') %>%
    ungroup() %>%
    arrange(max_cluster, desc(cluster), name)
  
  return(output)
  
}




# word cloud to represent clusters
plot_wordcloud <- function(cont_plot_df, cat_plot_df, select_cluster = 1, remove = c('PriorAnyAntibiotic', 'PriorCoamoxiclav', 'PriorBetaLactam', 'PriorAmoxicillin', 'InactiveBaseline')){
  
  df_1 <- cont_plot_df %>%
    filter(cluster == select_cluster) %>%
    select(name, diff, abs_diff) %>%
    # mutate(name = str_replace(name, "_", " ")) %>%
    mutate(name = str_replace(name, "Specialty_Other", "Other specialty")) %>%
    mutate(name = str_remove(name, "_1"),
           name = str_remove(name, "_0"),
           name = str_remove(name, "Specialty_"),
           name = str_replace(name, "_", " ")) %>%
    mutate(cleaned_name = case_when(grepl('AVPU', name) ~ name,
                                    grepl('ESBL', name) ~ name,
                                    grepl('Ecoli', name) ~ 'E. coli',
                                    grepl('cluster', name) ~ name,
                                    grepl('BMI', name) ~ name,
                                    grepl('SBP', name) ~ name,
                                    grepl('DBP', name) ~ name,
                                    name == 'alt' ~ 'ALT',
                                    grepl('crp', name) ~ 'CRP',
                                    grepl('egfr', name) ~ 'eGFR',
                                    grepl('mch', name) ~ 'MCH',
                                    grepl('mcv', name) ~ 'MCV',
                                    grepl('AgeAtAdmission', name) ~ 'Age',
                                    grepl('Charlson sum', name) ~ 'Charlson',
                                    grepl('Elixhauser sum', name) ~ 'Elixhauser',
                                    grepl('LinkedSex F', name) ~ 'Female',
                                    grepl('LinkedSex M', name) ~ 'Male',
                                    grepl('IMDPercentile', name) ~ 'IMD percentile',
                                    grepl('AdmissionDay', name) ~ make_clean_names(name, sep_out = ' ', case='title'),
                                    TRUE ~ make_clean_names(name, sep_out = ' ', case='sentence'))) %>%
    mutate(cleaned_name = str_replace(cleaned_name, "Admission Day of Week", "Admission on")) %>%
    mutate(cleaned_name = str_remove(cleaned_name, "_2")) %>%
    filter(name %notin% remove) %>%
    
    mutate(scaled_value = rescale(abs_diff, to = c(0, 1), from = range(abs_diff)))
  
  # rename variables for presentation
  
  df_2 <- cat_plot_df %>%
    filter(cluster == select_cluster) %>%
    select(name, diff) %>%
    mutate(name = str_replace(name, "Specialty_Other", "Other specialty")) %>%
    mutate(name = str_remove(name, "_1"),
           name = str_remove(name, "_0"),
           name = str_remove(name, "Specialty_"),
           name = str_replace(name, "_", " ")) %>%
    mutate(cleaned_name = case_when(grepl('AVPU', name) ~ name,
                                    grepl('ESBL', name) ~ name,
                                    grepl('Ecoli', name) ~ 'E. coli',
                                    grepl('cluster', name) ~ name,
                                    grepl('BMI', name) ~ name,
                                    grepl('SBP', name) ~ name,
                                    grepl('DBP', name) ~ name,
                                    name == 'alt' ~ 'ALT',
                                    grepl('crp', name) ~ 'CRP',
                                    grepl('egfr', name) ~ 'eGFR',
                                    grepl('mch', name) ~ 'MCH',
                                    grepl('mcv', name) ~ 'MCV',
                                    grepl('AgeAtAdmission', name) ~ 'Age',
                                    grepl('Charlson sum', name) ~ 'Charlson',
                                    grepl('Elixhauser sum', name) ~ 'Elixhauser',
                                    grepl('LinkedSex F', name) ~ 'Female',
                                    grepl('LinkedSex M', name) ~ 'Male',
                                    grepl('IMDPercentile', name) ~ 'IMD percentile',
                                    grepl('AdmissionDay', name) ~ make_clean_names(name, sep_out = ' ', case='title'),
                                    TRUE ~ make_clean_names(name, sep_out = ' ', case='sentence'))) %>%
    mutate(cleaned_name = str_replace(cleaned_name, "Admission Day of Week", "Admission on")) %>%
    mutate(cleaned_name = str_remove(cleaned_name, "_2")) %>%
    filter(name %notin% remove) %>%
    filter(diff >= 0) %>%            # remove redundant categorical variables where difference from total mean is negative (since represented by other variables e.g., less community-onset --> more hospital-onset) 
    mutate(abs_diff = abs(diff)) %>%
    mutate(scaled_value = rescale(abs_diff, to = c(0, 1), from = range(abs_diff)))
  
  df_stack <- rbind(df_1, df_2) %>%
    mutate(rescale_100 = rescale(scaled_value, to = c(0, 100), from = range(scaled_value))) %>%
    mutate(rescale_100 = round_half_up(rescale_100, digits = 0))
  
  p1 <- df_stack %>% 
    mutate(sentiment = case_when(diff >= 0 ~ 'positive',
                                 diff <0 ~ 'negative')) %>%
    mutate(sentiment = factor(sentiment, levels=c('positive', 'negative'))) %>%
    acast(cleaned_name ~ sentiment, value.var = 'rescale_100', fill=0) %>%
    comparison.cloud(colors = c("darkred", "deepskyblue4"),
                     max.words = 100,
                     rot.per=0,
                     scale=c(3, 0.7))
  
  return(p1)
}


get_db_rl_estimates <- function(treated_n) {

  # subset data
  df_sub <- df_numeric %>%
    filter(hate_crime_treatment %in% c(1, treated_n)) %>%
    mutate(treated = ifelse(hate_crime_treatment == treated_n, 1, 0))
  
  df_sub$chinese_origin <- factor(df_sub$chinese_origin)
  df_sub$democracy_cat <- factor(df_sub$democracy_cat)
  df_sub$china_threatend_high <- factor(df_sub$china_threatend_high) 
  df_sub$DEM <- factor(df_sub$DEM) 
  df_sub$GOP <- factor(df_sub$GOP) 
  df_sub$immigrant <- factor(df_sub$immigrant) 
  df_sub$college <- factor(df_sub$college) 
  
  # known propensity score
  df_sub <- df_sub %>%
    mutate(ps = rep(mean(treated), nrow(df_sub)))
  
  # configure HTE
  basic_config() %>%
    add_known_propensity_score("ps") %>% # known because this is an experiment 
    add_outcome_model("SL.glmnet") %>% # ensemble of machine learning models 
    # discrete
    add_moderator("Stratified", chinese_origin, democracy_cat, china_threatend_high, DEM, GOP, immigrant, college) %>%
    # continuous
    #add_moderator("KernelSmooth", china_threat_index) %>%
    add_vimp(sample_splitting = FALSE) -> hte_cfg
  
  # estimate HTE
  df_sub %>%
    attach_config(hte_cfg) %>%
    # cross validation, 10 folds 
    make_splits(caseid, .num_splits = 10) %>%
    produce_plugin_estimates(
      racial_linked_fate, # y 
      treated, # a 
      chinese_origin, # x1  
      democracy_cat, # x2 
      china_threatend_high, # x3 
      DEM, # x4 
      GOP, # x5
      immigrant, 
      college
    ) %>%
    construct_pseudo_outcomes(racial_linked_fate, treated) -> df_sub_splits 
  
  df_sub_splits %>%
    estimate_QoI(chinese_origin, democracy_cat, china_threatend_high, DEM, GOP, immigrant, college) -> results
  
  return(results)
}

create_dummies <- function(df_numeric) {
  
  df_copy <- df_numeric
  
  # Chinese origin 
  df_copy$chinese_origin <- ifelse(df_copy$nat_origin == 1, 1, 0)
  
  # Party ID
  df_copy$DEM <- ifelse(df_copy$pid3 == 1 | df_copy$pid7 %in% c(1, 3), 1, 0)
  
  df_copy$GOP <- ifelse(df_copy$pid3 == 3 | df_copy$pid7 %in% c(7, 5), 1, 0)
  
  # Immigrant
  df_copy$immigrant <- ifelse(df_copy$immigrant %in% c(1, 2), 1, 0)
  
  # College
  df_copy$college <- ifelse(df_copy$educ %in% c(3:6), 1, 0)
  
  # Male
  df_copy$female <- ifelse(df_copy$gender4 == 2, 1, 0)
  
  # Age
  df_copy$age <- 2024 - df_copy$birthyr
  
  return(df_copy)
  
}

get_db_el_estimates <- function(treated_n) {
  
  # subset data
  df_sub <- df_numeric %>%
    filter(hate_crime_treatment %in% c(1, treated_n)) %>%
    mutate(treated = ifelse(hate_crime_treatment == treated_n, 1, 0))
  
  df_sub$chinese_origin <- factor(df_sub$chinese_origin)
  df_sub$democracy_cat <- factor(df_sub$democracy_cat)
  df_sub$china_threatend_high <- factor(df_sub$china_threatend_high) 
  df_sub$DEM <- factor(df_sub$DEM) 
  df_sub$GOP <- factor(df_sub$GOP) 
  df_sub$immigrant <- factor(df_sub$immigrant) 
  df_sub$college <- factor(df_sub$college) 
  
  # known propensity score
  df_sub <- df_sub %>%
    mutate(ps = rep(mean(treated), nrow(df_sub)))
  
  # configure HTE
  basic_config() %>%
    add_known_propensity_score("ps") %>% # known because this is an experiment 
    add_outcome_model("SL.glmnet") %>% # ensemble of machine learning models 
    # discrete
    add_moderator("Stratified", chinese_origin, democracy_cat, china_threatend_high, DEM, GOP, immigrant, college) %>%
    # continuous
    #add_moderator("KernelSmooth", china_threat_index) %>%
    add_vimp(sample_splitting = FALSE) -> hte_cfg
  
  # estimate HTE
  df_sub %>%
    attach_config(hte_cfg) %>%
    # cross validation, 10 folds 
    make_splits(caseid, .num_splits = 10) %>%
    produce_plugin_estimates(
      ethnic_linked_fate, # y 
      treated, # a 
      chinese_origin, # x1  
      democracy_cat, # x2 
      china_threatend_high, # x3 
      DEM, # x4 
      GOP, # x5
      immigrant, 
      college
    ) %>%
    construct_pseudo_outcomes(ethnic_linked_fate, treated) -> df_sub_splits 
  
  df_sub_splits %>%
    estimate_QoI(chinese_origin, democracy_cat, china_threatend_high, DEM, GOP, immigrant, college) -> results
  
  return(results)
}
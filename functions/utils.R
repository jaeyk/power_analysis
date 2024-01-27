
get_int_coder <- function(df, var_name) {
  
  # percentage agreement 
  pct_agr_score <- df %>%
    irr::agree()
  
  # calculate Cohen's kappa
  kappa_score <- df %>%
    irr::kappa2() 

  out <- data.frame(
    agree = pct_agr_score$value,
    kappa = kappa_score$value
  )
  
  out <- out %>% 
    mutate(var = var_name) %>%
    select(var, everything())
  
  return(out)
}

ols <- function(df, dv) { 
  
  lm_robust(eval(as.symbol(dv)) ~ factor(hate_crime_treatment), data = df) %>%
    tidy() %>%
    filter(!str_detect(term, "Int")) %>%
    mutate(outcome = dv)
  
}

nest_ols <- function(df, dv) {
  
  df %>%
    nest() %>%
    mutate(ols = map(data, ~ols(., dv))) %>%
    unnest(ols) %>% 
    mutate(term = case_when(
      str_detect(term, "4") ~ "Hate crime + political representation", 
      str_detect(term, "3") ~ "Hate crime + China threat", 
      str_detect(term, "2") ~ "Hate crime"
    )) %>%
    mutate(outcome = case_when(
      str_detect(outcome, "racial") ~ "Racial linked fate", 
      str_detect(outcome, "ethnic") ~ "Ethnic linked fate"
    ))  
}

mean_no_na <- function(x) mean(x, na.rm = T)

std_no_na <- function(x) sd(x, na.rm = T)/sqrt(length(x))

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
  df_copy$college <- ifelse(df_copy$educ %in% c(5:6), 1, 0)
  
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

conText <- function(formula, data, text_var = 'text', pre_trained, transform = TRUE, transform_matrix, bootstrap = TRUE, num_bootstraps = 20, stratify_by = NULL, permute = TRUE, num_permutations = 100, getcontexts = TRUE, window = 6, valuetype = "fixed", case_insensitive = TRUE, hard_cut = FALSE, verbose = TRUE){
  
  # extract target word
  target <- as.character(formula[[2]])
  # extract covariates
  covariates <- attr(terms(formula), which = "term.labels")
  # if getcontexts is specified, get context around target word
  if(getcontexts){
    # subset data to where target is present (this step helps speed up the get_context function)
    if(valuetype == 'fixed') target_present <- grep(target, dplyr::pull(data, text_var), fixed = TRUE)
    if(valuetype != 'fixed') target_present <- grep(target, dplyr::pull(data, text_var), ignore.case = case_insensitive)
    data <- data[target_present,]
    # add document id variable (used for merging back with covariates)
    data <- data %>% dplyr::mutate(docname = paste0("text", 1:nrow(.)))
    # apply get_context function (see get_context documentation)
    context <- get_context(x = dplyr::pull(data, text_var), target = target, window = window, valuetype = valuetype, case_insensitive = case_insensitive, hard_cut = hard_cut, verbose = verbose)
    # merge with metadata
    context <- dplyr::left_join(context, data[,c("docname", covariates)], by = "docname") %>% dplyr::mutate('(Intercept)' = 1)
    # otherwise take the full text as the context to be embedded
  }else{
    # add intercept
    context <- data %>% dplyr::mutate('(Intercept)' = 1) %>% dplyr::mutate(context = dplyr::pull(data, text_var))
  }
  
  # embed context to get dependent variable (note, aggregate is set to FALSE as we want an embedding for each instance)
  embeds_out <- embed_target(context$context, pre_trained = pre_trained, transform_matrix = transform_matrix, transform = transform, aggregate = FALSE, verbose)
  Y <- embeds_out$target_embedding
  if(verbose) cat('total observations included in regression:', nrow(Y), '\n')
  
  # regressors
  X <- context[embeds_out$obs_included,c('(Intercept)', covariates)]
  # run full sample ols
  full_sample_out <- run_ols(Y = Y, X = X)
  # outputs
  beta_cofficients <- full_sample_out$betas
  norm_tibble <- dplyr::tibble(Coefficient = names(full_sample_out$normed_betas), Normed_Estimate = unname(full_sample_out$normed_betas))
  # -------------------
  # bootstrapping
  # -------------------
  if(bootstrap){
    if(verbose) cat('starting bootstrapping \n')
    # bootstrapped ols
    bootstrap_out <- replicate(num_bootstraps, bootstrap_ols(Y = Y, X = X, stratify_by = stratify_by), simplify = FALSE)
    # average betas
    betas <- lapply(bootstrap_out, '[[', 'betas')
    bs_betas <- Reduce("+",betas)/length(betas)
    # summary statistics for normed betas
    normed_betas <- lapply(bootstrap_out, '[[', 'normed_betas') %>% do.call(rbind,.)
    mean_normed_betas <- apply(normed_betas, 2, mean)
    stderror_normed_betas <- 1/sqrt(num_bootstraps) * apply(normed_betas, 2, sd)
    bs_normed_betas <- dplyr::tibble(Coefficient = names(mean_normed_betas), Normed_Estimate = unname(mean_normed_betas), Std.Error = unname(stderror_normed_betas))
    # output
    beta_cofficients <- bs_betas
    norm_tibble <- bs_normed_betas
    # notice
    if(verbose) cat('done with bootstrapping \n')
  }
  # -------------------
  # permutation
  # -------------------
  if(permute){
    if(verbose) cat('starting permutations \n')
    # permuted ols
    permute_out <- replicate(num_permutations, permute_ols(Y = Y, X = X), simplify = FALSE)
    # compute empirical p-value
    permuted_normed_betas <- lapply(permute_out, '[[', 'normed_betas') %>% do.call(rbind,.)
    empirical_pvalue <- sweep(permuted_normed_betas, 2, 1/full_sample_out$normed_betas, '*')
    empirical_pvalue <- apply(empirical_pvalue, 2, function(x) sum(x>1)/length(x))
    # bind with results
    norm_tibble <- cbind(norm_tibble, Empirical_Pvalue = unname(empirical_pvalue))
    if(verbose) cat('done with permutations \n')
  }
  # -------------------
  # output
  # -------------------
  return(list(betas = beta_cofficients, normed_betas = norm_tibble))
}
# -----------------------------
#
# SUB-FUNCTIONS
#
# -----------------------------
#' Permute OLS
#'
#' Estimate empirical p-value using permutated regression
#'
#' @param Y vector of regression model's dependent variable (embedded context)
#' @param X data.frame of model independent variables (covariates)
#'
#' @return list with two elements, `betas` = list of beta_cofficients (D dimensional vectors);
#' `normed_betas` = tibble with the norm of the non-intercept coefficients
#'
permute_ols <- function(Y = NULL, X = NULL){
  # shuffle Y
  Y_permuted <- Y[sample(1:nrow(Y), replace = FALSE),]
  # run ols
  ols_out <- run_ols(Y = Y_permuted, X = X)
  # output
  return(ols_out)
}
#' Bootstrap OLS
#'
#' Bootstrap model coefficients and standard errors
#'
#' @param Y vector of regression model's dependent variable (embedded context)
#' @param X data.frame of model independent variables (covariates)
#' @param stratify_by covariates to stratify when bootstrapping
#'
#' @return list with two elements, `betas` = list of beta_cofficients (D dimensional vectors);
#' `normed_betas` = tibble with the norm of the non-intercept coefficients
#'
bootstrap_ols <- function(Y = NULL, X = NULL, stratify_by = NULL){
  # label instances
  X_bs <- cbind(obs = 1:nrow(X), X)
  # sample observations with replacement
  if(is.null(stratify_by)) X_bs <- dplyr::sample_n(X_bs, size = nrow(X_bs), replace = TRUE)else{
    X_bs <- X_bs %>% dplyr::group_by_at(stratify_by) %>% dplyr::sample_n(size = dplyr::n(), replace = TRUE) %>% dplyr::ungroup()
  }
  # subset Y to sampled observations
  Y_bs <- Y[X_bs$obs,]
  # remove observation label
  X_bs <- X_bs[,-1]
  # run ols
  ols_out <- run_ols(Y = Y_bs, X = X_bs)
  # output
  return(ols_out)
}
#' Run OLS
#'
#' Bootstrap model coefficients and standard errors
#'
#' @param Y vector of regression model's dependent variable (embedded context)
#' @param X data.frame of model independent variables (covariates)
#'
#' @return list with two elements, `betas` = list of beta_cofficients (D dimensional vectors);
#' `normed_betas` = tibble with the norm of the non-intercept coefficients
#'
run_ols <- function(Y = NULL, X = NULL){
  # convert X to a matrix
  X_mat <- as.matrix(X, ncol = ncol(X))
  # compute OLS bets hats
  betas <- solve(t(X_mat)%*%X_mat)%*%t(X_mat)%*%Y
  # normed betas
  vars <- setdiff(colnames(X), "(Intercept)") # identify non-intercept covariates (norm of intercept is not all that informative)
  normed_betas <- apply(matrix(betas[vars,], nrow = length(vars)), 1, function(x) norm(matrix(x, nrow = 1))) %>% setNames(vars)
  # output
  return(list('betas' = betas, 'normed_betas' = normed_betas))
}

## for bootstrapping 95% confidence intervals; Borrowed from Nick Camp's code from Jaren, Nick, and my shared project

theta <- function(x, xdata, na.rm = T) {
  mean(xdata[x], na.rm = na.rm)
}

ci.low <- function(x, na.rm = T) {
  mean(x, na.rm = na.rm) - quantile(bootstrap::bootstrap(1:length(x), 1000, theta, x, na.rm = na.rm)$thetastar, .025, na.rm = na.rm)
}

ci.high <- function(x, na.rm = T) {
  quantile(bootstrap::bootstrap(1:length(x), 1000, theta, x, na.rm = na.rm)$thetastar, .975, na.rm = na.rm) - mean(x, na.rm = na.rm)
}

# preprocessing

get_word_count <- function(data, stem = TRUE) {
  
  tidy_df <- data %>%
    # Tokenize
    unnest_tokens("word", value) %>%
    # Remove stop words
    anti_join(get_stopwords(), by = "word")
  
  if (stem == TRUE) {
    
    # Stemming
    tidy_df <- tidy_df %>% mutate(stem = wordStem(word))
    
    df_words <- tidy_df %>%
      count(date, stem, sort = TRUE)
    
    total_words <- df_words %>%
      group_by(date) %>%
      summarize(total = sum(n))
    
    joined_words <- left_join(df_words, total_words)
    
    tf_idf <- joined_words %>%
      # TF-IDF
      bind_tf_idf(stem, date, n)
    
  } else {
    
    df_words <- tidy_df %>% count(date, word, sort = TRUE)
    
    total_words <- df_words %>%
      group_by(date) %>%
      summarize(total = sum(n))
    
    joined_words <- left_join(df_words, total_words)
    
    tf_idf <- joined_words %>%
      # TF-IDF
      bind_tf_idf(word, date, n)
    
  }
  
  return(tf_idf)
  
}

clean_text <- function(full_text) {
  
  vec <- tolower(full_text) %>%
    # Remove all non-alpha characters
    gsub("[^[:alpha:]]", " ", .) 
  
  vec <- tm::removeWords(vec, words = c(stopwords(source = "snowball"), "asian", "asians", "american", "americans"))
  
  return(vec)
  
}

visualize_diag <- function(sparse_matrix, many_models){
  
  k_result <- many_models %>%
    mutate(exclusivity = purrr::map(topic_model, exclusivity),
           semantic_coherence = purrr::map(topic_model, semanticCoherence, sparse_matrix))
  
  
  k_result %>%
    transmute(K,
              "Exclusivity" = map_dbl(exclusivity, mean),
              "Semantic coherence" = map_dbl(semantic_coherence, mean)) %>%
    pivot_longer(cols = c("Exclusivity", "Semantic coherence"),
                 names_to = "Metric",
                 values_to = "Value") %>%
    ggplot(aes(K, Value, color = Metric)) +
    geom_line(size = 1.5, show.legend = FALSE) +
    labs(x = "K (number of topics)",
         y = NULL) +   
    facet_wrap(~Metric, scales = "free_y") +
    theme_minimal()
  
}

pull_est <- function(res) {
  
  res %>%
    filter(estimand == "MCATE") %>%
    filter(!str_detect(level, "Non-G|Non-D")) %>%
    pull(estimate)
  
}
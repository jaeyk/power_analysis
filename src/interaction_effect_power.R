# import pkgs 
if (!require(pacman)) install.packages("pacman")

pacman::p_load(
  DeclareDesign,
  glue, # string interoperability
  here, # computational reproducibility
  purrr, # functional programming 
  patchwork, 
  plotly
) 

equity_designer <- function(sample_size, group_size_ratio = 0.5, interaction_effect_ratio = 0.5, effect_size, assignment = "complete") {
  
  # Model  
  pop <- fabricate(
    N = sample_size,
    X = draw_binary(N = N, prob = group_size_ratio),
    h_effect_size = effect_size*interaction_effect_ratio,
    potential_outcomes(
      # binomial outcome 
      Y ~ effect_size*Z + X + h_effect_size*Z*X + rnorm(N)
    ))
  
  
  design <- declare_model(data = pop) +
    
    # Inquiry 
    declare_inquiries(
      CATE_X0 = mean(Y_Z_1[X == 0]- Y_Z_0[X == 0]), 
      CATE_X1 = mean(Y_Z_1[X == 1]- Y_Z_0[X == 1]), 
      diff_in_CATEs = CATE_X1 - CATE_X0) +
    
    # Data strategy   
    declare_sampling(S = complete_rs(N)) # (1) unit # you can change the complete random sampling with the stratified random sampling 
  
  if (assignment == "complete") {
    design <- design + declare_assignment(Z = complete_ra(N))
  }
  
  if (assignment == "block") {
    design <- design + declare_assignment(Z = block_ra(blocks = X))
  }
  
  design <- design +
    declare_measurement(Y = reveal_outcomes(Y ~ Z)) + # (3) measurement
    
    # Answer strategy  +
    declare_estimator(Y ~ Z + X + Z * X, 
                      term = "Z:X", 
                      inquiry = "diff_in_CATEs")
  
  return(design)
}

sample_sizes <- seq(2000, 10000, by = 100)
effect_sizes <- c(c(0.15 + 0.18)/2) 
interaction_effect_ratios <- c(0.5, 1)

# using complete randomization
equity_designs <- expand_design(
  equity_designer,
  sample_size = sample_sizes,
  interaction_effect_ratio = interaction_effect_ratios,
  effect_size = effect_sizes,
  assignment = "complete"
)

e_diag_outs <- diagnose_design(equity_designs)

search_space <- expand.grid(sample_size = sample_sizes, 
            interaction_effect_ratio = interaction_effect_ratios, 
            effect_size = effect_sizes)

diags <- mutate(e_diag_outs$diagnosands_df, type = "Difference in CATEs")

write_csv(e_diag_outs$diagnosands_df, here("processed_data", "e_diag_out.csv"))

diag_df_diff <- diags %>% filter(inquiry %in% c("ATE", "diff_in_CATEs"))

diag_df_diff$sample_size <- search_space$sample_size
diag_df_diff$effect_size <- search_space$effect_size
diag_df_diff$interaction_effect_ratio <- search_space$interaction_effect_ratio

power_plot <- diag_df_diff %>%
    ggplot(aes(
    x = sample_size, y = power,
    ymax = power + 1.96 * `se(power)`,
    ymin = power - 1.96 * `se(power)`,
    fill = type
  )) +
  geom_point(col = "blue") +
  geom_line(alpha = 0.7) +
  geom_ribbon(alpha = 0.3) +
  labs(
    title = "Difference in CATEs",
    subtitle = "Results are based on simulations",
    x = "Sample size",
    y = "Statistical power",
    fill = "Estimand type"
  ) +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red") +
  facet_grid(~interaction_effect_ratio) +
  theme(legend.position = "none")

png(filename = here("outputs", "int_power.png"),
    height = 7, 
    width = 10, 
    unit = "in", 
    res = 1200)

power_plot

dev.off()
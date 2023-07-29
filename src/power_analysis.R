# install pkgs

if (!require(pacman)) install.packages("pacman")

pacman::p_load(
  # DeclareDesign stuff
  DeclareDesign,
  DesignLibrary,

  # tidyverse stuff
  tidyverse,
  glue,
  purrr,
  broom,
  patchwork,
  here,
  plotly
)

# Create a custom multi arm designer
# Reference: https://declaredesign.org/r/designlibrary/articles/multi_arm.html

my_multi_arm_designer <- function(sample_size = 1000) {
  
  N <- sample_size
  outcome_means <- c(0.7, 0.7 + 0.05, 0.7 + 0.02)
  outcome_sds <- c(0.01, 0.02, 0.02)
  sd_i <- 0.3 

  population <- declare_population(
    N = N,
    u_1 = rnorm(N, 0, outcome_sds[1L]),
    u_2 = rnorm(N, 0, outcome_sds[2L]),
    u_3 = rnorm(N, 0, outcome_sds[3L]),
    u = rnorm(N) * sd_i
  )

  potential_outcomes <- declare_potential_outcomes(
    formula = Y ~
      (outcome_means[1] + u_1) * (Z == "1") +
      (outcome_means[2] + u_2) * (Z == "2") +
      (outcome_means[3] + u_3) * (Z == "3") +
      u,
    conditions = c("1", "2", "3"),
    assignment_variables = Z
  )

  estimand <- declare_inquiries(
    ate_Y_2_1 = mean(Y_Z_2 - Y_Z_1),
    ate_Y_3_1 = mean(Y_Z_3 - Y_Z_1),
    ate_Y_3_2 = mean(Y_Z_3 - Y_Z_2)
  )

  assignment <- declare_assignment(Z = complete_ra(N,
    num_arms = 3,
    conditions = c("1", "2", "3")
  ))

  reveal_Y <- declare_reveal(assignment_variables = Z)

  estimator <- declare_estimator(handler = function(data) {
    estimates <- rbind.data.frame(
      ate_Y_2_1 = difference_in_means(
        formula = Y ~ Z,
        data = data,
        condition1 = "1", condition2 = "2"
      ),
      ate_Y_3_1 = difference_in_means(
        formula = Y ~ Z,
        data = data,
        condition1 = "1", condition2 = "3"
      ),
      ate_Y_3_2 = difference_in_means(
        formula = Y ~ Z,
        data = data,
        condition1 = "2", condition2 = "3"
      )
    )

    names(estimates)[names(estimates) == "N"] <- "N_DIM"

    estimates$estimator <- c(
      "DIM (Z_2 - Z_1)", "DIM (Z_3 - Z_1)",
      "DIM (Z_3 - Z_2)"
    )

    estimates$inquiry <- rownames(estimates)

    estimates$estimate <- estimates$coefficients

    estimates$term <- NULL

    return(estimates)
  })

  multi_arm_design <- population + potential_outcomes + assignment + reveal_Y + estimand + estimator

  return(multi_arm_design)
}

# Test the designer 

diagnose_design(my_multi_arm_designer(sample_size = 2000))

# Create the design sets

designs <- expand_design(
  my_multi_arm_designer, # multi arm designer
  sample_size = seq(1500, 5000, by = 50)
)

# Diagnose

diag_outs <- diagnose_design(designs, sims = 500)

# Visualize

power_plot <- diag_outs$diagnosands_df %>%
  ggplot(aes(
    x = sample_size, y = power)) +
  geom_line() +
  geom_ribbon(aes(ymax = power + 1.96 * `se(power)`,
                  ymin = power - 1.96 * `se(power)`),
              alpha = 0.2) +
  facet_wrap(~estimator) +
  labs(
    x = "Sample size",
    y = "Diagnosand: power"
  ) +
  geom_hline(yintercept = 0.8, col = "red", linetype = "dashed")

power_plot 

ggsave(here("outputs", "power_analysis.png"))

plotly::ggplotly(power_plot)

# Calculate the interaction effect 

calculate_interaction_effect_power <- function(interaction_effect_ratio) {
  
  main_power <- round(pnorm(2.8, mean = 1.96, sd = 1),2)
  
  interaction_power <- round(pnorm(2.8*interaction_effect_ratio*0.5, mean = 1.96, sd = 1),2)
  
  interaction_sample_size <- (1/(interaction_effect_ratio*0.5))**2
  
  out <- data.frame(
    interaction_effect_ratio = interaction_effect_ratio, 
    interaction_power = interaction_power,
    interaction_sample_size = interaction_sample_size
  )
  
  return(out)
  
}

calculate_interaction_effect_power(interaction_effect_ratio = 1); 
calculate_interaction_effect_power(interaction_effect_ratio = 0.5)

# install pkgs 

if (!require(pacman)) install.packages("pacman")

pacman::p_load(
  # DeclareDesign stuff
  DeclareDesign, 
  fabricatr, 
  randomizr, 
  estimatr, 
  DesignLibrary,
  
  # Power analysis stuff
  bmem, # bootstrap
  powerMediation,
  easypower, 
  pwr,
  
  # tidyverse stuff
  tidyverse, 
  glue, 
  purrr,
  broom,
  patchwork,
  here
  )

# the conventional threshold for the power analysis is alpha = .8, (so by definition beta = .2). alpha = false positive, beta = false negative 

# cohen'd = effect size, small = 0.2, medium = 0.5, large = 0.8

# the required n size for one arm
calculate_power <- function(delta_level) {

  out <- power.t.test(power = .8, 
               sig.level = .05, 
               delta = delta_level) 
  
  out <- out %>% tidy()
  
  out <- out %>% 
    mutate(one_arm_n = n/2)
  
  return(out)
  
}


# the required n size for two different conditions 
calculate_parallel_power <- function(delta) {

  combined <- data.frame(
    strong_treatment_n = calculate_power(delta_level = delta)$one_arm_n,
    weak_treatment_n = calculate_power(delta_level = delta/2)$one_arm_n)
  
  return(combined)
  
}

# run the power analysis over the numeric vector
deltas <- seq(0.1, 0.5, by = 0.01)
sums <- map_dfr(deltas, calculate_parallel_power)

sums$delta <- deltas

n2_plot <- sums %>%
  pivot_longer(cols = matches("_n"),
               names_to = "treatment_conditions",
               values_to = "n") %>%
  mutate(treatment_conditions = str_replace_all(treatment_conditions, "_treatment_n", "")) %>%
  ggplot(aes(x = delta, y = n, col = treatment_conditions)) +
  geom_point() +
  labs(x = "delta", 
       y = "sample size")

ratio_plot <- sums %>%
  mutate(ratio = weak_treatment_n/strong_treatment_n) %>%
  ggplot(aes(x = delta, y = ratio)) +
  geom_point() +
  labs(x = "delta", 
       y = "ratio btw the N required for strong and weak treatments")

n2_plot + ratio_plot

ggsave(here("outputs", "unlinking_linked_fate.png"))

# declare design 

## Step 1

simple_designer <- function(sample_size, effect_size) {
  declare_model(
    N = sample_size, 
    U = rnorm(N),
    potential_outcomes(Y ~ effect_size * Z + U)
  ) +
    declare_inquiry(PATE = mean(Y_Z_1 - Y_Z_0)) +
    declare_sampling(S = complete_rs(N, n = 50)) +
    declare_assignment(Z = complete_ra(N, prob = 0.5)) +
    declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
    declare_estimator(
      Y ~ Z, .method = difference_in_means, inquiry = "PATE"
    )
}

## Step 2

design <- simple_designer(sample_size = 100, effect_size = 0.25)

## Step 3

designs <- expand_design(
  simple_designer, 
  sample_size = c(100, 500, 1000), 
  effect_size = 0.25
)

diag_outs <- diagnose_design(designs)

diag_outs$diagnosands_df$power

# Three arm design

custom_multi_arm_design <- function(N) {
  
  three_arm <- multi_arm_designer(
      N = N, 
      m_arms = 3,
      outcome_means = c(0, 0.2, 0.1)
      )

  diag_outs <- diagnose_design(three_arm)

  power <- diag_outs$diagnosands_df$power
  
  out <- data.frame("power" = power,
                    "condition" =  c("Control", "Strong", "Weak")) %>%
    mutate(sample_size = N)
  
  return(out)
  
}

outs <- map_dfr(seq(1000, 5000, by = 200), custom_multi_arm_design)

outs %>%
  ggplot(aes(x = sample_size, y = power, col = condition)) +
  geom_point() +
  theme_light() +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")

ggsave(here("outputs", "multi_arm_design.png"))
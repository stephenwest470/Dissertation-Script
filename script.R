# Required libraries
library(tidyverse)
library(ggnewscale)

# Input data
data <- read.csv("input.csv", stringsAsFactors = FALSE)

# Population reference point
known_population_midpoint <- 800
known_population <- data$population[data$midpoint == known_population_midpoint]

# Proxy weights
proxy_1_weight <- 10
proxy_2_weight <- 30
proxy_3_weight <- proxy_1_weight
proxy_4_weight <- proxy_2_weight

# Prior mean and variance
past_prior_mean <- 0.005
past_prior_sigma <- 0.001
future_prior_mean <- -0.0005
future_prior_sigma <- 0.^01

# Credible intervals (z-score)
post_sigma_90CI <- 1.645
post_sigma_50CI <- 0.675

# Proxy sums
data <- data %>%
  mutate(
    proxy_sum = proxy_1 + proxy_2 + proxy_3 + proxy_4,
    weighted_proxy_1 = (proxy_1 / length) * proxy_1_weight,
    weighted_proxy_2 = (proxy_2 / length) * proxy_2_weight,
    weighted_proxy_3 = (proxy_3 / length) * proxy_3_weight,
    weighted_proxy_4 = (proxy_4 / length) * proxy_4_weight,
    weighted_proxy_sum = weighted_proxy_1 + weighted_proxy_2 +
      weighted_proxy_3 + weighted_proxy_4
  )

# Proxy ROC
data <- data[order(data$midpoint), ]
data$proxy_ROC <- NA_real_
for (i in seq_len(nrow(data) - 1)) {
  dt <- data$midpoint[i + 1] - data$midpoint[i]
  data$proxy_ROC[i] <- -log(data$weighted_proxy_sum[i] /
                              data$weighted_proxy_sum[i + 1]) / dt
}

# Data variance
proxy_var <- var(data$proxy_ROC, na.rm = TRUE)
data_var <- proxy_var / pmax(data$weighted_proxy_sum, 1)

# Posterior ROC
data <- data %>%
  mutate(
    prior_mean = ifelse(midpoint < known_population_midpoint,
                        past_prior_mean, future_prior_mean
    ),
    prior_sigma = ifelse(midpoint < known_population_midpoint,
                         past_prior_sigma, future_prior_sigma
    )
  ) %>%
  mutate(
    post_var = 1 / ((1 / data_var) + (1 / prior_sigma^2)),
    post_sigma_raw = sqrt(post_var),
    bayes_mean = ((proxy_ROC / data_var) +
                    (prior_mean / prior_sigma^2)) * post_var,
    bayes_lower_90CI = bayes_mean - post_sigma_90CI * post_sigma_raw,
    bayes_upper_90CI = bayes_mean + post_sigma_90CI * post_sigma_raw,
    bayes_lower_50CI = bayes_mean - post_sigma_50CI * post_sigma_raw,
    bayes_upper_50CI = bayes_mean + post_sigma_50CI * post_sigma_raw,
    post_ROC_mean = bayes_mean,
    post_ROC_lower_90CI = bayes_lower_90CI,
    post_ROC_upper_90CI = bayes_upper_90CI,
    post_ROC_lower_50CI = bayes_lower_50CI,
    post_ROC_upper_50CI = bayes_upper_50CI
  ) %>%
  select(-post_var, -post_sigma_raw,
         -bayes_mean, -bayes_lower_90CI, -bayes_upper_90CI,
         -bayes_lower_50CI, -bayes_upper_50CI
  )

# Posterior population
index <- which(data$midpoint == known_population_midpoint)
data[index, c(
  "post_pop_mean", "pop_lower_90CI", "pop_upper_90CI",
  "pop_lower_50CI", "pop_upper_50CI"
)] <- rep(known_population, 5)

for (i in seq(index - 1, 1)) {
  if (!is.na(data$post_ROC_mean[i]) &&
      !is.na(data$post_pop_mean[i + 1])) {
    dt <- data$midpoint[i + 1] - data$midpoint[i]
    data$post_pop_mean[i] <- data$post_pop_mean[i + 1] *
      exp(-data$post_ROC_mean[i] * dt)
    data$pop_lower_90CI[i] <- data$pop_lower_90CI[i + 1] *
      exp(-data$post_ROC_lower_90CI[i] * dt)
    data$pop_upper_90CI[i] <- data$pop_upper_90CI[i + 1] *
      exp(-data$post_ROC_upper_90CI[i] * dt)
    data$pop_lower_50CI[i] <- data$pop_lower_50CI[i + 1] *
      exp(-data$post_ROC_lower_50CI[i] * dt)
    data$pop_upper_50CI[i] <- data$pop_upper_50CI[i + 1] *
      exp(-data$post_ROC_upper_50CI[i] * dt)
  }
}

for (i in seq(index + 1, nrow(data))) {
  if (!is.na(data$post_ROC_mean[i - 1]) &&
      !is.na(data$post_pop_mean[i - 1])) {
    dt <- data$midpoint[i] - data$midpoint[i - 1]
    data$post_pop_mean[i] <- data$post_pop_mean[i - 1] *
      exp(data$post_ROC_mean[i - 1] * dt)
    data$pop_lower_90CI[i] <- data$pop_lower_90CI[i - 1] *
      exp(data$post_ROC_lower_90CI[i - 1] * dt)
    data$pop_upper_90CI[i] <- data$pop_upper_90CI[i - 1] *
      exp(data$post_ROC_upper_90CI[i - 1] * dt)
    data$pop_lower_50CI[i] <- data$pop_lower_50CI[i - 1] *
      exp(data$post_ROC_lower_50CI[i - 1] * dt)
    data$pop_upper_50CI[i] <- data$pop_upper_50CI[i - 1] *
      exp(data$post_ROC_upper_50CI[i - 1] * dt)
  }
}

# Interpolate prior population estimates
data$population_na <- NA_real_
known_pop_ref <- which(!is.na(data$population))
for (seg in seq_len(length(known_pop_ref) - 1)) {
  i <- known_pop_ref[seg]; j <- known_pop_ref[seg + 1]
  r <- log(data$population[j] / data$population[i]) /
    (data$midpoint[j] - data$midpoint[i])
  data$population_na[i] <- data$population[i]
  for (k in (i + 1):(j - 1)) {
    dt <- data$midpoint[k] - data$midpoint[k - 1]
    data$population_na[k] <- data$population_na[k - 1] * exp(r * dt)
  }
  data$population_na[j] <- data$population[j]
}

data$proxy_label <- "Concurrent Features"

# Export posterior as CSV
export_posterior <- function(df, file = "posterior_population.csv") {
  df_out <- df %>%
    select(midpoint,
           pop_upper_90CI,
           pop_upper_50CI,
           post_pop_mean,
           pop_lower_50CI,
           pop_lower_90CI
           )
  write.csv(df_out, file, row.names = FALSE)
}

export_posterior(data)

# Plot settings

plot <- ggplot(data, aes(x = midpoint)) +
  new_scale_fill() +
  geom_bar(aes(y = weighted_proxy_sum, fill = proxy_label),
           stat = "identity", color = "black", size = 0.5
  ) +
  scale_fill_manual(name = "Population Proxy",
                    values = c("Concurrent Features" = "white"),
                    guide = guide_legend(order = 3)) +
  new_scale_fill() +
  geom_ribbon(aes(ymin = pop_lower_90CI / 1000,
                  ymax = pop_upper_90CI / 1000,
                  fill = "90% CI"),
              alpha = 0.4, na.rm = TRUE
  ) +
  geom_ribbon(aes(ymin = pop_lower_50CI / 1000,
                  ymax = pop_upper_50CI / 1000,
                  fill = "50% CI"),
              alpha = 0.8, na.rm = TRUE
  ) +
  scale_fill_manual(name = "Credible Intervals",
                    values = c("50% CI" = "grey50",
                               "90% CI" = "grey70"),
                    guide = guide_legend(order = 2)) +
  new_scale_color() +
  geom_line(aes(y = post_pop_mean / 1000, color = "Posterior Mean"),
            size = 1.5, na.rm = TRUE
  ) +
  geom_line(aes(y = population_na / 1000, color = "Prior Interpolation"),
            linetype = "dotted", size = 1, na.rm = TRUE
  ) +
  geom_point(aes(y = population / 1000, color = "Prior Estimates"),
             size = 3, shape = 21, fill = "grey20", stroke = 1, na.rm = TRUE
  ) +
  scale_color_manual(name = "Population",
                     values = c("Posterior Mean" = "black",
                                "Prior Interpolation" = "grey20",
                                "Prior Estimates" = "grey20"),
                     guide = guide_legend(order = 1)) +
  scale_x_continuous(name = "Year",
                     breaks = seq(min(data$midpoint),
                                  max(data$midpoint), by = 100),
                     minor_breaks = seq(min(data$midpoint),
                                        max(data$midpoint), by = 50),
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Proxy Sum",
                     sec.axis = sec_axis(~ . * 1000,
                                         labels = scales::comma,
                                         name = "Population",
                                         breaks = seq(0, 500000, by = 100000)),
                     limits = c(0, 500), breaks = seq(0, 500, by = 100),
                     expand = c(0, 0)) +
  labs(
    title = "Population of Okayama Prefecture (950 BC–1200 AD)",
    caption = paste0(
      "Proxy weights: Pit = ", proxy_1_weight,
      ", Post = ", proxy_2_weight,
      ", Terrace = ", proxy_3_weight,
      ", Other = ", proxy_4_weight, "\n",
      "Pre-800 AD prior mean = ", sprintf("%.2f%%", past_prior_mean * 100),
      ", σ = ", sprintf("%.2f%%", past_prior_sigma * 100), "\n",
      "Post-800 AD prior mean = ", sprintf("%.2f%%", future_prior_mean * 100),
      ", σ = ", sprintf("%.2f%%", future_prior_sigma * 100)
    )
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "inside",
    legend.justification = c(0.01, 0.99),
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0)
  )

print(plot)

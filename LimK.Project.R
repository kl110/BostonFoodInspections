# DS 340H Final Project
# Modelling
# Kami Lim

library(lme4)
library(lmerTest)
library(MASS)
library(tidyverse)
library(performance)

# Load data
visit_df <- read.csv('/Users/kamilalim/Desktop/Wellesley/Semester_8/DS_340H/Project/boston_visit_df.csv')

cat("Rows:", nrow(visit_df), "\n")
cat("Mean violations per visit:", mean(visit_df$violations_this_visit), "\n")
cat("Variance:", var(visit_df$violations_this_visit), "\n")
# Everything looks the same as the dataset in python!

############# Find missing values information
# Overall missing values per column
colSums(is.na(visit_df))

# Missing income specifically - which tracts?
missing_income <- visit_df[is.na(visit_df$median_income), ]
cat("Visits missing income:", nrow(missing_income), "\n")
cat("Tracts with missing income:", n_distinct(missing_income$tract_num), "\n")
print(table(missing_income$tract_num))

# Missing tract assignment (fell outside tract boundaries in spatial join)
cat("Visits missing tract:", sum(is.na(visit_df$tract_num)), "\n")

# Cross-tab missing income by business type
cat("Missing income by business type:\n")
print(table(visit_df$analysis_group, is.na(visit_df$median_income)))




# Give a unique identifier to each business
visit_df$establishment_id <- paste(visit_df$businessname_clean, visit_df$address_clean, sep="_")

#### Make the analysis group into a factor
visit_df$analysis_group <- factor(
  visit_df$analysis_group,
  levels = c('Independent', 'National Chain', 'Brand Chain', 'Institutional/Owner Chain')
)

# Scale income to $10,000 units for interpretable coefficients
# A 1 unit changes will mean a 10k increase in median houshold income
# Makes output interpretable
visit_df$median_income_10k <- visit_df$median_income / 10000

cat("Reference category:", levels(visit_df$analysis_group)[1], "\n")
cat("Income range: $", min(visit_df$median_income, na.rm=TRUE), 
    "- $", max(visit_df$median_income, na.rm=TRUE), "\n")
print(table(visit_df$analysis_group))


####### MODEL 1: Fit the Poisson Model ##############
# (1 | ...) means we're estimating random itnercept only, not random slope.
  # MEans that every group gets its own base violation rate that shifts up and down
model_poisson_main <- glmer(
  violations_this_visit ~  # outcome variable is the number of violations recorded in 1 visit
    analysis_group + # Fixed effects we want to estimate coeffs for
    median_income_10k + # For every $10k increase in neighborhood income, how does the expected violation count change when analysis group is constant
    (1 | tract_num) + # Some tracts have higher or lower violation rates we can't account for
    (1 | tract_num:establishment_id), # NESTED establishment level random intercept nested within tract: some locations get more violations than others for random reasons
  data = visit_df,
  family = poisson(link = "log")
)

summary(model_poisson_main)


####### MODEL 2: Fit the Poisson Mixed Effects Interaction Model ##############
model_poisson_interaction <- glmer(
  violations_this_visit ~ 
    analysis_group * median_income_10k +  # * = main effects + all interactions
    (1 | tract_num) +
    (1 | tract_num:establishment_id),
  data = visit_df,
  family = poisson(link = "log")
)

summary(model_poisson_interaction)
# AIC: 83653.0
# Only brand chain x 


########## Likelihood Ratio test
anova(model_poisson_main, model_poisson_interaction)
# Chisq: 7.4872  DF: 3    P-value: 0.05789 
# No significant difference between models
# Interaction doesn't better explain


# Check for overdispersion
check_overdispersion(model_poisson_main)
# Poisson thing: poisson has mean and variance equal to each other
# Overdispersion, points are more varied 
# Poisson regression chapter on overdispersion 
# Statistical sleuth
# Check how to look at poisson coefficients
# Try to add an overdispersion argument = true argument in poisson
# Decide with overdispersion or decide to run poisson with overdispersion, allowed in glmer function
# Other option is to run negative binomial
# Ignoring spatial correlation

####### MODEL 3: Negative Binomial Model ##############
library(glmmTMB)

model_nb_main <- glmmTMB(
  violations_this_visit ~
    analysis_group +
    median_income_10k +
    (1 | tract_num) +
    (1 | tract_num:establishment_id),
  data = visit_df,
  family = nbinom2(link = "log")  # nbinom2 is the standard negative binomial
)

model_nb_interaction <- glmmTMB(
  violations_this_visit ~
    analysis_group * median_income_10k +
    (1 | tract_num) +
    (1 | tract_num:establishment_id),
  data = visit_df,
  family = nbinom2(link = "log")
)

summary(model_nb_main)
summary(model_nb_interaction)
model_performance(model_nb_interaction)
anova(model_nb_main, model_nb_interaction)
# Interaction model represents data better
# Statistically significant

# Finding results
# Extract coefficients
coef_summary <- summary(model_nb_interaction)$coefficients$cond

# Create interpretation table
interpret_df <- data.frame(
  Estimate = coef_summary[, "Estimate"],
  IRR = exp(coef_summary[, "Estimate"]),
  Percent_Change = (exp(coef_summary[, "Estimate"]) - 1) * 100,
  P_Value = coef_summary[, "Pr(>|z|)"]
)

print(round(interpret_df, 4))


##################
############################################# PLOTS TO VISUALIZE THE MODEL PREDICTIONS
# ============================================================
# LOAD VISUALIZATION LIBRARIES
# ============================================================
library(ggplot2)
library(dplyr)
library(broom.mixed)  # extracts tidy coefficients from mixed effects models

# ============================================================
# VISUALIZATION 2: COEFFICIENT PLOT (FOREST PLOT)
# ============================================================
# Visualization 1 is the updated interactive map in python code
# Extract coefficients from both NB models
# tidy() from broom.mixed gives us estimates, std errors, and confidence intervals
# exponentiate=TRUE converts log coefficients to Incidence Rate Ratios (IRR)
# conf.int=TRUE adds 95% confidence intervals
# effects="fixed" extracts only fixed effects, not random effects

coef_main <- tidy(model_nb_main, 
                  effects="fixed", 
                  exponentiate=TRUE, 
                  conf.int=TRUE) %>%
  mutate(model = "Main Effects")

coef_interaction <- tidy(model_nb_interaction, 
                         effects="fixed", 
                         exponentiate=TRUE, 
                         conf.int=TRUE) %>%
  mutate(model = "Interaction")

# Combine both models
coef_combined <- bind_rows(coef_main, coef_interaction)

# Clean up term names for readable labels
# The raw names from glmmTMB are verbose so we replace them
coef_combined <- coef_combined %>%
  filter(term != "(Intercept)") %>%  # intercept isn't substantively interesting
  mutate(term_clean = case_when(
    term == "analysis_groupNational Chain" ~ "National Chain",
    term == "analysis_groupBrand Chain" ~ "Brand Chain",
    term == "analysis_groupInstitutional/Owner Chain" ~ "Institutional/Owner Chain",
    term == "median_income_10k" ~ "Income ($10k)",
    term == "analysis_groupNational Chain:median_income_10k" ~ "National Chain × Income",
    term == "analysis_groupBrand Chain:median_income_10k" ~ "Brand Chain × Income",
    term == "analysis_groupInstitutional/Owner Chain:median_income_10k" ~ "Institutional/Owner × Income",
    TRUE ~ term
  )) %>%
  # Set order for y-axis display
  mutate(term_clean = factor(term_clean, levels = c(
    "Income ($10k)",
    "Institutional/Owner Chain",
    "Brand Chain",
    "National Chain",
    "Institutional/Owner × Income",
    "Brand Chain × Income",
    "National Chain × Income"
  )))

# Separate main effects and interaction terms for cleaner display
coef_main_only <- coef_combined %>%
  filter(model == "Main Effects")

coef_interaction_only <- coef_combined %>%
  filter(model == "Interaction")

# ---- Plot 2a: Main Effects Model ----
# Shows IRRs for business type and income
# IRR > 1 = more violations than Independent
# IRR < 1 = fewer violations than Independent
# Reference line at 1.0 = no difference from Independent

p_main <- ggplot(
  coef_main_only %>% filter(term_clean %in% c("National Chain", "Brand Chain", 
                                              "Institutional/Owner Chain", "Income ($10k)")),
  aes(x = estimate, y = term_clean)
) +
  # Point for the estimate
  geom_point(size = 3, color = "#2c7bb6") +
  # Horizontal line for 95% CI
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, color = "#2c7bb6", linewidth = 0.8) +
  # Reference line at IRR = 1 (no difference from Independent)
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  # Add IRR value labels next to each point
  geom_text(aes(label = sprintf("IRR = %.2f", estimate)), 
            hjust = -0.2, size = 3.5) +
  scale_x_continuous(
    limits = c(0.5, 1.6),
    name = "Incidence Rate Ratio (IRR)\n(relative to Independent businesses)"
  ) +
  labs(
    title = "Negative Binomial Mixed Effects Model: Main Effects",
    subtitle = "IRR < 1 = fewer violations than Independent | IRR > 1 = more violations than Independent",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.text.y = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

print(p_main)
ggsave('/Users/kamilalim/Desktop/Wellesley/Semester_8/DS_340H/Project/coef_plot_main.png',
       p_main, width=9, height=5, dpi=150)

# ---- Plot 2b: Interaction Model ----
# Shows all terms including interactions
# Separated into two panels: main effects and interaction terms

coef_int_plot <- coef_interaction_only %>%
  mutate(
    term_type = ifelse(grepl("×", term_clean), "Interaction Terms", "Main Effects"),
    significant = p.value < 0.05  # flag significant terms
  )

p_interaction <- ggplot(coef_int_plot, aes(x = estimate, y = term_clean, color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", estimate)), 
            hjust = -0.3, size = 3.2) +
  scale_color_manual(
    values = c("TRUE" = "#d7191c", "FALSE" = "#2c7bb6"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
    name = "Significance"
  ) +
  scale_x_continuous(
    limits = c(0.4, 1.8),
    name = "Incidence Rate Ratio (IRR)\n(relative to Independent businesses)"
  ) +
  # Split into two panels by term type
  facet_wrap(~term_type, scales = "free_y", ncol = 1) +
  labs(
    title = "Negative Binomial Mixed Effects Model: Interaction Model",
    subtitle = "Red = significant (p < 0.05) | Blue = not significant",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.text.y = element_text(size = 11),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

print(p_interaction)
ggsave('/Users/kamilalim/Desktop/Wellesley/Semester_8/DS_340H/Project/coef_plot_interaction.png',
       p_interaction, width=9, height=7, dpi=150)

# ============================================================
# VISUALIZATION 3: PREDICTED VIOLATIONS BY INCOME AND BUSINESS TYPE
# ============================================================
# This plot shows the model's predicted violation counts
# across the observed income range for each business type
# It makes the interaction tangible — you can see Brand Chain's
# steeper slope compared to the other groups

# Create a prediction grid
# We vary income across its observed range
# and compute predictions for each business type
income_range <- seq(
  min(visit_df$median_income_10k, na.rm=TRUE),
  max(visit_df$median_income_10k, na.rm=TRUE),
  length.out = 100  # 100 evenly spaced income values
)

# One row per combination of income value and business type
pred_grid <- expand.grid(
  median_income_10k = income_range,
  analysis_group = levels(visit_df$analysis_group)
)
pred_grid$analysis_group <- factor(pred_grid$analysis_group,
                                   levels = levels(visit_df$analysis_group))

# Generate predictions from the NB interaction model
# re.form = NA means we set all random effects to 0
# (predicting for an average establishment in an average tract)
# type = "response" back-transforms from log scale to violation counts
pred_grid$predicted <- predict(
  model_nb_interaction,
  newdata = pred_grid,
  re.form = NA,       # ignore random effects — population-level prediction
  type = "response"   # return predictions on count scale, not log scale
)

# Add confidence intervals using the bootstrap method
# This gives us uncertainty bands around the predicted lines
# nsim = 500 bootstrap samples for stable CIs
# (reduce to 100 if this is slow)
pred_ci <- predict(
  model_nb_interaction,
  newdata = pred_grid,
  re.form = NA,
  type = "response",
  se.fit = TRUE       # return standard errors for CI calculation
)
pred_grid$se <- pred_ci$se.fit
pred_grid$ci_low <- pred_grid$predicted - 1.96 * pred_grid$se
pred_grid$ci_high <- pred_grid$predicted + 1.96 * pred_grid$se

# Convert income back to dollars for readable x-axis
pred_grid$median_income_dollar <- pred_grid$median_income_10k * 10000

# Add observed data points (establishment-level means) for context
obs_points <- visit_df %>%
  group_by(establishment_id, analysis_group, median_income) %>%
  summarise(mean_vpv = mean(violations_this_visit), 
            n_visits = n(), 
            .groups = 'drop') %>%
  filter(n_visits >= 3)  # only show establishments with 3+ visits for reliability

# Define consistent colors for business types
group_colors <- c(
  "Independent" = "#2ca02c",
  "National Chain" = "#1f77b4",
  "Brand Chain" = "#ff7f0e",
  "Institutional/Owner Chain" = "#d62728"
)

p_predicted <- ggplot() +
  # Confidence band around predicted line
  geom_ribbon(
    data = pred_grid,
    aes(x = median_income_dollar, 
        ymin = ci_low, ymax = ci_high, 
        fill = analysis_group),
    alpha = 0.15  # semi-transparent so lines are visible
  ) +
  # Predicted lines — one per business type
  geom_line(
    data = pred_grid,
    aes(x = median_income_dollar, y = predicted, color = analysis_group),
    linewidth = 1.2
  ) +
  # Observed establishment means as jittered points
  # Alpha and size scaled by number of visits
  geom_point(
    data = obs_points,
    aes(x = median_income, y = mean_vpv, color = analysis_group),
    alpha = 0.2, size = 1.5
  ) +
  scale_color_manual(values = group_colors, name = "Business Type") +
  scale_fill_manual(values = group_colors, name = "Business Type") +
  scale_x_continuous(
    labels = scales::dollar_format(scale = 0.001, suffix = "k"),
    name = "Census Tract Median Household Income"
  ) +
  scale_y_continuous(
    name = "Predicted Violations Per Visit",
    limits = c(0, 12)  # adjust if predictions go higher
  ) +
  labs(
    title = "Predicted Violations Per Visit by Business Type and Neighborhood Income",
    subtitle = "Lines = model predictions | Shading = 95% CI | Points = observed establishment means (3+ visits)",
    caption = "Negative binomial mixed effects model, random effects set to population average"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    plot.caption = element_text(size = 9, color = "grey50"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_predicted)
ggsave('/Users/kamilalim/Desktop/Wellesley/Semester_8/DS_340H/Project/predicted_violations_plot.png',
       p_predicted, width=10, height=6, dpi=150)











##############################
# ============================================================
# SIMPLE VERSION: Main Effects Plot
# ============================================================
# Just show the two significant findings clearly

simple_coef <- data.frame(
  business_type = c("National Chain", "Brand Chain", 
                    "Institutional/Owner Chain", "Income\n(per $10k)"),
  irr = c(0.69, 1.02, 0.95, 0.99),
  ci_low = c(0.645, 0.949, 0.859, 0.984),
  ci_high = c(0.738, 1.094, 1.050, 0.999),
  significant = c(TRUE, FALSE, FALSE, TRUE)
)

simple_coef$business_type <- factor(simple_coef$business_type,
                                    levels = c("Income\n(per $10k)", "Institutional/Owner Chain", 
                                               "Brand Chain", "National Chain"))

ggplot(simple_coef, aes(x = irr, y = business_type, color = significant)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("TRUE" = "#d7191c", "FALSE" = "grey60"),
                     guide = "none") +
  annotate("text", x = 0.88, y = 4.4, label = "Fewer violations →", 
           size = 3.5, color = "grey40") +
  annotate("text", x = 1.12, y = 4.4, label = "← More violations", 
           size = 3.5, color = "grey40") +
  scale_x_continuous(limits = c(0.6, 1.2),
                     name = "Violation Rate Ratio\n(compared to Independent businesses)") +
  labs(
    title = "Which businesses have more or fewer violations?",
    subtitle = "Red = statistically significant finding | Grey = no significant difference",
    y = NULL,
    caption = "Reference group: Independent businesses | Adjusted for neighborhood income"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave('/Users/kamilalim/Desktop/Wellesley/Semester_8/DS_340H/Project/simple_coef_plot.png',
       width = 8, height = 4.5, dpi = 150)

# ============================================================
# SIMPLE VERSION: Predicted Violations Plot
# ============================================================
# Cleaner version without observed points — just the lines

ggplot(pred_grid, aes(x = median_income_dollar, y = predicted, color = analysis_group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = analysis_group),
              alpha = 0.1, color = NA) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(
    labels = scales::dollar_format(scale = 0.001, suffix = "k"),
    name = "Neighborhood Median Household Income"
  ) +
  scale_y_continuous(
    name = "Predicted Violations Per Inspection",
    limits = c(0, 7)
  ) +
  labs(
    #title = "Brand chain restaurants have more violations\nin lower-income neighborhoods",
    #subtitle = "Predicted violations per inspection visit by business type and neighborhood income",
    title = "Predicted violations per inspection visit by business type and neighborhood income",
    #caption = "Negative binomial mixed effects model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave('/Users/kamilalim/Desktop/Wellesley/Semester_8/DS_340H/Project/simple_predicted_plot.png',
       width = 9, height = 5.5, dpi = 150)

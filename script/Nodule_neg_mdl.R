
df_model$GM_label <- ifelse(df_model$GM_peak >=0.5 ,1,0)
vars <- c("Size_nodule","Nodule_CT" ,"Number_nodule", "GM_label","WBC_final", "CRP", "ESR", "age_adjusted")
cat("Before dropping NAs:", nrow(df_model), "\n")

df_gm <- df_model[complete.cases(df_model[, vars]), ] #114 , one has na for size

library(MASS)

# Fit negative binomial model
nb_model <- glm.nb(Number_nodule ~ GM_value_1+ CRP + ESR+age_adjusted + WBC_final, data = df_gm)

# Summary output
summary(nb_model)
exp(coef(nb_model))
exp(confint(nb_model))

library(ggplot2)
library(dplyr)

# Extract coefficients and CIs
coef_table <- summary(nb_model)$coefficients
irr <- exp(coef(nb_model))
ci <- exp(confint(nb_model))

# Create data frame for plotting (exclude intercept)
plot_data <- data.frame(
    Predictor = rownames(ci)[-1],  # Remove intercept
    IRR = irr[-1],
    Lower = ci[-1, 1],
    Upper = ci[-1, 2]
) %>%
    mutate(
        Predictor = recode(Predictor,
                           "GM_value_1" = "GM Value",
                           "CRP" = "CRP",
                           "ESR" = "ESR",
                           "age_adjusted" = "Age (adjusted)",
                           "WBC_final" = "WBC Count"),
        Predictor = factor(Predictor, levels = rev(Predictor))  # Reverse order for plot
    )

# Create the forest plot
ggplot(plot_data, aes(x = IRR, y = Predictor)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    labs(
        title = "Negative Binomial Regression: Incidence Rate Ratios",
        x = "Incidence Rate Ratio (95% CI)",
        y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 11),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )









#-------------------------------------------------------------------------------
si_model <- glm.nb(Size_nodule ~ GM_value_1 + CRP +ESR +age_adjusted + WBC_final, data = df_gm)

# Summary output
summary(si_model)
exp(coef(si_model))
exp(confint(si_model))

# Or use dispersion statistic
disp_stat <- sum(residuals(nb_model, type="pearson")^2) / nb_model$df.residual
# 0.9636665
disp_stat_s <- sum(residuals(si_model, type="pearson")^2) / si_model$df.residual
#0.5430842
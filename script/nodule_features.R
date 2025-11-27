
vars <- c("Size_nodule","Nodule_CT" ,"Number_nodule", "GM_label","WBC_final", "CRP", "ESR", "age_adjusted")
cat("Before dropping NAs:", nrow(df_model), "\n")

df_gm <- df_model[complete.cases(df_model[, vars]), ] #114 , one has na for size
library(MASS)
library(MASS)

# Fit negative binomial model
nb_model <- glm.nb(Number_nodule ~ GM_value_1+ CRP + ESR+age_adjusted + WBC_final, data = df_gm)

# Summary output
summary(nb_model)
exp(coef(nb_model))
exp(confint(nb_model))
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
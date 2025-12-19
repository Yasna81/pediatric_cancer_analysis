
library(logistf)
library(pROC)
library(caret)
library(dplyr)
library(ggplot2)
library(ResourceSelection)
library(car)
# --- 1. Clean data (drop rows with any missing predictor/outcome) ---
vars <- c("Nodule_CT", "WBC_final", "CRP", "ESR", "age_adjusted")
cat("Before dropping NAs:", nrow(df_model), "\n")

df_final <- df_model[complete.cases(df_model[, vars]), ] #115

df_final$nodule_label <- ifelse(df_final$Nodule_CT == "Yes", 1 , 0)
df_final$nodule_label <- factor(df_final$nodule_label, levels = c(0,1), labels = c("No","Yes"))
# --- 2. 5-fold STRATIFIED cross-validation (this is correct now) ---
set.seed(123)
folds <- createFolds(df_final$nodule_label, k = 5, list = TRUE, returnTrain = FALSE)

auc_cv <- sens_cv <- spec_cv <- f1_cv <- acc_cv <- numeric(5)

for(i in 1:5) {
    train <- df_final[-folds[[i]], ]
    test  <- df_final[ folds[[i]], ]
    
    model <- logistf(nodule_label ~ WBC_final + CRP + ESR + age_adjusted, data = train)
    
    probs <- predict(model, newdata = test, type = "response")
    
    # AUC
    roc_obj <- roc(test$nodule_label, probs, quiet = TRUE)
    auc_cv[i] <- auc(roc_obj)
    
    # Optimal threshold (Youden)
    optimal <- coords(roc_obj, "best", ret = "threshold")$threshold
    preds   <- factor(ifelse(probs > optimal, "Yes", "No"), levels = c("No","Yes"))
    
    cm <- confusionMatrix(preds, test$nodule_label, positive = "Yes")
    
    acc_cv[i]  <- cm$overall["Accuracy"]
    sens_cv[i] <- cm$byClass["Sensitivity"]
    spec_cv[i] <- cm$byClass["Specificity"]
    f1_cv[i]   <- cm$byClass["F1"]
}

cat("=== 5-FOLD STRATIFIED CV RESULTS ===\n")
cat("Mean AUC         :", round(mean(auc_cv), 3), "±", round(sd(auc_cv), 3), "\n")
cat("Mean Sensitivity :", round(mean(sens_cv), 3), "±", round(sd(sens_cv), 3), "\n")
cat("Mean Specificity :", round(mean(spec_cv), 3), "±", round(sd(spec_cv), 3), "\n")
cat("Mean F1-score    :", round(mean(f1_cv), 3), "±", round(sd(f1_cv), 3), "\n")
cat("Mean Accuracy    :", round(mean(acc_cv), 3), "±", round(sd(acc_cv), 3), "\n\n")


#> cat("=== 5-FOLD STRATIFIED CV RESULTS ===\n")
#=== 5-FOLD STRATIFIED CV RESULTS ===
#    > cat("Mean AUC         :", round(mean(auc_cv), 3), "±", round(sd(auc_cv), 3), "\n")
#Mean AUC         : 0.727 ± 0.077 
#> cat("Mean Sensitivity :", round(mean(sens_cv), 3), "±", round(sd(sens_cv), 3), "\n")
#Mean Sensitivity : 0.707 ± 0.192 
#> cat("Mean Specificity :", round(mean(spec_cv), 3), "±", round(sd(spec_cv), 3), "\n")
#Mean Specificity : 0.65 ± 0.185 
#> cat("Mean F1-score    :", round(mean(f1_cv), 3), "±", round(sd(f1_cv), 3), "\n")
#Mean F1-score    : 0.739 ± 0.158 
#> cat("Mean Accuracy    :", round(mean(acc_cv), 3), "±", round(sd(acc_cv), 3), "\n\n")
#Mean Accuracy    : 0.687 ± 0.167 

# --- 3. Final model on full clean data ---
final_model <- logistf(nodule_label ~ WBC_final + CRP + ESR + age_adjusted, data = df_final)

cat("=== FINAL MODEL ODDS RATIOS (95% CI) ===\n")
ORs <- exp(cbind(OR = coef(final_model), confint(final_model)))
print(round(ORs, 3))

# --- 4. Goodness-of-fit on final model (Hosmer-Lemeshow) ---
hl <- hoslem.test(as.numeric(df_final$nodule_label)-1, predict(final_model, type="response"), g=10)
cat("\nHosmer-Lemeshow goodness-of-fit p-value:", round(hl$p.value, 4), "\n")

# --- 5. VIF (on standard glm – Firth doesn’t give VIF) ---
glm_model <- glm(nodule_label ~ WBC_final + CRP + ESR + age_adjusted, data = df_final, family = binomial)
cat("\nVIF (tolerance check):\n")
print(round(vif(glm_model), 2))

summary(final_model)

#
df_final$cv_prob <- NA

for(i in 1:5) {
    train <- df_final[-folds[[i]], ]
    test  <- df_final[ folds[[i]], ]
    
    model <- logistf(nodule_label ~ WBC_final + CRP + ESR + age_adjusted, data = train)
    
    df_final$cv_prob[folds[[i]]] <- predict(model, newdata = test, type = "response")
}

library(pROC)

roc_cv <- roc(df_final$nodule_label, df_final$cv_prob)

plot(auc_cv, col = "blue", lwd = 3, main = "5-Fold Cross-Validated ROC Curve")
text(0.6, 0.2, paste("AUC =", round(auc(roc_cv), 3)), cex = 1.2)




# Create a data frame manually
coef_table <- data.frame(
    term = names(final_model$coefficients),
    estimate = final_model$coefficients,
    lower = final_model$ci.lower,
    upper = final_model$ci.upper,
    p_value = final_model$prob
)

# Remove the intercept if you don't want to plot it
coef_table <- coef_table[coef_table$term != "(Intercept)", ]

# Compute Odds Ratios
coef_table$OR <- exp(coef_table$estimate)
coef_table$CI_lower <- exp(coef_table$lower)
coef_table$CI_upper <- exp(coef_table$upper)

coef_table$term <-recode(coef_table$term,
                         "age_adjusted" = "Age (years)",
                         "WBC_final" = "WBC",
                         "CRP" = "CRP",
                         "ESR" = "ESR")
library(ggplot2)

ggplot(coef_table, aes(x = term, y = OR)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() + # horizontal layout
    labs(
        title = "Odds Ratios with 95% CI",
        x = "Predictor",
        y = "Odds Ratio (log scale)"
    ) +
    scale_y_log10() + # log scale for ORs
    theme_classic()


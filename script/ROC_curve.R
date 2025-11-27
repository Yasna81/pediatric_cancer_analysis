



# Empty plot
plot(NULL, xlim = c(1, 0), ylim = c(0, 1), 
     xlab = "1 - Specificity (False Positive Rate)", 
     ylab = "Sensitivity (True Positive Rate)",
     main = "ROC Curves from 5-fold Stratified Cross-Validation",
     asp = 1, las = 1)

# Colors + storage
colors <- rainbow(5)
auc_values <- numeric(5)

# Draw each fold
for(i in 1:5) {
    train <- df_final[-folds[[i]], ]
    test  <- df_final[folds[[i]], ]
    
    model <- logistf(nodule_label ~ WBC_final + CRP + ESR + age_adjusted, data = train)
    probs <- predict(model, newdata = test, type = "response")
    roc_obj <- roc(test$nodule_label, probs, quiet = TRUE)
    
    auc_values[i] <- round(as.numeric(auc(roc_obj)), 3)
    
    lines(roc_obj, col = colors[i], lwd = 2.5)
}

# Diagonal reference line
abline(a = 0, b = 1, lty = 2, col = "gray60", lwd = 1.5)

# Legend with AUC values only
legend("bottomright",
       legend = paste0("Fold ", 1:5, "  (AUC = ", auc_values, ")"),
       col = colors,
       lwd = 2.5,
       cex = 1.05,
       bty = "n",
       title = "5-fold CV results",
       title.adj = 0.2)

# Mean ± SD as prominent text (top-left or bottom-left — choose what you prefer)
mean_auc <- round(mean(auc_values), 3)
sd_auc   <- round(sd(auc_values), 3)

text(x = 0.98, y = 0.08, 
     labels = paste0("Mean AUC = ", mean_auc, " ± ", sd_auc),
     cex = 1.4, font = 2, adj = 1, col = "black")
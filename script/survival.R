#survival analysis 
df_8$last_follow <- as.Date("2002-10-29")
df_8$surve_time <- as.numeric(
    ifelse(df_8$Outcome == "Dead" ,
           df_8$Outcome_date - df_8$Admission_date,
           df_8$last_follow - df_8$Admission_date)
)
library(survival)
df_9 <- df_8[!is.na(df_8$surve_time) & !is.na(df_8$Outcome),]
df_9 <- df_9[df_9$surve_time >=0,]
df_9$outcom_num <- ifelse(df_9$Outcome == "Dead" ,1,0)
surv_obj <- Surv(time = df_9$surve_time, event = df_9$outcom_num)
km_fit <- survfit((surv_obj) ~ 1  , data = df_9 )
summary(km_fit)
sum(km_fit$n.risk)
sum(km_fit$n.event) #41
sum(km_fit$n.censor)#95 #one is dead before admission
plot(km_fit,
     main = "Kaplan-Meier Curve",
     xlab = "Time",
     ylab = "Survival Probability",
     lwd = 2,
     col = "blue")


summary(km_fit, times = c(7, 14,28,60,90,2660))

library(survminer)
ggsurvplot(km_fit,
           data=df_9,
           risk.table = TRUE,
           conf.int = TRUE,
           xlab = "Days since admission",
           ylab = "Survival probability",
           palette = "Dark2")
#different groups :
survdiff((surv_obj) ~ Nodule_CT, data = df_9,na.action = na.omit) #Chisq= 0.2  on 1 degrees of freedom, p= 0.6 
survdiff((surv_obj) ~ EORCT_MSG_criteria, data = df_9) #Chisq= 0.8  on 1 degrees of freedom, p= 0.3
survdiff((surv_obj) ~ Relapse , data = df_9)#Chisq= 0.2  on 1 degrees of freedom, p= 0.7
survdiff((surv_obj) ~ Atypic_positive_1 , data = df_9)#, 0.2 ,p value is 1



cox_model_2 <- coxph((surv_obj) ~ CRP + ESR + WBC_final + GM_peak + age_adjusted  , data = df_9 )
summary(cox_model_2)

cox_model_6 <- coxph((surv_obj) ~  Number_nodule+Size_nodule , data = df_9 )
summary(cox_model_6)

cox_model_ct <- coxph((surv_obj) ~ Nodule_CT + Halo_CT + GGO_CT + Cavity_CT , data = df_9)
summary(cox_model_ct)
library(car)
vif(cox_model_2) #no multicollinearity 





#  Load required packages
library(survival)
library(rms)
library(ggplot2)


#  Set up for rms package
dd <- datadist(df_9)
options(datadist = "dd")

#  Fit Cox model with restricted cubic spline on GM
# You can add other covariates if needed (e.g., + age + CRP)
fit <- cph(surv_obj ~ rcs(GM_peak, 4), data = df_9, x = TRUE, y = TRUE, surv = TRUE)

# Predict hazard ratios across GM range
pred <- Predict(fit, GM_peak, fun = exp)

#  Plot with ggplot2
ggplot(pred, aes(x = GM_peak, y = yhat)) +
    geom_line(color = "#2C3E50", linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#3498DB") +
    labs(
        title = "Hazard Ratio vs  maximum GM (Galactomannan)",
        x = "GM Value (ODI)",
        y = "Hazard Ratio"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 14)

#ESR


#  Fit Cox model with restricted cubic spline on GM
# You can add other covariates if needed (e.g., + age + CRP)
fit <- cph(surv_obj ~ rcs(ESR, 4), data = df_9, x = TRUE, y = TRUE, surv = TRUE)

# Predict hazard ratios across GM range
pred <- Predict(fit, ESR, fun = exp)

#  Plot with ggplot2
ggplot(pred, aes(x = ESR, y = yhat)) +
    geom_line(color = "#2C3E50", linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#3498DB") +
    labs(
        title = "Hazard Ratio vs  ESR (Erythrocyte Sedimentation Rate) ",
        x = "ESR Value (mm/hr)",
        y = "Hazard Ratio"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 14)

#crp


#  Fit Cox model with restricted cubic spline on GM
# You can add other covariates if needed (e.g., + age + CRP)
fit <- cph(surv_obj ~ rcs(CRP, 4), data = df_9, x = TRUE, y = TRUE, surv = TRUE)

# Predict hazard ratios across GM range
pred <- Predict(fit, CRP, fun = exp)

#  Plot with ggplot2
ggplot(pred, aes(x = CRP, y = yhat)) +
    geom_line(color = "#2C3E50", linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#3498DB") +
    labs(
        title = "Hazard Ratio vs  CRP (C-Reactive Protien)",
        x = "CRP Value (mg/L)",
        y = "Hazard Ratio"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 14)

#age 


#  Fit Cox model with restricted cubic spline on GM
# You can add other covariates if needed (e.g., + age + CRP)
fit <- cph(surv_obj ~ rcs(age_adjusted, 4), data = df_9, x = TRUE, y = TRUE, surv = TRUE)

# Predict hazard ratios across GM range
pred <- Predict(fit, age_adjusted , fun = exp)

#  Plot with ggplot2
ggplot(pred, aes(x = GM_peak, y = yhat)) +
    geom_line(color = "#2C3E50", linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#3498DB") +
    labs(
        title = "Hazard Ratio vs Age",
        x = "Age (years)",
        y = "Hazard Ratio"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 14)

#wbc


#  Fit Cox model with restricted cubic spline on GM
# You can add other covariates if needed (e.g., + age + CRP)
fit <- cph(surv_obj ~ rcs(WBC_final, 4), data = df_9, x = TRUE, y = TRUE, surv = TRUE)

# Predict hazard ratios across GM range
pred <- Predict(fit, WBC_final, fun = exp)

#  Plot with ggplot2
ggplot(pred, aes(x = WBC_1, y = yhat)) +
    geom_line(color = "#2C3E50", linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#3498DB") +
    labs(
        title = "Hazard Ratio vs  WBC (White Blood Cell Count)",
        x = "WBC Value (10^9/L)",
        y = "Hazard Ratio"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 14)
#general info of the cohort :
table(df_7$Fungal_infection_site)
table(df_7$Fungus_species)
table(df_7$Transplant)
table(df_7$Location_CT)
summary(df_7$age_adjusted)
age_table <- df_7 %>% group_by(Sex) %>%
    summarise("n of patients" = n(),
              "mean" = mean(age_adjusted),
              "median" = median(age_adjusted))
table(df_7$Cancer_Group)

#----------------------------------clinical-eda-----------------------------------
library(dplyr)
summary_EORTC_0 <- df_7 %>%
    group_by(EORCT_MSG_criteria) %>%
    summarise("number of patients"  = n() ,
              "median WBC" = median(WBC_final,na.rm = TRUE),
              "Q1 WBC" = quantile(WBC_final,probs = 0.25 ,na.rm = TRUE),
              "Q2 WBC" = quantile(WBC_final,probs = 0.75 ,na.rm = TRUE),
              "median ESR" = median(ESR,na.rm = TRUE),
              "Q1 ESR"= quantile(ESR,probs = 0.25 ,na.rm = TRUE),
              "Q2 ESR" = quantile(ESR,probs = 0.75 ,na.rm = TRUE),
              "median CRP" = median(CRP,na.rm = TRUE),
              "iqr 1-crp" = quantile(CRP,probs = 0.25 ,na.rm = TRUE),
              "iqr 2-crp" = quantile(CRP,probs = 0.75 ,na.rm = TRUE),
              "median GM"= median(GM_peak,na.rm = TRUE),
              "iqr 1" = quantile(GM_peak,probs = 0.25,na.rm= TRUE),
              "iqr 2" = quantile(GM_peak,probs = 0.75 ,na.rm = TRUE))

#stat test :
kruskal_wbc <-kruskal.test(WBC_final ~ EORCT_MSG_criteria, data = df_7) #p-value = 0.8849
kruskal_CRP <-kruskal.test(CRP ~ EORCT_MSG_criteria, data = df_7) #p-value = 0.6554
kruskal_ESR <-kruskal.test(ESR ~ EORCT_MSG_criteria, data = df_7) #p-value = 0.335
kruskal_GM <-kruskal.test(GM_peak ~EORCT_MSG_criteria, data = df_7)#p-value = 0.5743
#df_4 has possible cases>>
wilcox_df1 <- df_4 %>% filter(EORCT_MSG_criteria!="NA" & EORCT_MSG_criteria != "False positive" & WBC_final!= "NA")
library(dplyr)
summary_EORTC <- wilcox_df1 %>%
    group_by(EORCT_MSG_criteria) %>%
    summarise("number of patients"  = n() ,
              "median WBC" = median(WBC_final,na.rm = TRUE),
              "Q1 WBC" = quantile(WBC_final,probs = 0.25 ,na.rm = TRUE),
              "Q2 WBC" = quantile(WBC_final,probs = 0.75 ,na.rm = TRUE),
              "median ESR" = median(ESR,na.rm = TRUE),
              "Q1 ESR"= quantile(ESR,probs = 0.25 ,na.rm = TRUE),
              "Q2 ESR" = quantile(ESR,probs = 0.75 ,na.rm = TRUE),
              "median CRP" = median(CRP,na.rm = TRUE),
              "iqr 1-crp" = quantile(CRP,probs = 0.25 ,na.rm = TRUE),
              "iqr 2-crp" = quantile(CRP,probs = 0.75 ,na.rm = TRUE),
              "median GM"= median(GM_peak,na.rm = TRUE),
              "iqr 1" = quantile(GM_peak,probs = 0.25,na.rm= TRUE),
              "iqr 2" = quantile(GM_peak,probs = 0.75 ,na.rm = TRUE))

kruskal_wbc <-kruskal.test(WBC_final ~ EORCT_MSG_criteria, data = wilcox_df1) #p-value = 0.2722
kruskal_CRP <-kruskal.test(CRP ~ EORCT_MSG_criteria, data = wilcox_df1) #p-value = 0.142
kruskal_ESR <-kruskal.test(ESR ~ EORCT_MSG_criteria, data = wilcox_df1) #p-value = 0.05722
kruskal_GM <-kruskal.test(GM_peak ~EORCT_MSG_criteria, data = wilcox_df1)#p-value = 0.009633 :)))
#df_4 has possible cases>>

##
df_7$FN_episodes <- as.character(df_7$FN_episodes)
df_7$FN_episodes <- ifelse(tolower(df_7$FN_episodes) == "yes",1,df_4$FN_episodes)
df_7$FN_episodes <- as.numeric(df_7$FN_episodes)

summary_table_FN <- df_7 %>%
    group_by(FN) %>%
    summarise(
        "number of patients"  = n() ,
        "Median FN episodes" = median(FN_episodes,na.rm = TRUE),
        "IQR 1 FN episodes" = quantile(FN_episodes,probs = 0.25 ,na.rm = TRUE),
        "IQR 2 FN episodes" = quantile(FN_episodes,probs = 0.75 ,na.rm = TRUE),
        "Min FN episodes" = min(as.numeric(FN_episodes)),
        "Max FN episodes" = max(as.numeric(FN_episodes)),
        "median WBC" = median(WBC_final,na.rm = TRUE),
        "Q1 WBC" = quantile(WBC_final,probs = 0.25 ,na.rm = TRUE),
        "Q2 WBC" = quantile(WBC_final,probs = 0.75 ,na.rm = TRUE),
        "median ESR" = median(ESR,na.rm = TRUE),
        "Q1 ESR"= quantile(ESR,probs = 0.25 ,na.rm = TRUE),
        "Q2 ESR" = quantile(ESR,probs = 0.75 ,na.rm = TRUE),
        "median CRP" = median(CRP,na.rm = TRUE),
        "iqr 1-crp" = quantile(CRP,probs = 0.25 ,na.rm = TRUE),
        "iqr 2-crp" = quantile(CRP,probs = 0.75 ,na.rm = TRUE),
        "median GM"= median(GM_peak,na.rm = TRUE),
        "iqr 1" = quantile(GM_peak,probs = 0.25,na.rm= TRUE),
        "iqr 2" = quantile(GM_peak,probs = 0.75 ,na.rm = TRUE))


kruskal.test(WBC_final ~ FN,data= df_7) #p-value = 0.08
kruskal.test(CRP ~ FN,data= df_7) #P_value = 0.07
kruskal.test(ESR ~ FN,data= df_7) #p-value = 0.08
kruskal.test(GM_peak ~ FN , data= df_7) #p-value 0.3
##--------------------------------------

library(dplyr)
Nodule_table <- df_7 %>%
    filter(!is.na(Nodule_CT)) %>%
    group_by(Nodule_CT) %>%
    summarise(
        "number of patients"  = n() ,
        "median WBC" = median(WBC_final,na.rm = TRUE),
        "Q1 WBC" = quantile(WBC_final,probs = 0.25 ,na.rm = TRUE),
        "Q2 WBC" = quantile(WBC_final,probs = 0.75 ,na.rm = TRUE),
        "median ESR" = median(ESR,na.rm = TRUE),
        "Q1 ESR"= quantile(ESR,probs = 0.25 ,na.rm = TRUE),
        "Q2 ESR" = quantile(ESR,probs = 0.75 ,na.rm = TRUE),
        "median CRP" = median(CRP,na.rm = TRUE),
        "iqr 1-crp" = quantile(CRP,probs = 0.25 ,na.rm = TRUE),
        "iqr 2-crp" = quantile(CRP,probs = 0.75 ,na.rm = TRUE),
        "median GM"= median(GM_peak,na.rm = TRUE),
        "iqr 1" = quantile(GM_peak,probs = 0.25,na.rm= TRUE),
        "iqr 2" = quantile(GM_peak,probs = 0.75 ,na.rm = TRUE))
#26 patients have unknown status for nodule and were excluded.

wilcox_df <- df_7 %>% filter(!is.na(Nodule_CT))
wilcox.test(WBC_final ~ Nodule_CT, data = wilcox_df) #W = 1985, p-value = 0.8
wilcox.test(ESR ~ Nodule_CT, data = wilcox_df) #W = 1924.5, p-value = 0.133
wilcox.test(CRP ~ Nodule_CT, data = wilcox_df)#W = 1465.5, p-value = 0.29
wilcox.test(GM_peak ~ Nodule_CT, data= wilcox_df)#0.3164


#---
Halo_table <- df_7 %>%
    filter(!is.na(Halo_CT)) %>%
    group_by(Halo_CT) %>%
    summarise(
        "number of patients"  = n() ,
        "median WBC" = median(WBC_final,na.rm = TRUE),
        "Q1 WBC" = quantile(WBC_final,probs = 0.25 ,na.rm = TRUE),
        "Q2 WBC" = quantile(WBC_final,probs = 0.75 ,na.rm = TRUE),
        "median ESR" = median(ESR,na.rm = TRUE),
        "Q1 ESR"= quantile(ESR,probs = 0.25 ,na.rm = TRUE),
        "Q2 ESR" = quantile(ESR,probs = 0.75 ,na.rm = TRUE),
        "median CRP" = median(CRP,na.rm = TRUE),
        "iqr 1-crp" = quantile(CRP,probs = 0.25 ,na.rm = TRUE),
        "iqr 2-crp" = quantile(CRP,probs = 0.75 ,na.rm = TRUE),
        "median GM"= median(GM_peak,na.rm = TRUE),
        "iqr 1" = quantile(GM_peak,probs = 0.25,na.rm= TRUE),
        "iqr 2" = quantile(GM_peak,probs = 0.75 ,na.rm = TRUE))


wilcox_df_1 <- df_7 %>% filter(!is.na(Halo_CT))

wilcox.test(WBC_final ~ Halo_CT, data = wilcox_df_1) #0.55
wilcox.test(ESR ~ Halo_CT, data = wilcox_df_1) #0.4965
wilcox.test(CRP ~ Halo_CT, data = wilcox_df_1)#0.73
wilcox.test(GM_peak ~ Halo_CT, data= wilcox_df_1)#0.23
#28 patinets have unknown stats for halo
#_______________________________________________________________________________
#typic vs atypic 
atypic_cols <- c("GGO_CT", "Consolidation_CT","Pleuraleffusion_CT")
typic_cols <- c("Nodule_CT","Halo_CT","Reversehalo_CT","Cavity_CT")
df_7$Atypic_positive_1 <- apply(df_7[, atypic_cols], 1, function(x) {
    if (all(is.na(x))) {
        NA # keep unknown separate
    } else if (any(x == "Yes", na.rm = TRUE)) {
        "Yes"
    } else {
        "No"
    }
})

library(dplyr)
CT_table_at_1 <- df_7 %>%
    group_by(Atypic_positive_1) %>%
    summarise(
        "number of patients"  = n(),
        "median WBC" = median(WBC_final,na.rm = TRUE),
        "Q1 WBC" = quantile(WBC_final,probs = 0.25 ,na.rm = TRUE),
        "Q2 WBC" = quantile(WBC_final,probs = 0.75 ,na.rm = TRUE),
        "median ESR" = median(ESR,na.rm = TRUE),
        "Q1 ESR"= quantile(ESR,probs = 0.25 ,na.rm = TRUE),
        "Q2 ESR" = quantile(ESR,probs = 0.75 ,na.rm = TRUE),
        "median CRP" = median(CRP,na.rm = TRUE),
        "iqr 1-crp" = quantile(CRP,probs = 0.25 ,na.rm = TRUE),
        "iqr 2-crp" = quantile(CRP,probs = 0.75 ,na.rm = TRUE),
        "median GM"= median(GM_peak,na.rm = TRUE),
        "iqr 1" = quantile(GM_peak,probs = 0.25,na.rm= TRUE),
        "iqr 2" = quantile(GM_peak,probs = 0.75 ,na.rm = TRUE))

wilcox_df_2 <- df_7 %>% filter(!is.na(Atypic_positive_1))

wilcox.test(WBC_final ~ Atypic_positive_1, data = wilcox_df_2) #0.036
wilcox.test(ESR ~ Atypic_positive_1, data = wilcox_df_2) #0.6812
wilcox.test(CRP ~ Atypic_positive_1, data = wilcox_df_2)#0.4
wilcox.test(GM_peak ~ Atypic_positive_1, data= wilcox_df_2)#0.995
#26 patinets have unknown stats for halo




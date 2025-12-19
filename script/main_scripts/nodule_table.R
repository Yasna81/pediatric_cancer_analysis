#nodule size and number features 
f_n_table <- df_1 %>% filter(EORCT_MSG_criteria!="NA" & EORCT_MSG_criteria != "False positive" )

summary_nodule_f <- f_n_table %>%
    group_by(EORCT_MSG_criteria) %>%
    summarise("number of patients"  = n() ,
              "median n_n" = median(Number_nodule,na.rm = TRUE),
              "Q1 n" = quantile(Number_nodule,probs = 0.25 ,na.rm = TRUE),
              "Q2 n" = quantile(Number_nodule,probs = 0.75 ,na.rm = TRUE),
              "median s" = median(Size_nodule,na.rm = TRUE),
              "Q1 s"= quantile(Size_nodule,probs = 0.25 ,na.rm = TRUE),
              "Q2 s" = quantile(Size_nodule,probs = 0.75 ,na.rm = TRUE))
kruskal_n <-kruskal.test(Number_nodule ~ EORCT_MSG_criteria, data = f_n_table) #p-value = 0.335
kruskal_s <-kruskal.test(Size_nodule ~EORCT_MSG_criteria, data = f_n_table)#p-value = 0.5743
#fn 
#df_6 we are not going to consider wbc, but no possible
summary_table_FN <- df_6 %>%
    group_by(FN) %>%
    summarise(
        "number of patients"  = n() ,
        "Median FN episodes" = median(as.numeric(FN_episodes),na.rm = TRUE),
        "IQR 1 FN episodes" = quantile(as.numeric(FN_episodes),probs = 0.25 ,na.rm = TRUE),
        "IQR 2 FN episodes" = quantile(as.numeric(FN_episodes),probs = 0.75 ,na.rm = TRUE),
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

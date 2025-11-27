#importing data 
#install.packages("haven")
library("haven")
data <- read_sav("~/project/amir_project/data/Dr. Amanati_GM_23.02.1403.sav")
#spss to R labels may change
var_labels <- sapply(data,function(x) attr(x,"label"))
print(var_labels)
#converting labels to factors :
install.packages("dplyr")
library(dplyr)
df <- data %>% mutate(across(where(is.labelled),as_factor))
#subseting based on CT2 only NAs are ommited/ yes , nl and unclassicals were kept
df_1 <- df %>% filter(CT2 %in% c("Yes","NL","Unclassical CT") )
#df_1 has 169 rows. 114 yes/4 unclasical/51 Nl
#--------------------------------------analysis------------------------------------
#sensivity analysis :
df$lable <- ifelse(df$CT2 == "NA","out","in")
wilcox.test(WBC_2 ~ lable, data=df) #0.03
wilcox.test(CRP ~ lable, data = df) #1.198e-11
wilcox.test(ESR ~ lable, data = df) #5.55e-05
wilcox.test(GM_value_1 ~ lable , data = df) #1.607e-06

#table of sensivity analysis 
#standerd mean difference
install.packages("tableone")
#library(tableone) 
vars <- c("WBC_2","ESR","CRP","GM_value_1")
tab <- CreateTableOne(vars= vars ,strata = "lable",data = df , test= FALSE)
print(tab,smd = TRUE)
#absolute difference :
#install.packages("cobalt")
library(cobalt)
bal.tab(lable ~ WBC_2 + ESR + CRP + GM_value_1, data = df ,un = TRUE, disp.means =TRUE)
#table
library(dplyr)
summary_table <- df%>%
    group_by(lable) %>%
    summarise(
        "number of patients"  = n() ,
        "median WBC" = median(WBC_2,na.rm = TRUE),
        "Q1 WBC" = quantile(WBC_2,probs = 0.25 ,na.rm = TRUE),
        "Q2 WBC" = quantile(WBC_2,probs = 0.75 ,na.rm = TRUE),
        "median ESR" = median(ESR,na.rm = TRUE),
        "Q1 ESR"= quantile(ESR,probs = 0.25 ,na.rm = TRUE),
        "Q2 ESR" = quantile(ESR,probs = 0.75 ,na.rm = TRUE),
        "median CRP" = median(CRP,na.rm = TRUE),
        "iqr 1-crp" = quantile(CRP,probs = 0.25 ,na.rm = TRUE),
        "iqr 2-crp" = quantile(CRP,probs = 0.75 ,na.rm = TRUE),
        "median GM"= median(GM_value_1,na.rm = TRUE),
        "iqr 1" = quantile(GM_value_1,probs = 0.25,na.rm= TRUE),
        "iqr 2" = quantile(GM_value_1,probs = 0.75 ,na.rm = TRUE)
    )




#----------------------------------------------data-processing----------------------------------------
# we have wbc-1 and wbc-2 and avg value, so :
summary(df_1$WBC_1)
summary(df_1$WBC_2)
cor(df_1$WBC_1,df_1$WBC_2, use= "complete.obs")
cor(df_1$WBC_2,df_1$Avg_WBC, use= "complete.obs")
#no one is corrolated.
#lets check missingness in clinical values:
#install.packages("mice")
library(mice)
vars_to_check <- df[,c("WBC_1","WBC_2","ESR","CRP","GM_value_1")]
md.pattern(vars_to_check)
#p value 
install.packages("MissMech")
library(MissMech)
TestMCARNormality(df_1[,c("WBC_1","ESR","CRP","GM_value_1")])
#p value for all these is 0.003 means the missingnes in WBC_1 is not completely at random.
TestMCARNormality(df_1[,c("WBC_2","ESR","CRP","GM_value_1")])
# p value for all these clinical variables is 0.856: so MCAR(completly at random)
#imputation :
#for wbc we use wbc-2 as its the most clinically relevnt data at the time of infection, the missingnes is MCAR 
# for NA values in WBC-2 we use WBC-1 value.
df_1$WBC_used <- ifelse(!is.na(df_1$WBC_2),"WBC_2","WBC_1")
df_1$WBC_final <- ifelse(!is.na(df_1$WBC_2),df_1$WBC_2, df_1$WBC_1)
summary(df_1$WBC_final) #5 NAs
#other variables :
df_2 <- df_1
library(mice)
df_imp <- df_2[,c("ESR","CRP","GM_value_1")]
imp <- mice(df_imp,m=5,seed = 123)

df_comp <- complete(imp,1)
df_2$ESR <- df_comp$ESR
df_2$GM_value_1 <- df_comp$GM_value_1
df_2$CRP <- df_comp$CRP
#checking skewness :
hist(df_2$WBC_final_log) 
hist(df_2$CRP) #no log needed
hist(df_2$ESR) #no log needed
#for GM we need to aggregate
#log transform > all were skewed (will help the final model)
df_2$WBC_final_log <- log1p(df_2$WBC_final)
#df_2$CRP_log <- log1p(df_2$CRP)
#df_2$ESR_log <- log1p(df_2$ESR)
#GM values : we need the maxium value among all so :
library(dplyr)
df_3 <- df_2
gm_cols <- grep("^GM_value", names(df_3), value = TRUE)
df_3$GM_peak <- apply(df_3[, gm_cols], 1, function(x) {
    if(all(is.na(x))) NA else max(x, na.rm = TRUE)
})
hist(df_3$GM_peak_log)
df_3$GM_peak_log <- log1p(df_3$GM_peak)
#before moeling no need for logs 
#cancer diagnosis :
leukemias <- c("ALL", "AML", "T-cell ALL", "APL",  "Mix ALL+AML", "Mixed or undifferentiated leukemia","Mix ALL+AML")

lymphomas <- c("Non-Hodgkin Lymphoma", "Burkitt Lymphoma", "Hodgkin Lymphoma", "DLBCL", "LPD")

solid_tumors <- c(
    "Neuroblastoma", "Germ Cell Tumor", "Yolk Sak Tumor", "PNET", 
    "Hepatoblastomma", "Ependimoma", "Retinoblastoma", "Schwannoma", 
    "Nephroblastoma", "Wilmes Tumor", "Osteosarcoma", "Ewing Sarcoma", 
    "Rhabdomyosarcoma", "Soft tissue sarcoma", "Astrocytoma", 
    "Brain Tumor", "glioma"
)

non_malignant <- c(
    "Aplastic Anemia", "MDS", "Anemia", "Histiocytosis", "HLH",
    "CGD", "IgE syndrome", "PID", "LAD", 
    "Sickle cell Anemia", "Malignancy_rule outed"
)
library(dplyr)
df_4 <- df_3 %>%
    mutate(Cancer_Group = case_when(
        Diagnosis %in% leukemias ~ "Leukemia",
        Diagnosis %in% lymphomas ~ "Lymphoma",
        Diagnosis %in% solid_tumors ~ "Solid Tumor",
        Diagnosis %in% non_malignant ~ "Non-Malignant",
        TRUE ~ "Other"
    ))

#those who did not have nodules we change the nodule size and number of these cases from NA to 0.
df_4$Number_nodule[df_4$Nodule_CT == "No" ] <- 0
df_4$Size_nodule[df_4$Nodule_CT == "No" ] <- 0
# some patients have positive ct "yes" and a nodule size but the number of noudle is missing we assume that they have at least 1 nodule:
df_4$Number_nodule[is.na(df_4$Number_nodule) & df_4$Nodule_CT == "Yes" ] <- 1
#only one patient has a NA size of tumor. we leave it as it is .
#age for those under 1 year old 
df_4 <- df_4 %>% mutate(age_adjusted = if_else(is.na(Age),0.5,Age))
# fouces is on proven and probable.
df_5 <- df_4 %>% filter(EORCT_MSG_criteria %in% c("Proven","Probable")) # 12 patients with possible , 9 NA and 2 False positive.
# now we have 146 rows
#for a better view
df_6 <- df_5 %>% select(-starts_with("GM_value"), -starts_with("GM_date"),-starts_with("Name"))
# 5 patients do not have wbc1 or wbc2 or average
df_7 <- df_6 %>% filter(!(is.na(WBC_final)))
# now we have 141 patients left.
#saving cleaned data :
saveRDS(df_7,file ="cleaned-data-141.RDS")


# Load Packages and DataSets------------------------------------------------------------------
library(tidyverse)
library(mice)
library(gmodels)
library(caret)
library(survival) 
library(rms) 
library(patchwork)   # combining the plot
library(crayon)  # colorful text output
library(showtext)
font_add('Arial','/Library/Fonts/Arial.ttf')
showtext_auto() 
library(Cairo)

load("XXXXXXXXXXX.RData")  # Load case-control study data
load("XXXXXYYYYYYYYY.RData")
source("./funct.R")  # Load custom functions


# == Step 1 ==== Data Procession ==============================================================
#    !!!Data <- list()

      # Fill in disease diagnoses other than outcome variables:
tmp <- eczema_case_list$df_together_list$df_y 
for (col in colnames(tmp)){   # Factorize categorical variables before mice imputation
  tmp[, col] <- as.factor(tmp[, col])
}

Data$y_imp <-  mice(data = tmp %>% dplyr::select(-c('Child_withADbyDoctor')), 
                    method = rep('lasso.logreg', 5),  #5列
                    m = 1, maxit = 5, seed = 110)

Data$y_imp_i_ <- complete(Data$y_imp, fill = 1)

identical(rownames(complete(Data$y_imp, fill = 1)), rownames(complete(imp_eczema_list$imp_xBinary, fill = 1)))  # Interim data check

Data$Imp <- data.frame(complete(imp_eczema_list$imp_xBinary, fill = 1),
                          Data$y_imp_i_, 
                          complete(imp_eczema_list$imp_xMultiOrder, fill = 1),
                          complete(imp_eczema_list$imp_xMultiNonOrder, fill = 1),
                          complete(imp_eczema_list$imp_xContinue, fill = 1))

    # Case Group：
table(Data$Imp$Sexuality[Data$Imp_feast$event == '1']) 
prop.table(table(Data$Imp$Sexuality[Data$Imp_feast$event == '1']))

table(Data$Imp$Nationality[Data$Imp_feast$event == '1'])
prop.table(table(Data$Imp$Nationality[Data$Imp_feast$event == '1']))

summary(Data$Imp$realAgeToQuestionnaire_YMD_byM[Data$Imp_feast$event == '1']/12)
mean(Data$Imp$realAgeToQuestionnaire_YMD_byM[Data$Imp_feast$event == '1']/12)
sd(Data$Imp$realAgeToQuestionnaire_YMD_byM[Data$Imp_feast$event == '1']/12)


    # Control Group：
table(Data$Imp$Sexuality[Data$Imp_feast$event == '0']) 
prop.table(table(Data$Imp$Sexuality[Data$Imp_feast$event == '0']))

table(Data$Imp$Nationality[Data$Imp_feast$event == '0'])
prop.table(table(Data$Imp$Nationality[Data$Imp_feast$event == '0']))

summary(Data$Imp$realAgeToQuestionnaire_YMD_byM[Data$Imp_feast$event == '0']/12)
mean(Data$Imp$realAgeToQuestionnaire_YMD_byM[Data$Imp_feast$event == '0']/12)
sd(Data$Imp$realAgeToQuestionnaire_YMD_byM[Data$Imp_feast$event == '0']/12)



colnames(Data$Imp)
# dim(Data$Imp)

Data$Imp$QuestionaireCode <- eczema_case_list$caseControl_list$df_together$QuestionaireCode 

identical(Data$Imp$QuestionaireCode, eczema_case_list$caseControl_list$df_together$QuestionaireCode)  # Interim data check


Data$Imp <- data.frame(Data$Imp, eczema_case_list$caseControl_list$df_together[, c('match', 'type')]) %>% 
  dplyr::select(-c( 'Sexuality', 'Nationality', 'realAgeToQuestionnaire_YMD_byM', 
                   "height_birth", "weight_current", "height_current", 'bmi_birth', 'bmi_current' )) %>% 
  dplyr::mutate(birth_weight = factor(ifelse(weight_birth < 2.5, '2',  #1: Normal || 2: Underweight || 3: Overweight
                                             ifelse(weight_birth >4 , '3', '1')))) %>% 
  dplyr::select(!one_of('weight_birth')) %>% 
  dplyr::select(!one_of('QuestionaireCode')) %>% 
  dplyr::mutate(event = factor(ifelse(type == '1', '1', '0'))) %>% 
  dplyr::select(event, match, everything()) %>% 
  dplyr::select(!one_of('type'))


# Global Test---- check the effect ----------------------------------------------------------

get_CrossTable <- function(df_ = Data$test, outcome_ind_ = 1, outcome_nam_ = 'AD', factor_ind_ = 3:ncol(df_),
                           output_ = F, file_ = pastev(c(getwd(), '/ChisqTest.xlsx'), sep_ = '')){
  #df_: DF with variables to be tested (factor | multiple levels) and ending variables (factor | dichotomous)
        #! Note: Suggested format: col1: ending variable, col2: arbitrary (can be made up), col3-... : variables to be tested
  #outcome_ind_: column index of the ending variable ----e.g. Atopic Dermatitis (AD) / non-AD
  #factor_ind_: Column index of variable to be tested
  #output_: whether to export or not
  #file_: Path and name of the export Default: getwd()/ChisqTest.xlsx
  # pvalue may vary due to continuous correction or not, (two-sided test & n>40 & 2*2 list: p == 0.05: 3.84; p == 0.01: 10.83)
  
  # Basis: categorical chi-square:
  # Pearson chi-square Applies to all expected frequencies > 5 and not paired data ----- n ≥ 40 and T ≥ 5, using chi-square test
  # Yeta continuity correction Applies to any expected frequency less than 5 ------------n ≥ 40 and at least one 1 ≤ T < 5, using corrected chi-square test
  # Fisher exact probability method for any expected frequency less than 1 -------------n < 40 or at least 1 T < 1 , using Fisher exact probability method
  # Paired information using McNemar goodness-of-fit test
  
  # Trend chi-square:
  # Multicategorical ordered information, categorical chi-square continues trend chi-square after differences are found ----DescTools::CochranArmitageTest()

  get_chisqTable <- function(chisq_ret, col_, outcome_nam_){
    # chisq_ret: the return by CrossTable function
    dat_ <- data.frame(matrix(c(apply(ret_sav_[[col_]]$t, MARGIN = 1, FUN = function(x){paste(
      sum(x), '(', round(100*sum(x)/sum(ret_sav_[[col_]]$t), 2), '%)', sep = '')}), 
      paste0(ret_sav_[[col_]]$t[, '1'], '(', round(ret_sav_[[col_]][["prop.row"]][, '1']*100, 2), '%)'), rep(NA, nrow(ret_sav_[[col_]]$t)),
      paste0(ret_sav_[[col_]]$t, '(', round(ret_sav_[[col_]]$prop.col*100, 2), '%)')
    ), ncol = 5, nrow = nrow(ret_sav_[[col_]]$t), byrow = F, dimnames = list(row = paste(col_, rownames(ret_sav_[[col_]]$t), sep = '*'),
                                                                             col = c('各变量组别间样本数（%）', paste(outcome_nam_, '数（%）', sep = ''), 'Gap',
                                                                                     '对照组', paste(outcome_nam_, '组', sep = '')))))
    return(dat_)
  }
  
  
  library(gmodels)
  library(DescTools)
  
  ChisqT_ = list()
  ret_sav_ = list()
  fish_ = c('If no error, Fisher processing has been done, no need to pay attention to this result……')
  for (col_ in colnames(df_[, 3:ncol(df_)])){
    tmp_ <- CrossTable(table(df_[, col_], df_[, outcome_ind_]), chisq = T, expected = T)  
    
    if (all(tmp_$chisq$expected >= 5) & sum(tmp_$chisq$expected) >= 40){ 
      ret_sav_[[col_]] <- CrossTable(table(df_[, col_], df_[, outcome_ind_]), chisq = T,
                                     expected = T, fisher = F) 
      tmp2_ <- ret_sav_[[col_]]
      if (all(dim(tmp2_[["t"]]) <=2 )){
        ChisqT_[[col_]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
          dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(X2 = c(tmp2_[["chisq"]][["statistic"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(p.value = c(tmp2_[["chisq"]][["p.value"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(method = c(tmp2_[["chisq"]][["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::select(var, everything())
      }else if (any(dim(tmp2_[["t"]])>2)){
        ChisqT_[[col_]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
          dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(X2 = c(tmp2_$chisq$statistic, rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(p.value = c(tmp2_$chisq$p.value, rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(method = c(tmp2_[["chisq"]][["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::select(var, everything())
        
        
        ret_sav_[[paste(col_, '_trend_chi-square_test')]] <- CochranArmitageTest(table(df_[, col_], df_[, outcome_ind_]))
        tmp3_ <- ret_sav_[[paste(col_, '_trend_chi-square_test')]]
        ChisqT_[[paste(col_, '_trend_chi-square_test')]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
          dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(X2 = c(tmp3_[["statistic"]], 'Zvalue', rep(NA, nrow(tmp2_$t)-2))) %>% 
          dplyr::mutate(p.value = c(tmp3_[["p.value"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(method = c(tmp3_[["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::select(var, everything())
      }else{cat('please_check_code...................')}
      
    }else if(any(tmp_$chisq$expected < 5) & sum(tmp_$chisq$expected) >= 40 & all(tmp_$chisq$expected > 1)){
        
        ret_sav_[[col_]] <- CrossTable(table(df_[, col_], df_[, outcome_ind_]), chisq = T,
                                       expected = T, fisher = F) 
        
        tmp2_ <- ret_sav_[[col_]]
        if (all(dim(tmp2_[["t"]]) <=2 )){
          ChisqT_[[col_]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
            dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::mutate(X2 = c(tmp2_[["chisq.corr"]][["statistic"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::mutate(p.value = c(tmp2_[["chisq.corr"]][["p.value"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::mutate(method = c(tmp2_[["chisq.corr"]][["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::select(var, everything())
        }else if (any(dim(tmp2_[["t"]])>2)){
          ChisqT_[[col_]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
            dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::mutate(X2 = c(tmp2_$chisq$statistic, rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::mutate(p.value = c(tmp2_$chisq$p.value, rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::mutate(method = c(tmp2_[["chisq"]][["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::select(var, everything())
          
          
          ret_sav_[[paste(col_, '_trend_chi-square_test')]] <- CochranArmitageTest(table(df_[, col_], df_[, outcome_ind_]))
          tmp3_ <- ret_sav_[[paste(col_, '_trend_chi-square_test')]]
          ChisqT_[[paste(col_, '_trend_chi-square_test')]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
            dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::mutate(X2 = c(tmp3_[["statistic"]], 'Zvalue', rep(NA, nrow(tmp2_$t)-2))) %>% 
            dplyr::mutate(p.value = c(tmp3_[["p.value"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::mutate(method = c(tmp3_[["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
            dplyr::select(var, everything())
        }else{cat('please_check_code...................')}
    }else{ 
      fish_ <- c(fish_, col_)
      
      ret_sav_[[col_]] <- CrossTable(table(df_[, col_], df_[, outcome_ind_]), chisq = T,
                                     expected = T, fisher = T) 
      tmp2_ <- ret_sav_[[col_]]
      if (all(dim(tmp2_[["t"]]) <=2 )){
        ChisqT_[[col_]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
          dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(X2 = c(tmp2_[["chisq.corr"]][["statistic"]], 'Fisher method_report continuity-corrected chi-square value', rep(NA, nrow(tmp2_$t)-2))) %>% 
          dplyr::mutate(p.value = c(tmp2_[["fisher.ts"]][["p.value"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(method = c(tmp2_[["fisher.ts"]][["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::select(var, everything())
      }else if (any(dim(tmp2_[["t"]])>2)){
        ChisqT_[[col_]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
          dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(X2 = c(tmp2_[["chisq.corr"]][["statistic"]], 'Fisher method_report continuity-corrected chi-square value', rep(NA, nrow(tmp2_$t)-2))) %>% 
          dplyr::mutate(p.value = c(tmp2_[["fisher.ts"]][["p.value"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(method = c(tmp2_[["fisher.ts"]][["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::select(var, everything())
        
        
        ret_sav_[[paste(col_, '_trend_chi-square_test')]] <- CochranArmitageTest(table(df_[, col_], df_[, outcome_ind_]))
        tmp3_ <- ret_sav_[[paste(col_, '_trend_chi-square_test')]]
        ChisqT_[[paste(col_, '_trend_chi-square_test')]] <- get_chisqTable(chisq_ret = ret_sav_[[col_]], col_ = col_, outcome_nam_ = outcome_nam_) %>% 
          dplyr::mutate(var = c(col_, rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(X2 = c(tmp3_[["statistic"]], 'Zvalue', rep(NA, nrow(tmp2_$t)-2))) %>% 
          dplyr::mutate(p.value = c(tmp3_[["p.value"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::mutate(method = c(tmp3_[["method"]], rep(NA, nrow(tmp2_$t)-1))) %>% 
          dplyr::select(var, everything())
      }else{cat('please_check_code...................')}
    }
  }
  cat('Remember to randomly select a variable for table check........')
  togth_ <- do.call('rbind', ChisqT_)
  if(output_){openxlsx::write.xlsx(x = togth_, file = file_, rowNames = T)}
  return(list(ret_sav = ret_sav_, ChisqT = ChisqT_, togth = togth_))
}

#     !!!Ret <- list()

Ret$ChisqT <- get_CrossTable(Data$Imp, outcome_ind_ = 1, factor_ind_ = 3:ncol(Data$Imp), 
                       output_ = T, file_ = 'ChisqTest.xlsx')


ret <- CrossTable(table(Data$Imp[, 1], Data$Imp[, 'father_with_asthma']), chisq = T,
           expected = T, fisher = F) 



# == Step 2 ===== with The Objection of Early life exposure ================================
Data$Imp_feast <- Data$Imp

colnames(Data$Imp_feast)

for (col_i_ in c(84:99)){
  Data$Imp_feast[, col_i_] <- factor(ifelse(Data$Imp_feast[, col_i_] == '1' | Data$Imp_feast[, col_i_] == '4', '0', '1')) 
}

for (col_i_ in c(100:129)){ 
  Data$Imp_feast[, col_i_] <- factor(ifelse(Data$Imp_feast[, col_i_] == '1', '0', '1')) 
}
Data$Imp_feast[, 76] <- factor(ifelse(Data$Imp_feast[, 76] == '1', '0',
                                   ifelse(Data$Imp_feast[, 76] == '2' | Data$Imp_feast[, 76] == '3', '1', '2')))

Data$Imp_feast[, 81] <- factor(ifelse(Data$Imp_feast[, 81] == '0', '0',
                                   ifelse(Data$Imp_feast[, 81] == '1' | Data$Imp_feast[, 81] == '2', '1', '2')))


Data$Imp_feast <- Data$Imp_feast[, c(1, 2, 
                                       10,
                                       25:27, 
                                       28:30,
                                       51:54, 
                                       64:69,
                                       76,
                                       81, 
                                       84:95, 
                                       149:150,  
                                       153:155  
                                       )]

Ret$ChisqT_feast <- get_CrossTable(df_ = Data$Imp_feast, outcome_ind_ = 1, factor_ind_ = 3:ncol(Data$Imp_feast), 
                             output_ = T, file_ = 'ChisqTest_feast.xlsx')

        #Interim check：

df_total_logi <- Data$Imp_feast

table(df_total_logi$event); table(df_total_logi$Mode_of_birth); table(df_total_logi$Mode_of_birth, df_total_logi$event)
colnames(df_total_logi)

df_total_logi$virtualSurvivalTime <- ifelse(df_total_logi$event == '0', 1, 0)  # Generate virtual survival time variable
df_total_logi$event <- as.numeric(df_total_logi$event)

Ret$Log <- coxph(formula = Surv(virtualSurvivalTime, event) ~ father_with_asthma + father_with_allergic_rhinitis + 
                   father_with_AD + mother_with_allergic_rhinitis + mother_with_AD + 
                   Renovation of the dwelling\nbefore mp + Mold in the dwelling\nbefore mp + 
                   Renovation of the dwelling\nduring mp + Mold in the dwelling\nduring mp + Father smoking in the\ndwelling during mp + 
                   Renovation of the dwelling\nduring cfy + Mold in the dwelling\nduring cfy + 
                   Father smoking in the\ndwelling during cfy + Grow plants or have\npets during cfy + 
                   during_breastfeeding + antibiotic_therapy +
                   only_child + birth_weight, data = df_total_logi)

#     Ret$Comb$mi1 <- with(data = Data$Imp_feast, expr = coxph(formula))


Ret$Log_summary <- summary(Ret$Log)

Ret$Log_refine <- data.frame(cbind(Ret$Log_summary$coefficients, Ret$Log_summary$conf.int))

openxlsx::write.xlsx(x = Ret$Log_refine, file = 'Logi_feast.xlsx', rowNames = T)


# Additionally, Check the association between low birth weight and siblings ---------------------
df_further<- data.frame(list(
  event = ifelse(as.character(Data$Imp_feast$birth_weight) == '2', 1, 0),
  nothing = 'None',
  sibling = ifelse(as.character(Data$Imp_feast$only_child) == '0', 0, 1),
  senior = ifelse(as.character(Data$Imp_feast$only_child) == '1', 1, 0),
  monior = ifelse(as.character(Data$Imp_feast$only_child) == '1', 1, 0) 
  ))

get_CrossTable(df_ = df_further, outcome_ind_ = 1, factor_ind_ = 3:ncol(df_further), output_ = T, file_ = 'df_further.xlsx')


# Stratification analysis ------------------------------------------------------------------

    # data tidying
df_Strati <- Data$Imp_feast 
colnames(df_Strati)

df_Strati <- df_Strati %>% dplyr::mutate(
  parent_atopic = as.factor(ifelse(father_with_allergic_rhinitis == '1' | father_with_asthma == '1' |
                                  father_with_AD == '1' | mother_with_allergic_rhinitis == '1' | 
                                  mother_with_asthma == '1' | mother_with_AD == '1', 1, 0)),
  mom_atopic = as.factor(ifelse(mother_with_allergic_rhinitis == '1' | 
                                  mother_with_asthma == '1' | mother_with_AD == '1', 1, 0)),
  dad_atopic = as.factor(ifelse(father_with_allergic_rhinitis == '1' | father_with_asthma == '1' |
                                  father_with_AD == '1', 1, 0)),  # No further analysis on mom/dad, ---- sample size too small
  parent_AD = as.factor(ifelse(father_with_AD == '1' | mother_with_AD == '1', 1, 0)),
  parent_rhinitis_asthma = as.factor(ifelse(father_with_allergic_rhinitis == '1' | father_with_asthma == '1' |
                                            mother_with_allergic_rhinitis == '1' | mother_with_asthma == '1', 1, 0)),
  mom_rhinitis_asthma = as.factor(ifelse(mother_with_allergic_rhinitis == '1' | mother_with_asthma == '1', 1, 0)),
  dad_rhinitis_asthma = as.factor(ifelse(father_with_allergic_rhinitis == '1' | father_with_asthma == '1', 1, 0)),
  birth_way = as.factor(ifelse(Mode_of_birth == '1', 1, 0)),  # 1: Vaginal delivery; 0: Cesarean section
) %>% dplyr::select(event, match, parent_atopic, mom_atopic, dad_atopic, parent_AD,
                    parent_rhinitis_asthma, mom_rhinitis_asthma, dad_rhinitis_asthma, birth_way, everything())

attr(df_Strati$birth_way, 'label') <- '1：eutocia；0：Cesarean'


get_compositionDrop(df_Strati)


get_matchedCaseControlStrati_dimension <- function(df_, matched = 4, var_stra = 'mom_atopic', tar_limitValue = '1',
                                                   tar_limitValue_c = '0', tar_limitNum = 300){
  # downscale multi-match (1:4) cases into 1:1 matches by stratifying against stratification variables
  # df_: dataframe, requirements: first line event---binary&factor/c('0', '1'), second line match----number/(default 1:4 match under) ----- please ensure that stratification variables are c('0', '1')
  # matched: match ratio --- default 4
  # var_stra: the variable to be stratified
  # tar_limitValue: hierarchical value to prioritize for matching (tends to be the one with the smaller composition ratio)
  # tar_limitValue_c: the opposing value of the hierarchical value to be preferentially matched (tends to be the larger composition ratio)
  # tar_limitNum: number of hierarchical values (often with smaller composition ratios) to prioritize for matching
  ind_Strati <- c()
  yes_history <- c()
  
  match_record1 <- c()  # matched record for parent_atopic history
  match_record0 <- c()
  for (i in unique(df_$match)){
    tmp <- df_[df_$event == '1' & df_$match == i,]
    tmp2 <- df_[df_$event == '0' & df_$match == i,]
    for (m in 1:matched){
      if (length(yes_history[yes_history == tar_limitValue]) <= tar_limitNum){
        if (tmp2[m, var_stra] == tmp[, var_stra] & as.character(tmp[, var_stra]) == tar_limitValue){
          ind_Strati <- c(ind_Strati, rownames(tmp), rownames(tmp2[m, ]))
          cat(as.character(tmp[, var_stra]))
          yes_history <- c(yes_history, as.character(tmp[, var_stra]), as.character(tmp[, var_stra]))
          match_record1 <- c(match_record1, i)
          break
        }
      }
    }
  }
  
  match_rest <- unique(df_$match)[is.na(di_which(match_record1, unique(df_$match)))]  # the rest of match number
  for (i in match_rest){
    tmp <- df_[df_$event == '1' & df_$match == i,]
    tmp2 <- df_[df_$event == '0' & df_$match == i,]
    for (m in 1:matched){
      if (tmp2[m, var_stra] == tmp[, var_stra] & as.character(tmp[, var_stra]) == tar_limitValue_c){
        ind_Strati <- c(ind_Strati, rownames(tmp), rownames(tmp2[m, ]))
        cat(as.character(tmp[, var_stra]))
        yes_history <- c(yes_history, as.character(tmp[, var_stra]), as.character(tmp[, var_stra]))
        match_record0 <- c(match_record0, i)
        break
      }
    }
  }
  
  match_rest <- match_rest[is.na(di_which(match_record0, match_rest))]  # the rest of match number
  for (i in match_rest){
    tmp <- df_[df_$event == '1' & df_$match == i,]
    tmp2 <- df_[df_$event == '0' & df_$match == i,]
    for (m in 1:matched){
      if (tmp2[m, var_stra] == tmp[, var_stra] & as.character(tmp[, var_stra]) == tar_limitValue){
        ind_Strati <- c(ind_Strati, rownames(tmp), rownames(tmp2[m, ]))
        cat(as.character(tmp[, var_stra]))
        yes_history <- c(yes_history, as.character(tmp[, var_stra]), as.character(tmp[, var_stra]))
        match_record1 <- c(match_record1, i)
        break
      }
    }
  }
  return(ind_Strati)
}



Stra <- list(  # pool for stratification analysis
  dat = list(),
  ind = list(),
  desc = list(), 
  tab_count = list(),
  logi = list()
)

for (stra_p in c("parent_atopic", "mom_atopic", "dad_atopic", "parent_AD", "parent_rhinitis_asthma", "mom_rhinitis_asthma", "dad_rhinitis_asthma", "birth_way")){
  
  ind_Strati <- get_matchedCaseControlStrati_dimension(df_ = df_Strati, matched = 4, var_stra = stra_p,
                                                             tar_limitValue = '1', tar_limitValue_c = '0', tar_limitNum = 300)

  cat('\n', 'unique(df_[[stra_p]]):', unique(df_Strati[[stra_p]]))
  print(unique(df_Strati[[stra_p]]))

  length(unique(ind_Strati)) == length(ind_Strati)
  if(length(unique(ind_Strati)) == length(ind_Strati) ){
    cat('nOK, According to stratification variables', stra_p,'the extracted metadata indexes are all unique and non-repeating')
  }else{
    cat('!!ERROR:, Please check the hierarchical variables', stra_p,'procession of matching')
  }
  
  df_StratiOK <- df_Strati[ind_Strati,]
  
  tmp <- table(df_StratiOK[[stra_p]])
  names(tmp)
  
  if(length(unique(df_StratiOK$match)) == nrow(df_StratiOK)/2){
    cat('\nOK, According to stratification variables', stra_p,'the number of unique non-duplicates of the extracted data match is equal to1/2 * nrow()\n\n\n')
    for (stra_p_v in names(tmp)){
      Stra[['dat']][[paste(stra_p, '&', stra_p_v, sep = '')]] <- df_StratiOK[df_StratiOK[[stra_p]] == stra_p_v, ]
    }
    Stra[['ind']][[stra_p]] <- ind_Strati
    
  }else{
    cat('!!ERROR:, Please check the hierarchical variables', stra_p,'procession of matching\n')
  }
}


for (stra_p in names(Stra$dat)){
  dat <- Stra[["dat"]][[stra_p]]
  
  dat$virtualSurvivalTime <- ifelse(dat$event == '0', 1, 0)  # create virtual survival time variable
  dat$event <- ifelse(dat$event == '1', 1, 0)

  formula_ <- as.formula('Surv(virtualSurvivalTime, event) ~ during_breastfeeding + 
                         antibiotic_therapy + only_child + birth_weight')
  ret_dat <- coxph(formula = Surv(virtualSurvivalTime, event) ~ during_breastfeeding + 
                     antibiotic_therapy + only_child + birth_weight, data = dat)
  
  ret_dat_sum <- summary(ret_dat)
  
  ret_dat_refine <- data.frame(ret_dat_sum$coefficients, ret_dat_sum$conf.int[, c("lower .95", "upper .95")])
  or_string <- paste(round(ret_dat_refine$exp.coef., 2), '(', round(ret_dat_refine$lower..95, 2), ',',
                     round(ret_dat_refine$upper..95, 2), ')', 
                     ifelse(ret_dat_refine$Pr...z.. < 0.001, '***',
                            ifelse(ret_dat_refine$Pr...z.. < 0.01, '**',
                                   ifelse(ret_dat_refine$Pr...z.. < 0.05, '*', ''))), sep = '')
  
  Stra[['coef']][[stra_p]] <- as.data.frame(matrix(data = ret_dat_refine[, 'coef'], 
                                    ncol = 1, nrow = nrow(ret_dat_refine),
                                    dimnames = list(rownames = rownames(ret_dat_refine), colnames = stra_p)))
  Stra[['coef_se']][[stra_p]] <- as.data.frame(matrix(data = ret_dat_refine[, 'se.coef.'], 
                                    ncol = 1, nrow = nrow(ret_dat_refine),
                                    dimnames = list(rownames = rownames(ret_dat_refine), colnames = stra_p)))
  Stra[['p_value']][[stra_p]] <- as.data.frame(matrix(data = ret_dat_refine[, 'Pr...z..'], 
                                    ncol = 1, nrow = nrow(ret_dat_refine),
                                    dimnames = list(rownames = rownames(ret_dat_refine), colnames = stra_p)))
  Stra[['OR_CI']][[stra_p]] <- as.data.frame(matrix(data = or_string, 
                                    ncol = 1, nrow = nrow(ret_dat_refine),
                                    dimnames = list(rownames = rownames(ret_dat_refine), colnames = stra_p)))
  tab_nam <- c('sum')
  tab_count <- c(nrow(dat))
  tab_count_mer <- c(nrow(dat))
  for(var_ in c('during_breastfeeding', 'antibiotic_therapy', 'only_child', 'birth_weight')){
    tab_ <- table(dat[[var_]])
    tab_nam <- c(tab_nam, paste(var_, names(tab_), sep = '&'))
    tab_count <- c(tab_count, as.numeric(tab_))
    percen <- round(as.numeric(tab_)[1:(length(tab_)-1)]/nrow(dat), 4) 
    percen <- c(percen, round(1-sum(percen), 4))
    tab_count_mer <- c(tab_count_mer, paste0(as.numeric(tab_), 
                                             paste('(', percen*100, '%', ')', sep = '')))
  }
  Stra[["tab_count"]][[stra_p]] <- as.data.frame(matrix(data = tab_count_mer, 
                         ncol = 1, nrow = length(tab_nam),
                         dimnames = list(rownames = tab_nam, colnames = stra_p)))
  
  ret_dat_refine$merge <- or_string
  Stra[["logi"]][[stra_p]] <- ret_dat_refine
  
  # dat_desc <- data.frame(list(
  #   event = dat$event,
  #   nothing = 'None',
  #   dat[, c('during_breastfeeding', 'antibiotic_therapy', 'only_child', 'birth_weight')]
  # ))
  # desc_ret <- get_CrossTable(df_ = dat_desc, outcome_ind_ = 1, factor_ind_ = 3:ncol(dat_desc),
  #                            output_ = F)
  # Stra[["desc"]][[stra_p]] <- desc_ret$togth
}; Stra[["coef"]] <- do.call('cbind', Stra[["coef"]]); Stra[["coef_se"]] <- do.call('cbind', Stra[["coef_se"]]); Stra[["p_value"]] <- do.call(
  'cbind', Stra[["p_value"]]); Stra[["tab_count"]] <- do.call('cbind', Stra[["tab_count"]]); Stra[["OR_CI"]] <- do.call('cbind', Stra[["OR_CI"]])


# export the stratification data

stra_output <- Stra[c("tab_count", "coef", "coef_se", "p_value", "OR_CI")]
openxlsx::write.xlsx(x = stra_output, file = 'stra_output.xlsx', rowNames = T)

save(Stra, file = 'Stra.RData')

var_adjust <- as.data.frame(list(orignVar = rownames(Stra[["coef"]]),
        newVar = c("1-3 months of exclusive breastfeeding", "4+ months of exclusive breastfeeding",
            "1-2 times of antibiotic therapy during cfy", "3+ times of antibiotic therapy during cfy",
            "Have older siblings", "Have younger siblings", "Low birth weight (~2500 g)", "High birth weight (4500+ g)")))
View(var_adjust)



    # adjust the variable's names, and be ready for drawing
heat_ggplot <- list()
heatmap_dat <- list(
  atopic = Stra[["coef"]][, c("parent_atopic&0", "parent_atopic&1", "mom_atopic&0", "mom_atopic&1", "dad_atopic&0", "dad_atopic&1")],
  rhinitis_asthma = Stra[["coef"]][, c("parent_rhinitis_asthma&0", "parent_rhinitis_asthma&1", "mom_rhinitis_asthma&0",
                                       "mom_rhinitis_asthma&1", "dad_rhinitis_asthma&0", "dad_rhinitis_asthma&1")],
  others = Stra[["coef"]][, c("parent_AD&0", "parent_AD&1", "birth_way&0", "birth_way&1")]
); for (i in names(heatmap_dat)){heatmap_dat[[i]]['var'] <- var_adjust[, 'newVar']
  rownames(heatmap_dat[[i]]) <- 1:nrow(heatmap_dat[[i]])
  heatmap_dat[[i]] <- heatmap_dat[[i]]%>% reshape2:: melt(id.vars = c('var'), 
                      measure.vars = colnames(heatmap_dat[[i]])[-length(colnames(heatmap_dat[[i]]))],
                      variable.name = "class_regl", value.name = 'B_V') %>% 
    dplyr::mutate(var = factor(var, levels = rev(var_adjust[, 'newVar']), labels = rev(var_adjust[, 'newVar'])))
  
  # Plot the heatmap
  heat_ggplot[[i]] <- ggplot(heatmap_dat[[i]], 
                          aes(x = class_regl,y = var, fill = B_V)) +
    geom_raster() +
    scale_fill_gradientn(colors = c('#40c057', 'white', "#e74c3c"), 
                         values = c(0, 0.15, .35, .45, .5, 1)) + 
    labs(x = NULL, y = NULL) + 
    scale_x_discrete(position = "top") +
    theme_classic() +
    theme_default
}

heat_ggplot_output <- heat_ggplot[[1]]/heat_ggplot[[2]]/heat_ggplot[[3]] + 
  plot_annotation(tag_levels = "A")

ggsave(filename = 'Stratification_heatmap.pdf', plot = heat_ggplot_output, 
       width = 45, height = 45, units = 'cm')




# Association between low birth weight and breastfeeding ---------------------------------

df_lowerBirthWeight <- df_Strati[, c('low_birthWeight', 'during_breastfeeding')]

df_lowerBirthWeight <-  df_lowerBirthWeight %>% dplyr::mutate(
  breastfeed = ifelse(during_breastfeeding == '0', 1, 0) 
) %>% dplyr::select(breastfeed, everything())

table(df_lowerBirthWeight$low_birthWeight)

table(df_lowerBirthWeight$low_birthWeight, df_lowerBirthWeight$breastfeed)

CrossTable(table(df_lowerBirthWeight$low_birthWeight, df_lowerBirthWeight$during_breastfeeding), chisq = T)

CrossTable(table(df_lowerBirthWeight$low_birthWeight, df_lowerBirthWeight$breastfeed), chisq = T)



colnames(df_lowerBirthWeight)

df_lowerBirthWeight$virtualSurvivalTime <- ifelse(df_lowerBirthWeight$breastfeed == '0', 1, 0)  # create the virtual survival time
df_lowerBirthWeight$breastfeed <- as.numeric(df_lowerBirthWeight$breastfeed)

breastfeed_log <- coxph(formula = Surv(virtualSurvivalTime, breastfeed) ~ low_birthWeight,
                        data = df_lowerBirthWeight)

#     Ret$Comb$mi1 <- with(data = Data$Imp_feast, expr = coxph(formula))

breastfeed_summary <- summary(breastfeed_log)
breastfeed_refine <- data.frame(cbind(breastfeed_summary$coefficients, breastfeed_summary$conf.int))
openxlsx::write.xlsx(x = breastfeed_refine, file = 'breastfeed_feast.xlsx', rowNames = T)





# === Step 3 ==== Construct the Machine Learning Models ------------------------------------

var_standart <- data.frame(list(
  var = c('event', 'match', 'father_with_asthma', 'father_with_allergic_rhinitis',  'father_with_AD',
              'mother_with_allergic_rhinitis', 'mother_with_AD', 'Renovation of the dwelling\nbefore mp',
              'Mold in the dwelling\nbefore mp', 'Renovation of the dwelling\nduring mp', 'Mold in the dwelling\nduring mp', 
              'Father smoking in the\ndwelling during mp', 'Renovation of the dwelling\nduring cfy', 'Mold in the dwelling\nduring cfy', 
              'Father smoking in the\ndwelling during cfy', 'Grow plants or have\npets during cfy', 'during_breastfeeding',
              'antibiotic_therapy', 'only_child', 'birth_weight'),
  var_eng = c('event', 'match', 'Asthma of dad', 'AR of dad', 'AD of dad', 'AR of mom', 'AD of mom',
              'Renovation of the dwelling\nbefore mp', 'Mold in the dwelling\nbefore mp', 'Renovation of the dwelling\nduring mp',
              'Mold in the dwelling\nduring mp', 'Father smoking in the\ndwelling during mp', 'Renovation of the dwelling\nduring cfy',
              'Mold in the dwelling\nduring cfy', 'Father smoking in the\ndwelling during cfy', 'Grow plants or have\npets during cfy', 
              'Months of exclusive\nbreastfeeding', 'Times of antibiotic\ntherapy during cfy', 'Siblings', 'birth weight'),
  var_eng_strict = c('event', 'match', 'Asthma_of_dad', 'AR_of_dad', 'AD_of_dad', 'AR_of_mom', 'AD_of_mom',
              'Renovation_of_the_dwelling_before_mp', 'Mold_in_the_dwelling_before_mp', 'Renovation_of_the_dwelling_during_mp',
              'Mold_in_the_dwelling_during_mp', 'Father_smoking_in_the_dwelling_during_mp', 'Renovation_of_the_dwelling_during_cfy',
              'Mold_in_the_dwelling_during_cfy', 'Father_smoking_in_the_dwelling_during_cfy', 'Grow_plants_or_have_pets_during_cfy',
              'Months_of_exclusive_breastfeeding', 'Times_of_antibiotic_therapy_during_cfy', 'Siblings', 'birth_weight')
))

df_model <- Data[['Imp_feast']][, var_standart[, 'var']]
rownames(df_model) <- 1:nrow(df_model); colnames(df_model) <- var_standart[match(colnames(df_model), var_standart[, 'var']), 'var_eng_strict']
write.table(label(df_model), file = 'df_model_label.txt')  # Export the dafult labels.
for (col in colnames(df_model)){attr(df_model[, col], 'label') <- NULL}  # down set to total data format.

df_model <- df_strToNum(df_model)
df_model <- df_model %>% dplyr::mutate(
  YongerSiblings = ifelse(Siblings == 2, 1, 0),
  OlderSiblings = ifelse(Siblings == 1, 1, 0),
  LowBirthWeight = ifelse(birth_weight == 2, 1, 0),  
  ParentalAllergyHistory = apply(df_model[, c("Asthma_of_dad", "AR_of_dad", "AD_of_dad", "AR_of_mom", "AD_of_mom")],
                                 MARGIN = 1, FUN = function(x){ifelse(sum(x) > 0, 'yes', 'no')})
) %>% dplyr::select(-c('Siblings')) %>% dplyr::select('event', 'match', 'ParentalAllergyHistory', everything())



#  Data deletion --- 25--01--16----------------------------------------------------------

df_test <- df_model %>% dplyr::mutate(ParentalAllergyHistory = ifelse(ParentalAllergyHistory == 'yes', 1, 0)) %>% 
  dplyr::select(-c("ParentalAllergyHistory", "birth_weight", "YongerSiblings"))

#priority_v <- list(Times_of_antibiotic_therapy_during_cfy = 1,
#                   Months_of_exclusive_breastfeeding = 2,
#                   OlderSiblings = 1, 
#                   Renovation_of_the_dwelling_during_mp = 1,
#                   LowBirthWeight = 2)
#     library(Boruta); boruta_ml <- Boruta(as.formula('event ~ .'), data=df_test, doTrace=1)
#     ml_vars <- c('event', 'match', getSelectedAttributes(boruta_ml, withTentative = F), names(priority_v))
#     df_test <- df_test[, ml_vars[!duplicated(ml_vars)]]

df_cutDownDiff <- NULL
for (match_i in unique(df_test[, 'match'])){
  dat <- df_test[df_test[, 'match'] == match_i, ]
  ref_v <- dat[dat[, 'event'] == 1, ]
  match_m <- t(as.data.frame(apply(dat, MARGIN = 1, FUN = function(x){x == ref_v})));colnames(match_m) <- colnames(dat)
  match_m2 <- list()
  for (col in colnames(match_m)){
    match_m2[[col]] <- sapply(match_m[, col], FUN = function(y){ifelse(y == F & col %in% names(priority_v), priority_v[[col]],
                                          ifelse(y == F & !(col %in% names(priority_v)), 1, 0))})
  };match_m2 <- do.call('cbind', match_m2)
  dat_cutDwon <- dat[dat[, 'event'] == 1 | rownames(dat) == names(which.max(apply(match_m2, MARGIN = 1, FUN = function(x){sum(x)}))), ]
  if (is.null(df_cutDownDiff)){
    df_cutDownDiff <- dat_cutDwon
  }else{df_cutDownDiff <- rbind(df_cutDownDiff, dat_cutDwon)}
}


df_test_uni <- df_cutDownDiff %>% dplyr::mutate(
  event = ifelse(event == 1, 'AD', 'Control')
);for (col in colnames(df_test_uni)){df_test_uni[, col] <- as.factor(df_test_uni[, col])}

univFeatureSelection_ret <- get_UniMultivariateAnalysis(df_ = df_test_uni, outcome_ind_ = 1, outcome_nam_ = 'AD',
                            factor_ind_ = 3:ncol(df_test_uni),
                            perform__ = list(chisq = T, text = T, logistic = T),
                            output_ = list(chisq = T, text = F, logistic = F), family = 'binomial',
                            file_ = 'univFeatureSelet_卡方_250119')

univ_ret <- univFeatureSelection_ret[["univariateRet"]][!is.na(univFeatureSelection_ret[["univariateRet"]][, "p.value"]), ]
df_test_uni2 <- df_test_uni[, unique(c('event', 'match', univ_ret[, 'var'][univ_ret[, "p.value"] < 0.3]))] %>% 
  dplyr::mutate(event = ifelse(event == 'AD', 1, 0)) %>% df_strToNum()

# The Structure Functions of All Models==========================================================================================

#logistic Regression---------------------------------------------------------------------

getLR_model <- function(model_dat_, formula_, optParam = NA, VarMeasure = F){
  fit_ <- glm(formula_, family = binomial, data = model_dat_)
  # VarMeasure: cnduct the variable evaluation
  
  log_coeff <- summary(fit_)$coefficients %>% as.data.frame()
  
  Var_measure <- list()
  for (col in colnames(log_coeff)){
    Var_measure[[col]] <- as.data.frame(
      matrix(c(log_coeff[, col][2:nrow(log_coeff)]), nrow = 1, ncol = nrow(log_coeff) - 1, byrow = T,
             dimnames = list(rownames = c('col'), colnames = rownames(log_coeff[2:nrow(log_coeff), ]))))
  }
  if (VarMeasure){
    return(Var_measure)
  }else{
    return(fit_)
  }
}

getLR_Pred <- function(model_, pred_X){
  # document: get the predictive value of logistic regression model
  # model_: input model
  # pred_X: predictive variables
  pred_ <- predict(model_, newdata = pred_X, type='response')
  return(pred_)
}

#Integration model of bagging ------bootstrap aggregating bootstrap aggregating method (similar|with put back)----Random Forest-----------------------------------------------------
library(randomForest)  #random forests
modelLookup('rf')

getRF_optParam <- function(Tune.dat_, formula_, fitFuncs, predFuncs){
  grid_ <- expand.grid(
    ntree = round(seq(50, 400, length.out = 5), 0),  # 5
    mtry = round(seq(1, 7, length.out = 6), 0),  # 6
    nodesize = round(seq(2, 12, length.out = 4), 0)  # 3
  )
  
  Ret_ <- list()  # three parameters
  for (para.dat_p in names(Tune.dat_)){
    para.dat_ <- Tune.dat_[[para.dat_p]]
    X_ <- para.dat_[['X_']]
    Y_ <- as.factor(para.dat_[['Y_']])
    
    ret_ <- c()
    for (i in 1:nrow(grid_)){
      optParam <- list(ntree = grid_[i, 'ntree'], mtry = grid_[i, 'mtry'], nodesize = grid_[i, 'nodesize'])
      
      samp_i <- c(1:nrow(X_))[1:round(nrow(X_)*.7)]
      samp_i2 <- c(1:nrow(X_))[-samp_i]
      
      fit2_ <- fitFuncs(model_dat_ = data.frame(event = Y_[samp_i], X_[samp_i, ]),
                        formula_ = formula_, optParam = optParam) #-------------------model 2
      pred2_.test <- predFuncs(model_ = fit2_, pred_X = X_[samp_i2, ])
      auroc.test2_ <- pROC::roc(ifelse(Y_[samp_i2] == '1', 1, 0), as.numeric(pred2_.test),
                                direction = c('<'),  # control < case
      )
      if (auroc.test2_[["auc"]] < 0.5){
        auroc.test2_ <- pROC::roc(ifelse(Y_[samp_i2] == '1', 1, 0), as.numeric(pred2_.test),
                                  direction = c('>')  # control > case
        )
      }
      cat('Looking for optimal Auroc parameters(超参数·：', paste(optParam),')：', i, '/', nrow(grid_), '----', auroc.test2_[["auc"]][1], '\n')
      ret_ <- c(ret_, auroc.test2_[["auc"]][1])
    }
    Ret_[[para.dat_p]] <- ret_
  }
  
  tmp_ <- do.call('cbind', Ret_) %>% as.data.frame()
  write.csv(tmp_, '_tmp_plot.csv')
  tmp_ <- read.csv('_tmp_plot.csv', row.names = 1)  # first column is rownames
  tmp_2 <- apply(tmp_[,c('one', 'two', 'thr')], 1, function(x){mean(x)})
  
  cat('AUROC 最优时，参数为')
  opt_i <- which.max(tmp_2) %>% print()
  tmp_2[which.max(tmp_2)] %>% print()
  grid_[which.max(tmp_2), ] %>% print()
  
  return(list(ntree = grid_[opt_i, 'ntree'], mtry = grid_[opt_i, 'mtry'], nodesize = grid_[opt_i, 'nodesize']))
}


getRF_model <- function(model_dat_, formula_, optParam = NA, VarMeasure = F){
  # VarMeasure: Whether or not to perform variable evaluation
  
  model_dat_2 <- model_dat_
  if (VarMeasure){
    if("maternalFeeding" %in% colnames(model_dat_2)){
      model_dat_2 <- model_dat_2 %>% dplyr::mutate(
        maternalFeeding_2 = as.factor(ifelse(maternalFeeding == '2', 1, 0)),
        maternalFeeding_3 = as.factor(ifelse(maternalFeeding == '3', 1, 0)),
        maternalFeeding_4 = as.factor(ifelse(maternalFeeding == '4', 1, 0)),
      ) %>% dplyr::select(-c('maternalFeeding'))
    }
    if ("antibiotic_therapy" %in% colnames(model_dat_2)){
      model_dat_2 <- model_dat_2 %>% dplyr::mutate(
        antibiotic_therapy_2 = as.factor(ifelse(antibiotic_therapy == '2', 1, 0)),
        antibiotic_therapy_3 = as.factor(ifelse(antibiotic_therapy == '3', 1, 0))
      ) %>% dplyr::select(-c('antibiotic_therapy'))
    }
    if ("momEducation" %in% colnames(model_dat_2)){
      model_dat_2 <- model_dat_2 %>% dplyr::mutate(
        momEducation_2 = as.factor(ifelse(momEducation == '2', 1, 0)),
        momEducation_3 = as.factor(ifelse(momEducation == '3', 1, 0))
      ) %>% dplyr::select(-c('momEducation'))
    }
  }
  
  fit_tune_ <- randomForest(formula_, data = model_dat_2, 
                            ntree = optParam[['ntree']], mtry = optParam[['mtry']], 
                            nodesize = optParam[['nodesize']])
  
  if (VarMeasure){ return(list(MeanDecreaseGini = as.data.frame(t(fit_tune_[["importance"]])))) }else{
    return(fit_tune_)
  }
}

getRF_Pred <- function(model_, pred_X){
  # document: Get the predicted VALUE of the model on the new data, return vector (probability)
  # model_: incoming model
  # pred_X: predictor variable
  pred_ <- predict(model_, newdata = pred_X, type = 'prob')
  pred_ <- pred_[, 2]
  return(pred_)
}

#Integration model of boosting evolution (like for like | no put back)----XGBoost---------------------------------------------------
library(xgboost)  # gradient boosting
xgb.info <- getModelInfo('xgbTree')$xgbTree
modelLookup('xgbTree')


getXGB_optParam <- function(Tune.dat_, formula_, fitFuncs, predFuncs){
  grid_ <- expand.grid(
    nrounds = round(seq(150, 750, length.out = 3), 0),  # the biggest number of boosting rounds
    max_depth = round(seq(2, 8, length.out = 3), 0),  #the biggest depth of the tree
    eta = round(seq(.05, .75, length.out = 3), 3),  #learning rate
    gamma = round(seq(1, 10, length.out = 3), 2),  #the least loss reduction required to make a further partition on a leaf node of the tree
    colsample_bytree = round(seq(.3, 0.95, length.out = 3), 2),  
    min_child_weight = 2,  # the bigger, the conservative model
    subsample = 0.5###round(seq(.3, .75, length.out = 3), 2)  # proportion of training data to use in each boosting round
  )
  
  
  Ret_ <- list()  # three times tuning
  for (para.dat_p in names(Tune.dat_)){
    para.dat_ <- Tune.dat_[[para.dat_p]]
    X_ <- para.dat_[['X_']]
    Y_ <- para.dat_[['Y_']]
    
    ret_ <- c()
    for (i in 1:nrow(grid_)){
      
      optParam <- list(param = list(objective = 'binary:logistic', booster = 'gbtree', eval_metric = 'error',
                                    eta = grid_[i, 'eta'],
                                    max_depth = grid_[i, 'max_depth'],
                                    subsample = grid_[i, 'subsample'],
                                    colsample_bytree = grid_[i, 'colsample_bytree'],
                                    gamma = grid_[i, 'gamma']
      ),
      nrounds = grid_[i, 'nrounds'])
      
      samp_i <- c(1:nrow(X_))[1:round(nrow(X_)*.7)]
      samp_i2 <- c(1:nrow(X_))[-samp_i]
      
      fit2_ <- fitFuncs(model_dat_ = data.frame(event = Y_[samp_i], X_[samp_i, ]),
                        formula_ = formula_, optParam = optParam) #-------------------model 2
      pred2_.test <- predFuncs(model_ = fit2_, pred_X = X_[samp_i2, ])
      auroc.test2_ <- pROC::roc(ifelse(Y_[samp_i2] == '1', 1, 0), as.numeric(pred2_.test),
                                direction = c('<')  # control < case
      )
      if (auroc.test2_[["auc"]] < 0.5){
        auroc.test2_ <- pROC::roc(ifelse(Y_[samp_i2] == '1', 1, 0), as.numeric(pred2_.test),
                                  direction = c('>')  # control > case
        )
      }
      cat('Looking for optimal Auroc parameters：', i, '/', nrow(grid_), '----', auroc.test2_[["auc"]][1], '\n')
      ret_ <- c(ret_, auroc.test2_[["auc"]][1])
      
      if (i%%50 == 49){Sys.sleep(5)}
    }
    Ret_[[para.dat_p]] <- ret_
  }
  
  tmp_ <- do.call('cbind', Ret_) %>% as.data.frame()
  write.csv(tmp_, '_tmp_plot.csv')
  tmp_ <- read.csv('_tmp_plot.csv', row.names = 1)  # the first column is index
  tmp_2 <- apply(tmp_[,c('one', 'two', 'thr')], 1, function(x){mean(x)})
  
  cat('AUROC is optimal when the parameters are')
  opt_i <- which.max(tmp_2) %>% print()
  tmp_2[which.max(tmp_2)] %>% print()
  grid_[which.max(tmp_2), ] %>% print()
  
  return(list(param = list(objective = 'binary:logistic', booster = 'gbtree', eval_metric = 'error',
                           eta = grid_[opt_i, 'eta'],
                           max_depth = grid_[opt_i, 'max_depth'],
                           subsample = grid_[opt_i, 'subsample'],
                           colsample_bytree = grid_[opt_i, 'colsample_bytree'],
                           gamma = grid_[opt_i, 'gamma']
  ),
  nrounds = grid_[opt_i, 'nrounds']))
}



getXGB_model <- function(model_dat_, formula_, optParam = NA, VarMeasure = F){
  # VarMeasure: conduct the variable evaluation or not
  outcome <- as.character(formula_)[2]
  train.x_ <- as.matrix(df_strToNum(model_dat_[, -which(colnames(model_dat_) == outcome)]))
  train.y_ <- ifelse(model_dat_[[outcome]] == '0', 0, 1 )
  trainMat_ <- xgb.DMatrix(data = train.x_, 
                           label = train.y_)
  fit_ <- xgb.train(params = optParam[['param']], data = trainMat_, nrounds = optParam[['nrounds']])
  if (VarMeasure){ return(NULL) }else{ return(fit_) }
}

getXGB_Pred <- function(model_, pred_X){
  # document: get the prediction of the model on the new data, return vector (probability)
  # model_: input model
  # pred_X: predict variables
  train.x_ <- as.matrix(df_strToNum(pred_X))
  pred_ <- predict(model_, train.x_, type = 'prob')
  return(pred_)
}


# end====================================================================================================

# Fit and train the model ---------------------------------------------------------------------

ML_Ret <- list()

for (m in c('lr', 'rf', 'xgb')){  # Around 100 mins
  time1 <- Sys.time()
  if (m == 'lr'){fitFuncs <- getLR_model; predFuncs = getLR_Pred
  grid_ <- expand.grid(
    dummy1 = 1, dummy2 = 1
  )
  }else if (m == 'rf'){fitFuncs <- getRF_model; predFuncs = getRF_Pred
  grid_ <- expand.grid(
    ntree = round(seq(10, 1000, length.out = 15), 0),  # 5
    mtry = round(seq(1, 7, length.out = 6), 0),  # 6
    nodesize = round(seq(2, 12, length.out = 8), 0)  # 3
  )
  }else if (m == 'xgb'){fitFuncs <- getXGB_model; predFuncs = getXGB_Pred
  grid_ <- expand.grid(
    nrounds = round(seq(10, 750, length.out = 8), 0),
    max_depth = round(seq(1, 8, length.out = 4), 0),
    eta = round(seq(.05, .75, length.out = 4), 3),
    gamma = round(seq(1, 10, length.out = 3), 2),
    colsample_bytree = round(seq(.3, 0.95, length.out = 3), 2),  
    min_child_weight = 2,
    subsample = 0.5
  )}
  
  X_ <- df_test_uni2[, 3:ncol(df_test_uni2)]
  Y_ <- as.factor(df_test_uni2[, 'event'])
  
  ret_ <- c()
  ind_param <- createDataPartition(df_test_uni2[['event']], times = 5, p = 0.7, list = FALSE)
  for (ind_param_ in colnames(ind_param)){
    ind_test <- ind_param[, ind_param_]
    
    for (i in 1:nrow(grid_)){
      optParam <- if (T){
        if (m == 'lr'){NA
        }else if(m == 'rf'){list(ntree = grid_[i, 'ntree'], mtry = grid_[i, 'mtry'], nodesize = grid_[i, 'nodesize'])
        }else if(m == 'xgb'){list(param = list(objective = 'binary:logistic', booster = 'gbtree', eval_metric = 'error',
                                               eta = grid_[i, 'eta'],
                                               max_depth = grid_[i, 'max_depth'],
                                               subsample = grid_[i, 'subsample'],
                                               colsample_bytree = grid_[i, 'colsample_bytree'],
                                               gamma = grid_[i, 'gamma']
        ), nrounds = grid_[i, 'nrounds'])}
      }
      
      fit2_ <- fitFuncs(model_dat_ = data.frame(event = Y_[ind_test], X_[ind_test, ]),
                        formula_ = formula_, optParam = optParam) #-------------------model 2
      pred2_.test <- predFuncs(model_ = fit2_, pred_X = X_[-ind_test, ])
      auroc.test2_ <- pROC::roc(ifelse(Y_[-ind_test] == '1', 1, 0), as.numeric(pred2_.test),
                                direction = c('<'),  # control < case
      )
      if (auroc.test2_[["auc"]] < 0.5){
        auroc.test2_ <- pROC::roc(ifelse(Y_[-ind_test] == '1', 1, 0), as.numeric(pred2_.test),
                                  direction = c('>')  # control > case
        )
      }
      cat('Looking for the optimal Auroc parameter (hyperparameter-：', paste(optParam),')：', i, '/', nrow(grid_), '----', auroc.test2_[["auc"]][1], '\n')
      ret_[[ind_param_]] <- c(ret_[[ind_param_]], auroc.test2_[["auc"]][1])
    }
    
    time2 <- Sys.time()
    if ((round((as.numeric(time2) - as.numeric(time1))/60, 0)+1) %% 15 == 0){Sys.sleep(60)}
  }
  ret_2 <- do.call('cbind', ret_) %>% as.data.frame()
  ret_2 <- ret_2 %>% dplyr::mutate(
    average = apply(ret_2, MARGIN = 1, FUN = function(x){mean(x, na.rm = T)})
  )
  ML_Ret[[m]][['ret_']] <- ret_2
  ML_Ret[[m]][['optimaParam']] <- if (T){
    i_opt <- which.max(ret_2[, 'average'])
    if (m == 'lr'){NA
    }else if(m == 'rf'){list(ntree = grid_[i_opt, 'ntree'], mtry = grid_[i_opt, 'mtry'], nodesize = grid_[i_opt, 'nodesize'])
    }else if(m == 'xgb'){list(param = list(objective = 'binary:logistic', booster = 'gbtree', eval_metric = 'error',
                                           eta = grid_[i_opt, 'eta'],
                                           max_depth = grid_[i_opt, 'max_depth'],
                                           subsample = grid_[i_opt, 'subsample'],
                                           colsample_bytree = grid_[i_opt, 'colsample_bytree'],
                                           gamma = grid_[i_opt, 'gamma']
    ), nrounds = grid_[i, 'nrounds'])}
  }
}


for (m in c('lr', 'rf', 'xgb')){
  if (m == 'lr'){fitFuncs <- getLR_model; predFuncs = getLR_Pred; optParam <- NA
  }else if (m == 'rf'){fitFuncs <- getRF_model; predFuncs = getRF_Pred; optParam <- ML_Ret[["rf"]][["optimaParam"]]
  }else if (m == 'xgb'){fitFuncs <- getXGB_model; predFuncs = getXGB_Pred; optParam <- ML_Ret[["xgb"]][["optimaParam"]]}
  
  X_ <- df_test_uni2[, 3:ncol(df_test_uni2)]
  Y_ <- as.factor(df_test_uni2[, 'event'])
  
  ret_ <- list()
  ind_part <- createDataPartition(df_test_uni2[['event']], times = 100, p = 0.7, list = FALSE)
  for (ind_p in colnames(ind_part)){
    ind_test <- ind_part[, ind_p]
    formula_ <- as.formula('event ~ .')
    
    fit2_ <- fitFuncs(model_dat_ = data.frame(event = Y_[ind_test], X_[ind_test, ]),
                      formula_ = formula_, optParam = optParam) #-------------------model 2
    pred2_.train <- predFuncs(model_ = fit2_, pred_X = X_[ind_test, ])
    auroc.train2_ <- pROC::roc(ifelse(Y_[ind_test] == '1', 1, 0), as.numeric(pred2_.train),
                              direction = c('<'),  # control < case
    ); if (auroc.train2_[["auc"]] < 0.5){
      auroc.train2_ <- pROC::roc(ifelse(Y_[-ind_test] == '1', 1, 0), as.numeric(pred2_.test),
                                direction = c('>')  # control > case
      )}
    
    pred2_.test <- predFuncs(model_ = fit2_, pred_X = X_[-ind_test, ])
    auroc.test2_ <- pROC::roc(ifelse(Y_[-ind_test] == '1', 1, 0), as.numeric(pred2_.test),
                              direction = c('<'),  # control < case
    ); if (auroc.test2_[["auc"]] < 0.5){
      auroc.test2_ <- pROC::roc(ifelse(Y_[-ind_test] == '1', 1, 0), as.numeric(pred2_.test),
                                direction = c('>')  # control > case
      )}
    
    pred_X <- c(seq(0, .01, length.out = 200), seq(.01, .99, length.out = 98), seq(.99, 1, length.out = 200))
    train.pred <- get_nlrPred(
      dat = data.frame(list(spe = auroc.train2_$specificities, sen = auroc.train2_$sensitivities)), pred_X = pred_X)
    test.pred <- get_nlrPred(
      dat = data.frame(list(spe = auroc.test2_$specificities, sen = auroc.test2_$sensitivities)), pred_X = pred_X)

    cat('auroc为：', auroc.test2_[["auc"]][1], '\n')
    ret_[['roc']][['partition']] <- c(ret_[['roc']][['partition']], rep(ind_p, 2))
    ret_[['roc']][['class']] <- c(ret_[['roc']][['class']], 'train', 'test')
    ret_[['roc']][['auroc']] <- c(ret_[['roc']][['auroc']], auroc.train2_[["auc"]][1], auroc.test2_[["auc"]][1])
    ret_[['rocPlot']][['partition']] <- c(ret_[['rocPlot']][['partition']], rep(ind_p, 996))
    ret_[['rocPlot']][['class']] <- c(ret_[['rocPlot']][['class']], rep('train', 498), rep('test', 498))
    ret_[['rocPlot']][['spe']] <- c(ret_[['rocPlot']][['auroc']], pred_X, pred_X)
    ret_[['rocPlot']][['sen']] <- c(ret_[['rocPlot']][['auroc']], train.pred[["pred"]], test.pred[["pred"]])
    ret_[['rocPlot']][['R2']] <- c(ret_[['rocPlot']][['auroc']], rep(train.pred[["R2"]], 498), rep(test.pred[["R2"]], 498))
  }
  ML_Ret[[m]][['auroc']] <- do.call('cbind', ret_[['roc']]) %>% as.data.frame()
  ML_Ret[[m]][['rocPlot']] <- do.call('cbind', ret_[['rocPlot']]) %>% as.data.frame()
}

sapply(c('lr', 'rf', 'xgb'), FUN = function(x){
  mean(as.numeric(ML_Ret[[x]][["auroc"]][ML_Ret[[x]][["auroc"]][, 'class'] == 'test', 'auroc']))})

#     save(ML_Ret, file = 'ML_Ret.RData')    # ----------------------------25 --01--22

export_parameter_ml <- as.data.frame(list(
  models = names(ML_Ret),
  parameters = c('NA', pastev(paste0(names(ML_Ret[["rf"]][["optimaParam"]]), ': ',
                         as.character(ML_Ret[["rf"]][["optimaParam"]])), sep_ = ', '),
                 pastev(c('nrounds: 750', paste0(names(ML_Ret[["xgb"]][["optimaParam"]][["param"]]), ': ',
                         as.character(ML_Ret[["xgb"]][["optimaParam"]][["param"]]))), sep_ = ', '))
)); write.xlsx(export_parameter_ml, file = 'export_parameter_ml.xlsx')


# Model---ROC Plotting----------------------------------------------------
for (m in names(ML_Ret)){ML_Ret[[m]][['rocPlot']][, 'model'] <- m}
df_rocPlot_merg <- do.call('rbind', list(lr = ML_Ret[['lr']][["rocPlot"]], rf = ML_Ret[['rf']][["rocPlot"]],
                                     xgb = ML_Ret[['xgb']][["rocPlot"]])) %>% as.data.frame() %>% 
  dplyr::mutate(sen = as.numeric(sen), spe = as.numeric(spe))

roc_plot_df <- aggregate(x = df_rocPlot_merg[, 'sen'],
                         by = list(df_rocPlot_merg[, 'spe'], df_rocPlot_merg[, 'model'], df_rocPlot_merg[, 'class']),
                         FUN = function(x){c(mean(x, na.rm = T), sd(x, na.rm = T),
                                             sd(x, na.rm = T)/sqrt(length(x)))}) %>% as.data.frame()
write.csv(roc_plot_df, '_tmp_plot.csv')
roc_plot_df <- read.csv('_tmp_plot.csv', row.names = 1)
colnames(roc_plot_df) <- c('spe', 'model', 'class', 'sen_adjust', 'sd', 'se')
roc_plot_df <- roc_plot_df %>% dplyr::mutate(
  model = factor(model, levels = c('lr', 'rf', 'xgb'), labels = c('LR', 'RF', 'XGBoost')),
  class = factor(class, levels = c('train', 'test'), labels = c('train', 'test'))
)


par(mar = c(2, 2, 2, 2))
gg_RocPlot <- ggplot(data = roc_plot_df, aes(x = 1-spe, y = sen_adjust, colour = model)) +
  geom_line(size = .7, alpha = 0.85) + 
  scale_color_manual(values = c('#9775fa', '#12b886', '#fa5252', '#4dabf7', '#ff922b')) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', size = .7, color = '#adb5bd') +
  xlab('1 - Specificities') + 
  scale_y_continuous(expand=c(0,0),
                     breaks = c(.25, .50, .75, 1)) +
  scale_x_continuous(expand=c(.002,0)) +
  ylab('Sensitivities') +
  facet_grid(. ~ class) +
  facet_grid(. ~ class, scales = 'free', space = 'free_x') +
  theme(
    axis.text.x = element_text(family = 'Arial'
                               #, face = 'italic'
                               , colour = 'black'
                               , size = rel(1.5) 
    ),
    axis.text.y = element_text(family = 'Arial'
                               #, face = 'italic'
                               , colour = 'black'
                               , size = rel(1.5)),
    axis.title = element_text(family = 'Arial'
                              #, face = 'italic'
                              , colour = 'black'
                              , size = rel(1.5)),
    axis.ticks = element_line(color = "black", size = 0.6,lineend = 2),
    axis.line = element_line(linetype = 1,color = "black",size = 0.8, lineend = 'square'),
    strip.text = element_text(face = 'bold', size = rel(1.15)), 
    panel.spacing.x = unit(.4, 'cm'), 
    legend.title = element_text(
      face = 'italic',
      family = 'Arial',
      colour = 'black',
      size = rel(1.1)
    ),
    legend.text = element_text(
      face = 'italic',
      family = 'Arial',
      colour = 'black',
      size = rel(0.9)
    ),
    legend.position = c(0.92, 0.3)
  )

print(gg_RocPlot)
ggsave(filename = pastev(c('train_test_ROC', '.pdf'), sep_ = ''), plot = gg_RocPlot, 
       width = 8 * 2.54, height = 4.5*2.54, units = 'cm')


# Calculating the difference between groups
for (m in names(ML_Ret)){ML_Ret[[m]][['auroc']][, 'model'] <- m}
df_ROC_merg <- do.call('rbind', list(lr = ML_Ret[['lr']][["auroc"]], rf = ML_Ret[['rf']][["auroc"]],
                                     xgb = ML_Ret[['xgb']][["auroc"]]))
df_ROC_merg_train <- list(
  lr_rf = as.data.frame(list(group = as.factor(c(rep('lr', 100), rep('rf', 100))),
                             values = as.numeric(c(df_ROC_merg[df_ROC_merg[, 'class'] == 'train' & df_ROC_merg[, 'model'] == 'lr', 'auroc'],
                                                   df_ROC_merg[df_ROC_merg[, 'class'] == 'train' & df_ROC_merg[, 'model'] == 'rf', 'auroc'])))),
  lr_xgb = as.data.frame(list(group = as.factor(c(rep('lr', 100), rep('xgb', 100))),
                              values = as.numeric(c(df_ROC_merg[df_ROC_merg[, 'class'] == 'train' & df_ROC_merg[, 'model'] == 'lr', 'auroc'],
                                                    df_ROC_merg[df_ROC_merg[, 'class'] == 'train' & df_ROC_merg[, 'model'] == 'xgb', 'auroc'])))),
  rf_xgb = as.data.frame(list(group = as.factor(c(rep('rf', 100), rep('xgb', 100))),
                              values = as.numeric(c(df_ROC_merg[df_ROC_merg[, 'class'] == 'train' & df_ROC_merg[, 'model'] == 'rf', 'auroc'],
                                                    df_ROC_merg[df_ROC_merg[, 'class'] == 'train' & df_ROC_merg[, 'model'] == 'xgb', 'auroc']))))
)

roc_train_hypothesisTest <- get_hypothesisTesting(object = df_ROC_merg_train, reverse = F,  # For performing independent samples t-tests
                                                 ignoreNormal = T,
                                                 appoximately_normal = F,
                                                 adj_method = 'bonferroni',
                                                 exact_ = T,
                                                 maxGroup = 3, 
                                                 filename_ = paste0('ignoreNormal__roc_train__', strsplit(as.character(Sys.time()), ' ', fixed = T)[[1]][1])
)


df_ROC_merg_test <- list(
  lr_rf = as.data.frame(list(group = as.factor(c(rep('lr', 100), rep('rf', 100))),
        values = as.numeric(c(df_ROC_merg[df_ROC_merg[, 'class'] == 'test' & df_ROC_merg[, 'model'] == 'lr', 'auroc'],
                  df_ROC_merg[df_ROC_merg[, 'class'] == 'test' & df_ROC_merg[, 'model'] == 'rf', 'auroc'])))),
  lr_xgb = as.data.frame(list(group = as.factor(c(rep('lr', 100), rep('xgb', 100))),
        values = as.numeric(c(df_ROC_merg[df_ROC_merg[, 'class'] == 'test' & df_ROC_merg[, 'model'] == 'lr', 'auroc'],
                  df_ROC_merg[df_ROC_merg[, 'class'] == 'test' & df_ROC_merg[, 'model'] == 'xgb', 'auroc'])))),
  rf_xgb = as.data.frame(list(group = as.factor(c(rep('rf', 100), rep('xgb', 100))),
        values = as.numeric(c(df_ROC_merg[df_ROC_merg[, 'class'] == 'test' & df_ROC_merg[, 'model'] == 'rf', 'auroc'],
                  df_ROC_merg[df_ROC_merg[, 'class'] == 'test' & df_ROC_merg[, 'model'] == 'xgb', 'auroc']))))
)

roc_test_hypothesisTest <- get_hypothesisTesting(object = df_ROC_merg_test, reverse = F,  # For performing independent samples t-tests
                                   ignoreNormal = T, 
                                   appoximately_normal = F,
                                   adj_method = 'bonferroni',
                                   exact_ = T,
                                   maxGroup = 3,
                                   filename_ = paste0('ignoreNormal__roc_test__', strsplit(as.character(Sys.time()), ' ', fixed = T)[[1]][1])
)


# Feature importance of model（XGB & Random Forest）--------------------------------------------------

    # SHAP Value
library(shapviz)  # Build SHAP interpretable machine learning
shap_X <- df_test_uni2[, 3:ncol(df_test_uni2)]
shap_Y <- as.factor(df_test_uni2[, 'event'])
shap_ind <- createDataPartition(df_test_uni2[['event']], times = 1, p = 0.7, list = FALSE)
shap_fit <- getXGB_model(model_dat_ = data.frame(event = shap_Y[shap_ind], shap_X[shap_ind, ]),
               formula_ = as.formula('event ~ .'), optParam = ML_Ret[["xgb"]][["optimaParam"]])
pred_.shap <- getXGB_Pred(model_ = shap_fit, pred_X = shap_X[-shap_ind, ])
auroc.shap <- pROC::roc(ifelse(shap_Y[-shap_ind] == '1', 1, 0), as.numeric(pred_.shap),
                           direction = c('<')  # control < case
                        ); auroc.shap[["auc"]]

xgb.shap <- shapviz(object = shap_fit, X_pred = as.matrix(df_strToNum(shap_X)))

pdf('xgb_shap.pdf', width = 8.5, height = 4.5) 
#tiff(pastev(c('shap_',i, '.tif'), sep_ = ''), width = 12, height = 9, units = 'cm', res=300, compression = 'lzw')
sv_importance(xgb.shap, kind = "both"
              , alpha = 0.2, bar_width = .7
              , fill = '#8ce99a'
              , bee_width = 0.4
) +
  scale_color_gradient(low = '#5c7cfa', high = '#fa5252') + 
  theme_classic() +
  theme(
    axis.text.x = element_text(family = 'Arial'
                               #, face = 'italic'
                               , colour = 'black'
                               , size = rel(1.2) 
    ),
    axis.text.y = element_text(family = 'Arial'
                               #, face = 'italic'
                               , colour = '#495057'
                               , size = rel(1.2)
    ),
    axis.title = element_text(family = 'Arial' 
                              #, face = 'italic'
                              , colour = 'black'
                              , size = rel(1.3)),
    axis.ticks = element_line(color = "black", size = 0.5,lineend = 2),
    axis.line = element_line(linetype = 1,color = "black",size = 0.6, lineend = 'square'),
    legend.text = element_text(
      face = 'italic',
      family = 'Arial',
      colour = 'black',
      size = rel(0.9)
    )
  )
dev.off()


# Random Forest：

rf_fit <- getRF_model(model_dat_ = data.frame(event = shap_Y[shap_ind], shap_X[shap_ind, ]),
                         formula_ = as.formula('event ~ .'), optParam = ML_Ret[["rf"]][["optimaParam"]])
pred_.rf <- getRF_Pred(model_ = rf_fit, pred_X = shap_X[-shap_ind, ])
auroc.rf <- pROC::roc(ifelse(shap_Y[-shap_ind] == '1', 1, 0), as.numeric(pred_.shap),
                        direction = c('<')  # control < case
); auroc.rf[["auc"]]
rf_importance_df <- as.data.frame(list(var = rownames(rf_fit[["importance"]]), values = as.numeric(rf_fit[["importance"]])))


gg_rf <- ggplot(data = rf_importance_df, aes(x = values, y = reorder(var, values))) +
  geom_point(size = 3, color = '#da77f2') + 
  theme_bw() +
  ylab(NULL) +
  xlab('Importance of variables\n(MeanDecreaseGini)') +
  theme_classic() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey60', linetype = 'dashed'),
    axis.text.x = element_text(family = 'Arial'
                               #, face = 'italic'
                               , colour = 'black'
                               , size = rel(1.5) 
    ),
    axis.text.y = element_text(family = 'Arial'
                               #, face = 'italic'
                               , colour = 'black'
                               , size = rel(1.5)),
    axis.title = element_text(family = 'Arial'  
                              #, face = 'italic'
                              , colour = 'black'
                              , size = rel(1.5)),
    axis.ticks = element_line(color = "black", size = 0.6,lineend = 2),
    axis.line = element_line(linetype = 1,color = "black",size = 0.8, lineend = 'square')
  )
ggsave(filename = 'gg_rf.pdf', plot = gg_rf, width = 8 * 2.54, height = 4.5*2.54, units = 'cm')

# export the data using in python ----------------------------------------------------------------------------

df_ml_python <- df_test_uni2[df_test_uni2[, 'match'] %in% c(1:550), ]
write.xlsx(df_ml_python, file = 'df_ml_python.xlsx')

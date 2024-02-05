library(haven)
library(naniar)
library(tidyverse)
rm(list=ls())

###------------------------------------------------------------------------

rescale <- function(ref_base, value) {
  ref_value<-summary(value)
  v_rescaled <- as.data.frame(value) %>%
    mutate(Qn = ifelse(value>=ref_value[[5]],ref_value[[5]],
                       ifelse(value>=ref_value[[3]],ref_value[[3]],
                              ifelse(value>=ref_value[[2]],ref_value[[2]],
                                     ref_value[[1]]))),
           Qn1 = ifelse(value>=ref_value[[5]],ref_value[[6]],
                        ifelse(value>=ref_value[[3]],ref_value[[5]],
                               ifelse(value>=ref_value[[2]],ref_value[[3]],
                                      ref_value[[2]]))),
           rQn = ifelse(value>=ref_value[[5]],ref_base[[5]],
                        ifelse(value>=ref_value[[3]],ref_base[[3]],
                               ifelse(value>=ref_value[[2]],ref_base[[2]],
                                      ref_base[[1]]))),
           rQn1 = ifelse(value>=ref_value[[5]],ref_base[[6]],
                         ifelse(value>=ref_value[[3]],ref_base[[5]],
                                ifelse(value>=ref_value[[2]],ref_base[[3]],
                                       ref_base[[2]]))),
           Y = ifelse(Qn1!=Qn,(rQn1-rQn)/(Qn1-Qn),0),
           rescaled = ifelse(Y==0,
                             (rQn+rQn1)/2,
                             rQn1 + (Y*(value-Qn1)))) %>%
    
    select(rescaled)
  return(unlist(v_rescaled))
}

###------------------------------------------------------------------------

demog_clean<- read_sas("Raw Data/demog_clean.sas7bdat")   %>%
  select(-SUBJECT, -TIME_SINCE_STROKE_UNITS) %>%
  mutate(DIABETES = ifelse(DIABETES %in% c("YES","NO"),
                           DIABETES,NA),
         APHASIA_BASE = ifelse(APHASIA_BASE=="",NA,APHASIA_BASE))

eq5d_comp <- read_sas("Raw Data/eq5d_comp_revised.sas7bdat")  %>%
  rename(TIME = ASSESSMENT_TIME) %>%
  mutate(TIME_cat = ifelse(TIME<=28, "0-4 WEEKS",
                           ifelse(TIME<=90,"5 WEEKS - 3 MONTHS",
                                  ifelse(TIME<=180, "4-6 MONTHS",
                                         ">6 MONTHS"))),
         SELF_CARE = ifelse(SELF_CARE==9,NA,SELF_CARE),
         USUAL_ACTIVITIES = ifelse(USUAL_ACTIVITIES==9,NA,USUAL_ACTIVITIES),
         ANXIETY = ifelse(ANXIETY==9,NA,ANXIETY), 
         ANXIETY_bin = ANXIETY>1,
         VAS = ifelse(VAS>100,NA,VAS)) %>%
  select(-SUBJECT, -EQ5D_TIME) %>%
  filter(!is.na(PAIN)) 
eq5d_3<-eq5d_comp %>% filter(EQ5D=="EQ5D-3L") %>% select(-EQ5D) %>%
  rename(MOBILITY3 = MOBILITY)
# all measures on 1-3 scale
eq5d_5<-eq5d_comp %>% filter(EQ5D=="EQ5D-5L") %>% select(-EQ5D)  %>%
  rename(MOBILITY5 = MOBILITY)
# all measures on 1-5 scale
# EQ5D_TIME (etc) is time in trial, and ASSESSMENT_TIME is time since stroke onset
#  VAS is the overall perception of health

times<-unique(eq5d_comp$TIME_cat)

temp <- read_sas("Raw Data/mobility_comp_update.sas7bdat")  %>%
  rename(TIME = ASSESSMENT_TIME,
         MOBILITY = BI_MOBILITY) %>%
  select(-SUBJECT,-BI_MOBILITY_TIME,-BI_TOTAL_TIME,-BI_TOT_ASSESSMENT_TIME,OTT) %>%
  filter(TIME != "") %>%
  mutate(MOBILITY = ifelse(MOBILITY==99 | TRIAL=="BOTULS",NA,
                           ifelse(MOBILITY<4,MOBILITY*5,
                                  MOBILITY)),
         MOB_TRANSFORMED = ifelse(MOBILITY==0,1,
                                  ifelse(MOBILITY==10,2,
                                         3)),
         BI_TOTAL = ifelse(BI=="0-20",
                           BI_TOTAL*5,
                           BI_TOTAL),
         TIME_cat = ifelse(TIME<=28, "0-4 WEEKS",
                           ifelse(TIME<=90,"5 WEEKS - 3 MONTHS",
                                  ifelse(TIME<=180, "4-6 MONTHS",
                                         ">6 MONTHS")))) 

BI_total<-temp %>% select(TRIAL, UNIQUE_ID,BI_TOTAL, TIME,TIME_cat) %>% filter(!is.na(BI_TOTAL))
mobility_comp<-temp %>% select(TRIAL, UNIQUE_ID,MOBILITY, MOB_TRANSFORMED,TIME,TIME_cat) %>% filter(!is.na(MOBILITY))
rm(temp)

independence <- read_sas("Raw Data/independence_comp.sas7bdat") %>%
  select(-SUBJECT) %>%
  filter(INDEPENDENCE!=9) %>%
  rename(TIME = ASSESSMENT_TIME) %>%
  mutate(TIME_cat = ifelse(TIME<=28, "0-4 WEEKS",
                           ifelse(TIME<=90,"5 WEEKS - 3 MONTHS",
                                  ifelse(TIME<=180, "4-6 MONTHS",
                                         ">6 MONTHS"))))


pain_scale <- read_sas("Raw Data/pain_scale.sas7bdat")  %>%
  rename(TIME = ASSESSMENT_TIME) %>%
  select(-SUBJECT,-PAIN0_10_TIME,-OTT)  %>%
  mutate(PAIN0_10 = ifelse(PAIN0_10==99,NA,PAIN0_10),
         TIME_cat = ifelse(TIME<=28, "0-4 WEEKS",
                           ifelse(TIME<=90,"5 WEEKS - 3 MONTHS",
                                  ifelse(TIME<=180, "4-6 MONTHS",
                                         ">6 MONTHS"))))
# measure on 1-10 scale 
# OTT is Onset to treatment time, basically a surrogate marker for time from stroke to enrolment in the study.


sf36_21_paindomain<- read_sas("Raw Data/sf36_21_paindomain.sas7bdat")  %>%
  mutate(sf_pain = (P_BODY-1)) %>%
  select(UNIQUE_ID,sf_pain,TRIAL) 
  
sf36_emotion<- read_sas("Raw Data/sf36_emotion.sas7bdat") %>%
  select(UNIQUE_ID,sf_emotion,TRIAL)

sf36_original_scores<- read_sas("Raw Data/sf36_original_scores.sas7bdat") %>%
  select(-c(SUBJECT,BASELINE,GH_COMP,GH_NORMAL,P_BODY,P_INTERFERE,
            EE_NERV, EE_DUMPS, EE_CALM,EE_BLUE, EE_HAPPY, 
            SF36_TIME, ASSESSMENT_TIME)) %>%
  filter(!is.na(GH_NOW))

sf_physical_functioning<-sf36_original_scores %>%
  select(c(UNIQUE_ID,TRIAL, ACT_VIG:ACT_BATH)) %>%
  gather("measure","score",-UNIQUE_ID) %>%
  filter(!is.na(score)) %>%
  mutate(score_clean = ifelse(score==1,0,
                              ifelse(score==2,50,
                                     100))) %>%
  group_by(UNIQUE_ID) %>%
  summarise(physical_functioning = mean(score_clean))

sf_physical_health<-sf36_original_scores %>%
  select(c(UNIQUE_ID,TRIAL, PH_TIME:PH_DIFF)) %>%
  gather("measure","score",-UNIQUE_ID) %>%
  filter(!is.na(score)) %>%
  mutate(score_clean = ifelse(score==1,0,100)) %>%
  group_by(UNIQUE_ID) %>%
  summarise(physical_health = mean(score_clean))

sf_emotional_health<-sf36_original_scores %>%
  select(c(UNIQUE_ID,TRIAL, EH_TIME:EH_CARE)) %>%
  gather("measure","score",-UNIQUE_ID) %>%
  filter(!is.na(score)) %>%
  mutate(score_clean = ifelse(score==1,0,100)) %>%
  group_by(UNIQUE_ID) %>%
  summarise(emotional_health = mean(score_clean))

sf_social_functioning<-sf36_original_scores %>%
  select(UNIQUE_ID,TRIAL, SA_INTERFERE,SOCIAL_ACT)  %>%
  mutate(SA_INTERFERE = ifelse(SA_INTERFERE==1,100,
                            ifelse(SA_INTERFERE==2,75,
                                   ifelse(SA_INTERFERE==3,50,
                                          ifelse(SA_INTERFERE==4,25,
                                                 0)))),
         SOCIAL_ACT = ifelse(SOCIAL_ACT==1,0,
                            ifelse(SOCIAL_ACT==2,25,
                                   ifelse(SOCIAL_ACT==3,50,
                                          ifelse(SOCIAL_ACT==4,75,
                                                 100))))) %>%
  gather("measure","score",-UNIQUE_ID, -TRIAL) %>%
  filter(!is.na(score)) %>%
  group_by(UNIQUE_ID) %>%
  summarise(social_functioning = mean(score))

sf_general_health<-sf36_original_scores %>%
  select(UNIQUE_ID,TRIAL, GH_NOW, GH_ILL,GH_WORSE,GH_EXCELLENT) %>%
  mutate(GH_NOW = ifelse(GH_NOW==1,100,
                            ifelse(GH_NOW==2,75,
                                   ifelse(GH_NOW==3,50,
                                          ifelse(GH_NOW==4,25,
                                                 0)))),
         GH_EXCELLENT = ifelse(GH_EXCELLENT==1,100,
                                  ifelse(GH_EXCELLENT==2,75,
                                         ifelse(GH_EXCELLENT==3,50,
                                                ifelse(GH_EXCELLENT==4,25,
                                                       0)))),
         GH_ILL = ifelse(GH_ILL==1,0,
                            ifelse(GH_ILL==2,25,
                                   ifelse(GH_ILL==3,50,
                                          ifelse(GH_ILL==4,75,
                                                 100)))),
         GH_WORSE = ifelse(GH_WORSE==1,0,
                            ifelse(GH_WORSE==2,25,
                                   ifelse(GH_WORSE==3,50,
                                          ifelse(GH_WORSE==4,75,
                                                 100))))) %>%
  gather("measure","score",-UNIQUE_ID,-TRIAL) %>%
  filter(!is.na(score)) %>%
  group_by(UNIQUE_ID) %>%
  summarise(general_health = mean(score))

sf_fatigue<-sf36_original_scores %>%
  select(UNIQUE_ID,TRIAL, EE_PEP,EE_ENERGY,EE_WORN,EE_TIRED) %>%
  mutate(EE_PEP = ifelse(EE_PEP==1,100,
                          ifelse(EE_PEP==2,80,
                                 ifelse(EE_PEP==3,60,
                                        ifelse(EE_PEP==4,40,
                                               ifelse(EE_PEP==5,20,
                                                      0))))),
         EE_ENERGY = ifelse(EE_ENERGY==1,100,
                               ifelse(EE_ENERGY==2,80,
                                      ifelse(EE_ENERGY==3,60,
                                             ifelse(EE_ENERGY==4,40,
                                                    ifelse(EE_ENERGY==5,20,
                                                           0))))),
         EE_WORN = ifelse(EE_WORN==1,0,
                             ifelse(EE_WORN==2,20,
                                    ifelse(EE_WORN==3,40,
                                           ifelse(EE_WORN==4,60,
                                                  ifelse(EE_WORN==5,80,
                                                         100))))),
         EE_TIRED = ifelse(EE_TIRED==1,0,
                              ifelse(EE_TIRED==2,20,
                                     ifelse(EE_TIRED==3,40,
                                            ifelse(EE_TIRED==4,60,
                                                   ifelse(EE_TIRED==5,80,
                                                          100)))))) %>%
  gather("measure","score",-UNIQUE_ID,-TRIAL) %>%
  filter(!is.na(score)) %>%
  group_by(UNIQUE_ID) %>%
  summarise(fatigue = mean(score))

sf36<-full_join(sf_emotional_health,sf_general_health)
sf36<-full_join(sf36,sf_physical_functioning)
sf36<-full_join(sf36,sf_physical_health)
sf36<-full_join(sf36,sf_social_functioning)
sf36<-full_join(sf36,sf_fatigue)
sf36<-full_join(sf36,sf36_21_paindomain)
sf36<-full_join(sf36,sf36_emotion)
rm(sf_emotional_health,sf_general_health,sf_physical_functioning,sf_physical_health,
   sf_social_functioning,sf_fatigue,sf36_original_scores,sf36_21_paindomain,sf36_emotion)

##############################################################################
###  Rescaling Pain
##############################################################################

ref_base<-summary(unlist(eq5d_3 %>% select(PAIN)))

pain_all<- pain_scale %>%
  rename(PAIN = PAIN0_10) %>%
  mutate(Measure = "Pain_scale",
         PAIN_rescaled = rescale(ref_base,PAIN)) %>%
  bind_rows(eq5d_3 %>% select(UNIQUE_ID,PAIN,TIME,TIME_cat,TRIAL,RESPONDENT) %>% 
              mutate(Measure = "EQ5D_3", 
                     PAIN_rescaled=PAIN)) %>%
  bind_rows(eq5d_5 %>% select(UNIQUE_ID,PAIN,TIME,TIME_cat,TRIAL) %>%  
              mutate(Measure = "EQ5D_5",
                     PAIN_rescaled = rescale(ref_base,PAIN))) %>%
  bind_rows(sf36 %>% select(UNIQUE_ID,sf_pain,TRIAL) %>%
              rename(PAIN = sf_pain) %>%
              mutate(Measure = "SF36", 
                     TIME=90, 
                     TIME_cat = "5 WEEKS - 3 MONTHS",
                     PAIN_rescaled = rescale(ref_base,PAIN))) %>%
  filter(!is.na(PAIN)) %>%
  mutate(RESP_proxy = ifelse(Measure=="EQ5D_3",
                             ifelse(is.na(RESPONDENT) | RESPONDENT=="","UNKNOWN",
                                    ifelse(RESPONDENT %in% c("SUBJECT","PATIENT"),
                                           "SELF",
                                           "PROXY")),
                             NA))

population<-unique(pain_all$UNIQUE_ID) 
rm(ref_base)

BI_total <- BI_total %>% filter(UNIQUE_ID %in% population)
demog_clean <- demog_clean %>% filter(UNIQUE_ID %in% population)
independence <- independence %>% filter(UNIQUE_ID %in% population)
mobility_comp <- mobility_comp %>% filter(UNIQUE_ID %in% population)
sf36 <- sf36 %>% filter(UNIQUE_ID %in% population)


##############################################################################
###  Rescaling Mobility
##############################################################################

ref_base<-summary(unlist(mobility_comp %>% filter(!is.na(MOB_TRANSFORMED)) %>% select(MOB_TRANSFORMED)))

mobility_all<-mobility_comp %>% 
  mutate(Measure = "BI",
         Mobility_rescaled = MOB_TRANSFORMED) %>%
  bind_rows(eq5d_3 %>% select(UNIQUE_ID,TRIAL, MOBILITY3,TIME,TIME_cat) %>% 
              rename(MOBILITY = MOBILITY3) %>% 
              mutate(Measure = "EQ5D_3", 
                     MOBILITY_rev = ifelse(MOBILITY==1,3,
                                           ifelse(MOBILITY==3,1,
                                                  MOBILITY)), 
                     Mobility_rescaled = rescale(ref_base,MOBILITY_rev))) %>%
  bind_rows(eq5d_5 %>% select(UNIQUE_ID,MOBILITY5,TIME,TIME_cat,TRIAL) %>%  
              rename(MOBILITY = MOBILITY5) %>% 
              mutate(Measure = "EQ5D_5",
                     MOBILITY_rev = ifelse(MOBILITY==5,1,
                                           ifelse(MOBILITY==1,5,
                                                  ifelse(MOBILITY==2,4,
                                                         ifelse(MOBILITY==4,2,
                                                                MOBILITY)))),
                     Mobility_rescaled = rescale(ref_base,MOBILITY_rev))) %>%
  filter(!is.na(MOBILITY)) 
rm(ref_base,mobility_comp)

##############################################################################
###  Final Set Up
##############################################################################

pain_cor<-pain_all %>% 
  group_by(UNIQUE_ID,TIME_cat,Measure) %>%
  summarise(PAIN = median(PAIN)) %>% 
  spread(Measure, PAIN) %>%
  ungroup %>% 
  left_join(independence %>% 
              select(-TIME, -TRIAL) %>% 
              group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(INDEPENDENCE = median(INDEPENDENCE))) %>%
  left_join(BI_total %>% 
              select(-TIME, -TRIAL) %>% 
              group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(BI_TOTAL = median(BI_TOTAL))) %>%
  left_join(mobility_all %>%
              select(UNIQUE_ID, TIME_cat, Mobility_rescaled) %>%
              group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(Mobility_rescaled = median(Mobility_rescaled),
                        Mobility_flag = 1*(sum(Mobility_rescaled!=3)>0))) %>%
  left_join(demog_clean)

rm(list=setdiff(ls(), c("pain_cor", "eq5d_3", "eq5d_5", "pain_all", "pain_scale","rescale")))

##############################################################################
###  Demographics (Table 1)
############################################################################## 

demog<-pain_cor %>%
  group_by(UNIQUE_ID) %>%
  distinct(AGE, SEX, BNIH, APHASIA_BASE, BEST_LANGUAGE_BASE) 
summary(demog$AGE)
table(demog$SEX)
table(demog$APHASIA_BASE, useNA = "ifany")
table(demog$BEST_LANGUAGE_BASE)
summary(demog$BNIH)

tapply(pain_all$PAIN, pain_all$RESP_proxy, summary)
wilcox.test(unlist(pain_all %>% filter(RESP_proxy=="PROXY" & Measure=="EQ5D_3") %>% select(PAIN)), 
            unlist(pain_all %>% filter(RESP_proxy=="SELF" & Measure=="EQ5D_3") %>% select(PAIN)), 
            alternative = "two.sided")
rm(pain_all)

##############################################################################
###  Validation of pain transformations in studies with both EQ5D-3 and 0-10 pain
##############################################################################

ref_base<-summary(unlist(eq5d_3 %>% select(PAIN)))
pain_validation<-pain_scale %>% 
  mutate(PAIN_rescaled = rescale(ref_base,PAIN0_10)) %>%
  select(-TRIAL,-TIME_cat) %>%
  left_join(eq5d_3 %>% select(UNIQUE_ID,PAIN,TIME)) %>%
  filter(!is.na(PAIN) & !is.na(PAIN_rescaled))
summary(pain_validation$PAIN_rescaled-pain_validation$PAIN)
cor(pain_validation$PAIN_rescaled, pain_validation$PAIN, method="pearson")
cor(pain_validation$PAIN_rescaled, pain_validation$PAIN, method="spearman")
lm<-lm(PAIN_rescaled ~ PAIN, pain_validation)
round(summary(lm)$r.squared,3)

pain_validation2<-pain_scale %>% 
  mutate(PAIN_rescaled = rescale(ref_base,PAIN0_10)) %>%
  select(-TRIAL,-TIME_cat) %>%
  left_join(eq5d_5 %>% select(UNIQUE_ID,PAIN,TIME)) %>%
  mutate(PAIN_rescaled_5 = rescale(ref_base,PAIN)) %>%
  filter(!is.na(PAIN_rescaled_5) & !is.na(PAIN_rescaled))
summary(pain_validation2$PAIN_rescaled-pain_validation2$PAIN_rescaled_5)
cor(pain_validation2$PAIN_rescaled_5, pain_validation2$PAIN, method="pearson")
cor(pain_validation2$PAIN_rescaled_5, pain_validation2$PAIN, method="spearman")
lm<-lm(PAIN_rescaled_5 ~ PAIN, pain_validation2)
round(summary(lm)$r.squared,3)

##############################################################################
###  Correlation between pain measures
############################################################################## 

temp<-pain_cor %>% filter(!is.na(Pain_scale) & !is.na(EQ5D_3)) 
nrow(temp)
round(cor.test(temp$EQ5D_3,temp$Pain_scale, method = "pearson")$estimate,3)
round(cor.test(temp$EQ5D_3,temp$Pain_scale, method = "pearson")$p.value,3)
lm<-lm(EQ5D_3 ~ Pain_scale, temp)
round(sqrt(summary(lm)$r.squared),3)
round(summary(lm)$r.squared,3)
round(cor(temp$EQ5D_3,temp$Pain_scale, method = "spearman"),3)
lm<-lm(EQ5D_3 ~ SEX + AGE + Pain_scale, temp)
round(summary(lm)$r.squared,3)

temp<-pain_cor %>% filter(!is.na(Pain_scale) & !is.na(EQ5D_5)) 
nrow(temp)
round(cor.test(temp$EQ5D_5,temp$Pain_scale, method = "pearson")$estimate,3)
round(cor.test(temp$EQ5D_5,temp$Pain_scale, method = "pearson")$p.value,3)
lm<-lm(EQ5D_5 ~ Pain_scale, temp)
round(sqrt(summary(lm)$r.squared),3)
round(summary(lm)$r.squared,3)
round(cor(temp$EQ5D_5,temp$Pain_scale, method = "spearman"),3)
lm<-lm(EQ5D_5 ~ SEX + AGE + Pain_scale, temp)
round(summary(lm)$r.squared,3)

##############################################################################
### Table 2:  Correlation between pain and independence
##############################################################################

temp<-pain_cor %>% filter(!is.na(SF36) & !is.na(INDEPENDENCE))
nrow(temp)
round(cor.test(temp$SF36,temp$INDEPENDENCE)$estimate,3)
round((cor.test(temp$SF36,temp$INDEPENDENCE)$estimate)^2,3)
round(cor.test(temp$SF36,temp$INDEPENDENCE)$p.value,3)
lm<-lm(SF36 ~ SEX + AGE + INDEPENDENCE, temp)
round(summary(lm)$r.squared,3)
round(cor.test(temp$SF36,temp$INDEPENDENCE,method='spearman')$estimate,3)
round(cor.test(temp$SF36,temp$INDEPENDENCE,method='spearman')$p.value,3)


temp<-pain_cor %>% filter(!is.na(EQ5D_3) & !is.na(INDEPENDENCE))
nrow(temp)
round(cor.test(temp$EQ5D_3,temp$INDEPENDENCE)$estimate,3)
round((cor.test(temp$EQ5D_3,temp$INDEPENDENCE)$estimate)^2,3)
round(cor.test(temp$EQ5D_3,temp$INDEPENDENCE)$p.value,3)
lm<-lm(EQ5D_3 ~ SEX + AGE + INDEPENDENCE, temp)
round(summary(lm)$r.squared,3)
round(cor.test(temp$EQ5D_3,temp$INDEPENDENCE,method='spearman')$estimate,3)
round(cor.test(temp$EQ5D_3,temp$INDEPENDENCE,method='spearman')$p.value,3)

temp<-pain_cor %>% filter(!is.na(SF36) & !is.na(BI_TOTAL))
nrow(temp)
round(cor.test(temp$SF36,temp$BI_TOTAL)$estimate,3)
round((cor.test(temp$SF36,temp$BI_TOTAL)$estimate)^2,3)
round(cor.test(temp$SF36,temp$BI_TOTAL)$p.value,3)
lm<-lm(SF36 ~ SEX + AGE + BI_TOTAL, temp)
round(summary(lm)$r.squared,3)
round(cor.test(temp$SF36,temp$BI_TOTAL,method='spearman')$estimate,3)
round(cor.test(temp$SF36,temp$BI_TOTAL,method='spearman')$p.value,3)

temp<-pain_cor %>% filter(!is.na(EQ5D_3) & !is.na(BI_TOTAL))
nrow(temp)
round(cor.test(temp$EQ5D_3,temp$BI_TOTAL)$estimate,3)
round((cor.test(temp$EQ5D_3,temp$BI_TOTAL)$estimate)^2,3)
round(cor.test(temp$EQ5D_3,temp$BI_TOTAL)$p.value,3)
lm<-lm(EQ5D_3 ~ SEX + AGE + BI_TOTAL, temp)
round(summary(lm)$r.squared,3)
round(cor.test(temp$EQ5D_3,temp$BI_TOTAL,method='spearman')$estimate,3)
round(cor.test(temp$EQ5D_3,temp$BI_TOTAL,method='spearman')$p.value,3)

temp<-pain_cor %>% filter(!is.na(EQ5D_5) & !is.na(BI_TOTAL))
nrow(temp)
round(cor.test(temp$EQ5D_5,temp$BI_TOTAL)$estimate,3)
round((cor.test(temp$EQ5D_5,temp$BI_TOTAL)$estimate)^2,3)
round(cor.test(temp$EQ5D_5,temp$BI_TOTAL)$p.value,3)
lm<-lm(EQ5D_5 ~ SEX + AGE + BI_TOTAL, temp)
round(summary(lm)$r.squared,3)
round(cor.test(temp$EQ5D_5,temp$BI_TOTAL,method='spearman')$estimate,3)
round(cor.test(temp$EQ5D_5,temp$BI_TOTAL,method='spearman')$p.value,3)

temp<-pain_cor %>% filter(!is.na(Pain_scale) & !is.na(BI_TOTAL))
nrow(temp)
round(cor.test(temp$Pain_scale,temp$BI_TOTAL)$estimate,3)
round((cor.test(temp$Pain_scale,temp$BI_TOTAL)$estimate)^2,3)
round(cor.test(temp$Pain_scale,temp$BI_TOTAL)$p.value,3)
lm<-lm(Pain_scale ~ SEX + AGE + BI_TOTAL, temp)
round(summary(lm)$r.squared,3)
round(cor.test(temp$Pain_scale,temp$BI_TOTAL,method='spearman')$estimate,3)
round(cor.test(temp$Pain_scale,temp$BI_TOTAL,method='spearman')$p.value,3)


##############################################################################
###  EQ5D Scores
############################################################################## 

eq5d_3_cor<-pain_cor %>%
  filter(!is.na(EQ5D_3) & !is.na(Pain_scale)) %>%
  left_join(eq5d_3 %>%
              select(UNIQUE_ID, TIME_cat, MOBILITY3,SELF_CARE, USUAL_ACTIVITIES,VAS,ANXIETY) %>%
              group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(MOBILITY3 = median(MOBILITY3),
                        SELF_CARE = median(SELF_CARE),
                        USUAL_ACTIVITIES = median(USUAL_ACTIVITIES),
                        VAS = median(VAS),
                        ANXIETY = median(ANXIETY))) 

round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + ANXIETY, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + USUAL_ACTIVITIES, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + SELF_CARE, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3 + ANXIETY, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3 + USUAL_ACTIVITIES, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3 + SELF_CARE, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + ANXIETY + USUAL_ACTIVITIES, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + ANXIETY + SELF_CARE, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + USUAL_ACTIVITIES + SELF_CARE, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3 + USUAL_ACTIVITIES + SELF_CARE, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3 + ANXIETY + SELF_CARE, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3 + ANXIETY + USUAL_ACTIVITIES, eq5d_3_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3 + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, eq5d_3_cor))$r.squared),3)

round(sqrt(summary(lm(EQ5D_3 ~ Pain_scale + MOBILITY3 + ANXIETY + USUAL_ACTIVITIES + SELF_CARE + VAS, eq5d_3_cor))$r.squared),3)


eq5d_5_cor<-pain_cor %>%
  filter(!is.na(EQ5D_5) & !is.na(Pain_scale)) %>%
  left_join(eq5d_5 %>%
              select(UNIQUE_ID, TIME_cat, MOBILITY5,SELF_CARE, USUAL_ACTIVITIES,VAS,ANXIETY) %>%
              group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(MOBILITY5 = median(MOBILITY5),
                        SELF_CARE = median(SELF_CARE),
                        USUAL_ACTIVITIES = median(USUAL_ACTIVITIES),
                        VAS = median(VAS),
                        ANXIETY = median(ANXIETY))) 

round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + ANXIETY, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + USUAL_ACTIVITIES, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + SELF_CARE, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5 + ANXIETY, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5 + USUAL_ACTIVITIES, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5 + SELF_CARE, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + ANXIETY + USUAL_ACTIVITIES, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + ANXIETY + SELF_CARE, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + USUAL_ACTIVITIES + SELF_CARE, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5 + USUAL_ACTIVITIES + SELF_CARE, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5 + ANXIETY + SELF_CARE, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5 + ANXIETY + USUAL_ACTIVITIES, eq5d_5_cor))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5 + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, eq5d_5_cor))$r.squared),3)

round(sqrt(summary(lm(EQ5D_5 ~ Pain_scale + MOBILITY5 + ANXIETY + USUAL_ACTIVITIES + SELF_CARE + VAS, eq5d_5_cor))$r.squared),3)

#############################################################################
###  Composite Pain Measure
############################################################################## 

round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY3)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + ANXIETY)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + USUAL_ACTIVITIES)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY3 + ANXIETY)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY3 + USUAL_ACTIVITIES)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY3 + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + ANXIETY + USUAL_ACTIVITIES)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + ANXIETY + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                              mutate(Composite = Pain_scale + USUAL_ACTIVITIES + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + ANXIETY + USUAL_ACTIVITIES + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY3 + USUAL_ACTIVITIES + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY3 + ANXIETY + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY3 + ANXIETY + USUAL_ACTIVITIES)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_3 ~ Composite, eq5d_3_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY3 + ANXIETY + USUAL_ACTIVITIES + SELF_CARE)))$r.squared),3)



round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY5)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + ANXIETY)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + USUAL_ACTIVITIES)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY5 + ANXIETY)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY5 + USUAL_ACTIVITIES)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY5 + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + ANXIETY + USUAL_ACTIVITIES)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + ANXIETY + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + USUAL_ACTIVITIES + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + ANXIETY + USUAL_ACTIVITIES + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY5 + USUAL_ACTIVITIES + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY5 + ANXIETY + SELF_CARE)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY5 + ANXIETY + USUAL_ACTIVITIES)))$r.squared),3)
round(sqrt(summary(lm(EQ5D_5 ~ Composite, eq5d_5_cor %>% 
                        mutate(Composite = Pain_scale + MOBILITY5 + ANXIETY + USUAL_ACTIVITIES + SELF_CARE)))$r.squared),3)
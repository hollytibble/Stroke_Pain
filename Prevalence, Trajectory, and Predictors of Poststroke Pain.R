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
  mutate(TIME_cat = ifelse(ASSESSMENT_TIME<=28, "0-4 WEEKS",
                           ifelse(ASSESSMENT_TIME<=90,"5 WEEKS - 3 MONTHS",
                                  ifelse(ASSESSMENT_TIME<=180, "4-6 MONTHS",
                                         ">6 MONTHS")))) %>%
  group_by(UNIQUE_ID,TIME_cat) %>%
  summarise(INDEPENDENCE_med = median(INDEPENDENCE, na.rm=T),
            w_independence = max(INDEPENDENCE, na.rm=T)) 


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
  mutate(sf_pain = P_BODY) %>%
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
              mutate(Measure = "EQ5D_3L", 
                     PAIN_rescaled=PAIN)) %>%
  bind_rows(eq5d_5 %>% select(UNIQUE_ID,PAIN,TIME,TIME_cat,TRIAL) %>%  
              mutate(Measure = "EQ5D_5L",
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
###  Final Set Up Bits
##############################################################################

pain_graph<- pain_all %>%
  mutate(TIME_cat = ifelse(TIME<=180,TIME_cat,
                           ifelse(TIME<=365, "6 - 12 MONTHS",
                                  ifelse(TIME<=730,"1 - 2 YEARS",
                                         ">2 YEARS"))),
         TRIAL_c = ifelse(TRIAL %in% c("SAINT I", "SAINT II","ECASSII","CLEAR 3"),
                          "ACUTE",
                          ifelse(TRIAL %in% c("RATULS","BOTULS"),
                                 "CHRONIC",
                                 "MIXED")))
pain_graph$TIME_cat<-factor(pain_graph$TIME_cat, 
                            levels = c("0-4 WEEKS","5 WEEKS - 3 MONTHS",
                                       "4-6 MONTHS","6 - 12 MONTHS",
                                       "1 - 2 YEARS", ">2 YEARS"))
pain_graph$Anon<-paste0("Study ",as.numeric(as.factor(pain_graph$TRIAL)))
times2<-setdiff(unique(pain_graph$TIME_cat),NA)
p_measures<-unique(pain_graph$Measure)

demographics <- demog_clean %>%
  left_join(mobility_all %>% group_by(UNIQUE_ID) %>% 
              summarise(immobility = min(Mobility_rescaled,na.rm=T))) %>%
  mutate(TRIAL_c = ifelse(TRIAL %in% c("SAINT I", "SAINT II","ECASSII","CLEAR 3"),
                          "ACUTE",
                          ifelse(TRIAL %in% c("RATULS","BOTULS"),
                                 "CHRONIC",
                                 "MIXED")))

demog_clean<-demog_clean %>%
  mutate(APHASIA_BASE = ifelse(APHASIA_BASE=="",NA,APHASIA_BASE),
         TRIAL_c2 = ifelse(TRIAL %in% c("SAINT I", "SAINT II","ECASSII","CLEAR 3"),
                           "ACUTE",
                           "NON-ACUTE"),
         DIABETES = factor(ifelse(is.na(DIABETES),"MISSING",DIABETES),
                           levels=c("NO","MISSING","YES")),
         APHASIA_BASE = factor(ifelse(is.na(APHASIA_BASE),"MISSING",APHASIA_BASE),
                               levels=c("NO","MISSING","YES")),
         BNIH_c = factor(ifelse(is.na(BNIH),"MISSING",
                                ifelse(BNIH<=10,"0-10",
                                       ifelse(BNIH<=20,"11-20",
                                              ifelse(BNIH<=30,"21-30",
                                                     "31-40")))),
                         levels=c("0-10","11-20","21-30", "31-40","MISSING")))



##############################################################################
###  Validate pain transformations in studies with both EQ5D-3 and 0-10 pain
##############################################################################

ref_base<-summary(unlist(eq5d_3 %>% select(PAIN)))
pain_validation<-pain_scale %>% 
  mutate(PAIN_rescaled = rescale(ref_base,PAIN0_10)) %>%
  select(-TRIAL,-TIME_cat) %>%
  left_join(eq5d_3 %>% select(UNIQUE_ID,PAIN,TIME)) %>%
  filter(!is.na(PAIN) & !is.na(PAIN_rescaled))
length(unique(pain_validation$UNIQUE_ID))
summary(pain_validation$PAIN_rescaled-pain_validation$PAIN)
cor(pain_validation$PAIN_rescaled,pain_validation$PAIN,
    method="spearman",use="pairwise.complete.obs")
cor(pain_validation$PAIN0_10,pain_validation$PAIN,
    method="spearman",use="pairwise.complete.obs")

pain_validation2<-pain_scale %>% 
  mutate(PAIN_rescaled = rescale(ref_base,PAIN0_10)) %>%
  select(-TRIAL,-TIME_cat) %>%
  left_join(eq5d_5 %>% select(UNIQUE_ID,PAIN,TIME)) %>%
  mutate(PAIN_rescaled_5 = rescale(ref_base,PAIN)) %>%
  filter(!is.na(PAIN_rescaled_5) & !is.na(PAIN_rescaled))
length(unique(pain_validation2$UNIQUE_ID))
summary(pain_validation2$PAIN_rescaled-pain_validation2$PAIN_rescaled_5)
cor(pain_validation2$PAIN_rescaled,pain_validation2$PAIN_rescaled_5,
    method="spearman",use="pairwise.complete.obs")
cor(pain_validation2$PAIN,pain_validation2$PAIN0_10,
    method="spearman",use="pairwise.complete.obs")
rm(pain_validation,pain_validation2)

##############################################################################
###  Table 1: Baseline characteristics
##############################################################################

summary(demographics$AGE)
table(demographics$SEX)
round(table(demographics$SEX)*100/nrow(demographics),1)
summary(demographics$TIME_SINCE_STROKE)
table(demographics$STROKE_TYPE)
summary(demographics$BNIH)
round(sum(!is.na(demographics$BNIH))*100/nrow(demographics),1)
table(demographics$immobility)
round(table(demographics$immobility)*100/sum(!is.na(demographics$immobility)),1)
table(demographics$DIABETES)
round(table(demographics$DIABETES)*100/sum(!is.na(demographics$DIABETES)),1)
round(sum(!is.na(demographics$DIABETES))*100/nrow(demographics),1)
table(demographics$APHASIA_BASE)
round(table(demographics$APHASIA_BASE)*100/sum(!is.na(demographics$APHASIA_BASE)),1)
sum(!is.na(demographics$APHASIA_BASE))
round(sum(!is.na(demographics$APHASIA_BASE))*100/nrow(demographics),1)

table(demographics$TRIAL_c)
rm(demographics)

##############################################################################
###  Pain Assessments
##############################################################################

nrow(pain_all %>% filter(Measure=="Pain_scale") %>% distinct(UNIQUE_ID))

table(pain_all %>% filter(Measure=="EQ5D_3") %>% select(RESP_proxy))

tapply(pain_all$PAIN, pain_all$RESP_proxy, summary)
ggplot(pain_all %>% filter(Measure=="EQ5D_3")) + geom_boxplot(aes(x=RESP_proxy, y = PAIN))
wilcox.test(unlist(pain_all %>% filter(RESP_proxy=="PROXY" & Measure=="EQ5D_3") %>% select(PAIN)), 
            unlist(pain_all %>% filter(RESP_proxy=="SELF" & Measure=="EQ5D_3") %>% select(PAIN)), 
            alternative = "two.sided")



table(pain_all$RESPONDENT)
pain_all %>% filter(Measure=="EQ5D_3") %>% 
  group_by(RESP_proxy,TIME_cat) %>%
  summarise(x = paste0(median(PAIN)," (",
                       quantile(PAIN,prob=0.25),"-",
                       quantile(PAIN,prob=0.75),") [",
                       n(),"]")) %>%
  spread(TIME_cat,x)

##############################################################################
###  Figure 1: Smoothed estimates of standardised pain scores compared to each 
# pain assessment tool score (up to 2 years post-stroke)
##############################################################################

ggplot(bind_rows(pain_graph %>% select(TIME,PAIN,Measure) %>%
                   mutate(Measure = ifelse(Measure=="Pain_scale", "Numeric Pain Scale [0-10]",
                                           ifelse(Measure=="EQ5D_3L", "EQ5D_3L [1-3]",
                                                  "EQ5D_5L [1-5]"))),
                 pain_graph %>% select(TIME,PAIN_rescaled) %>%
                   rename(PAIN = PAIN_rescaled) %>% mutate(Measure="Standardised [1-3]")) %>%
         filter(TIME<=365.25*2)) +
  geom_smooth(aes(x=TIME,y=PAIN, linetype=Measure), color="black") +
  theme_bw() + xlab("Days since Stroke") +ylab("Pain Score") +
  labs(color = "Measure [scale]") + 
  coord_cartesian(ylim =c(0,5)) +
  theme(text = element_text(size=15),
        legend.key.width=unit(1,"cm")) 


##############################################################################
###  Table 2: Pain measurements by time period and scale
##############################################################################

# by timepoint and measure
for (time in times2) {
  for (measure in p_measures) {
    print(time)
    print(measure)
    print(nrow(pain_graph %>% filter(!is.na(PAIN) & 
                                       TIME_cat==time & Measure==measure) 
               %>% distinct(UNIQUE_ID)))
    print(round(nrow(pain_graph %>% 
                       filter(!is.na(PAIN) &  TIME_cat==time & Measure==measure)  %>% 
                       distinct(UNIQUE_ID))*100/length(unique(pain_graph$UNIQUE_ID)),1))
  }
}

# by timepoint, any measure
for (time in times2) {
  print(time)
  print(nrow(pain_graph %>% filter(!is.na(PAIN) & 
                                     TIME_cat==time) 
             %>% distinct(UNIQUE_ID)))
  print(round(nrow(pain_graph %>% 
                     filter(!is.na(PAIN) &  TIME_cat==time)  %>% 
                     distinct(UNIQUE_ID))*100/length(unique(pain_graph$UNIQUE_ID)),1))
}

# by measure, any timepoint
for (measure in p_measures) {
  print(measure)
  print(nrow(pain_graph %>% filter(!is.na(PAIN) & Measure==measure) 
             %>% distinct(UNIQUE_ID)))
  print(round(nrow(pain_graph %>% 
                     filter(!is.na(PAIN) &  Measure==measure)  %>% 
                     distinct(UNIQUE_ID))*100/length(unique(pain_graph$UNIQUE_ID)),1))
}


##############################################################################
###  Table 3A: proportion of people with some pain over time
##############################################################################

for (time in times2) {
  for (measure in p_measures) {
    print(time)
    print(measure)
    min<-ifelse(measure %in% c("Pain_Scale"),0,1)
    temp<-pain_graph %>% filter(!is.na(PAIN) & 
                                  TIME_cat==time & 
                                  Measure==measure)
    if(nrow(temp)>0){
      temp<-temp %>%
        group_by(UNIQUE_ID) %>% 
        summarise(PAIN_person = max(PAIN)>min)
      print(sum(temp$PAIN_person))
      print(length(unique(temp$UNIQUE_ID)))
      print(sum(temp$PAIN_person)*100/length(unique(temp$UNIQUE_ID)))
      print(binom.test(sum(temp$PAIN_person),length(unique(temp$UNIQUE_ID)))$conf.int[1:2])
    }
  }
}

for (time in times2) {
  print(time)
  temp<-pain_graph %>% filter(!is.na(PAIN_rescaled) & TIME_cat==time) %>%
    group_by(UNIQUE_ID) %>% summarise(PAIN_person = max(PAIN_rescaled)>1)
  print(sum(temp$PAIN_person))
  print(length(unique(temp$UNIQUE_ID)))
  print(sum(temp$PAIN_person)*100/length(unique(temp$UNIQUE_ID)))
  print(binom.test(sum(temp$PAIN_person),length(unique(temp$UNIQUE_ID)))$conf.int[1:2])
}

##############################################################################
###  Table 3B: proportion of people with EXTREME pain over time
##############################################################################

for (time in times2) {
  for (measure in p_measures) {
    print(time)
    print(measure)
    max<-ifelse(measure=="EQ5D_3L",3,
                ifelse(measure=="EQ5D_5L",5,
                       ifelse(measure=="SF36",5,
                              10)))
    temp<-pain_graph %>% filter(!is.na(PAIN) & 
                                  TIME_cat==time & 
                                  Measure==measure)
    if(nrow(temp)>0){
      temp<-temp %>%
        group_by(UNIQUE_ID) %>% 
        summarise(PAIN_person = max(PAIN)==max)
      print(sum(temp$PAIN_person))
      print(length(unique(temp$UNIQUE_ID)))
      print(round(sum(temp$PAIN_person)*100/length(unique(temp$UNIQUE_ID)),1))
      print(100*round(binom.test(sum(temp$PAIN_person),
                                 length(unique(temp$UNIQUE_ID)))$conf.int[1:2],3))
    }
  }
}

for (time in times2) {
  print(time)
  temp<-pain_graph %>% filter(!is.na(PAIN_rescaled) & TIME_cat==time) %>%
    group_by(UNIQUE_ID) %>% summarise(PAIN_person = max(PAIN_rescaled)==3)
  print(sum(temp$PAIN_person))
  print(length(unique(temp$UNIQUE_ID)))
  print(round(sum(temp$PAIN_person)*100/length(unique(temp$UNIQUE_ID)),1))
  print(100*round(binom.test(sum(temp$PAIN_person),
                             length(unique(temp$UNIQUE_ID)))$conf.int[1:2],3))
}


##############################################################################
###  Figure 2: Pain over time in trials with multiple measurements per person
##############################################################################

temp<-pain_graph %>% 
   filter(TIME<=365.25) %>%
   add_count(UNIQUE_ID, TRIAL_c) %>%
   filter(n>1) %>%
  mutate(temp = ifelse(TRIAL_c=="ACUTE", "Acute", "Non-Acute"), 
         Anon = paste0(Anon,": ", temp))
length(unique(temp$UNIQUE_ID))
 
ggplot(temp)  +
  geom_smooth(aes(x=TIME,y=PAIN_rescaled, linetype=Anon),color="black") +
  theme_bw() + xlab("Days since Stroke") +ylab("Standardised Pain Score") +
  labs(linetype = "Anonymised \nTrial") +
  theme(text = element_text(size=15),
        legend.key.width=unit(1,"cm")) 

rm(temp)

##############################################################################
###  Table 4: Unadjusted associations between the presence of extreme pain and 
###       participant characteristics
##############################################################################

table5<-pain_all %>% 
  group_by(UNIQUE_ID,TIME_cat)  %>%
  summarise(Pain_flag = 1*(sum(PAIN_rescaled==3)>0)) %>%
  left_join(independence) %>%
  left_join(mobility_all %>% group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(w_mobility   = min(Mobility_rescaled, na.rm=T),
                        Mobility_flag = 1*(sum(Mobility_rescaled!=3)>0))) %>%
  left_join(demog_clean) %>%
  mutate(APHASIA_BASE = ifelse(APHASIA_BASE=="",NA,APHASIA_BASE))

for(time in unique(pain_all$TIME_cat)) {
  print(time)
  temp<-table5 %>% filter(TIME_cat==time)
  print(nrow(temp))
  print(length(unique(temp$TRIAL)))
  print(tapply(temp$AGE,temp$Pain_flag,summary))
  print(wilcox.test(temp$AGE~temp$Pain_flag))
  print(table(temp$SEX,temp$Pain_flag))
  print(prop.table(table(temp$SEX,temp$Pain_flag),2)*100)
  print(chisq.test(temp$SEX,temp$Pain_flag))
  print(tapply(temp$BNIH,temp$Pain_flag,summary))
  print(wilcox.test(temp$BNIH~temp$Pain_flag))
  print(table(temp$APHASIA_BASE,temp$Pain_flag))
  print(prop.table(table(temp$APHASIA_BASE,temp$Pain_flag),2)*100)
  print(chisq.test(temp$APHASIA_BASE,temp$Pain_flag))
  if(time!="0-4 WEEKS") {   # no diabetes measurements for trials with data this early
    print(table(temp$DIABETES,temp$Pain_flag))
    print(prop.table(table(temp$DIABETES,temp$Pain_flag),2)*100)
    print(chisq.test(temp$DIABETES,temp$Pain_flag))
  }
  print(tapply(temp$w_independence,temp$Pain_flag,summary))
  print(wilcox.test(temp$w_independence~temp$Pain_flag))
  print(tapply(temp$w_mobility,temp$Pain_flag,summary))
  print(wilcox.test(temp$w_mobility~temp$Pain_flag))
}



##############################################################################
###  Table 5: Adjusted Logistic Regression: Presence of pain, stratified
###       by timepoints
##############################################################################

table6<- pain_all %>% 
  group_by(UNIQUE_ID,TIME_cat)  %>%
  summarise(Pain_flag = 1*(sum(PAIN_rescaled>1)>0)) %>%
  left_join(independence) %>%
  left_join(mobility_all %>% group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(Mobility_flag = 1*(sum(Mobility_rescaled!=3)>0))) %>%
  left_join(demog_clean) %>%
  ungroup() %>%
  select(UNIQUE_ID, TIME_cat, TRIAL_c2, TRIAL, Pain_flag,
         AGE, SEX, BNIH, DIABETES, Mobility_flag, APHASIA_BASE, 
         INDEPENDENCE_med)

nrow(table6 %>% filter(TIME_cat=="0-4 WEEKS" & TRIAL_c2!="CHRONIC"))
nrow(table6 %>% filter(TIME_cat=="0-4 WEEKS" & TRIAL_c2=="CHRONIC"))
nrow(distinct(table6 %>% filter(TIME_cat=="0-4 WEEKS") %>% select(TRIAL)))


temp<-table6 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" & TRIAL_c2=="ACUTE" & complete.cases(.))
nrow(distinct(temp %>% select(TRIAL)))
nrow(temp)
table(temp$SEX)
table(temp$Mobility_flag)
table(temp$DIABETES)
table(temp$APHASIA_BASE)
glm<-glm(Pain_flag ~ AGE + SEX + BNIH + DIABETES + Mobility_flag + APHASIA_BASE + INDEPENDENCE_med,
         data=temp, family = "binomial")
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)




round(prop.table(table(table6 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                           TRIAL_c2=="NON-ACUTE") %>% 
                         select(DIABETES), useNA = "ifany"))*100)
round(prop.table(table(table6 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                           TRIAL_c2=="NON-ACUTE") %>% 
                         select(INDEPENDENCE_med), useNA = "ifany"))*100)
temp<-table6 %>% select(-DIABETES, -INDEPENDENCE_med) %>%
  filter(TIME_cat=="5 WEEKS - 3 MONTHS" & 
           TRIAL_c2=="NON-ACUTE" & 
           APHASIA_BASE!="MISSING" &
           complete.cases(.))
nrow(distinct(temp %>% select(TRIAL)))
nrow(temp)
table(temp$SEX)
table(temp$Mobility_flag)
table(temp$APHASIA_BASE)
glm<-glm(Pain_flag ~ AGE + SEX + APHASIA_BASE + BNIH + Mobility_flag,
         data=temp,family = "binomial")
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)



round(prop.table(table(table6 %>% filter(TIME_cat=="4-6 MONTHS" &
                                           TRIAL_c2=="ACUTE") %>% 
                         select(DIABETES)))*100)
temp<-table6 %>% filter(TIME_cat=="4-6 MONTHS" & TRIAL_c2=="ACUTE" & complete.cases(.))
nrow(distinct(temp %>% select(TRIAL)))
nrow(temp)
table(temp$SEX)
table(temp$Mobility_flag)
table(temp$APHASIA_BASE)
glm<-glm(Pain_flag ~ AGE + SEX + BNIH + Mobility_flag + APHASIA_BASE + INDEPENDENCE_med,
         data= temp, family = "binomial")
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)





round(prop.table(table(table6 %>% filter(TIME_cat=="4-6 MONTHS" &
                                           TRIAL_c2=="NON-ACUTE") %>% 
                         select(DIABETES), useNA = "ifany"))*100)
round(prop.table(table(is.na(table6 %>% filter(TIME_cat=="4-6 MONTHS" &
                                                 TRIAL_c2=="NON-ACUTE") %>% 
                               select(INDEPENDENCE_med))))*100)

temp<-table6 %>% select(-DIABETES,-INDEPENDENCE_med) %>%
  filter(TIME_cat=="4-6 MONTHS" & TRIAL_c2=="NON-ACUTE" & complete.cases(.))
nrow(distinct(temp %>% select(TRIAL)))
nrow(temp)
table(temp$SEX)
table(temp$Mobility_flag)
table(temp$APHASIA_BASE)
glm<-glm(Pain_flag ~ AGE + SEX + APHASIA_BASE + Mobility_flag +BNIH,
         data=temp, family = "binomial")
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)


round(prop.table(table(table6 %>% filter(TIME_cat==">6 MONTHS" & 
                                           TRIAL_c2=="ACUTE") %>% 
                         select(DIABETES)=="MISSING"))*100)
temp<-table6 %>% select(-DIABETES) %>%
  filter(TIME_cat==">6 MONTHS" & TRIAL_c2=="ACUTE" & complete.cases(.))
nrow(distinct(temp %>% select(TRIAL)))
nrow(temp)
table(temp$SEX)
table(temp$Mobility_flag)
table(temp$APHASIA_BASE)
glm<-glm(Pain_flag ~ AGE + SEX + BNIH + Mobility_flag + APHASIA_BASE + INDEPENDENCE_med,
         data=temp, family = "binomial")
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)


round(prop.table(table(table6 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                           TRIAL_c2=="NON-ACUTE") %>% 
                         select(DIABETES)=="MISSING"))*100)
round(prop.table(table(is.na(table6 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                                 TRIAL_c2=="NON-ACUTE") %>% 
                               select(INDEPENDENCE_med))))*100)
temp<-table6 %>% select(-DIABETES,-INDEPENDENCE_med) %>%
  filter(TIME_cat==">6 MONTHS" & TRIAL_c2=="NON-ACUTE" & complete.cases(.))
nrow(distinct(temp %>% select(TRIAL)))
nrow(temp)
table(temp$SEX)
table(temp$Mobility_flag)
table(temp$APHASIA_BASE)
glm<-glm(Pain_flag ~ AGE + SEX + APHASIA_BASE + BNIH+ Mobility_flag,
         data=temp, family = "binomial")
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)


##############################################################################
###  Table 6: Adjusted Logistic Regression: Presence of extreme pain, stratified
###       by timepoints
##############################################################################

table7<- pain_all %>% 
  group_by(UNIQUE_ID,TIME_cat)  %>%
  summarise(Pain_flag = 1*(sum(PAIN_rescaled==3)>0)) %>%
  left_join(independence) %>%
  left_join(mobility_all %>% group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(Mobility_flag = 1*(sum(Mobility_rescaled!=3)>0))) %>%
  left_join(demog_clean) %>%
  ungroup() 
  
glm<-glm(Pain_flag ~ AGE + SEX + BNIH + DIABETES + Mobility_flag + APHASIA_BASE + INDEPENDENCE_med,
         data=table7 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" & 
                                    TRIAL_c2=="ACUTE"), 
         family = "binomial")
length(glm$residuals)
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)

round(prop.table(table(table7 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                             TRIAL_c2=="NON-ACUTE") %>% 
                         select(DIABETES), useNA = "ifany"))*100)
round(prop.table(table(table7 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                             TRIAL_c2=="NON-ACUTE") %>% 
                         select(INDEPENDENCE_med), useNA = "ifany"))*100)
glm<-glm(Pain_flag ~ AGE + SEX + APHASIA_BASE + BNIH + Mobility_flag,
         data=table7 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                    TRIAL_c2=="NON-ACUTE" &
                                    APHASIA_BASE!="MISSING"),
         family = "binomial")
length(glm$residuals)
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)

round(prop.table(table(table7 %>% filter(TIME_cat=="4-6 MONTHS" & 
                                             TRIAL_c2=="ACUTE") %>% 
                         select(DIABETES), useNA = "ifany"))*100)
glm<-glm(Pain_flag ~ AGE + SEX + BNIH + APHASIA_BASE + INDEPENDENCE_med + Mobility_flag,
         data=table7 %>% filter(TIME_cat=="4-6 MONTHS" & 
                                    TRIAL_c2=="ACUTE"), 
         family = "binomial")
length(glm$residuals)
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)

round(prop.table(table(table7 %>% filter(TIME_cat=="4-6 MONTHS" &
                                             TRIAL_c2=="NON-ACUTE") %>% 
                         select(DIABETES)))*100)
round(prop.table(table(is.na(table7 %>% filter(TIME_cat=="4-6 MONTHS" &
                                                   TRIAL_c2=="NON-ACUTE") %>% 
                               select(INDEPENDENCE_med))))*100)
glm<-glm(Pain_flag ~ AGE + SEX + APHASIA_BASE + Mobility_flag +BNIH,
         data=table7 %>% filter(TIME_cat=="4-6 MONTHS" & 
                                    TRIAL_c2=="NON-ACUTE"), 
         family = "binomial")
length(glm$residuals)
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)


round(prop.table(table(is.na(table7 %>% filter(TIME_cat==">6 MONTHS" & 
                                                   TRIAL_c2=="ACUTE") %>% 
                               select(DIABETES))))*100)
glm<-glm(Pain_flag ~ AGE + SEX + BNIH + Mobility_flag + APHASIA_BASE + INDEPENDENCE_med,
         data=table7 %>% filter(TIME_cat==">6 MONTHS" & 
                                    TRIAL_c2=="ACUTE"), 
         family = "binomial")
length(glm$residuals)
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)


round(prop.table(table(table7 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                             TRIAL_c2=="NON-ACUTE") %>% 
                         select(DIABETES)=="MISSING"))*100)
round(prop.table(table(is.na(table7 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                                   TRIAL_c2=="NON-ACUTE") %>% 
                               select(INDEPENDENCE_med))))*100)
glm<-glm(Pain_flag ~ AGE + SEX + APHASIA_BASE + BNIH+ Mobility_flag,
         data=table7 %>% filter(TIME_cat==">6 MONTHS" & 
                                    TRIAL_c2=="NON-ACUTE"), 
         family = "binomial")
length(glm$residuals)
round(exp(summary(glm)$coefficients[,1]),3)
round(exp(confint.default(glm)),3)
round(summary(glm)$coefficients[,4],3)




##############################################################################
###  Table 7: Adjusted Linear Regression: Amount of pain, stratified by 
###     timepoints
############################################################################## 

table8<-pain_all %>% 
  group_by(UNIQUE_ID,TIME_cat) %>%
  summarise(pain_med = median(PAIN)) %>% 
  left_join(independence) %>%
  left_join(mobility_all %>% group_by(UNIQUE_ID,TIME_cat) %>%
              summarise(Mobility_flag = 1*(sum(Mobility_rescaled!=3)>0))) %>%
  left_join(demog_clean) %>%
  ungroup() 

lm<-lm(pain_med ~ AGE + SEX + BNIH + DIABETES + Mobility_flag + APHASIA_BASE + INDEPENDENCE_med,
       data=table8 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" & 
                                TRIAL_c2=="ACUTE"))
length(lm$residuals)
round(exp(summary(lm)$coefficients[,1]),3)
round(exp(confint.default(lm)),3)
round(summary(lm)$coefficients[,4],3)

round(prop.table(table(table8 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                           TRIAL_c2=="NON-ACUTE") %>% 
                         select(DIABETES)=="MISSING", useNA = "ifany"))*100)
round(prop.table(table(is.na(table8 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" &
                                                 TRIAL_c2=="NON-ACUTE") %>% 
                               select(INDEPENDENCE_med)), useNA = "ifany"))*100)
lm<-lm(pain_med ~ AGE + SEX + APHASIA_BASE +BNIH + Mobility_flag,
       data=table8 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS" & 
                                TRIAL_c2=="NON-ACUTE" & 
                                APHASIA_BASE!="MISSING"))
length(lm$residuals)
round(exp(summary(lm)$coefficients[,1]),3)
round(exp(confint.default(lm)),3)
round(summary(lm)$coefficients[,4],3)

lm<-lm(pain_med ~ AGE + SEX + BNIH + Mobility_flag +APHASIA_BASE+ INDEPENDENCE_med,
       data=table8 %>% filter(TIME_cat=="4-6 MONTHS" & 
                                TRIAL_c2=="ACUTE"))
length(lm$residuals)
round(exp(summary(lm)$coefficients[,1]),3)
round(exp(confint.default(lm)),3)
round(summary(lm)$coefficients[,4],3)

round(prop.table(table(table8 %>% filter(TIME_cat=="4-6 MONTHS" &
                                           TRIAL_c2=="NON-ACUTE") %>% 
                         select(DIABETES)=="MISSING", useNA = "ifany"))*100)
round(prop.table(table(is.na(table8 %>% filter(TIME_cat=="4-6 MONTHS" &
                                                 TRIAL_c2=="NON-ACUTE") %>% 
                               select(INDEPENDENCE_med)), useNA = "ifany"))*100)
lm<-lm(pain_med ~ AGE + SEX + BNIH + Mobility_flag +APHASIA_BASE,
       data=table8 %>% filter(TIME_cat=="4-6 MONTHS" & 
                                TRIAL_c2=="NON-ACUTE"))
length(lm$residuals)
round(exp(summary(lm)$coefficients[,1]),3)
round(exp(confint.default(lm)),3)
round(summary(lm)$coefficients[,4],3)

lm<-lm(pain_med ~ AGE + SEX + BNIH + Mobility_flag +APHASIA_BASE + INDEPENDENCE_med,
       data=table8 %>% filter(TIME_cat==">6 MONTHS" & 
                                TRIAL_c2=="ACUTE"))
length(lm$residuals)
round(exp(summary(lm)$coefficients[,1]),3)
round(exp(confint.default(lm)),3)
round(summary(lm)$coefficients[,4],3)

lm<-lm(pain_med ~ AGE + SEX + BNIH + APHASIA_BASE + Mobility_flag,
       data=table8 %>% filter(TIME_cat==">6 MONTHS" & 
                                TRIAL_c2=="NON-ACUTE"))
length(lm$residuals)
round(exp(summary(lm)$coefficients[,1]),3)
round(exp(confint.default(lm)),3)
round(summary(lm)$coefficients[,4],3)

rm(table8)


##############################################################################
###  Table 8: Linear Regression of amount of pain (on EQ5D 3-level and 5-level) 
# adjusted for EQ-5D domains
############################################################################## 

table9<-eq5d_3 %>%
  group_by(UNIQUE_ID,TIME_cat) %>%
  summarise(PAIN = median(PAIN),
            MOBILITY3 = median(MOBILITY3),
            SELF_CARE = median(SELF_CARE),
            USUAL_ACTIVITIES = median(USUAL_ACTIVITIES),
            VAS = median(VAS),
            ANXIETY = median(ANXIETY)) %>%
  inner_join(demog_clean) %>%
  filter(APHASIA_BASE!="MISSING") %>%
  ungroup 

# EQ5D_3 5 weeks to 3 months
table9 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS") %>% 
  select(BNIH, VAS, ANXIETY, APHASIA_BASE,
         USUAL_ACTIVITIES,SELF_CARE,MOBILITY3) %>% 
  miss_var_summary()
lm<-lm(PAIN ~ MOBILITY3 +  BNIH + AGE + SEX + APHASIA_BASE + 
         +VAS + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, 
       data=table9 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS"))
length(lm$residuals)
nrow(distinct(table9 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS") %>% select(TRIAL)))
round(exp(summary(lm)$coefficients[,1]),3)
round(summary(lm)$coefficients[,4],3)
paste0("(",round(exp(confint.default(lm)),3)[,1],", ",
       round(exp(confint.default(lm)),3)[,2],")")

# EQ5D_3 4 -6 months
table9 %>% filter(TIME_cat=="4-6 MONTHS") %>% 
  select(BNIH, VAS, ANXIETY, APHASIA_BASE,
         USUAL_ACTIVITIES,SELF_CARE,MOBILITY3) %>% 
  miss_var_summary()
lm<-lm(PAIN ~ MOBILITY3 + AGE + SEX + APHASIA_BASE + 
         + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, 
       data=table9 %>% filter(TIME_cat=="4-6 MONTHS"))
length(lm$residuals)
nrow(distinct(table9 %>% filter(TIME_cat=="4-6 MONTHS") %>% select(TRIAL)))
round(exp(summary(lm)$coefficients[,1]),3)
round(summary(lm)$coefficients[,4],3)
paste0("(",round(exp(confint.default(lm)),3)[,1],", ",
       round(exp(confint.default(lm)),3)[,2],")")

table9 %>% filter(TIME_cat==">6 MONTHS") %>% 
  select(BNIH, VAS, ANXIETY, APHASIA_BASE,
         USUAL_ACTIVITIES,SELF_CARE,MOBILITY3) %>% 
  miss_var_summary()
lm<-lm(PAIN ~ MOBILITY3 + AGE + SEX + APHASIA_BASE + 
          + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, 
       data=table9 %>% filter(TIME_cat==">6 MONTHS"))
length(lm$residuals)
nrow(distinct(table9 %>% filter(TIME_cat==">6 MONTHS") %>% select(TRIAL)))
round(exp(summary(lm)$coefficients[,1]),3)
round(summary(lm)$coefficients[,4],3)
paste0("(",round(exp(confint.default(lm)),3)[,1],", ",
       round(exp(confint.default(lm)),3)[,2],")")


table9_2<-eq5d_5 %>%
  group_by(UNIQUE_ID,TIME_cat) %>%
  summarise(PAIN = median(PAIN),
            MOBILITY5 = median(MOBILITY5),
            SELF_CARE = median(SELF_CARE),
            USUAL_ACTIVITIES = median(USUAL_ACTIVITIES),
            VAS = median(VAS),
            ANXIETY = median(ANXIETY)) %>%
  inner_join(demog_clean) %>%
  filter(APHASIA_BASE!="MISSING") %>%
  ungroup 

lm<-lm(PAIN ~ MOBILITY5 +  BNIH + AGE + SEX + APHASIA_BASE + 
         +VAS + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, 
       data=table9_2 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS"))
length(lm$residuals)
nrow(distinct(table9_2 %>% filter(TIME_cat=="5 WEEKS - 3 MONTHS") %>% select(TRIAL)))

table9_2 %>% filter(TIME_cat=="4-6 MONTHS") %>% 
  select(BNIH, VAS, ANXIETY, APHASIA_BASE,
         USUAL_ACTIVITIES,SELF_CARE,MOBILITY5) %>% 
  miss_var_summary()
lm<-lm(PAIN ~ MOBILITY5 + BNIH + AGE + SEX + APHASIA_BASE + 
         + VAS  + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, 
       data=table9_2 %>% filter(TIME_cat=="4-6 MONTHS"))
length(lm$residuals)
nrow(distinct(table9_2 %>% filter(TIME_cat=="4-6 MONTHS") %>% select(TRIAL)))
round(exp(summary(lm)$coefficients[,1]),3)
round(summary(lm)$coefficients[,4],3)
paste0("(",round(exp(confint.default(lm)),3)[,1],", ",
       round(exp(confint.default(lm)),3)[,2],")")

table9_2 %>% filter(TIME_cat==">6 MONTHS") %>% 
  select(BNIH, VAS, ANXIETY, APHASIA_BASE,
         USUAL_ACTIVITIES,SELF_CARE,MOBILITY5) %>% 
  miss_var_summary()
lm<-lm(PAIN ~ MOBILITY5 + BNIH+ AGE + SEX + APHASIA_BASE + 
         +VAS + ANXIETY + USUAL_ACTIVITIES + SELF_CARE, 
       data=table9_2 %>% filter(TIME_cat==">6 MONTHS"))
length(lm$residuals)
nrow(distinct(table9_2 %>% filter(TIME_cat==">6 MONTHS") %>% select(TRIAL)))
round(exp(summary(lm)$coefficients[,1]),3)
round(summary(lm)$coefficients[,4],3)
paste0("(",round(exp(confint.default(lm)),3)[,1],", ",
       round(exp(confint.default(lm)),3)[,2],")")

##############################################################################
###  Figure 3: Relationship between pain and fatigue (SF36)
############################################################################## 

sum(!is.na(sf36$sf_pain) & !is.na(sf36$fatigue))

ggplot(sf36) +
  geom_boxplot(aes(x=as.factor(sf_pain),y=fatigue)) +
  xlab("SF36 Pain") +
  ylab("SF36 Fatigue") +
  theme_bw()

cor.test(sf36$sf_pain,sf36$fatigue,method="spearman",use="pairwise.complete.obs")



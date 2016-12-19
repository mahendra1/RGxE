################################################################################
#Program Name : RGxE                                                                                     
#Version: 1.1
#Purpose: It computes various stability analysis statistics for Genotype x 
#         Environment Interaction data
#         Overall ANOVA with different combination of random and fixed effects; 
#         Regression coefficients - slope and Deviation from regression; 
#         Shukla Sigma; Shukla SSquare; Wricke's Ecovalence; Kang's YS; 
#         Location statistics: Genotype F ratio across location and environment;
#         correlation of individual location with average location; 
#         Descriptive statistics on genotype, location, environment, year and 
#         replication;  
#         Cluster analysis of location
#libraries and packages used: dplyr, tidyr, sqldf, broom, agricolae 
#                             (stability.par), cluster, lme4, afex 
#Developers: Mahendra Dia,Todd C. Wehner and Consuelo Arellano
#            NC State University, NC, USA
#Last Modified: Feburary 24, 2016
################################################################################

#Install packages

install.packages("tidyr")
install.packages("dplyr")
install.packages ("sqldf")
install.packages("lme4")
install.packages("afex")
install.packages("broom")
install.packages("agricolae")


#call libraries

library (tidyr)
library (dplyr)
library (sqldf)
library(lme4)
library(afex)
library(broom)
library(agricolae)
library(cluster)
library(grDevices)

#list current objects

ls()

#to remove list of all objects from workspace

rm(list=ls())

#to identify current directory

getwd()

#set current directory
###user NEED to replace input data file path###

setwd("E:/PhD Research Work/PhD Articles/Articles for Publication/BLUE_BLUP_Prediction_Interval_SAS_R")

#check current directory

getwd()

########################################################################
####################        For Windows user         ###################
########################################################################

######################################
##########  Option 1  ################
######################################

#import input data
###user NEED to replace input data file name###

#tempa <-read.csv("RGxEInputData2_2016_02_15.csv", header = TRUE) 

######################################
##########  Option 2  ################
######################################

#opens pop up window to select file
#file.name <- file.choose()

#read csv file
#tempa <-read.csv(file.name)

########################################################################
####################          For Mac user           ###################
########################################################################

file.name <- file.choose()
file.name

#copy path from console and paste here
file.name <- "/Users/toddwehner/Desktop/RCodeCurrent/GxeR14Data4.csv"

#create object path for where the data will go at the end
# I append the file name with "_out" so I know its data coming out of R
out.file <- "/Users/toddwehner/Desktop/RCodeCurrent/GxeR14Out4.csv"
tempa <-read.csv(file.name)

########################################################################
######           Define trait that need to be analyzed           #######
########################################################################

#user need to define trait name for analysis
#example: MKMGHA = Marketable mega-gram per hectare

tempa <- tempa %>% rename(Trait = MKMGHA) %>%
          select(YR, LC, RP, CLT, Trait)

#view top 6 rows of input data

head(tempa)

#view bottom 6 rows of input data

tail(tempa)

#Get dimension of data

dim(tempa)

#Get data type

class(tempa)

#Get structure of data (character vs. numeric vs. matrix vs. vector vs. factor)
#make sure: numeric value should be numeric

str(tempa)

#validate the variable types

tempa$YR <- as.factor(tempa$YR)
tempa$RP <- as.factor(tempa$RP)
tempa$LC <- as.factor(tempa$LC)
tempa$CLT <- as.factor(tempa$CLT)
tempa$Trait <- as.numeric(tempa$Trait)


########################################################################
##                ANOVA: Compute analysis of variance                 ##
########################################################################

#Generate unique id for replication for anova

tempa$RPid<-as.factor(paste(tempa$YR, tempa$LC, tempa$RP, sep="."))

########################################################################
###         ANOVA Case 1: CLT, YR, LC and RP - All Random            ###
########################################################################

#full model

fit.f1<-lmer(Trait~ 1 + (1|YR)  + (1|LC) + (1|CLT) + (1|YR:LC) + (1|YR:CLT) + 
               (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)
#model summary

summary1 <- summary(fit.f1)

#variance of random factors

variance <- as.data.frame(summary1$varcor)

#drop rownames

rownames(variance) <- NULL
variance1 <- variance %>% select (-var1, -var2) %>% 
  rename(sov=grp, Variance=vcov, stddev=sdcor)

#Type 3  test of hypothesis
#Type III Wald chisquare tests

anova(fit.f1, type="III")

#Type 1  test of hypothesis

anova(fit.f1, type="marginal", test="F")

#model fitness

anovacase1 <- plot(fit.f1, 
                   main="Model fitness Case 1: CLT, YR, LC and RP - All Random",
                   xlab="Predicated Value", ylab="Residual")

#LRT - likelihood ratio test for computing significance of random effect 
#null model for YR

fit.f1y<-lmer(Trait~ 1 + (1|LC) + (1|CLT) + (1|YR:LC) + (1|YR:CLT) + 
                (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance

anova1y <-anova(fit.f1,fit.f1y)

#convert anova into data frame

anova1y <- data.frame(anova1y)

#convert rownames into column

anova1y$name <- rownames(anova1y)

#drop rownames

rownames(anova1y) <- NULL
anova1y <- anova1y %>% filter(name=="fit.f1") %>%
  mutate(sov="YR") %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))

#null model for LC

fit.f1l<-lmer(Trait~ 1 + (1|YR)  + (1|CLT) + (1|YR:LC) + (1|YR:CLT) + 
                (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance

anova1l <-anova(fit.f1,fit.f1l)

#convert anova into data frame

anova1l <- data.frame(anova1l)

#convert rownames into column

anova1l$name <- rownames(anova1l)

#drop rownames

rownames(anova1l) <- NULL
anova1l <- anova1l %>% filter(name=="fit.f1") %>%
  mutate(sov="LC") %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))

#null model for CLT

fit.f1c<-lmer(Trait~ 1 + (1|YR)  + (1|LC) + (1|YR:LC) + (1|YR:CLT) + 
                (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance

anova1c <-anova(fit.f1,fit.f1c)

#convert anova into data frame

anova1c <- data.frame(anova1c)

#convert rownames into column

anova1c$name <- rownames(anova1c)

#drop rownames

rownames(anova1c) <- NULL
anova1c <- anova1c %>% filter(name=="fit.f1") %>%
  mutate(sov="CLT") %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))

#null model for YR:LC
fit.f1yl<-lmer(Trait~ 1 + (1|YR)  + (1|LC) + (1|CLT) + (1|YR:CLT) + 
                 (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance

anova1yl <-anova(fit.f1,fit.f1yl)

#convert anova into data frame

anova1yl <- data.frame(anova1yl)

#convert rownames into column

anova1yl$name <- rownames(anova1yl)

#drop rownames

rownames(anova1yl) <- NULL
anova1yl <- anova1yl %>% filter(name=="fit.f1") %>%
  mutate(sov="YR:LC") %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))

#null model for YR:CLT

fit.f1yc<-lmer(Trait~ 1 + (1|YR)  + (1|LC) + (1|CLT) + (1|YR:LC) +  
                 (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance

anova1yc <-anova(fit.f1,fit.f1yc)

#convert anova into data frame

anova1yc <- data.frame(anova1yc)

#convert rownames into column

anova1yc$name <- rownames(anova1yc)

#drop rownames

rownames(anova1yc) <- NULL
anova1yc <- anova1yc %>% filter(name=="fit.f1") %>%
  mutate(sov="YR:CLT") %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))

#null model for LC:CLT

fit.f1lc<-lmer(Trait~ 1 + (1|YR)  + (1|LC) + (1|CLT) + (1|YR:LC) + 
                 (1|YR:CLT) + 
                 (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance

anova1lc <-anova(fit.f1,fit.f1lc)

#convert anova into data frame

anova1lc <- data.frame(anova1lc)

#convert rownames into column

anova1lc$name <- rownames(anova1lc)

#drop rownames

rownames(anova1lc) <- NULL
anova1lc <- anova1lc %>% filter(name=="fit.f1") %>%
  mutate(sov="LC:CLT") %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))

#null model for YR:LC:CLT

fit.f1ylc<-lmer(Trait~ 1 + (1|YR)  + (1|LC) + (1|CLT) + (1|YR:LC) + 
                  (1|YR:CLT) + 
                  (1|LC:CLT) +   (1|RPid), data=tempa)

#level of significance

anova1ylc <-anova(fit.f1,fit.f1ylc)

#convert anova into data frame

anova1ylc <- data.frame(anova1ylc)

#convert rownames into column

anova1ylc$name <- rownames(anova1ylc)

#drop rownames

rownames(anova1ylc) <- NULL
anova1ylc <- anova1ylc %>% filter(name=="fit.f1") %>%
  mutate(sov="YR:LC:CLT") %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))

#null model for RP

fit.f1r<-lmer(Trait~ 1 + (1|YR)  + (1|LC) + (1|CLT) + (1|YR:LC) + 
                (1|YR:CLT) + 
                (1|LC:CLT) + (1|YR:LC:CLT) , data=tempa)

#level of significance

anova1r <-anova(fit.f1,fit.f1r)

#convert anova into data frame

anova1r <- data.frame(anova1r)

#convert rownames into column

anova1r$name <- rownames(anova1r)

#drop rownames

rownames(anova1r) <- NULL
anova1r <- anova1r %>% filter(name=="fit.f1") %>%
  mutate(sov="RPid") %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))

#Merge anova and level of significance

anova1 <- bind_rows(anova1y, anova1l)%>% bind_rows(anova1c)%>%
  bind_rows(anova1yl)%>% bind_rows(anova1yc)%>%
  bind_rows(anova1r)%>% bind_rows(anova1lc)%>%bind_rows(anova1ylc)
anova1 <- as.data.frame(anova1)

#Merge final output

anova_randall <- variance1%>% left_join(anova1 , by ="sov")
anova_randall$Pr_Chisq[anova_randall$stddev == 0] <- NA


#Compute BLUP for CLT, CLT across location with prediction interval
#BLUP - Best linear unbiased predictor

randeffect1 <- ranef(fit.f1)

#BLUP for clt
BLUP_CLT <- as.data.frame(randeffect1$CLT)

#convert rownames into column

BLUP_CLT$genotype <- rownames(BLUP_CLT)

#drop rownames

rownames(BLUP_CLT) <- NULL

#rename variable name

BLUP_CLT <- BLUP_CLT %>% select(genotype,Blup = starts_with("(Intercept)"))

#select top 1 row

BLUP_CLT1 <- BLUP_CLT %>% filter(row_number()==1)

#select row # >1

BLUP_CLT2 <- BLUP_CLT %>% filter(row_number()>1)

#compute BLUP value

BLUP_CLT3 <- BLUP_CLT2 %>% mutate(Blup = Blup + BLUP_CLT1$Blup)

#final output for BLUP for cultivars

BLUP_CLT4 <- BLUP_CLT1 %>% bind_rows(BLUP_CLT3)
BLUP_CLT4 <- as.data.frame(BLUP_CLT4)


########################################################################
##         ANOVA Case 2: CLT, YR and LC - Fixed; RP - Random          ##
########################################################################

#full model

fit.f2<-lmer(Trait~ YR*LC*CLT + (1|RPid), data=tempa)

#model summary

summary2 <- summary(fit.f2)

##variance of random factors

variance2 <- as.data.frame(summary2$varcor)

#drop rownames

rownames(variance2) <- NULL
variance2 <- variance2 %>% select (-var1, -var2) %>% 
  rename(sov=grp, Mean_Sq=vcov, stddev=sdcor)

#Type 3  test of hypothesis

anova.f2t3 <- as.data.frame(anova(fit.f2, type="III"))

#convert rownames into column

anova.f2t3$name <- rownames(anova.f2t3)

#drop rownames

rownames(anova.f2t3) <- NULL

#Type 3  test of hypothesis

anova.f2t2 <- anova(fit.f2, type="marginal", test="F")

#model fitness

anovacase2 <- plot(fit.f2, 
                   main="Model fitness Case 2: CLT, YR and LC - Fixed; RP - Random",
                   xlab="Predicated Value", ylab="Residual")

#level of significance: "KR" is implemented corresponding to  
#the Kenward-Rogers approximation for degrees of freedom

anova_f2t3 <- mixed(Trait~ YR*LC*CLT + (1|RPid), data=tempa)

anova_f2t3 <- as.data.frame(anova_f2t3$anova_table)

#convert rownames into column

anova_f2t3$name <- rownames(anova_f2t3)

#drop rownames

rownames(anova_f2t3) <- NULL

anova_f2t3 <- anova_f2t3 %>% select(name, Prob_F = starts_with("Pr(>F)"))

#final output for case 2: CLT, YR and LC - Fixed

anova_f2 <- anova.f2t3 %>% left_join(anova_f2t3, by = "name")%>%
  rename(sov = name)%>%
  select(sov, Df,Sum_Sq=starts_with("Sum Sq"), 
         Mean_Sq = starts_with("Mean Sq"),
         F_value=starts_with("F value"), Prob_F)

anova_f2 <- anova_f2 %>% bind_rows(variance2)
anova_f2a <- as.data.frame(anova_f2) # print for output

########################################################################
##         ANOVA CASE 3: CLT - Fixed; YR, LC and RP - Random          ##
########################################################################

fit.f3<-lmer(Trait~ 1 + (1|YR)  + (1|LC) + CLT + (1|YR:LC) + 
               (1|YR:CLT) + 
               (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#model summary

summary3 <- summary(fit.f3)

##variance of random factors

variance3 <- as.data.frame(summary3$varcor)

#drop rownames

rownames(variance3) <- NULL
variance3 <- variance3 %>% select (-var1, -var2) %>% 
  rename(sov=grp, Mean_Sq=vcov, stddev=sdcor)

#Type 3 test of hypothesis

anova.f3t3 <- as.data.frame(anova(fit.f3, type="III"))

#convert rownames into column

anova.f3t3$name <- rownames(anova.f3t3)

#drop rownames

rownames(anova.f3t3) <- NULL

#Type 3  test of hypothesis

anova.f3t2 <- anova(fit.f3, type="marginal", test="F")

#model fitness

anovacase3 <- plot(fit.f3, 
                   main="Model fitness Case 3: CLT - Fixed; YR, LC and RP - Random",
                   xlab="Predicated Value", ylab="Residual")

#p value for fixed effects

anova_f3t3 <- mixed(Trait~ 1 + (1|YR)  + (1|LC) + CLT + (1|YR:LC) + 
                      (1|YR:CLT) + (1|LC:CLT) + (1|YR:LC:CLT) +  
                      (1|RPid), data=tempa)

anova_f3t3 <- as.data.frame(anova_f3t3$anova_table)

#convert rownames into column

anova_f3t3$name <- rownames(anova_f3t3)

#drop rownames

rownames(anova_f3t3) <- NULL

anova_f3t3 <- anova_f3t3 %>% select(name, Prob_F = starts_with("Pr(>F)"))

#final output of case 3 for fixed effect - CLT 

anova_f3 <- anova.f3t3 %>% left_join(anova_f3t3, by = "name")%>%
  rename(sov = name)%>%
  select(sov, Df,Sum_Sq=starts_with("Sum Sq"), 
         Mean_Sq = starts_with("Mean Sq"),
         F_value=starts_with("F value"), Prob_F)

#create function for Likelihood ratio test
#a=outputdatasetname; example-anova1r
#b=full model name; example-fit.f1
#c=reduced model name; example-fit.f1r
#d=effect name; example- "RPid", NOTE: call it in quotation

anova_lrt <- function (a,b,c,d){
  #level of significance
  a <-anova(b,c)
  #convert anova into data frame
  a <- data.frame(a)
  #convert rownames into column
  a$name <- rownames(a)
  # drop rownames
  rownames(a) <- NULL
  a <- a %>% filter(name=="b") %>%
    mutate(sov=d) %>% select(sov, Pr_Chisq = starts_with("Pr..Chisq."))
  # return the result
  return(a)
}

#level of significance for random effects
#null model for CLT

fit.f3c <- lmer(Trait~ 1 + (1|YR) + (1|LC) + (1|YR:LC) + (1|YR:LC) +
                  (1|YR:CLT) + (1|LC:CLT) + (1|YR:LC:CLT) +  
                  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova3c <- anova_lrt(anova3c,fit.f3,fit.f3c,"CLT")

#null model for LC

fit.f3l <- lmer(Trait~ 1 + CLT + (1|YR) + (1|YR:LC) + (1|YR:CLT) + 
                  (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt
anova3l <- anova_lrt(anova3l,fit.f3,fit.f3l,"LC")

#null model for YR

fit.f3y <- lmer(Trait~ 1 + CLT + (1|LC) + (1|YR:LC) + (1|YR:CLT) + 
                  (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova3y <- anova_lrt(anova3y,fit.f3,fit.f3y,"YR")

#null model for YR:LC

fit.f3yl <- lmer(Trait~ 1 + CLT + (1|LC) + (1|YR) + (1|YR:CLT) + 
                   (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova3yl <- anova_lrt(anova3yl,fit.f3,fit.f3yl,"YR:LC")

#null model for YR:CLT

fit.f3yc <- lmer(Trait~ 1 + CLT + (1|LC) + (1|YR) + (1|YR:LC) +  
                   (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova3yc <- anova_lrt(anova3yc,fit.f3,fit.f3yc,"YR:CLT")

#null model for LC:CLT

fit.f3lc <- lmer(Trait~ 1 + CLT + (1|LC) + (1|YR) + (1|YR:LC) + 
                   (1|YR:CLT) + 
                   (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova3lc <- anova_lrt(anova3lc,fit.f3,fit.f3lc,"LC:CLT")

#null model for YR:LC:CLT

fit.f3ylc <- lmer(Trait~ 1 + CLT + (1|LC) + (1|YR) + (1|YR:LC) + 
                    (1|YR:CLT) + 
                    (1|LC:CLT) + (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova3ylc <- anova_lrt(anova3ylc,fit.f3,fit.f3ylc,"YR:LC:CLT")

#null model for RPid or replication

fit.f3r <- lmer(Trait~ 1 + CLT + (1|LC) + (1|YR) + (1|YR:LC) + 
                  (1|YR:CLT) + 
                  (1|LC:CLT) + (1|YR:LC:CLT), data=tempa)

#level of significance
#call function anova_lrt

anova3r <- anova_lrt(anova3r,fit.f3,fit.f3r,"RPid")

#Merge anova and level of significance

anova3 <- bind_rows(anova3y, anova3l)%>% bind_rows(anova3yl)%>% 
  bind_rows(anova3yc)%>% bind_rows(anova3r)%>% 
  bind_rows(anova3lc)%>%bind_rows(anova3ylc)
anova3 <- as.data.frame(anova3)

#Merge final output

anova_cfix <- variance3%>% left_join(anova3 , by ="sov")
anova_cfix$Pr_Chisq[anova_cfix$stddev == 0] <- NA

########################################################################
##         ANOVA Case 4: LC - Fixed; YR, CLT and RP - Random          ##
########################################################################

fit.f4 <- lmer(Trait~ 1 +  LC + (1|YR) + (1|CLT) + (1|YR:LC) + 
                 (1|YR:CLT) + 
                 (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#model summary

summary4 <- summary(fit.f4)

##variance of random factors

variance4 <- as.data.frame(summary4$varcor)

# drop rownames

rownames(variance4) <- NULL
variance4 <- variance4 %>% select (-var1, -var2) %>% 
  rename(sov=grp, Mean_Sq=vcov, stddev=sdcor)

#Type 3 test of hypothesis

anova.f4t3 <- as.data.frame(anova(fit.f4, type="III"))

#convert rownames into column

anova.f4t3$name <- rownames(anova.f4t3)

#drop rownames

rownames(anova.f4t3) <- NULL

#Type 3 test of hypothesis

anova.f4t2 <- anova(fit.f3, type="marginal", test="F")

#model fitness

anovacase4 <- plot(fit.f4, 
                   main="Model fitness Case 4: LC - Fixed; YR, CLT and RP - Random",
                   xlab="Predicated Value", ylab="Residual")

#level of significance for fixed effects 

anova_f4t3 <- mixed(Trait~ 1 + LC + (1|YR)  + (1|CLT) + (1|YR:LC) + 
                      (1|YR:CLT) + 
                      (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

anova_f4t3 <- as.data.frame(anova_f4t3$anova_table)

#convert rownames into column

anova_f4t3$name <- rownames(anova_f4t3)

#drop rownames

rownames(anova_f4t3) <- NULL

anova_f4t3 <- anova_f4t3 %>% select(name, Prob_F = starts_with("Pr(>F)"))

#final output of case 4 for fixed effect - LC 

anova_f4 <- anova.f4t3 %>% left_join(anova_f4t3, by = "name")%>%
  rename(sov = name)%>%
  select(sov, Df,Sum_Sq=starts_with("Sum Sq"), 
         Mean_Sq = starts_with("Mean Sq"),
         F_value=starts_with("F value"), Prob_F)

#level of significance of random effect
#null model for CLT

fit.f4c <- lmer(Trait~ 1 + LC + (1|YR) +  (1|YR:LC) + (1|YR:LC) + 
                  (1|YR:CLT)+  
                  (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova4c <- anova_lrt(anova4c,fit.f4,fit.f4c,"CLT")

#null model for LC

fit.f4l <- lmer(Trait~ 1 + (1|CLT) + (1|YR) + (1|YR:LC) + (1|YR:CLT) + 
                  (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova4l <- anova_lrt(anova4l,fit.f4,fit.f4l,"LC")

#null model for YR

fit.f4y <- lmer(Trait~ 1 + LC + (1|CLT) + (1|YR:LC) + (1|YR:CLT) + 
                  (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova4y <- anova_lrt(anova4y,fit.f4,fit.f4y,"YR")

#null model for YR:LC

fit.f4yl <- lmer(Trait~ 1 + LC + (1|CLT) + (1|YR) + (1|YR:CLT) + 
                   (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova4yl <- anova_lrt(anova4yl,fit.f4,fit.f4yl,"YR:LC")

#null model for YR:CLT

fit.f4yc <- lmer(Trait~ 1 + LC + (1|CLT) + (1|YR) + (1|YR:LC) +  
                   (1|LC:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova4yc <- anova_lrt(anova4yc,fit.f4,fit.f4yc,"YR:CLT")

#null model for LC:CLT

fit.f4lc <- lmer(Trait~ 1 + LC + (1|CLT) + (1|YR) + (1|YR:LC) + 
                   (1|YR:CLT) + 
                   (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova4lc <- anova_lrt(anova4lc,fit.f4,fit.f4lc,"LC:CLT")

#null model for YR:LC:CLT

fit.f4ylc <- lmer(Trait~ 1 + LC + (1|CLT) + (1|YR) + (1|YR:LC) + 
                    (1|YR:CLT) + 
                    (1|LC:CLT) + (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova4ylc <- anova_lrt(anova4ylc,fit.f4,fit.f4ylc,"YR:LC:CLT")

#null model for RPid or replication

fit.f4r <- lmer(Trait~ 1 + LC + (1|CLT) + (1|YR) + (1|YR:LC) + 
                  (1|YR:CLT) + 
                  (1|LC:CLT) + (1|YR:LC:CLT), data=tempa)

#level of significance
#call function anova_lrt

anova4r <- anova_lrt(anova4r,fit.f4,fit.f4r,"RPid")

#Merge anova and level of significance

anova4 <- bind_rows(anova4c, anova4y)%>% bind_rows(anova4yl)%>% 
  bind_rows(anova4yc)%>% bind_rows(anova4r)%>% 
  bind_rows(anova4lc)%>%bind_rows(anova4ylc)
anova4 <- as.data.frame(anova4)

#Merge final output

anova_lfix <- variance4%>% left_join(anova4 , by ="sov")
anova_lfix$Pr_Chisq[anova_lfix$stddev == 0] <- NA

########################################################################
##      ANOVA Case 5: CLT and LC - Fixed; YR, and RP - Random         ##
########################################################################

fit.f5 <- lmer(Trait~ 1 +  CLT + LC + LC:CLT + (1|YR) + (1|YR:LC) + 
                 (1|YR:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#model summary

summary5 <- summary(fit.f5)

##variance of random factors

variance5 <- as.data.frame(summary5$varcor)

#drop rownames

rownames(variance5) <- NULL
variance5 <- variance5 %>% select (-var1, -var2) %>% 
  rename(sov=grp, Mean_Sq=vcov, stddev=sdcor)

#Type 3 test of hypothesis

anova.f5t3 <- as.data.frame(anova(fit.f5, type="III"))

#convert rownames into column

anova.f5t3$name <- rownames(anova.f5t3)

#drop rownames

rownames(anova.f5t3) <- NULL

#Type 3 test of hypothesis

anova.f5t2 <- anova(fit.f3, type="marginal", test="F")

#model fitness

anovacase5 <- plot(fit.f5, 
                   main="Model fitness Case 5: CLT and LC - Fixed; YR and RP - Random",
                   xlab="Predicated Value", ylab="Residual")

#level of significance for fixed effects

anova_f5t3 <- mixed(Trait~ 1 +  CLT + LC + LC:CLT + (1|YR) +
                      (1|YR:LC) + (1|YR:CLT) + (1|YR:LC:CLT) +  
                      (1|RPid), data=tempa)

anova_f5t3 <- as.data.frame(anova_f5t3$anova_table)

#convert rownames into column

anova_f5t3$name <- rownames(anova_f5t3)

#drop rownames

rownames(anova_f5t3) <- NULL

anova_f5t3 <- anova_f5t3 %>% select(name, Prob_F = starts_with("Pr(>F)"))

#final output of case 5 for fixed effect - CLT and LC 

anova_f5 <- anova.f5t3 %>% left_join(anova_f5t3, by = "name")%>%
  rename(sov = name)%>%
  select(sov, Df,Sum_Sq=starts_with("Sum Sq"), 
         Mean_Sq = starts_with("Mean Sq"),
         F_value=starts_with("F value"), Prob_F)

#level of significance for random effects
#null model for CLT

fit.f5c <- lmer(Trait~ 1 +  LC + LC:CLT + (1|YR) + (1|YR:LC) + 
                  (1|YR:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova5c <- anova_lrt(anova5c,fit.f5,fit.f5c,"CLT")

#null model for LC

fit.f5l <- lmer(Trait~ 1 +  CLT + LC:CLT + (1|YR) + (1|YR:LC) + 
                  (1|YR:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova5l <- anova_lrt(anova5l,fit.f5,fit.f5l,"LC")

#null model for YR

fit.f5y <- lmer(Trait~ 1 +  CLT + LC + LC:CLT + (1|YR:LC) + 
                  (1|YR:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova5y <- anova_lrt(anova5y,fit.f5,fit.f5y,"YR")

#null model for YR:LC

fit.f5yl <- lmer(Trait~ 1 +  CLT + LC + LC:CLT + (1|YR) +  
                   (1|YR:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova5yl <- anova_lrt(anova5yl,fit.f5,fit.f5yl,"YR:LC")

#null model for YR:CLT

fit.f5yc <- lmer(Trait~ 1 +  CLT + LC + LC:CLT + (1|YR) + (1|YR:LC) + 
                   (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova5yc <- anova_lrt(anova5yc,fit.f5,fit.f5yc,"YR:CLT")

#null model for LC:CLT

fit.f5lc <- lmer(Trait~ 1 +  CLT + LC + (1|YR) + (1|YR:LC) + 
                   (1|YR:CLT) + (1|YR:LC:CLT) +  (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova5lc <- anova_lrt(anova5lc,fit.f5,fit.f5lc,"LC:CLT")

#null model for YR:LC:CLT

fit.f5ylc <- lmer(Trait~ 1 +  CLT + LC + LC:CLT + (1|YR) + (1|YR:LC) + 
                    (1|YR:CLT) + (1|RPid), data=tempa)

#level of significance
#call function anova_lrt

anova5ylc <- anova_lrt(anova5ylc,fit.f5,fit.f5ylc,"YR:LC:CLT")

#null model for RPid or replication

fit.f5r <- lmer(Trait~ 1 +  CLT + LC + LC:CLT + (1|YR) + (1|YR:LC) + 
                  (1|YR:CLT) + (1|YR:LC:CLT) , data=tempa)

#level of significance
#call function anova_lrt

anova5r <- anova_lrt(anova5r,fit.f5,fit.f5r,"RPid")

#Merge anova and level of significance

anova5 <- bind_rows(anova5y, anova5yl)%>% 
  bind_rows(anova5yc)%>% bind_rows(anova5r)%>% 
  bind_rows(anova5ylc)
anova5 <- as.data.frame(anova5)

#Merge final output

anova_clfix <- variance5%>% left_join(anova5 , by ="sov")
anova_clfix$Pr_Chisq[anova_clfix$stddev == 0] <- NA

########################################################################
#   Compute Mean and CV of genotype, location, year, rep, environment  #
########################################################################

#Compute environment -  Location by year combination

tempa2 <- tempa %>%
  mutate (ENV = paste(LC,YR, sep='-')) %>%
  #remove missing records
  na.omit() 

head(tempa2)

#validate data

tempa2$YR <- as.factor(tempa2$YR)
tempa2$RP <- as.factor(tempa2$RP)
tempa2$LC <- as.factor(tempa2$LC)
tempa2$CLT <- as.factor(tempa2$CLT)
tempa2$ENV <- as.factor(tempa2$ENV)
tempa2$Trait <- as.numeric(tempa2$Trait)

#compute descritptive statistics
#trait mean over genotype and environment

mean_ge <- tempa2 %>%
  group_by (CLT, ENV) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_ge1 <- mean_ge %>% 
  spread (ENV, Trait) #transpose using library tidyr

mean_ge2 <- as.data.frame(mean_ge1)

#trait mean over genotype and year

mean_gy <- tempa2 %>%
  group_by (CLT, YR) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_gy1 <- mean_gy %>% 
  spread (YR, Trait) #transpose using library tidyr

mean_gy2 <- as.data.frame(mean_gy1)

#trait sum over genotype and year

sum_gy <- tempa2 %>%
  group_by (CLT, YR) %>%
  summarize (Trait = sum(Trait,na.rm=FALSE))

sum_gy1 <- sum_gy %>% 
  spread (YR, Trait) #transpose using library tidyr

sum_gy2 <- as.data.frame(sum_gy1)

#trait standard deviation (sd) over genotype and year

sd_gy <- tempa2 %>%
  group_by (CLT, YR) %>%
  summarize (Trait = sd(Trait,na.rm=FALSE))

sd_gy1 <- sd_gy %>% 
  spread (YR, Trait) #transpose using library tidyr

sd_gy2 <- as.data.frame(sd_gy1)

#trait mean over genotype and location

mean_gl <- tempa2 %>%
  group_by (CLT, LC) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_gl1 <- mean_gl %>% 
  spread (LC, Trait) #transpose using library tidyr

mean_gl2 <- as.data.frame(mean_gl1)

#trait sum over genotype and location

sum_gl <- tempa2 %>%
  group_by (CLT, LC) %>%
  summarize (Trait = sum(Trait,na.rm=FALSE))

sum_gl1 <- sum_gl %>% 
  spread (LC, Trait) #transpose using library tidyr

sum_gl2 <- as.data.frame(sum_gl1)

#trait standard deviation (sd) over genotype and location

sd_gl <- tempa2 %>%
  group_by (CLT, LC) %>%
  summarize (Trait = sd(Trait,na.rm=FALSE))

sd_gl1 <- sd_gl %>% 
  spread (LC, Trait) #transpose using library tidyr

sd_gl2 <- as.data.frame(sd_gl1)

#trait mean over genotype, location and year

mean_gly <- tempa2 %>%
  group_by (CLT, LC, YR) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_gly1 <- mean_gly %>% 
  spread (LC, Trait) #transpose using library tidyr

mean_gly2 <- as.data.frame(mean_gly1)

#trait sum over genotype, location and year

sum_gly <- tempa2 %>%
  group_by (CLT, LC, YR) %>%
  summarize (Trait = sum(Trait,na.rm=FALSE))

sum_gly1 <- sum_gly %>% 
  spread (LC, Trait) #transpose using library tidyr

sum_gly2 <- as.data.frame(sum_gly1)

#trait standard deviation (sd) over genotype, location and year

sd_gly <- tempa2 %>%
  group_by (CLT, LC, YR) %>%
  summarize (Trait = sd(Trait,na.rm=FALSE))

sd_gly1 <- sd_gly %>% 
  spread (LC, Trait) #transpose using library tidyr

sd_gly2 <- as.data.frame(sd_gly1)

#trait mean over genotype, location and rep

mean_glr <- tempa2 %>%
  group_by (CLT, LC, RP) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_glr1 <- mean_glr %>% 
  spread (LC, Trait) #transpose using library tidyr

mean_glr2 <- as.data.frame(mean_glr1)

#trait mean over location and year

mean_ly <- tempa2 %>%
  group_by (LC, YR) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_ly1 <- mean_ly %>% 
  spread (YR, Trait) #transpose using library tidyr

mean_ly2 <- as.data.frame(mean_ly1)

#trait mean over location and replication

mean_lr <- tempa2 %>%
  group_by (LC, RP) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_lr1 <- mean_lr %>% 
  spread (RP, Trait) #transpose using library tidyr

mean_lr2 <- as.data.frame(mean_lr1)

#trait mean over location

mean_l <- tempa2 %>%
  group_by (LC ) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_l1 <- as.data.frame(mean_l)

#trait mean over year

mean_y <- tempa2 %>%
  group_by (YR ) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_y1 <- as.data.frame(mean_y)

#trait mean over genotype

mean_g <- tempa2 %>%
  group_by (CLT ) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_g1 <- as.data.frame(mean_g)

#trait sum over genotype

sum_g <- tempa2 %>%
  group_by (CLT ) %>%
  summarize (Trait = sum(Trait,na.rm=FALSE))

sum_g1 <- as.data.frame(sum_g)

#trait standard deviation (sd) over genotype

sd_g <- tempa2 %>%
  group_by (CLT ) %>%
  summarize (Trait = sd(Trait,na.rm=FALSE))

sd_g1 <- as.data.frame(sd_g)

#trait mean over environment

mean_e <- tempa2 %>%
  group_by (ENV) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))

mean_e1 <- as.data.frame(mean_e)

#trait CV over genotype and location
#CV = (standard deviation /mean) *100 

cv_gl <- tempa2 %>%
  group_by (CLT, LC) %>%
  summarise (Trait_m = mean(Trait,na.rm=FALSE) ,
             Trait_s = sd(Trait,na.rm=FALSE)) %>%
  mutate (Trait = (Trait_s/Trait_m)*100) %>% # CV  
  select (-Trait_m, -Trait_s)

cv_gl1 <- cv_gl %>% 
  spread (LC, Trait) #transpose using library tidyr

cv_gl2 <- as.data.frame(cv_gl1)

########################################################################
##    Compute univariate stability statistics - regression analysis   ##
########################################################################

#Compute regression (slope) and deviation from regression
#compute environmental index

dsterm <- tempa2 %>% 
  group_by (ENV, RP, YR, LC) %>% 
  summarize (ENVTrait = mean(Trait,na.rm=FALSE))

dst02 <- tempa2 %>%
  left_join(dsterm, by=c("ENV", "RP")) %>% #Left join on multiple columns
  arrange (CLT) %>%
  rename (YR= YR.x, LC = LC.x )

#fit model

fit_model <- dst02 %>% 
  group_by(CLT) %>% #group regression analysis by cultivar 
  do (model=lm(Trait~ENVTrait + ENV + RP, data=.))

#parameter estimates

paramlm <- as.data.frame(fit_model %>% tidy(model))
glancelm <- as.data.frame(fit_model %>% glance(model))
augmentlm <- as.data.frame(fit_model %>% augment(model))

outmsed <- lapply(fit_model$model, anova) #anova output
outmsed2 <- as.data.frame(do.call(rbind, outmsed)) #convert list into data.frame

#convert rownames into column

outmsed2$SOV <- rownames(outmsed2)

# drop rownames

rownames(outmsed2) <- NULL

#remove numeric values from string of rownames using function gsub

outmsed2 <- outmsed2 %>% mutate(SOV = gsub("\\d+","",SOV)) 

#extract unique cultivar name and merge to outmsed2 dataset

genotypes <- dst02 %>% select(CLT) %>% distinct (CLT) %>% arrange(CLT)

#Stack 4 times to match number of rows with outmsed2 dataset

genotypes1 <- genotypes %>% bind_rows(genotypes) %>% 
  bind_rows(genotypes) %>% bind_rows(genotypes) %>% arrange(CLT)

#attach list of cultivars to outmsed2

outmsed3 <- as.data.frame(outmsed2 %>% bind_cols(genotypes1))

#transpose outmsed3

outmsed4 <- outmsed3 %>% 
  select (CLT, SOV, MS = starts_with("Mean")) %>%  #rename variables
  filter (SOV != "RP") 

#transpose MS values

MSDS <- outmsed4 %>% 
  spread (SOV, MS) %>% #transpose using library tidyr
  arrange (CLT)

#Transpose degress of freedoms for F-test

FDS3 <- outmsed3 %>%
  filter (SOV != "RP")%>%
  select (CLT, SOV, Df) %>%
  spread (SOV, Df) %>% #transpose using library tidyr
  rename (DF_ENVTrait = ENVTrait, DF_Residuals = Residuals, 
          DF_ENV = ENV)

#Subset parameters - paramlm dataset

REGCOEFGS <- paramlm %>% 
  filter (term == "ENVTrait") %>%
  select (-statistic, -p.value)

#Test and level of significance of regression and deviation from regression
#Merge MSDS, FDS3, REGCOEFGS

slope <- MSDS %>% inner_join (REGCOEFGS, by = "CLT") %>% 
  inner_join (FDS3, by = "CLT") %>%
  rename (MSE = Residuals, LREGMS=ENVTrait, DEVLMS = ENV, 
          BI= estimate, STDERR = std.error)

#test significance levels

slope1 <- slope %>%
  mutate (T_H01 = (BI-1)/STDERR , #Null Hypothesis for slope = 1
          PT_H01 = 2*pt(-abs(T_H01), DF_Residuals),
          F_DEVREG=DEVLMS/MSE, #NULL HYPOTHESIS: PREDICTED-ACTUAL = 0 
          PF_HO0= 1-pf(F_DEVREG, DF_ENV, DF_Residuals)) 

#add legend for level of significance

slope2 <- slope1 %>% 
  mutate (SIG_SLOPE = ifelse(PT_H01 <= 0.001, "***",
                             ifelse(PT_H01 <= 0.01, "**",
                                    ifelse(PT_H01 <= 0.05, "*","")))) %>%
  mutate (SIG_DEVREG = ifelse(PF_HO0 <= 0.001, "***",
                              ifelse(PF_HO0 <= 0.01, "**",
                                     ifelse(PF_HO0 <= 0.05, "*",""))))

#final regression output

options(digits=5)
univariate2 <- slope2 %>%
  mutate (SLOPE = paste(BI,SIG_SLOPE, sep=""), 
          DEVREG = paste(DEVLMS,SIG_DEVREG, sep="") ) %>%
  select (CLT,SLOPE, DEVREG )

########################################################################
#   Compute univariate stability statistics - shukla, ecovalence, YS   #
########################################################################

#Compute Shukla, Wricke Ecovalense, Kangs YS

repno <- tempa2 %>% 
  summarise (total_rep = n_distinct(RP)) #count total number of rep

dstgl <- tempa2 %>% 
  group_by (CLT, LC) %>%
  #Summarize genotype performance across locations 
  summarize (Trait = mean(Trait,na.rm=FALSE)) 

dstgl1 <- dstgl %>% 
  spread (LC, Trait) #transpose values

#convert into data frame so that row containing structure information is deleted

dstgl2 <- as.data.frame(dstgl1) 

#create rownames

rownames(dstgl2) <- dstgl2[ ,1]
shukla <- dstgl2[,-1]

#compute MS error term

tempa3 <- glm(Trait ~ LC + YR + LC:YR + RP %in% (LC:YR) + CLT + CLT:LC + 
                CLT:YR + CLT:LC:YR, family = gaussian , data= tempa2 )

#model summary

summary1 <- summary.glm(tempa3)

#Error SS

error_ss <- as.data.frame(summary1$deviance)
error_ss1 <- error_ss %>%
  #rename variable
  select (Deviance = starts_with("summary"), everything())

#Error DF

error_df <- as.data.frame(summary1$df.residual)
error_df1 <- error_df %>%
  #rename variable
  select (Df = starts_with("summary"), everything())

#MS of error

mse1 <- as.data.frame(error_ss1/error_df1)
mse <- mse1 %>%
  #rename variable
  rename (MS = Deviance)

# MSError is used populated from ANOVA

univariate1a <- stability.par(shukla, rep= repno$total_rep , MSerror=mse$MS, 
                              alpha=0.1, main="Genotype")

#pool results into individual columns

univariate1b <- univariate1a$statistics

#create column genotype from rownames

univariate1b$genotype <- rownames(univariate1b) 
rownames(univariate1b) = NULL #remove rownames
names(univariate1b)[3] <- "significane_sigma"   # rename duplicate name dot  
names(univariate1b)[5] <- "significane_s2"   # rename duplicate name dot 
names(univariate1b)[2] <- "sigma"   # rename
names(univariate1b)[4] <- "ssquare"   # rename

univariate1c <- univariate1a$stability

#create column genotype from rownames

univariate1c$genotype <- rownames(univariate1c) 
rownames(univariate1c) = NULL #remove rownames
names(univariate1c)[8] <- "legend"   # rename variable ... to legend

#Merge 

univariate1d <- univariate1b %>%
  inner_join (univariate1c , by = "genotype") %>% 
  # deselect all columns between Yield and Stab.rating
  select (-Yield: -Stab.rating) %>% 
  # arrange the column order for final output
  select (CLT=genotype, Mean, sigma, 
          significane_sigma, ssquare, 
          significane_s2, Ecovalence,YSi, legend) 


#Final stability statistics
#Merge Univariate2 and Univariate1d

univariate <- univariate2 %>%
  inner_join(univariate1d, by = "CLT") %>%
  mutate (SIGMA=paste(sigma,significane_sigma, ""),
          SIGMA_SQUARE=paste(ssquare,significane_s2, ""),
          YS_Kang =paste(YSi,legend, "")) %>%
  select (Genotype = CLT, Mean, SLOPE, DEVREG, SIGMA, 
          SIGMA_SQUARE,Ecovalence, YS_Kang)

########################################################################
##        Compute location statistics - genotype F ratio across       ##
##        location and environment; location correlation              ## 
########################################################################

#Location values
#F-value of genotype across location
#fit model

fit_modellc <- tempa2 %>% 
  group_by(LC) %>% #group regression analysis by location
  do (model1=lm(Trait~CLT + YR + CLT:YR + RP%in%YR , data=.))

#parameter estimates

paramlmlc <- as.data.frame(fit_modellc %>% tidy(model1))
glancelmlc <- as.data.frame(fit_modellc %>% glance(model1))
augmentlmlc <- as.data.frame(fit_modellc %>% augment(model1))

outmsedlc <- lapply(fit_modellc$model1, anova) #anova output

#convert list into data.frame

outmsedlc2 <- as.data.frame(do.call(rbind, outmsedlc)) 

#convert rownames into column

outmsedlc2$SOV <- rownames(outmsedlc2)

#drop rownames

rownames(outmsedlc2) <- NULL

#remove numeric values from string of rownames using function gsub

outmsedlc2 <- outmsedlc2 %>% mutate(SOV = gsub("\\d+","",SOV)) 

#extract unique location name and merge to outmsedlc2 dataset

location <- dst02 %>% select(LC) %>% distinct (LC) %>% arrange(LC)

#Extract gentype F value for each location

locvalue <- outmsedlc2 %>% 
  filter(SOV == "CLT") %>%
  select (FRatioGenotype = starts_with("F value")) %>%
  bind_cols (location) %>% select (LC, FRatioGenotype)

locvalue <- as.data.frame (locvalue)

#Correlation between location and average location for each genotype 
#compute genotype mean at each location 

glcmean1 <- tempa2 %>% 
  group_by (CLT, LC) %>%
  summarize (glcmean = mean(Trait,na.rm=FALSE)) %>%
  as.data.frame(select (CLT, LC, glcmean))


#compute genotype mean across all location -average location

gmean1 <- tempa2 %>% 
  group_by (CLT ) %>%
  summarize (gmean = mean(Trait,na.rm=FALSE)) %>%
  as.data.frame(select (CLT, gmean))

#merge location mean with average location for each genotype

lgmean <- glcmean1 %>%
  left_join(gmean1, by="CLT") %>%
  arrange(LC) %>% select (-CLT)

#compute correlation with level of significance

lcgcorr <- lgmean %>%
  group_by(LC) %>%
  do(tidy(cor.test(.$glcmean, .$gmean)))

lcgcorr1 <- lcgcorr %>%
  select (LC, Corr_Value = starts_with ("estimate"), 
          Pvalue = starts_with("p.value"))

#post process correlation value

lcgcorr2 <- lcgcorr1 %>%
  mutate (SIG_CORR = ifelse(Pvalue <= 0.001, "***",
                            ifelse(Pvalue <= 0.01, "**",
                                   ifelse(Pvalue <= 0.05, "*",""))))

#concatenate p value symbol with correlation value

lcgcorr3 <- lcgcorr2 %>%
  mutate (LocCorrelation=paste(Corr_Value,SIG_CORR, sep="")) %>%
  select (LC, LocCorrelation)

lcgcorr3 <- as.data.frame(lcgcorr3)

#Final location value table for output
#compute location mean

Locmean <- tempa2 %>%
  group_by (LC ) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))%>%
  select (LC, Mean = starts_with("Trait"))

Locmean <- as.data.frame(Locmean)

#merge all location value outputs for print

LocationValue <- Locmean %>%
  inner_join (locvalue, by = "LC") %>% 
  inner_join(lcgcorr3, by = "LC") %>%
  rename (Location = LC)

########################################################################
###              F-value of genotype across environmen                ##
########################################################################

#fit model

fit_modelen <- tempa2 %>% 
  group_by(ENV) %>% #group regression analysis by location
  do (model2=lm(Trait~CLT + RP , data=.))

#parameter estimates

paramlmen <- as.data.frame(fit_modelen %>% tidy(model2))
glancelmen <- as.data.frame(fit_modelen %>% glance(model2))
augmentlmen <- as.data.frame(fit_modelen %>% augment(model2))

outmseden <- lapply(fit_modelen$model2, anova) #anova output

#convert list into data.frame

outmseden2 <- as.data.frame(do.call(rbind, outmseden)) 

#convert rownames into column

outmseden2$SOV <- rownames(outmseden2)

# drop rownames

rownames(outmseden2) <- NULL

#remove numeric values from string of rownames using function gsub

outmseden2 <- outmseden2 %>% mutate(SOV = gsub("\\d+","",SOV)) 

#extract unique environment name and merge to outmseden2 dataset

environment <- dst02 %>% select(ENV) %>% distinct (ENV) %>% arrange(ENV)

#Extract gentype F value for each location

locvalue2 <- outmseden2 %>% 
  filter(SOV == "CLT") %>%
  select (FRatioGenotype = starts_with("F value")) %>%
  bind_cols (environment) %>% select (ENV, FRatioGenotype)

locvalue2 <- as.data.frame (locvalue2)

########################################################################
###               Compute cluster analysis of location               ###
########################################################################

#location cluster analysis 
#Euclidean distance
#Ward Hierarchical Clustering

#trait mean over location

mean_l <- tempa2 %>%
  group_by (LC ) %>%
  summarize (Trait = mean(Trait,na.rm=FALSE))
mean_l1 <- as.data.frame(mean_l)

clusterdata <- mean_l1 %>% select (Trait) 
clusterdata <- na.omit(clusterdata)
distance <- dist(clusterdata, method = "euclidean") # distance matrix
hcluster <- hclust(d=distance, method="ward.D")
locationcluster <- plot(hcluster, labels=mean_l1$LC) # display dendogram



########################################################################
####                        Clear console                           ####
########################################################################

#clear console
cat("\014")

########################################################################
###                  Export output in .CSV or .TXT                   ###
###                  Print final output in console                   ###
########################################################################

#user can turn on and off .csv or .txt file by commenting/uncommenting 
#codes in below 2 lines

#sink(file="RGxEOutput.csv", append=FALSE, split=TRUE)
sink(file="RGxEOutput.txt", append=FALSE, split=TRUE)
""
""
"#############      Print date and time     #################"
""
print(Sys.time())
""
"########################################################################"
"###           Section 1: Identify level of significance of           ###"
"###              different effects, variances and BLUP               ###"
"########################################################################"
""
""
"########################################################################"
"###                ANOVA - ANALYSIS OF VARIANCE                     ####"
"########################################################################"
""
"########################################################################"
"##            ANOVA Case 1: CLT, YR, LC and RP - All Random           ##"
"########################################################################"
""
"##Print ANOVA"
"#Note: P values are generated using LRT-Likelihood Ratio Test via "
"#model comparison and anova"
""
print(anova_randall, row.names = FALSE)
""
"##Print BLUP-Best linear unbiased predictor value for random genotypes"
""
print(BLUP_CLT4, row.names = FALSE)
""
"########################################################################"
"##         ANOVA Case 2: CLT, YR and LC - Fixed; RP - Random          ##"
"########################################################################"
""
"##Print ANOVA"
"#Note: P values are computed for F ratio using KR method "
"#Kenward-Rogers approximation for degrees of freedom"
""
print(anova_f2a, row.names = FALSE)
""
"########################################################################"
"##          ANOVA Case 3: CLT - Fixed; YR, LC and RP - Random         ##"
"########################################################################"
""
"##Print ANOVA"
"#Note: For fixed effects P values are computed for F ratio using " 
"#KR method Kenward-Rogers approximation for degrees of freedom"
""
print(anova_f3, row.names = FALSE)
""
"#For random effect P values are generated using LRT-Likelihood Ratio Test "
"#via model comparison and anova"
""
print(anova_cfix, row.names = FALSE)
""
"########################################################################"
"##         ANOVA Case 4: LC - Fixed; YR, CLT and RP - Random          ##"
"########################################################################"
""
"##Print ANOVA"
"#Note: For fixed effects P values are computed for F ratio using"  
"#KR method Kenward-Rogers approximation for degrees of freedom"
""
print(anova_f4, row.names = FALSE)
""
"#For random effect P values are generated using LRT-Likelihood Ratio Test"  
"#via model comparison and anova"
""
print(anova_lfix, row.names = FALSE)
""
"########################################################################"
"##        ANOVA Case 5: CLT and LC - Fixed; YR, and RP - Random       ##"
"########################################################################"
""
"##Print ANOVA"
"#Note: For fixed effects P values are computed for F ratio using " 
"#KR method Kenward-Rogers approximation for degrees of freedom"
""
print(anova_f5, row.names = FALSE)
""
"#For random effect P values are generated using LRT-Likelihood Ratio Test"  
"#via model comparison and anova"
""
print(anova_clfix, row.names = FALSE)
""
""
"########################################################################"
"###       Section 2: Identify general statistcs (mean and CV) of     ###"
"###       genotype, location, year, replication and environment      ###"
"########################################################################"
""
""
"########################################################################"
"###                     Descriptive Statistics                       ###"
"########################################################################"
""
"#Descriptive Statistics - Means and CV"
""
"#trait mean across genotype and environment(location x year combination)"
""
print(mean_ge2, row.names = FALSE)
""
"#trait mean across genotype and years"
""
print(mean_gy2, row.names = FALSE)
""
"#trait mean across genotype and location"
""
print(mean_gl2, row.names = FALSE)
""
"#trait mean across genotype, location and year"
""
print(mean_gly2, row.names = FALSE)
""
"#trait mean across genotype, location and replication"
""
print(mean_glr2, row.names = FALSE)
""
"#trait mean across location and year"
""
print(mean_ly2, row.names = FALSE)
""
"#trait mean across location and replication"
""
print(mean_lr2, row.names = FALSE)
""
"#trait mean across location"
""
print(mean_l1, row.names = FALSE)
""
"#trait mean across year"
""
print(mean_y1, row.names = FALSE)
""
"#trait mean across genotype"
""
print(mean_g1, row.names = FALSE)
""
"#trait mean across environment"
""
print(mean_e1, row.names = FALSE)
""
"#trait coefficient of variation (cv) across gentoype and location"
""
print(cv_gl2, row.names = FALSE)
""
"trait sum across genotype and location"
""
print(sum_gl2, row.names = FALSE)
""
"trait sum across genotype, location and year"
""
print(sum_gly2, row.names = FALSE)
""
"trait sum across genotype and year"
""
print(sum_gy2, row.names = FALSE)
""
"trait sum across genotype"
""
print(sum_g1, row.names = FALSE)
""
"trait standard deviation (sd) across genotype and location"
""
print(sd_gl2, row.names = FALSE)
""
"trait standard deviation (sd) across genotype, location and year"
""
print(sd_gly2, row.names = FALSE)
""
"trait standard deviation (sd) across genotype and year"
""
print(sd_gy2, row.names = FALSE)
""
"trait standard deviation (sd) across genotype"
""
print(sd_g1, row.names = FALSE)
""
""
"########################################################################"
"###             Section 3: Choose stable and best genotype           ###"
"########################################################################"
""
""
"########################################################################"
"###                      Stability Statistics                        ###"
"########################################################################"
""
"#Univariate statistics - mean, slope, deviation from regression," 
"#Shukla,Ecovalence, Kang"
""
print(univariate, row.names = FALSE)
""
""
"########################################################################"
"##    Section 4: Choose discriminative and representative location    ##"
"########################################################################"
""
""
"########################################################################"
"###       location statistics - genotype F ratio across location     ###"
"###                and environment; location correlation             ###"
"########################################################################"
""
"#location value"
"# location mean, genotype F ratio across location," 
"#correlation of location with average location performace"
""
print(LocationValue, row.names = FALSE)
""
"#location value-genotype F ratio across location,"
""
print(locvalue2, row.names = FALSE)
sink()

########################################################################
###             ANOVA model fitness plots : Case 1 to 5              ###
########################################################################

#Print model fitness of Case I to Case V
#pdf file is sent to same folder where input file is located

#Model fitness Case 1: CLT, YR, LC and RP - All Random
#Model fitness Case 2: CLT, YR and LC - Fixed; RP - Random
#Model fitness Case 3: CLT - Fixed; YR, LC and RP - Random
#Model fitness Case 4: LC - Fixed; YR, CLT and RP - Random
#Model fitness Case 5: CLT and LC - Fixed; YR and RP - Random

pdf("modelfitness.pdf")
anovacase1
anovacase2
anovacase3
anovacase4
anovacase5
dev.off()

########################################################################
###            Cluster analysis of location performance              ###
########################################################################

#print dendogram for location
#pdf file is sent to same folder where input file is located
dev.print(pdf, 'locationcluster.pdf')
dev.off()

#End



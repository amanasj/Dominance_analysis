library(tidyverse)
library(openxlsx)
library(readxl)
library(lme4)
library(stats)
library(effects)
library(gtools)
library(patchwork)
library(flexplot)
library(interactions)
library(sjPlot)
library(performance)
library(emmeans)
library(kableExtra)
library(caret)
library(domir)
library(MuMIn)
library(plotrix)


data0 <- suppressMessages(read_excel(("../data/INPUT_SPREADSHEET_HERE.xlsx"), 
                                            sheet = 1, 
                                            col_names = T))



data <- data0[!is.na(data0$ID), ]

data[c(1:3)] <- lapply(data[c(1:3)], as.factor)
data[c(6:51)] <- lapply(data[c(6:51)], as.numeric)


df_OD <- data[ , !grepl( "OS" , names( data ) ) ]
df_OS <- data[ , !grepl( "OD" , names( data ) ) ]


df_OD <- df_OD %>% filter(df_OD$tx!="OD")
df_OS <- df_OS %>% filter(df_OS$tx!="OS")

#names(df_OD)
#names(df_OS)

colnames(df_OD) <-  sub(" OD", "", colnames(df_OD))
colnames(df_OS) <-  sub(" OS", "", colnames(df_OS))

df_untx <- rbind.data.frame(df_OD, df_OS)


df <- df_untx

levels(df$time)
df <- df %>% filter(time!="Day 1" & time!="Week 1")

df$time <- factor(df$time, levels = c("Baseline",
                                      "Month 1", "Month 3", "Month 6",
                                      "Month 9", "Month 12", "Month 18",
                                      "Month 24", "Month 30", "Month 36",
                                      "Month 40", "Month 42", "Month 48",
                                      "Month 54", "Month 60"))

levels(df$time)

names(df) 

#### remove any covariates with 2 or fewer observations per patient
df <- df[-c(10:17)]




### plot select predictor (change x variable as needed)
ggplot(aes(x=LLVA, y=BAF_area), data=df) +
  geom_point() +
  geom_smooth(aes(group=ID), method = "lm") +
  #geom_line(aes(group=ID), data=df[!is.na(df$MP),]) +
  facet_wrap(~ID, scales = "free_x")




df$time <- as.character(df$time)
df$time[df$time=="Baseline"] <- 0
df$time[df$time=="Month 1"] <- 1
df$time[df$time=="Month 3"] <- 3
df$time[df$time=="Month 6"] <- 6
df$time[df$time=="Month 9"] <- 9
df$time[df$time=="Month 12"] <- 12
df$time[df$time=="Month 18"] <- 18
df$time[df$time=="Month 24"] <- 24
df$time[df$time=="Month 30"] <- 30
df$time[df$time=="Month 36"] <- 36
df$time[df$time=="Month 40"] <- 40
df$time[df$time=="Month 42"] <- 42
df$time[df$time=="Month 48"] <- 48
df$time[df$time=="Month 54"] <- 54
df$time[df$time=="Month 60"] <- 60
df$time <- as.numeric(df$time)

names(df)

library(lubridate)
df$dob <- ymd(df$dob)
df$baseline_date <- ymd(df$baseline_date)
df$age_at_bl <- interval(start= df$dob, end=df$baseline_date)/duration(n=1, unit="years")


df_stats <- df
df_stats <- df_stats[c(1:5, 20,21, 6:17)]

df[c(6:20)] <- scale(df[c(6:20)])

#summary(df_stats)

################################################
############## Model selection #################
################################################


################################################################################
##### Lasso regression for variable selection - doesn't deal with missing data 
##### so either have to remove or impute. Lets use glmulti instead.. see below
################################################################################
# d <- df[c(1,3:5,7,26)]
# d <- na.omit(d)
# 
# library(glmmLasso)
# lasso <- glmmLasso(fix=BAF_area ~ time+LLVA+BCVA+MP,
#           rnd = list(ID=~1), 
#           data=d, 
#           lambda=10, 
#           family = gaussian(link="identity"),
#           switch.NR=FALSE, final.re=FALSE, control = list())
# 
# summary(lasso)
#################################################################################






########################################################################
########################################################################
# glmulti lets you find which variables are relevant to include in a model. 
# Domin then tells you the relative importance of each of those included variables
# N.B. DO NOT USE DOMIN FOR VARIABLE SELECTION.. ONLY VARIABLE IMPORTANCE
# Conversely do not use glmulti for variable importance (even though it gives 
# this as an output) since I don't think it doesn't deal with multicollinearity.
#######################################################################
library(glmulti)
library(flextable)
glmer.glmulti <- function(formula, data, random = "", ...){
      glmer(paste(deparse(formula),random), 
            data = data, ...)
  }

names(df)

## try glmulti plackage to compare all possible models
mixed_model <- glmulti(BAF_area ~ Pelli_Rob+time+LLVA+BCVA+MP+time+
                                  camb_contrast_SF1_mean+
                                  camb_contrast_SF2_mean+
                                  camb_contrast_SF4_mean+
                                  camb_contrast_SF8_mean+
                                  camb_contrast_SF10_mean+
                                  cct_Protan+cct_Deutan+cct_Tritan,
                       data   = df,
                       random  = "+(1|ID)",
                       crit   = aicc,       # AICC corrected AIC for small samples
                       level  = 1,          # 2 with interactions, 1 without
                       method = "g",        # "d", or "h", or "g" - use genetic (g) algorithm
                       family=gaussian(),
                       marginality = F,
                       confsetsize = 100,
                       fitfunc = glmer.glmulti   # Type of model (LM, GLM etc.)
)

print(mixed_model)
plot(mixed_model)
weightable(mixed_model)[1:10,] %>%
  regulartable() %>%       # beautifying tables
  autofit()

#tiff(filename = "plot1.tiff", width=2000, height=1500, res = 300)
plot(mixed_model, type = "s", cex.names=0.5)
#dev.off()



#################################################
######### Fit model and check model assumptions
#################################################
lmm1 <- lme4::lmer(BAF_area ~ MP+
                              LLVA+
                              BCVA+
                              time+
                              camb_contrast_SF1_mean+
                              camb_contrast_SF4_mean+
                              camb_contrast_SF8_mean+
                              camb_contrast_SF10_mean+
                              cct_Protan + cct_Deutan + cct_Tritan +
                              Pelli_Rob +
                              (1|ID), data = df, REML = T)

lmm1@call$formula   

sjPlot::tab_model(lmm1,
                  show.reflvl = F,
                  show.intercept = T,
                  p.style = "numeric_stars", digits.p = 2)

#performance::r2_nakagawa(lmm1)
#compare_performance(lmm1,lmm2, rank = T)

### check validity of model using performance package
model_performance(lmm1)   ## high marginal and conditional R^2 - model fits the data very well
#MuMIn::r.squaredGLMM(lmm1)
checks <- check_model(lmm1)
checks

plot(check_normality(lmm1), effects=c("fixed"), type="qq")
#plot(check_heteroscedasticity(lmm1))
#check_heteroscedasticity(lmm1)   
plot(check_outliers(lmm1))
#plot(check_residuals(lmm1))

qqnorm(residuals(lmm1))
qqline(residuals(lmm1), col = "maroon4", lwd = 2)
hist(residuals(lmm1)) 

par(pty="s")
vif <- plot(check_model(lmm1, check = "vif"))
qq <- plot(check_model(lmm1, check = "qq"))
norm <- plot(check_model(lmm1, check = "normality"))
reqq <- plot(check_model(lmm1, check = "reqq"))
lin <- plot(check_model(lmm1, check = "linearity"))  
out <- plot(check_model(lmm1, check = "outliers"))  



patchwork <- (qq | norm)/(out | reqq) 
#tiff(filename="model_check.tiff",width=4000,height=2000,compression="lzw",res=350)
patchwork
#dev.off()





#############################################
###### HOW TO DEAL WITH HETEROSCEDASTICITY
#############################################
## see https://stats.stackexchange.com/questions/91872/alternatives-to-one-way-anova-for-heteroskedastic-data/91881#91881
## for how to deal with heteroscedasticity
## A rule of thumb is that linear models are fairly robust to heterogeneity of variance so long as 
## the maximum variance is no more than 4× greater than the minimum variance, so we'll find that 
## ratio as well
#apply(df, 2, function(x){ var(x, na.rm=T) })
#var(df$BCVA, na.rm=T)/var(df$MP, na.rm=T)      ## <3X so ok. Besides was homoscedastic anyway.

#================================#================================#=============================================
### try a bayesian model and compare to above
#=============================================
lmm2 <- brms::brm(BAF_area ~ MP+
           LLVA+
           BCVA+
           time+
           camb_contrast_SF1_mean+
           camb_contrast_SF4_mean+
           camb_contrast_SF8_mean+
           camb_contrast_SF10_mean+
           cct_Protan + cct_Deutan + cct_Tritan +
           Pelli_Rob +
           (1|ID), data = df)


sjPlot::tab_model(lmm1, lmm2,
                  show.reflvl = F,
                  show.intercept = T,
                  p.style = "numeric_stars", digits.p = 2,
                  file = "model_comparison.doc")



### estimates are almost identical. Go with lmm1

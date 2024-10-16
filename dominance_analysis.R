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
theme_set(theme_minimal())


data0 <- suppressMessages(read_excel(("../data/INSERT_SPREADSHEET_HERE.xlsx"), 
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




### plot
ggplot(aes(x=Pelli_Rob, y=BAF_area), data=df) +
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

#summary(df_stats)

df <- df[-c(2,4,5,18,19,21)]

df[c(2:15)] <- scale(df[c(2:15)])

names(df)






####################################################################################################
#########################################  Dominance analysis  #####################################
####################################################################################################

    domin <-     domir::domin(BAF_area ~
                    MP+
                    LLVA+
                    BCVA+
                    time+
                    camb_contrast_SF1_mean+
                    camb_contrast_SF4_mean+
                    camb_contrast_SF8_mean+
                    camb_contrast_SF10_mean+
                    cct_Protan + cct_Deutan + cct_Tritan +
                    Pelli_Rob,
      lmer,
      #fitstat = list(performance::r2, "R2"),
      #fitstat = list(\(x) list(r2.nagelkerke = as.numeric(performance::r2(x))), "r2.nagelkerke"),
      #fitstat = list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"),
      fitstat = list(\(x) list(R2c = MuMIn::r.squaredGLMM(x)[[2]]), "R2c"),
      #fitstat = list(function(x) list(aic = extractAIC(x)[[2]]), "aic"),
      data = df,
      complete = T,
      #reverse = TRUE,
      #sets = list(
      #             c("time","LLVA","BCVA","MP")),
      #             c("LLVA","BCVA","MP"),
      #             c("LLVA","time","MP"),
      #             c("time","MP"),
      #             c("LLVA","MP")),
      #             c("BCVA","MP"),
      #             c("BCVA","LLVA")
      #            ),
      consmodel = "(1|ID)"
      )

    



###################################################################################################
############ Take a sample of patients with replacement in order to Bootstrap the domin results
###################################################################################################
set.seed(321)
domin_list <- list()

resamp_func <- function(dat, 
                        nboot, 
                        cluster = c("ID"), 
                        replace = TRUE){
  
  for (i in 1:nboot){ 
    
    tryCatch({
      
  # boot expects an indices argument but the sampling happens
  # via sample() as in the original source of the function
  
  # sample the clustering factor (i.e. ID)
  cls <- sample(unique(dat[[cluster[1]]]), replace=replace)
  
  # subset on the sampled clustering factors
  sub <- lapply(cls, function(b) subset(dat, dat[[cluster[1]]]==b))
  
  # join and return samples
  sub <- do.call(rbind, sub)
  
  # UNCOMMENT HERE TO SEE SAMPLED SUBJECTS 
   #print(sub)
  
   domin <-     domir::domin(BAF_area ~ 
                                MP+
                                LLVA+
                                BCVA+
                                time+
                                camb_contrast_SF1_mean+
                                camb_contrast_SF4_mean+
                                camb_contrast_SF8_mean+
                                camb_contrast_SF10_mean+
                                cct_Protan + cct_Deutan + cct_Tritan +
                                Pelli_Rob,
                                #lmer, 
                                reg =  function(y) lmer(formula = y, data = sub, REML = TRUE,
                                                        control = lmerControl(optimizer ="Nelder_Mead")),
                                fitstat = list(\(x) list(R2c = MuMIn::r.squaredGLMM(x)[[2]]), "R2c"), 
                                #fitstat = list(function(x) list(aic = extractAIC(x)[[2]]), "aic"),
                                #reverse = TRUE,
                                #data = sub, 
                                complete = T,
                                consmodel = "(1|ID)"
   )
   
    domin_df <- cbind.data.frame(domin$General_Dominance, 
                                 domin$Standardized, 
                                 domin$Ranks, 
                                 domin$Fit_Statistic_Overall,
                                 domin$Fit_Statistic_Constant_Model
                                 )
    domin_list[i] <- list(domin_df)
    
    }, error=function(e){cat("ERROR in da house:",conditionMessage(e), "\n")})
 }
  domin_list
} 


domin_list <- resamp_func(df, nboot=1500)



#####################################################
#####################################################
####  reformat bootstrapped results to get CI's
####
domin_list_2 <- domin_list[lengths(domin_list) != 0]
boot_df <- do.call(rbind, Map(cbind, bootstrap = seq_along(domin_list_2), domin_list_2))
boot_df$test <- row.names(boot_df)
rownames(boot_df) <- NULL
boot_df <- boot_df[c(1,7,2:6)]
boot_df$test <- gsub("[0-9]+$", "", boot_df$test)
sum(boot_df$`domin$Standardized`)

##look at a single bootstrapped result in more detail
domin_example <- domin_list_2[[997]]
sum(domin_example$`domin$Standardized`) ##always adds up to 1.. use this fit statistic

#############################################
###### create table for saving to docx
library(flextable)
domin_example <- cbind(Functional_test = rownames(domin_example), domin_example)
rownames(domin_example) <- 1:nrow(domin_example)
domin_example %>%
  flextable() %>%
  save_as_docx(path = "domin_table.docx")
##############################################


### function to find the mean and 95% CI for each test from bootstrapped samples
boot_stats <- function(testtype){
    boot_df_single <- boot_df %>% filter(boot_df$test == testtype)
    hist(boot_df_single$`domin$Standardized`)
    boot_df_single <- boot_df_single[order(boot_df_single$`domin$Standardized`),]
    
    boot_df_single <- boot_df_single %>% filter(boot_df_single$`domin$Standardized` > 0)
    mean_single <- mean(boot_df_single$`domin$Standardized`)
    lower_bound_single <- boot_df_single$`domin$Standardized`[round(0.025*nrow(boot_df_single))]
    upper_bound_single <- boot_df_single$`domin$Standardized`[round(0.975*nrow(boot_df_single))]
    mean_overall_fit <- mean(boot_df_single$`domin$Fit_Statistic_Overall`[1])
    
    
    single_stats <- c(test=boot_df_single$test[1], mean=round(mean_single, digits = 3), 
                   "2.5%"=round(lower_bound_single,digits=3), "97.5%"=round(upper_bound_single, digits=3),
                   overall_fit=round(mean_overall_fit, digits=3)
    )
    single_stats
    
}

stats_MP <- boot_stats(testtype="MP")
stats_time <- boot_stats(testtype="time")
stats_pellirob <- boot_stats(testtype="Pelli_Rob")
stats_csf1 <- boot_stats(testtype="camb_contrast_SF1_mean")
stats_csf4 <- boot_stats(testtype="camb_contrast_SF4_mean")
stats_csf8 <- boot_stats(testtype="camb_contrast_SF8_mean")
stats_csf10 <- boot_stats(testtype="camb_contrast_SF10_mean")
stats_protan <- boot_stats(testtype="cct_Protan")
stats_deutan <- boot_stats(testtype="cct_Deutan")
stats_tritan <- boot_stats(testtype="cct_Tritan")
stats_llva <- boot_stats(testtype="LLVA")
stats_bcva <- boot_stats(testtype="BCVA")



Final_DF <- rbind.data.frame(stats_MP,stats_time,stats_pellirob,stats_csf1,
                  stats_csf4,stats_csf8,stats_csf10,stats_protan,
                  stats_deutan,stats_tritan,stats_llva,stats_bcva)
colnames(Final_DF) <- c("Functional Test","Variable Dominance", 
                        "2.5%", "97.5%", "Overall model fit")
Final_DF[c(2:5)] <- lapply(Final_DF[c(2:5)], as.numeric)


Final_DF[1, 1] <- "Microperimetry"
Final_DF[2, 1] <- "Time"
Final_DF[3, 1] <- "Pelli-Robson"
Final_DF[4, 1] <- "CSF-SF1"
Final_DF[5, 1] <- "CSF-SF4"
Final_DF[6, 1] <- "CSF-SF8"
Final_DF[7, 1] <- "CSF-SF10"
Final_DF[8, 1] <- "CCT-Protan"
Final_DF[9, 1] <- "CCT-Deutan"
Final_DF[10, 1] <- "CCT-Tritan"
Final_DF[11, 1] <- "LLVA"
Final_DF[12, 1] <- "BCVA"

mean_overall_fit <- round(mean(Final_DF$`Overall model fit`), digits = 3)

Final_DF <- Final_DF %>% arrange(desc(`Variable Dominance`))


#####################################################




library(ggtext)
ggplot(data = Final_DF, aes(y=reorder(`Functional Test`, `Variable Dominance`), 
                            x=`Variable Dominance`)) +
  geom_bar(stat = "identity",
           aes(fill = `Functional Test`),
           position = position_dodge(), 
           show.legend = FALSE) +
  scale_fill_manual(values = c("lightblue4", "mistyrose4", "peachpuff3","wheat2",
                               "lightblue4", "mistyrose4", "peachpuff3","wheat2",
                               "lightblue4", "mistyrose4", "peachpuff3","wheat2"),
                    aesthetics = "fill") +
  geom_point(size=1, alpha=0.5) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.2, size=0.8, alpha=0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.6)) +
  theme(axis.title=element_text(size=26)) +
  theme(axis.title.x = element_text(vjust = -2)) +
  theme(axis.text.x = element_text(size=18, face="bold", color = "black"),
        axis.text.y = element_text(size=18, face="bold", color = "black")) +
  xlab(bquote(atop('Variable Dominance', ~~(R^'2' ~'contribution')))) +
  ylab("Functional Test") +
  geom_richtext(x=0.22, y=5.5, size=9, fill = c("cornsilk"), colour="darkblue",
                label = "<span>Longitudinal model: <br>
                <i>R<sup>2</sup><sub>(cond)</i></sub> = 0.996</span>") +
  coord_cartesian(xlim = c(0.01,0.295))
#dev.off()



######################################################################################################

####################################################
############## statistics/demographics #############
####################################################
df_baseline <- df_stats %>% filter(time==0)
stats <- sapply(df_baseline, summary)

library(psych)
table <- describe(df_baseline)

table <- table %>% select(mean, sd, median, min, max)
table <- table[-c(1,2,3,4,5), ]

row.names(table) <- c("BAF area (mm^2)","Age (yrs)","BCVA (letters)",
                      "LLVA (letters)","Microperimetry (dB)","Pelli-Robson (letters)",
                      "CSF-SF1","CSF-SF2","CSF-SF4",
                      "CSF-SF8","CSF-SF10",
                      "CCT-Protan (u'v' units)","CCT-Deutan (u'v' units)","CCT-Tritan (u'v' units)")
table$functional_measure <- rownames(table)
table <- table[c(6,1:5)]
rownames(table) <- NULL
colnames(table) <- c("Functional measure","Mean","Std deviation","Median","Min","Max")


table_print <- table %>%
  mutate_if(is.numeric, format, digits=2) %>%
  kbl(caption = "") %>%
  kable_classic(full_width = T, html_font = "Cambria", font_size=12) %>%
  kable_styling("striped") %>%
  kable_styling(full_width = F)
table_print
table_print %>% save_kable(file = "baseline_statistics.png", density = 300) 



d1 <- df %>%
  group_by(ID) %>%
  drop_na(BCVA) %>%
  summarize(BCVA = n()) 
d2 <- df %>%
  group_by(ID) %>%
  drop_na(LLVA) %>%
  summarize(LLVA = n()) 
d3 <- df %>%
  group_by(ID) %>%
  drop_na(MP) %>%
  summarize(MP = n()) 
d4 <- df %>%
  group_by(ID) %>%
  drop_na(camb_contrast_SF1_mean) %>%
  summarize(camb_contrast_SF1_mean = n()) 
d5 <- df %>%
  group_by(ID) %>%
  drop_na(cct_Protan) %>%
  summarize(cct_Protan = n()) 
d6 <- df %>%
  group_by(ID) %>%
  drop_na(BAF_area) %>%
  summarize(BAF_area = n()) 
d_tot <- cbind.data.frame(d1, d2[c(2)], d3[c(2)], d4[c(2)], d5[c(2)], d6[c(2)])

d_tot$ID <- as.character(d_tot$ID)
d_tot[1, 1] <- "A"
d_tot[2, 1] <- "B"
d_tot[3, 1] <- "C"
d_tot[4, 1] <- "D"
d_tot[5, 1] <- "E"
d_tot[6, 1] <- "F"
d_tot[7, 1] <- "G"
d_tot[8, 1] <- "H"
d_tot[9, 1] <- "I"
d_tot[10, 1] <- "J"
d_tot[11, 1] <- "K"
d_tot[12, 1] <- "L"


colnames(d_tot) <- c("Patient","BCVA", "LLVA", "Microperimetry", "CSF", "CCT", "BAF")
d_tot <- d_tot[c(1,7,2:6)]

table2_print <- d_tot %>%
  mutate_if(is.numeric, format, digits=2) %>%
  kbl(caption = "") %>%
  kable_classic(full_width = T, html_font = "Cambria", font_size=12) %>%  
  add_header_above(header = c(" " = 1, "no. of visits" = 6)) %>%
  kable_styling("striped") %>%
  kable_styling(full_width = F)

table2_print
#table2_print %>% save_kable(file = "no_of_visits.png", density = 300) 
####################################################







# Zoo microsatellite project script - plots and models for question 3, correlation between F estimates from ROH and microsats
# Olivia Tatro
# last modified 3/25/25

# goal: compare Cassidy's ROH inbreeding estimates with the coancestry and genhet estimates
# ROH F measures:
  # 1. FROH
  # 2. F
  # 3. MLH
# microsat F measures:
  # 1. Lynch Rt (coancestry sheet)
  # 2. Internal relatedness (genhet sheet)
  # 3. PHt (genhet sheet)
# we want to compare FROH and F with LynchRt and IR
# also compare MLH with PHt
# Reminder to filter ROH samples with F_MISS <= 0.3




# set wd, load libraries, load data ---------------------------------------

setwd("/Users/oliviatatro/Downloads/NEC inbreeding depression thesis")
getwd()

library(tidyverse)
library(cowplot)

# ROH inbreeding estimates for each individual
roh <- read.csv("data/ROH_inbreeding.csv")
# inbreeding estimates from coancestry
coancestry <- read.delim("data/InbreedingEstimates.txt", sep = "", header = TRUE)
# inbreeding estimates from genhet
genhet <- read.csv("genhet/genhet_output.csv")
genhet <- select(genhet, -1)





# data wrangling ----------------------------------------------------------

# selecting columns we care about in each df
roh <- select(roh, c(FID, F_MISS, FROH, `F`, MLH))
coancestry <- select(coancestry, c(ID, LynchRt))
genhet <- select(genhet, c(sampleid, PHt, IR))

# filter out ROH samples with F_MISS over 0.3
roh_filtered <- filter(roh, F_MISS <= 0.3)

# merge dfs together
all <- merge(coancestry, roh_filtered, by.x = "ID", by.y = "FID", all.x = TRUE)
all <- merge(all, genhet, by.x = "ID", by.y = "sampleid")







# FROH vs LynchRt ---------------------------------------------------------
# create scatterplot
ggplot(all, aes(x = LynchRt, y = FROH)) +
  geom_point() +
  geom_smooth(method = "lm", col = "cornflowerblue") +
  theme_classic()

# fit a linear model
m <- lm(FROH ~ LynchRt, all)
summary(m)





# FROH vs IR --------------------------------------------------------------
# create scatterplot
froh_ir <- 
  ggplot(all, aes(x = IR, y = FROH)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "blue") +
  xlab("Internal Relatedness") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 14)) +
  annotate("text", x = -0.1, y = 0.65, label = "p < 0.0001", col = "red", size = 5)
ggsave("plots/nrem/froh_ir.png")

# fit a linear model
m <- lm(FROH ~ IR, all)
summary(m)

filter(all, !is.na(FROH) & !is.na(IR))

# this is to look at if the points in the low cluster are missing loci
low_froh <- filter(all, FROH < 0.3)
  # 400625, 400682, Q18028 are not missing any loci
  # 400544, 400563, 400725, 400730, Q18029, Q19051 are missing 1 locus










# F vs LynchRt ------------------------------------------------------------
# create scatterplot
ggplot(all, aes(x = LynchRt, y = `F`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

# fit a linear model
m <- lm(`F` ~ LynchRt, all)
summary(m)






# F vs IR -----------------------------------------------------------------
# create scatterplot
ggplot(all, aes(x = IR, y = `F`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

# fit a linear model
m <- lm(`F` ~ IR, all)
summary(m)







# MLH vs PHt --------------------------------------------------------------
# create scatterplot
mlh_pht <-
  ggplot(all, aes(x = PHt, y = MLH)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "blue") +
  xlab("Proportion heterozygosity") +
  ylab("Multilocus heterozygosity") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 14)) +
  annotate("text", x = 0.4, y = 0.1, label = "p < 0.0001", col = "red", size = 5)

# fit a linear model
m <- lm(MLH ~ PHt, all)
summary(m)


plot_grid(froh_ir, mlh_pht, nrow = 1, ncol = 2, labels = c("A", "B"))
ggsave("plots/q3_figs.png", height = 4, width = 8)






















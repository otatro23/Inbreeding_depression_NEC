# Zoo microsatellite project script - compare parent relatedness to reproductive success - Question 1B
# Olivia Tatro
# last modified 4/12/25 

# goal: we want to compare the relatedness of dams and sires and inbreeding to their reproductive success metrics



# load libraries, set working directory, load data ------------------------
setwd("/Users/oliviatatro/Downloads/NEC inbreeding depression thesis")
getwd()

library(tidyverse)
library(cowplot)
library(gamlss)
library(VGAM)

# relatedness from coancestry
relatedness <- read.delim("data/RelatednessEstimates.txt", sep = ",")
# offspring metadata with their fitness and parent IDs
offspring <- read.csv("data/metadata_offspring_updated.csv")
# reproductive success of each parent pair
parent_pair_metadata <- read.csv("data/parent_pair_metadata.csv")
# inbreeding estimates from coancestry
coancestry <- read.delim("data/InbreedingEstimates.txt", sep = "", header = TRUE)
# inbreeding estimates from genhet
genhet <- read.csv("genhet/genhet_output.csv")
genhet <- select(genhet, -1) #remove index


# data wrangling ----------------------------------------------------------

# create successful_breeding_attempt column take out unneccessary columns
pair_success <- parent_pair_metadata %>% 
  select(DamID, SireID, offspring_produced, prop_weaned_litter, offspring_weaned) %>% 
  mutate(successful_breeding_attempt = offspring_produced > 0)

# clean up RelatednessEstimates to only have sires in one column and dams in another ------------------------------------------------------------------
t1 <- relatedness %>% # makes a df of rows in the relatedness comparisons that have a sire and dam in offspring_metadata
  filter(ID1 %in% offspring$SireID & ID2 %in% offspring$DamID) %>% 
  mutate(SireID = ID1, DamID = ID2)
t2 <- relatedness %>% # same as above but switches the columns in the relatedness sheet because they aren't separated by dam and sire
  filter(ID2 %in% offspring$SireID & ID1 %in% offspring$DamID) %>% 
  mutate(SireID = ID2, DamID = ID1)
relatedness<- rbind(t1, t2) # adds the two dfs together

# add QGT to pair success sheet
pair_success <- merge(pair_success, relatedness, by = c("SireID", "DamID"), all = TRUE)
pair_success <- pair_success %>% 
  select(c(1:6, 16)) # keeps all of the cols from parent_pair and queller-goodnight relatedness

# add genhet dam inbreeding to pair_success sheet
pair_success <- merge(pair_success, genhet, by.x = "DamID", by.y = "sampleid", all = TRUE) %>% 
  rename(dam_PHt = PHt, dam_Hs_obs = Hs_obs, dam_Hs_exp = Hs_exp, dam_IR = IR, dam_HL = HL)
# add genhet sire inbreeding to pair_success sheet
pair_success <- merge(pair_success, genhet, by.x = "SireID", by.y = "sampleid", all = TRUE) %>% 
  rename(sire_PHt = PHt, sire_Hs_obs = Hs_obs, sire_Hs_exp = Hs_exp, sire_IR = IR, sire_HL = HL)

# add coancestry dam LynchRt to pair_success sheet
pair_success <- merge(pair_success, coancestry, by.x = "DamID", by.y = "ID", all = TRUE) %>% 
  rename(dam_LynchRt = LynchRt) 
# add coancestry sire LynchRt to pair_success sheet
pair_success <- merge(pair_success, coancestry, by.x = "SireID", by.y = "ID", all = TRUE) %>% 
  rename(sire_LynchRt = LynchRt) 

# take out extra columns from coancestry inbreeding
pair_success <- select(pair_success, c(1:16, 18, 22))

# take out extra rows from merge
pair_success <- filter(pair_success, !is.na(offspring_produced))

# create a column for added dam/sire inbreeding
pair_success <- mutate(pair_success, add_IR = (dam_IR + sire_IR))

# create a column for multiplied dam/sire inbreeding
pair_success <- mutate(pair_success, mult_IR = (dam_IR * sire_IR))

# create a column, birth_0_1 to do a log regression plot for successful breeding attempt
pair_success <- mutate(pair_success, birth_0_1 = ifelse(successful_breeding_attempt, 1, 0))



#pdf("Q1B_plots_models_ignore_repeats.pdf",height=8,width=11)




# comparing QUELLER GOODNIGHT RELATEDNESS to various measures of success ----------------------------------

# summary stats
filter(pair_success, !is.na(QGT)) %>% 
  #group_by(birth_0_1) %>% 
  summarize(mean = mean(QGT),
            stdev = sd(QGT))


# mean litter size vs queller goodnight relatedness
unfiltered_qgt_litter <- 
  filter(pair_success, offspring_produced > 0) %>% 
  ggplot(., aes(x = QGT, y = offspring_produced)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = FALSE, col = "blue") +
  xlab("Queller-Goodnight Relatedness") +
  ylab("Litter Size") +
  theme_classic() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10), 
    ) +
  annotate("text", x = 0.4, y = 7.5, label = "p = 0.36", col = "black", size = 5)

# Linear regression model
model <- lm(offspring_produced ~ QGT, data = filter(pair_success, offspring_produced > 0))
summary(model) # p-value: 0.374; 0.0044 without the QGT outliers
par(mfrow = c(2,2))
plot(model)

# poisson regression model 
model <- glm(offspring_produced ~ QGT, data = filter(pair_success, offspring_produced > 0))
summary(model)


#filtered_qgt_litter <- 
  filter(pair_success, QGT < 0.5, offspring_produced > 0) %>% 
  ggplot(., aes(x = QGT, y = offspring_produced)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) +
  xlab("Queller-Goodnight Relatedness") +
  ylab("Litter Size") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 10)) +
  annotate("text", x = 0, y = 6.5, label = "p = 0.97", col = "maroon", size = 5)

# Linear regression model
model <- lm(offspring_produced ~ QGT, data = filter(pair_success, QGT < 0.5, offspring_produced > 0))
summary(model) # p-value: 0.374; 0.0044 without the QGT outliers

# poisson regression model 
model <- glm(offspring_produced ~ QGT, data = filter(pair_success, QGT < 0.5, offspring_produced > 0), family = "poisson")
summary(model)


# proportion of kits weaned per litter vs Queller goodnight relatedness
par(mfrow = c(1,1))
#filtered_qgt_wean <- # filtered
  filter(pair_success, QGT < 0.5) %>% 
  ggplot(., aes(x = QGT, y = prop_weaned_litter)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Queller-Goodnight Relatedness") +
  ylab("Proportion of kits weaned per litter") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 10)) +
  scale_x_continuous(limits = c(NA, 0.11)) +
  annotate("text", x = 0.05, y = 0.85, label = "p = 0.036", col = "maroon", size = 5)

# linear model
model <- lm(prop_weaned_litter ~ QGT, data = filter(pair_success, QGT < 0.5))
summary(model)

# beta inflated regression model 
model <- gamlss(prop_weaned_litter ~ QGT, family = BEINF, data = na.omit(filter(pair_success, QGT < 0.5)))
summary(model)
plot(model)


unfiltered_qgt_wean <- 
  ggplot(pair_success, aes(x = QGT, y = prop_weaned_litter)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Pair Queller-Goodnight Relatedness") +
  ylab("Proportion of kits weaned per litter") +
  theme_classic()  +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text = element_text(size = 10)) +
  annotate("text", x = 0.45, y = 0.85, label = "p = 0.010", col = "red", size = 5)

#Linear regression model 
model <- lm(prop_weaned_litter ~ QGT, data = pair_success)
summary(model) # p-value: 0.0101; 0.03636 without outliers
par(mfrow = c(2,2))
plot(model)


# logistic regression (trials = offspring_produced, successes = offspring weaned)
model <- glm(cbind(offspring_weaned, offspring_produced) ~ QGT, family = binomial, data = pair_success)
summary(model)

# plotting the logistic regression model 
pair_success$QGT_wean_predicted <- predict(model, newdata = pair_success, type = "response")
unfiltered_qgt_wean <- 
  ggplot(pair_success, aes(x = QGT, y = prop_weaned_litter))+
  geom_point()+
  geom_smooth(data = pair_success, aes(y = QGT_wean_predicted), col = "blue",  size = 1) +
  xlab("Pairwise Queller-Goodnight Relatedness") +
  ylab("Kit weaning rate") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 10)) +
  annotate("text", x = 0.5, y = 0.9, label = "p = 0.010", col = "red", size = 5)



# plots queller goodnight relatedness for pairings that were and were not successful
par(mfrow = c(1,1))
#qgt_birth <- 
  filter(pair_success, !is.na(successful_breeding_attempt), QGT < 0.5) %>% 
  ggplot(., aes(x = QGT, y = successful_breeding_attempt)) +
  geom_boxplot() +
  ylab(NULL) +
  xlab("Queller-Goodnight Relatedness") +
  #geom_jitter() +
  scale_y_discrete(labels = c("TRUE" = "Birth", "FALSE" = "No birth")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 14)) +
  annotate("text", x = 0.25, y = 2.25, label = "p = 0.0013", col = "maroon", size = 5)


# ttest (is there a difference in relatedness between parents in pairings that produced and did not produce at least 1 kit)

model <- lm(QGT ~ successful_breeding_attempt, data = filter(pair_success, QGT < 0.5))
summary(model) # p-value: 0.6975; 0.0013 when you remove QGT outliers
par(mfrow = c(2,2))
plot(model)


# log regression plot of birth vs qgt unfiltered
unfiltered_qgt_birth <-
  ggplot(pair_success, aes(x = QGT, y = birth_0_1)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1) +
  #geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, col = "blue") +
  xlab("Queller-Goodnight Relatedness") +
  ylab("Probability of birth") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10)) +
  annotate("text", x = 0.4, y = 0.9, label = "p = 0.69", col = "black", size = 5)


model <- glm(birth_0_1 ~ QGT, data = pair_success, family = binomial(link = "logit"))
summary(model) #p=0.686 # 0.010 after taking out outliers
plot(model)


# log regression plot of birth vs qgt filtered qgt < 0.5
filtered_qgt_birth <- 
  filter(pair_success, QGT < 0.5) %>% 
  ggplot(., aes(x = QGT, y = birth_0_1)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1) +
  #geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, color = "blue") +
  xlab("Pairwise Queller-Goodnight Relatedness") +
  ylab("Probability of birth") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 10)) +
  annotate("text", x = 0.25, y = 0.9, label = "p = 0.010", col = "red", size = 5)

model <- glm(birth_0_1 ~ QGT, data = filter(pair_success, QGT < 0.5), family = binomial(link = "logit"))
summary(model) #p=0.010 
plot(model)
ggsave("plots/URC/filtered_qgt_birth.png", height = 5, width = 5)

# for NR 713
plot_grid(unfiltered_log_regression, filtered_log_regression, nrow = 1, ncol = 2, labels = c("A","B"))
ggsave("../NR_713/HW/NR_713_HW_5/birth_qgt.png", width = 6, height = 3)

# plots for NREM 
plot_grid(unfiltered_log_regression, unfiltered_qgt_litter, unfiltered_qgt_wean, nrow = 1, ncol = 3)
ggsave("plots/nrem/relatedness.png", width = 12, height = 4)

plot_grid(filtered_log_regression, filtered_qgt_litter, filtered_qgt_wean, nrow = 1, ncol = 3)
ggsave("plots/nrem/filtered_relatedness.png", width = 12, height = 4)











# compare DAM INTERNAL RELATEDNESS to measures of reproductive success -----------------

# mean litter size vs dam IR
damIR_litter <- 
  filter(pair_success, offspring_produced > 0) %>% 
  ggplot(., aes(x = dam_IR, y = offspring_produced)) +
  geom_point(size = 1) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = FALSE, col = 'blue') +
  xlab("Dam Internal Relatedness") +
  ylab("Litter size") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10)) +
  annotate("text", x = 0.5, y = 7.5, label = "p = 0.49", col = "black", size = 5)


# Linear regression model 
model <- lm(offspring_produced ~ dam_IR, data = filter(pair_success, offspring_produced > 0))
summary(model) #p-value: 0.33
par(mfrow = c(2,2))
plot(model)

# poisson regression model 
model <- glm(offspring_produced ~ dam_IR, data = filter(pair_success, offspring_produced > 0), family = "poisson")
summary(model) # pval = 0.49 

# zero inflated poisson regression model 
zip_model <- vglm(offspring_produced ~ dam_IR, family = zipoisson, data = pair_success)
summary(zip_model)

# proportion weaned per litter vs dam IR
par(mfrow = c(1,1))
par(mfrow = c(1,1))
#damIR_wean <- 
  ggplot(pair_success, aes(x = dam_IR, y = prop_weaned_litter)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Dam Internal Relatedness") +
  ylab("Kit Weaning Rate") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)) +
  annotate("text", x = 0.5, y = 0.9, label = "p = 0.36", col = "black", size = 5)

#Linear regression model 
model <- lm(prop_weaned_litter ~ dam_IR, data = pair_success)
summary(model) # p-value: 0.3618
par(mfrow = c(2,2))
plot(model)


# logistic regression (trials = offspring_produced, successes = offspring weaned)
model <- glm(cbind(offspring_weaned, offspring_produced) ~ dam_IR, family = binomial, data = pair_success)
summary(model)

# plotting the logistic regression model 
pair_success$dam_wean_predicted <- predict(model, newdata = pair_success, type = "response")
damIR_wean <- 
  ggplot(pair_success, aes(x = dam_IR, y = prop_weaned_litter))+
  geom_point()+
  geom_smooth(data = pair_success, aes(y = dam_wean_predicted), col = "blue",  size = 1) +
  xlab("Dam internal relatedness") +
  ylab("Kit weaning rate") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)) +
  annotate("text", x = 0.5, y = 0.9, label = "p = 0.38", col = "black", size = 5)




# plots dam IR for breeding attempts that were and were not successful
par(mfrow = c(1,1))
#damIR_birth <- 
  filter(pair_success, !is.na(successful_breeding_attempt)) %>% 
  ggplot(., aes(x = dam_IR, y = successful_breeding_attempt)) +
  geom_boxplot() +
  ylab(NULL) +
  xlab("Dam Internal Relatedness") +
  #geom_jitter() +
  scale_y_discrete(labels = c("TRUE" = "Birth", "FALSE" = "No birth")) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 14)) +
  annotate("text", x = 0.45, y = 2.25, label = "p = 0.0095", col = "maroon", size = 5)

# ttest comparing dam IR for pairings that had at least 1 kit vs those who did not
model <- lm(dam_IR ~ successful_breeding_attempt, data = pair_success)
model <- lm(successful_breeding_attempt ~ dam_IR, data = pair_success)

summary(model) # p-value: 0.009537
par(mfrow = c(2,2))
plot(model)


#logistic regression for birth vs damIR
damIR_birth <- 
  ggplot(pair_success, aes(x = dam_IR, y = birth_0_1)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1) +
  #geom_point(size = 1) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, col = "blue") +
  xlab("Dam Internal Relatedness") +
  ylab("Probability of birth") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10), 
    axis.title = element_text(size = 14)) +
  annotate("text", x = 0.45, y = 0.9, label = "p = 0.012", col = "red", size = 5)

nrow(filter(pair_success, !is.na(dam_IR)))

# calculate means for each group 
filter(pair_success, !is.na(dam_IR)) %>% 
  group_by(birth_0_1) %>% 
  summarize(mean = mean(dam_IR),
            stdev = sd(dam_IR))

plogis(1.042-0.09031055*(-2.8978)) - plogis(1.042+0.58558559*(-2.8978))

ggsave("../NR_713/HW/NR_713_HW_5/damIR_log_regression.png", width = 4, height = 4)

model <- glm(birth_0_1 ~ dam_IR, data = pair_success, family = "binomial")
summary(model) #p=0.01
plot(model)





plot_grid(damIR_birth, damIR_litter, damIR_wean, nrow = 1, ncol = 3)
ggsave("plots/nrem/dam_IR.png", width = 12, height = 4)









# compare SIRE INTERNAL RELATEDNESS to measures of reproductive success -----------------

# mean litter size vs sire IR
par(mfrow = c(1,1))
filter(pair_success, offspring_produced > 0) %>% 
ggplot(., aes(x = sire_IR, y = offspring_produced)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) +
  xlab("Sire Internal Relatedness") +
  ylab("Litter size") +
  theme_classic()

#Linear regression model 
model <- lm(sire_IR ~ offspring_produced, data = filter(pair_success, offspring_produced > 0))
summary(model) # p-value: 0.074
par(mfrow = c(2,2))
plot(model)

model <- glm(offspring_produced ~ sire_IR, data = filter(pair_success, offspring_produced > 0), family = "poisson")
summary(model) # p-value: 0.177

# proportion weaned per litter vs sire IR
par(mfrow = c(1,1))
ggplot(pair_success, aes(x = sire_IR, y = prop_weaned_litter)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Sire Internal Relatedness") +
  ylab("Proportion of kits weaned per litter") +
  theme_classic()

#Linear regression model 
model <- lm(sire_IR ~ prop_weaned_litter, data = pair_success)
summary(model) # p-value: 0.5183
par(mfrow = c(2,2))
plot(model)

# beta inflated regression model 
model <- gamlss(prop_weaned_litter ~ sire_IR, family = BEINF, data = na.omit(pair_success))
summary(model) # p-value: 0.23


# plots sire IR for breeding attempts that were and were not successful
par(mfrow = c(1,1))
filter(pair_success, !is.na(successful_breeding_attempt)) %>% 
  ggplot(., aes(x = successful_breeding_attempt, y = sire_IR)) +
  geom_boxplot() +
  geom_smooth(method = "lm") +
  ylab("Sire Internal Relatedness") +
  xlab(NULL) +
  geom_jitter() +
  scale_x_discrete(labels = c("TRUE" = "Produced at least 1 kit", "FALSE" = "Did not produce any kits")) +
  theme_classic()

# ttest comparing sire IR for pairings that had at least 1 kit vs those who did not
model <- lm(sire_IR ~ successful_breeding_attempt, data = pair_success)
summary(model) # p-value: 0.7105
par(mfrow = c(2,2))
plot(model)

# logistic regression plot 
ggplot(pair_success, aes(x = sire_IR, y = birth_0_1)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1) +
  #geom_point(size = 1) +
  stat_smooth(method = "glm", method.args = list(family = "binomial")) +
  xlab("Sire Internal Relatedness") +
  ylab("Probability of birth") +
  theme_classic()

# logistic regression model 
model <- glm(birth_0_1 ~ sire_IR, data = pair_success, family = "binomial")
summary(model) # pval = 0.706








# compare ADDED IR to measures of reproductive success -----------------

# mean litter size vs added IR
par(mfrow = c(1,1))
ggplot(pair_success, aes(x = add_IR, y = offspring_produced)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("added Internal Relatedness") +
  ylab("Litter size") +
  theme_classic()

#Linear regression model 
model <- lm(add_IR ~ offspring_produced, data = pair_success)
summary(model) # p-value: 0.3853
par(mfrow = c(2,2))
plot(model)


# proportion weaned per litter vs added IR
par(mfrow = c(1,1))
ggplot(pair_success, aes(x = add_IR, y = prop_weaned_litter)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("added Internal Relatedness") +
  ylab("Proportion of kits weaned per litter") +
  theme_classic()

#Linear regression model 
model <- lm(add_IR ~ prop_weaned_litter, data = pair_success)
summary(model) # p-value: 0.3046
par(mfrow = c(2,2))
plot(model)


# plots added IR for breeding attempts that were and were not successful
par(mfrow = c(1,1))
filter(pair_success, !is.na(successful_breeding_attempt)) %>% 
  ggplot(., aes(x = successful_breeding_attempt, y = add_IR)) +
  geom_boxplot() +
  geom_smooth(method = "lm") +
  ylab("added Internal Relatedness") +
  xlab(NULL) +
  geom_jitter() +
  scale_x_discrete(labels = c("TRUE" = "Produced at least 1 kit", "FALSE" = "Did not produce any kits")) +
  theme_classic()

# ttest comparing added IR for pairings that had at least 1 kit vs those who did not
model <- lm(add_IR ~ successful_breeding_attempt, data = pair_success)
summary(model) # p-value: 0.06649
par(mfrow = c(2,2))
plot(model)

# log regression plot
ggplot(pair_success, aes(x = sire_IR, y = birth_0_1)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1) +
  #geom_point(size = 1) +
  stat_smooth(method = "glm", method.args = list(family = "binomial")) +
  xlab("Sire Internal Relatedness") +
  ylab("Probability of birth") +
  theme_classic()

# log regression model 
model <- glm(birth_0_1 ~ sire_IR, pair_success, family = "binomial")
summary(model)










# compare MULTIPLIED IR to measures of reproductive success -----------------

# mean litter size vs multiplied IR
par(mfrow = c(1,1))
ggplot(pair_success, aes(x = mult_IR, y = offspring_produced)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Multiplied Internal Relatedness") +
  ylab("Litter size") +
  theme_classic()

#Linear regression model 
model <- lm(mult_IR ~ offspring_produced, data = pair_success)
summary(model) # p-value: 0.2137
par(mfrow = c(2,2))
plot(model)


# proportion weaned per litter vs multiplied IR
par(mfrow = c(1,1))
ggplot(pair_success, aes(x = mult_IR, y = prop_weaned_litter)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Multiplied Internal Relatedness") +
  ylab("Proportion of kits weaned per litter") +
  theme_classic()

#Linear regression model 
model <- lm(mult_IR ~ prop_weaned_litter, data = pair_success)
summary(model) # p-value: 0.3987
par(mfrow = c(2,2))
plot(model)


# plots multiplied IR for breeding attempts that were and were not successful
par(mfrow = c(1,1))
filter(pair_success, !is.na(successful_breeding_attempt)) %>% 
  ggplot(., aes(x = successful_breeding_attempt, y = mult_IR)) +
  geom_boxplot() +
  geom_smooth(method = "lm") +
  ylab("Multiplied Internal Relatedness") +
  xlab(NULL) +
  geom_jitter() +
  scale_x_discrete(labels = c("TRUE" = "Produced at least 1 kit", "FALSE" = "Did not produce any kits")) +
  theme_classic()

# ttest comparing multiplied IR for pairings that had at least 1 kit vs those who did not
model <- lm(sire_IR*dam_IR ~ successful_breeding_attempt, data = pair_success)
summary(model) # p-value: 0.02719
par(mfrow = c(2,2))
plot(model)







# save --------------------------------------------------------------------
write_csv(pair_success, "data/pair_success.csv")

dev.off()



# summary stats -----------------------------------------------------------
pair_success %>%
  #filter(QGT < 0.5) %>% 
  summarize(mean_dam_IR = mean(dam_IR, na.rm = TRUE), sd_dam_IR = sd(dam_IR, na.rm = TRUE),
            mean_sire_IR = mean(sire_IR, na.rm = TRUE), sd_sire_IR = sd(sire_IR, na.rm = TRUE),
            mean_qgt = mean(QGT, na.rm = TRUE), sd_QGT = sd(QGT, na.rm = TRUE))

pair_success %>%
  #filter(QGT <0.5) %>% 
  #group_by(birth_0_1) %>% 
  summarize(mean_dam_IR = mean(dam_IR, na.rm = TRUE), sd_dam_IR = sd(dam_IR, na.rm = TRUE),
            mean_sire_IR = mean(sire_IR, na.rm = TRUE), sd_sire_IR = sd(sire_IR, na.rm = TRUE),
            mean_qgt = mean(QGT, na.rm = TRUE), sd_QGT = sd(QGT, na.rm = TRUE))








# Assemble figure for URC -------------------------------------------------
dist_birth <- readRDS("plots/URC/dist_birth.rds")
dist_litter <- readRDS("plots/URC/dist_litter.rds")
dist_wean <- readRDS("plots/URC/dist_wean.rds")

plot_grid(damIR_birth, unfiltered_qgt_birth, dist_birth,
          damIR_litter, unfiltered_qgt_litter, dist_litter,
          damIR_wean, unfiltered_qgt_wean, dist_wean,
          nrow = 3, ncol = 3, align = "hv", axis = "tblr",
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
ggsave("plots/URC/FIGURE.png", width = 12, height = 12)







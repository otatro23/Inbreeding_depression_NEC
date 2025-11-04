# Zoo microsatellite project script - models for spatial effect on reproductive success - Question 2
# Olivia Tatro
# last modified 4/12/25

# goal: plot and model the relationship between geographic distance and focus areas and reproductive success

# set wd, load libraries, load data ---------------------------------------

setwd("/Users/oliviatatro/Downloads/NEC inbreeding depression thesis")
getwd()

library(tidyverse)
library(fields)
library(cowplot)

parent_pair_metadata <- read.csv("data/parent_pair_metadata.csv")
founder_metadata <- read.csv("data/NEC_Breeding_Records_Body_Weights(Animal Location Data).csv")
#unique(sort(founder_metadata$patch))




# data wrangling ----------------------------------------------------------

# clean pair success dataframe
pair_success <- parent_pair_metadata %>% 
  select(DamID, SireID, offspring_produced, prop_weaned_litter, offspring_weaned) %>% # select only necessar columns
  mutate(successful_breeding_attempt = offspring_produced > 0) %>% # create column for if there was a successful breeding attempt
  filter(!is.na(DamID) & !is.na(SireID) & DamID != "" & SireID != "") # filter out blank sire and dam IDs


# add location data to pair success data 
# prepares sires from founder_metadata to be added to the pair success sheet
  # note: this also contains dams in the sireID column, but they will be filtered out when it is merged
sires <- select(founder_metadata, c(2,4,9,10,14)) %>% # keeps ID, focus area, coords, adn patch
  rename(sire_focus_area = focus_area,
         sire_patch = patch,
         sire_x = Trap.coordinate..x.,
         sire_y = Trap.coordinate..y.,
         SireID = Accession..)

# prepares dams from founder_metadata to be added to the pair success sheet
  # note: this also contains sires in the DamID column, but they will be filtered out when it is merged
dams <- select(founder_metadata, c(2,4,9,10,14)) %>% # keeps ID, focus area, and coords, and patch
  rename(dam_focus_area = focus_area,
         dam_patch = patch,
         dam_x = Trap.coordinate..x.,
         dam_y = Trap.coordinate..y.,
         DamID = Accession..)


# adds sire location data based on SireID and dam location data based on DamID
location_success_data <- merge(pair_success, sires, by = "SireID")
location_success_data <- merge(location_success_data, dams, by = "DamID")


# create a column for matching focus area and patch
location_success_data <- location_success_data %>% 
  mutate(matching_focus_area = sire_focus_area == dam_focus_area) %>% 
  mutate(matching_patch = sire_patch == dam_patch)
#unique(location_success_data$dam_focus_area)

# calculate distances
location_success_data <- location_success_data %>%
  filter(!is.na(sire_x) & !is.na(dam_x)) %>% # filter out pairs that don't have coords for both parents
  mutate(dist = NA)

# go through each row and calculate distance between dam/sire coords and add that to the "dist" column
for(i in 1:nrow(location_success_data)){
  sire_coords <- as.matrix(location_success_data[i,8:9])
  dam_coords <- as.matrix(location_success_data[i,12:13])
  location_success_data[i,17] <- as.numeric(rdist.earth(sire_coords, dam_coords, miles = FALSE) )
}


# create a column, birth_0_1 to do a log regression plot for successful breeding attempt
location_success_data <- mutate(location_success_data, birth_0_1 = ifelse(successful_breeding_attempt, 1, 0))









# compare FOCUS AREAS to measures of reproductive success -----------------

# litter size vs matching_focus_area 
filter(location_success_data, offspring_produced > 1) %>% 
ggplot(., aes(x = matching_focus_area, y = offspring_produced)) +
  geom_boxplot() +
  geom_jitter(size = 0.5, col = "cornflowerblue") +
  xlab(NULL) +
  ylab("Litter size") +
  scale_x_discrete(labels = c("TRUE" = "Same focus area", "FALSE" = "Different focus areas")) +
  theme_classic()

# ttest with groups: same focus area and different focus area
model <- lm(offspring_produced ~ matching_focus_area, filter(location_success_data, offspring_produced > 0))
summary(model) # p = 0.95
par(mfrow = c(2,2)) # create 4 spaces for plots
plot(model) # plot assumptions


# proportion weaned vs matching_focus_area 
par(mfrow=c(1,1)) # reset to have 1 space for plot
  ggplot(location_success_data, aes(x = matching_focus_area, y = prop_weaned_litter)) +
  geom_boxplot() +
  geom_jitter(size = 0.5, col = "cornflowerblue") +
  xlab(NULL) +
  ylab("Proportion of kits weaned") +
  scale_x_discrete(labels = c("TRUE" = "Same focus area", "FALSE" = "Different focus areas")) +
  theme_classic()

# ttest with groups: same focus area and different focus area
model <- lm(prop_weaned_litter ~ matching_focus_area, location_success_data)
summary(model) # p = 0.1868
par(mfrow = c(2,2))# create 4 spaces for plots
plot(model)# plot assumptions

location_success_data %>% 
  filter(!is.na(prop_weaned_litter)) %>% 
  group_by(matching_focus_area) %>% 
  summarize(mean = mean(prop_weaned_litter))

# successful_breeding_attempt vs matching_focus_area 
# boxplot
par(mfrow=c(1,1)) # reset to have 1 space for plot
ggplot(location_success_data, aes(x = interaction(matching_focus_area, successful_breeding_attempt))) +
  geom_bar() +
  theme_classic() 

# ttest with groups: same focus area and different focus area
model <- lm(successful_breeding_attempt ~ matching_focus_area, location_success_data)
summary(model) # p = 0.2996
par(mfrow = c(2,2)) # create 4 spaces for plots
plot(model) # plot assumptions







# compare PATCH to measures of reproductive success -----------------

# litter size vs matching patch
filter(location_success_data, offspring_produced > 0) %>% 
ggplot(., aes(x = matching_patch, y = offspring_produced)) +
  geom_boxplot() +
  geom_jitter(size = 0.5, col = "cornflowerblue") +
  xlab(NULL) +
  ylab("Litter size") +
  theme_classic() 

# ttest with groups: close (<5 km) and far
model <- lm(offspring_produced ~ matching_patch, filter(location_success_data, offspring_produced > 0))
summary(model) # p = 0.48
par(mfrow = c(2,2)) # create 4 spaces for plots
plot(model) # plot assumptions


# proportion weaned vs matching patch 
par(mfrow=c(1,1)) # reset to have 1 space for plot
ggplot(location_success_data, aes(x = matching_patch, y = prop_weaned_litter)) +
  geom_boxplot() +
  geom_jitter(size = 0.5, col = "cornflowerblue") +
  xlab(NULL) +
  ylab("Proportion of kits weaned") +
  theme_classic()

# ttest with groups: same focus area and different focus area
model <- lm(prop_weaned_litter ~ matching_patch, location_success_data)
summary(model) # p = 0.76
par(mfrow = c(2,2))# create 4 spaces for plots
plot(model)# plot assumptions


# successful_breeding_attempt vs matching patch
# boxplot
par(mfrow=c(1,1)) # reset to have 1 space for plot
ggplot(location_success_data, aes(x = interaction(matching_patch, successful_breeding_attempt))) +
  geom_bar() +
  theme_classic() 

# ttest with groups: same focus area and different focus area
model <- lm(successful_breeding_attempt ~ matching_patch, location_success_data)
summary(model) # p = 0.17
par(mfrow = c(2,2)) # create 4 spaces for plots
plot(model) # plot assumptions








# compare DISTANCE between trap locations to measures of reproductive success -----------------

# litter size vs distance 
par(mfrow=c(1,1)) # reset to have 1 space for plot
dist_litter <- 
  filter(location_success_data, offspring_produced > 0) %>% 
  ggplot(., aes(x = dist, y = offspring_produced)) +
  geom_point(size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = FALSE, col = "blue") +
  xlab("Distance between trap sites (km)") +
  ylab("Litter size") +
  scale_y_continuous(limits = c(0,8)) +
  theme(axis.text = element_blank()) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10)
    ) +
  annotate("text", x = 250, y = 7.5, label = "p = 0.98", col = "black", size = 5)

# save plot object to transfer to Q1B script 
saveRDS(dist_litter, file = "plots/URC/dist_litter.rds")
#ggsave("plots/litter_size_vs_distance.png")

# linear regression model
model <- glm(offspring_produced ~ dist, filter(location_success_data, offspring_produced > 0), family = "poisson")
summary(model) # p = 0.98
par(mfrow = c(2,2)) # create 4 spaces for plots
plot(model) # plot assumptions

# poisson regression model 
model <- glm(offspring_produced ~ dist, data = filter(location_success_data, offspring_produced > 0), family = "poisson")
summary(model)


# proportion weaned vs distance
par(mfrow=c(1,1)) # reset to have 1 space for plot
dist_wean <-
  ggplot(location_success_data, aes(x = dist, y = prop_weaned_litter)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Distance between trap sites (km)") +
  ylab("Proportion of kits weaned per litter") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text = element_text(size = 10))+
  annotate("text", x = 250, y = 1, label = "p = 0.6", col = "black", size = 5)
saveRDS(dist_wean, file = "plots/URC/dist_wean.rds")

# linear regression model
model <- lm(prop_weaned_litter ~ dist, location_success_data)
summary(model) # p = 0.6052
par(mfrow = c(2,2))# create 4 spaces for plots
plot(model)# plot assumptions



# logistic regression (trials = offspring_produced, successes = offspring weaned)
model <- glm(cbind(offspring_weaned, offspring_produced) ~ dist, family = binomial, data = location_success_data)
summary(model) # p = 0.88

# plotting the logistic regression model 
location_success_data$dist_wean_predicted <- predict(model, newdata = location_success_data, type = "response")
dist_wean <- 
  ggplot(location_success_data, aes(x = dist, y = prop_weaned_litter))+
  geom_point()+
  geom_smooth(data = location_success_data, aes(y = dist_wean_predicted), col = "blue",  size = 1) +
  xlab("Distance between trap sites (km)") +
  ylab("Kit weaning rate") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),,
    axis.text = element_text(size = 10)) +
  annotate("text", x = 250, y = 0.9, label = "p = 0.88", col = "black", size = 5)
saveRDS(dist_wean, file = "plots/URC/dist_wean.rds")




# successful_breeding_attempts vs distance 
# boxplot
par(mfrow=c(1,1)) # reset to have 1 space for plot
#dist_birth <- 
  ggplot(location_success_data, aes(x = dist, y = successful_breeding_attempt)) +
  geom_boxplot() +
  #geom_jitter(size = 0.5, col = "cornflowerblue") +
  ylab(NULL) +
  xlab("Distance between trap sites (km)") +
  scale_y_discrete(labels = c("FALSE" = "No birth", "TRUE" = "Birth")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 14)) +
  annotate("text", x = 250, y = 2.5, label = "p = 0.019", col = "maroon", size = 5)


# ttest with groups: same focus area and different focus area
model <- glm(birth_0_1 ~ dist, location_success_data, family = "binomial")
summary(model) # p = 0.019
par(mfrow = c(2,2)) # create 4 spaces for plots
plot(model) # plot assumptions


#logistic regression for birth vs distance
par(mfrow=c(1,1))
dist_birth <- 
  ggplot(location_success_data, aes(x = dist, y = birth_0_1)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, col = "blue") +
  theme_classic() +
  ylab("Probability of birth") +
  xlab("Distance between trap sites (km)") +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10)
    ) +
  annotate("text", x = 250, y = 0.9, label = "p = 0.018", col = "red", size = 5)
saveRDS(dist_birth, file = "plots/URC/dist_birth.rds")

model <- glm(birth_0_1 ~ dist, data = location_success_data)
summary(model) # p = 0.0178
plot(model)

#probability of birth at 365.6 km
plogis(0.4768041+365.6*(0.0006614))-0.4768041






plot_grid(birth_dist_log, dist_litter, dist_wean, nrow = 1, ncol = 3)
ggsave("plots/nrem/distance.png", width = 12, height = 4)








# compare Distance under 5km to measures of reproductive success -----------------
location_success_data <- mutate(location_success_data, dist_under_5 = (dist <= 5))

# litter size vs dist_under_5 
filter(location_success_data, offspring_produced > 1) %>% 
  ggplot(., aes(x = dist_under_5, y = offspring_produced)) +
  geom_boxplot() +
  geom_jitter(size = 0.5, col = "cornflowerblue") +
  xlab(NULL) +
  ylab("Litter size") +
  scale_x_discrete(labels = c("TRUE" = "dist < 5", "FALSE" = "dist > 5")) +
  theme_classic()

# ttest 
model <- lm(offspring_produced ~ dist_under_5, filter(location_success_data, offspring_produced > 0))
summary(model) # p = 0.94
par(mfrow = c(2,2)) # create 4 spaces for plots
plot(model) # plot assumptions


# proportion weaned vs dist_under_5
par(mfrow=c(1,1)) # reset to have 1 space for plot
ggplot(location_success_data, aes(x = dist_under_5, y = prop_weaned_litter)) +
  geom_boxplot() +
  geom_jitter(size = 0.5,  col = "cornflowerblue") +
  xlab(NULL) +
  ylab("Proportion of kits weaned") +
  scale_x_discrete(labels = c("TRUE" = "dist < 5", "FALSE" = "dist > 5")) +
  theme_classic()

# ttest 
model <- lm(prop_weaned_litter ~ dist_under_5, location_success_data)
summary(model) # p = 0.13
par(mfrow = c(2,2))# create 4 spaces for plots
plot(model)# plot assumptions


# successful_breeding_attempt vs dist_under_5 
# boxplot
par(mfrow=c(1,1)) # reset to have 1 space for plot
ggplot(location_success_data, aes(x = interaction(dist_under_5, successful_breeding_attempt))) +
  geom_bar() +
  theme_classic() 

# ttest with groups: same focus area and different focus area
model <- lm(successful_breeding_attempt ~ dist_under_5, location_success_data)
summary(model) # p = 0.3112
par(mfrow = c(2,2)) # create 4 spaces for plots
plot(model) # plot assumptions







# save --------------------------------------------------------------------
write_csv(location_success_data, "data/location_success_data.csv")



# summary stats  ----------------------------------------------------------
location_success_data %>% 
  #group_by(birth_0_1) %>% 
  summarize(mean_dist = mean(dist), sd = sd(dist))


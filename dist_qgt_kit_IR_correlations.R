# Zoo microsatellite project script - models for kit inbreeding ~ fitness - Question 1A
# Olivia Tatro
# last modified 4/12/25

# goal: create lms and glms for successfully_weaned vs each inbreeding metric

# set wd, load libraries, load data ---------------------------------------

setwd("/Users/oliviatatro/Downloads/NEC inbreeding depression thesis")

library(tidyverse)

# has pairs and distance
pair_dist <- read.csv("data/location_success_data.csv")
# has pairs and relatedness
pair_relatedness <- read.csv("data/pair_success.csv")
# offspring metadata
offspring <- read.csv("data/metadata_offspring_updated.csv")

pair <- merge(pair_dist, pair_relatedness, by = c("SireID", "DamID"))

filtered <- filter(pair, QGT <= 0.2)


# correlation between distance and relatedness
ggplot(pair, aes(x = dist, y = QGT))+
  geom_point() +
  geom_smooth(method = "lm", col = "blue", se = FALSE) +
  xlab("Distance between trap sites (km)")+
  ylab("Pair Queller Goodnight Relatedness")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)) +
  annotate("text", x = 275, y = 0.6, label = "p = 0.0001", col = "red", size = 5)



ggsave("plots/URC/dist_qgt.png", width = 5, height = 5)

model <- lm(QGT ~ dist, data = pair)
summary(model)
plot(model)







# correlation between parent relatedness and kit IR 

kits_parents <- merge(pair_relatedness, offspring, by = c("SireID", "DamID"))
ggplot(kits_parents, aes(x = QGT, y = IR))+
  geom_point() +
  xlab("Parent Queller Goodnight Relatedness")+
  ylab("Kit Internal Relatedness") +
  theme_bw() +
  geom_smooth(method = "lm")
ggsave("plots/kitIR_parentQGT.png")





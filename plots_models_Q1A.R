# Zoo microsatellite project script - compare kit inbreeding and fitness - Question 1A
# Olivia Tatro
# last modified 3/11/25

# we want to compare measures of inbreeding of kits with their success
  # also compare this with parent success

# load libraries, set working directory, load data ------------------------
setwd("/Users/oliviatatro/Downloads/NEC inbreeding depression thesis")
getwd()

library(tidyverse)

# inbreeding estimates from coancestry
coancestry <- read.delim("data/InbreedingEstimates.txt", sep = "", header = TRUE)
# offspring metadata file with kit success and parent IDs
offspring <- read.csv("data/metadata_offspring_updated.csv")
# inbreeding estimates from genhet
genhet <- read.csv("genhet/genhet_output.csv")

# merge offspring metadata and inbreeding estimates by ID--------------------
# since I saved the df from this script, we don't need to do this again, but I am leaving the code in 
# offspring_coancestry <- merge(offspring, coancestry, by = "ID")
# instead do this so that the rest of the code works:
offspring_coancestry <- offspring


#pdf("Q1A_plots.pdf",height=8,width=11)

# Weaning success vs lynch ritland inbreeding ---------------------------------------
offspring_coancestry %>% 
  filter(!is.na(Successfully.Weaned.) & "-" != Successfully.Weaned.) %>% 
  ggplot(., aes(x = `Successfully.Weaned.`, y = LynchRt)) +
  geom_boxplot() +
  ylab("Lynch-Ritland Inbreeding Coefficient") +
  #geom_jitter(size = 0.75) +
  scale_x_discrete(labels = c("Yes" = "Weaned", "No" = "Not weaned")) + #changes boxplot labels
  xlab(NULL) +  # removes x-axis label 
  theme_light()
  

# Died_day_0 vs lynch ritland inbreeding ----------------------------------
# compare kit inbreeding and if they died at day 0 
offspring_coancestry <- offspring_coancestry %>% 
  mutate(died_day_0 = Age.of.Death == 0)

# There are only a few points for kits that died on day 0, so this metric doesn't make sense to use
temp <- offspring_coancestry %>% 
  filter(died_day_0)
nrow(temp) # 6 

offspring_coancestry %>% 
  filter(!is.na(died_day_0)) %>% 
  ggplot(., aes(y = LynchRt, x = died_day_0)) +
  geom_boxplot() +
  ylab("Kit inbreeding") +
  geom_jitter()


# genhet inbreeding comparisons -----------------------------------
# genhet variables:
  # PHt = proportion of heterozygous loci 
    # #heterozygous loci / #genotypes loci
  # Hs_obs = standardized heterozygosity based on the mean observe heterozygosity
    # PHt / mean observed heterozygosity of typed loci
  # Hs_exp = tandardized heterozygosity based on the mean expected heterozygosity
    # PHt / mean observed heterozygosity of typed loci
  # IR = internal relatedness
    # sharing of rare alleles counts more that sharing common alleles
    # can overestimate homozygosity for individuals with rare alleles
  # HL = homozygosity by locus 
    # considers the contribution of each locus instead of each allele
    # weight given to each locus depends on its variability 

# add all the genhet measures to the offspring metadata
genhet <- rename(genhet, ID = sampleid)
# don't need the following line because I saved the data sheet from this script 
#offspring_genhet <- merge(offspring, genhet, by = "ID")
# do this instead so that the rest of the code works
offspring_genhet <- offspring

# compare PHt and if the kit successfully weaned
offspring_genhet %>% 
  filter(!is.na(Successfully.Weaned.) & "-" != Successfully.Weaned.) %>% 
  ggplot(., aes(x = Successfully.Weaned., y = PHt)) +
  geom_boxplot() +
  #geom_jitter(size = 0.75) +
  scale_x_discrete(labels = c("Yes" = "Weaned", "No" = "Not weaned")) + #changes boxplot labels
  xlab(NULL) +  # removes x-axis label 
  theme_light() 

# compare Hs_obs and if the kit successfully weaned
offspring_genhet %>% 
  filter(!is.na(Successfully.Weaned.) & "-" != Successfully.Weaned.) %>% 
  ggplot(., aes(x = Successfully.Weaned., y = Hs_obs)) +
  geom_boxplot() +
  #geom_jitter(size = 0.75) +
  scale_x_discrete(labels = c("Yes" = "Weaned", "No" = "Not weaned")) + #changes boxplot labels
  xlab(NULL) +  # removes x-axis label 
  theme_light()

# compare Hs_exp and if the kit successfully weaned
offspring_genhet %>% 
  filter(!is.na(Successfully.Weaned.) & "-" != Successfully.Weaned.) %>% 
  ggplot(., aes(x = Successfully.Weaned., y = Hs_exp)) +
  geom_boxplot() +
  #geom_jitter(size = 0.75) +
  scale_x_discrete(labels = c("Yes" = "Weaned", "No" = "Not weaned")) + #changes boxplot labels
  xlab(NULL) +  # removes x-axis label 
  theme_light()


# compare IR and if the kit successfully weaned
offspring_genhet %>% 
  filter(!is.na(Successfully.Weaned.) & "-" != Successfully.Weaned.) %>% 
  ggplot(., aes(x = Successfully.Weaned., y = IR)) +
  geom_boxplot() +
  ylab("Internal Relatedness") +
  #geom_jitter(size = 0.75) +
  scale_x_discrete(labels = c("Yes" = "Weaned", "No" = "Unweaned")) + #changes boxplot labels 
  xlab(NULL) +  # removes x-axis label 
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),   
    axis.text.y = element_text(size = 10)
  )
# logistic regression for weaning vs IR 
offspring_genhet <- filter(offspring_genhet, !is.na(Successfully.Weaned.) & "-" != Successfully.Weaned.) %>% 
  mutate(weaned_01 = ifelse(Successfully.Weaned. == "Yes", 1, 0))
ggplot(offspring_genhet, aes(x = IR, y = weaned_01)) +
  geom_jitter(height = 0.015, size = 0.5) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, col = "blue") +
  xlab("Kit Internal Relatedness")+
  ylab("Probability of weaning") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 10)) +
  annotate("text", x = 0.5, y = 0.9, label = "p = 0.45", col = "black", size = 5)

model <- glm(weaned_01 ~ IR, offspring_genhet, family = "binomial")
summary(model)
ggsave("plots/URC/kitIR_wean.png", width = 3, height = 3)

filter(offspring_genhet, !is.na(IR)) %>% 
  group_by(Successfully.Weaned.) %>% 
  summarise(meanIR = mean(IR),
            sdIR = sd(IR),
            minIR = min(IR), 
            maxIR = max(IR))

#ggsave("plots/IR_vs_weaning.png", width = 5, height = 5)

# compare HL (homozygosity by locus) and if the kit successfully weaned
offspring_genhet %>% 
  filter(!is.na(Successfully.Weaned.) & "-" != Successfully.Weaned.) %>% 
  ggplot(., aes(x = Successfully.Weaned., y = HL)) +
  geom_boxplot() +
  ylab("Homozygosity by locus") +
  #geom_jitter(size = 0.75) +
  scale_x_discrete(labels = c("Yes" = "Weaned", "No" = "Not weaned")) + #changes boxplot labels
  xlab(NULL) +  # removes x-axis label 
  theme_light()


# save --------------------------------------------------------------------
#offspring_metadata_updated <- merge(offspring_coancestry, offspring_genhet)
#offspring_metadata_updated <- merge(offspring_metadata_updated, offspring, all = TRUE)

#offspring_metadata_updated <- offspring_metadata_updated %>% 
  #filter(!is.na(ID)) %>% 
  #select(1:15, 17:26, 29, 32, 35:39) #filters out extra accession ID and extra columns from coancestry and genhet

#write_csv(offspring_metadata_updated, "data/metadata_offspring_updated.csv")

#dev.off()


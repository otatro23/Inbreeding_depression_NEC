# Zoo microsatellite project script - calculate how many kits were released each year
# Olivia Tatro
# last modified 4/12/25

# load libraries, set working directory, load data ------------------------
setwd("/Users/oliviatatro/Downloads/NEC inbreeding depression thesis")
getwd()

library(tidyverse)

# offspring metadata file with kit IDs, year, and weaning status
offspring_metadata <- read.csv("data/metadata_offspring_updated.csv")


# clean offspring df -----------------------------------------------------------

# keeps kit ID and weaning status and finds the year of birth from DOB
offspring <- offspring_metadata %>% 
  select(c(1,8)) %>% 
  mutate(date = as.Date(Date.of.Release, format = "%m/%d/%y")) %>% #converts DOB character to a date
  mutate(release_year = format(date, "%Y")) %>% # takes just the year out of the date column
  select(1,4) %>% # keeps kit ID, weaning status, and year of birth
  filter(!is.na(release_year) & release_year != 2002) 

ggplot(offspring, aes(x = release_year)) +
  geom_bar(fill = "cornflowerblue") +
  theme_classic() +
  ylab("Number of Kits Released") +
  xlab("Year") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 8), 
    axis.title.y = element_text(size = 16),   
    axis.text.y = element_text(size = 10)
  )
ggsave("plots/kits_released_per_year.png")

#release <- NA
#for(i in 2011:2024){
  #print(nrows(filter(offspring, release_year == i)))
#}
  




# create goal/kits born/kits released barplot for nrem-----------------------------
kits_nrem <- offspring_metadata %>% 
  select(c(ID, Date.of.Release))

born <- nrow(kits_nrem)
released <- nrow(filter(kits_nrem, !is.na(Date.of.Release)))

to_plot <- data.frame(
  kit_cat = c("Born", "Released"), 
  count = c(born/14, released/14)
)

ggplot(to_plot, aes(x = kit_cat, y = count)) +
  geom_bar(stat = "identity", fill = "cornflowerblue") +
  geom_hline(yintercept = 250, linetype="dashed", size =1, col = "maroon") +
  theme_classic() +
  ylab("Mean number of kits per year") +
  xlab(NULL) +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),   
    axis.text.y = element_text(size = 10)
  )

ggsave("plots/nrem/born_released.png", width=3, height=4)


genotyped <- data.frame(
  cat = c("Dams", "Sires", "Weaned Kits", "Unweaned Kits"), 
  count = c(27, 13, 117, 50)
)

ggplot(genotyped, aes(x = cat, y = count)) +
  geom_bar(stat = "identity", fill = "cornflowerblue")+
  geom_text(aes(label = count), vjust = -0.5, size = 4)+
  theme_classic() +
  xlab(NULL) +
  ylab("Genotyped samples") +
  scale_y_continuous(limits = c(0, 125)) +
  theme(
  axis.text.x = element_text(angle = 30,size = 14, hjust = 1), 
  axis.title.y = element_text(size = 16),   
  axis.text.y = element_text(size = 10)
)

ggsave("plots/nrem/genotype_counts.png", width = 4, height = 4)  
  


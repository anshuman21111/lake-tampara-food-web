
library(tidyverse)
library(ggpubr)

### Read in the data ####

species_metad <- read.csv2("data/species_data_frame.csv")

species_level <- read.csv("data/processed/species_metrics.csv")%>%
  mutate(web = factor(web,
                      levels = c("Premonsoon","Monsoon","Postmonsoon"))) %>%
  left_join(species_metad, by = join_by("species"=="SPECIES"))

network_level <- read.csv("data/processed/network_metrics.csv")[c(1,3,2),] %>%
  mutate(links_per_species = links/size)

#### Visualize the data ####

#### By web #### 

# Size and connectance
ggplot(data = network_level, aes(x = web, group = 1)) +
  geom_line(aes(y = size)) +
  geom_line(aes(y = C*1000), linetype = 2) +
  scale_y_continuous(
    #First axis 
    name = "Number of species (solid line)",
    # Second axis
    sec.axis = sec_axis(~./1000, name = "Connectance (dashed line)")
  ) +
    theme_pubr()

# Size and links per species
ggplot(data = network_level, aes(x = web, group = 1)) +
  geom_line(aes(y = size)) +
  geom_line(aes(y = links_per_species*7), linetype = 2) +
  scale_y_continuous(
    #First axis 
    name = "Number of species (solid line)",
    # Second axis
    sec.axis = sec_axis(~./7, name = "Links per species (dashed line)")
  ) +
  theme_pubr()

#### By species ####

# Trophic level by season - looks like it increases postmonsoon

ggdensity(data = species_level,
          x = "TL",
          fill = "web", 
          facet.by = "guild") 

ggdensity(data = species_level,
          x = "degree",
          fill = "web",
          facet.by = "guild")

ggdensity(data = species_level,
          x = "omnivory",
          fill = "web",
          facet.by = "guild")

#### Change by species across seasons ####
species_level <- species_level %>%
  group_by(species) %>%
  mutate(TL_change = TL - mean(TL),
         degree_change = degree - mean(degree),
         omniv_change = omnivory - mean(omnivory))


ggplot(data = species_level) + 
  geom_line(aes(x = web, y = TL_change, group = species)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  labs(title = "Change in trophic level across seasons") +
  theme_pubr()

ggplot(data = species_level) + 
  geom_line(aes(x = web, y = degree_change, group = species)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  labs(title = "Change in degree across seasons") +
  theme_pubr()


ggplot(data = species_level) + 
  geom_line(aes(x = web, y = omniv_change, group = species)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  labs(title = "Change in omnivory across seasons") +
  theme_pubr()

#### By guild ####

guild_level <- species_level %>%
  pivot_longer(cols = c("TL","degree","omnivory"), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(guild = factor(guild)) %>%
  group_by(guild, web, metric) %>%
  summarize(median = median(value),
            low = quantile(value, 0.25),
            high = quantile(value, 0.75)) %>%
  mutate(xaxis = as.numeric(web)*100 + as.numeric(guild))

ggplot(data = guild_level %>% filter(metric == "TL")) + 
  geom_line(aes(x = xaxis, y = median, color = guild, group = guild)) + 
  geom_point(aes(x = xaxis, y = median, color = guild, group = guild)) + 
  geom_errorbar(aes(x = xaxis, ymin = low, ymax = high, color = guild, group = guild), width = 0.2) + 
  labs(title = "Median trophic level across seasons by guild") +
  scale_x_continuous(breaks=c(104,204,304),
                   labels=c("Premonsoon", "Monsoon", "Postmonsoon")) +
  theme_pubr()


ggplot(data = guild_level %>% filter(metric == "degree")) + 
  geom_line(aes(x = xaxis, y = median, color = guild, group = guild)) + 
  geom_point(aes(x = xaxis, y = median, color = guild, group = guild)) + 
  geom_errorbar(aes(x = xaxis, ymin = low, ymax = high, color = guild, group = guild), width = 0.2) + 
  labs(title = "Median trophic level across seasons by guild") +
  scale_x_continuous(breaks=c(104,204,304),
                     labels=c("Premonsoon", "Monsoon", "Postmonsoon")) +
  theme_pubr()

ggplot(data = guild_level %>% filter(metric == "omnivory")) + 
  geom_line(aes(x = xaxis, y = median, color = guild, group = guild)) + 
  geom_point(aes(x = xaxis, y = median, color = guild, group = guild)) + 
  geom_errorbar(aes(x = xaxis, ymin = low, ymax = high, color = guild, group = guild), width = 0.2) + 
  labs(title = "Median trophic level across seasons by guild") +
  scale_x_continuous(breaks=c(104,204,304),
                     labels=c("Premonsoon", "Monsoon", "Postmonsoon")) +
  theme_pubr()

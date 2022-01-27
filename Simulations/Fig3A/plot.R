library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)
library(dplyr)

# set theme
my_theme <- theme_light(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(size = 20)) + 
  theme(rect = element_rect(fill = "transparent")) +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))

# read data
path <- "regulationFactor.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")

# scale data
coeff <- 10/0.1
data$`regulation factor` <- 100 * data$`regulation factor`
data$`damage at end` <- coeff * data$`damage at end`

# prepare data
meltedData = melt(data, id = c("regulation factor", "phase", "damage formation", 
                               "damage repair","status"))
meltedData %>%
  group_by(`regulation factor`, phase, status, variable) %>%
  summarize(mean = mean(value, na.rm = TRUE),
            q5 = quantile(value, 0.05, na.rm = TRUE),
            q95 = quantile(value, 0.95, na.rm = TRUE)) -> groupedData

# plot (rls, time and damage in phases) 
fig = ggplot() +
  geom_ribbon(data = subset(groupedData, `regulation factor` >= 5.0),
              aes(x = `regulation factor`, ymin = q5, ymax = q95, 
                  fill = variable), alpha = 0.4) +
  geom_line(data = subset(groupedData, `regulation factor` >= 5.0), 
            aes(x = `regulation factor`,
                y = mean, color = variable), alpha = 0.5, size = 1) +
  geom_ribbon(data = subset(groupedData, `regulation factor` <= 5.0),
              aes(x = `regulation factor`, ymin = q5, ymax = q95, 
                  fill = variable), alpha = 0.6) +
  geom_line(data = subset(groupedData, `regulation factor` <= 5.0), 
            aes(x = `regulation factor`,
                y = mean, color= variable), size = 1) +
  geom_vline(xintercept = 5.0, color = "grey", size = 1) +
  facet_grid( ~ phase) +
  xlab(expression("regulation factor [10"^"-2"~"]")) +
  ylab("") +
  xlim(c(0, 8)) +
  scale_fill_manual(values = c("#4C996B", "#708FA3", "#252525"),
                    labels = c("time [h]", "divisions",
                               "fraction of damage at end")) +
  scale_color_manual(values = c("#4C996B", "#708FA3", "#252525"),
                     labels = c("time [h]", "divisions",
                                "fraction of damage at end")) +
  scale_y_continuous(name = "", 
                     sec.axis = sec_axis(~./coeff, name = "")) +
  my_theme
ggsave(file = "phases.svg", plot = fig, height = 5 , width = 10)
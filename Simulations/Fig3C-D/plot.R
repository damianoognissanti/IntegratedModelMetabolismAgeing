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
path <- "ngamEffects.txt" 
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
data <- subset(data, status != "alive")

# oxygen is taken up, so change sign
data$`norm. O2 uptake (end)` <- -1 * data$`norm. O2 uptake (end)` 

# get damage at end of life from phase 3 (unless it is not reached)
subset1 <- subset(data, phase == "phase 3")
for (i in 1 : nrow(subset1)) {
  if (subset1$time[i] == 0.0) {
    tmpData <- data[data$NGAM == subset1$NGAM[i] &
                      data$`damage formation` == subset1$`damage formation`[i] &
                      data$`damage repair` == subset1$`damage repair`[i] &
                      data$phase == "phase 2", ] 
    subset1$damage[i] <- tmpData$damage 
  }
}

# prepare data
subset1 <- subset1[c(1, 5)]
subset1 %>%
  group_by(`NGAM`) %>%
  summarize(mean = mean(damage, na.rm = TRUE),
            q5 = quantile(damage, 0.05, na.rm = TRUE),
            q95 = quantile(damage, 0.95, na.rm = TRUE)) -> groupedData

# plot (damage at end of life)
fig = ggplot(groupedData,  aes(x = NGAM)) +
  my_theme +
  geom_hline(yintercept = 0.46, size = 1, color= "grey") +
  geom_ribbon(aes(x = `NGAM`, ymin = q5, ymax = q95), 
              fill = "#252525", alpha = 0.4) +
  geom_line(aes(x = `NGAM`, y = mean), color = "#252525", 
            alpha = 0.5, size = 1) +
  xlab(expression("NGAM"["max"]~"[mmol (gDW h)"^"-1"~"]")) +
  ylab("fraction of damaged proteins D at end of life")
ggsave(file = "damage.svg", plot = fig)
  
 # prepare data
meltedData <- melt(data[c(1, 2, 7:10, 13)],
                  id = c("NGAM", "phase", "status"))
meltedData <- subset(meltedData, phase !="phase 3")
meltedData %>%
  group_by(`NGAM`, phase, status, variable) %>%
  summarize(mean = mean(value, na.rm = TRUE),
            q5 = quantile(value, 0.05, na.rm = TRUE),
            q95 = quantile(value, 0.95, na.rm = TRUE)) -> groupedData

# plot (exchange fluxes in beginning and end of phase 2) 
fig = ggplot() +
  my_theme +
  geom_ribbon(data = subset(groupedData, `NGAM` >= 5.0),
              aes(x = `NGAM`, ymin = q5, ymax = q95, 
                  fill = variable), alpha = 0.4) +
  geom_line(data = subset(groupedData, `NGAM` >= 5.0), 
            aes(x = `NGAM`,
                y = mean, color = variable), alpha = 0.5, size = 1) +
  geom_ribbon(data = subset(groupedData, `NGAM` <= 5.0),
              aes(x = `NGAM`, ymin = q5, ymax = q95, 
                  fill = variable), alpha = 0.6) +
  geom_line(data = subset(groupedData, `NGAM` <= 5.0), 
            aes(x = `NGAM`,
                y = mean, color= variable), size  = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", size = 1, color = "black") + 
  facet_grid( ~ phase) +
  xlab(expression("NGAM"["max"]~"[mmol (gDW h)"^"-1"~"]")) +
  ylab("normalised production/uptake") +
  xlim(c(0, 1.0)) +
  scale_fill_manual(values = c("#AA4E39", "#552F00", "#D4A66A", "#AA7839"), 
                    labels = c("ethanol", "O2", "CO2", "acetate")) +
  scale_color_manual(values = c("#AA4E39", "#552F00", "#D4A66A", "#AA7839"), 
                    labels = c("ethanol", "O2", "CO2", "acetate")) +
  theme(aspect.ratio = 1)
ggsave(file = "exchange.svg", plot = fig, height = 5 , width = 10)

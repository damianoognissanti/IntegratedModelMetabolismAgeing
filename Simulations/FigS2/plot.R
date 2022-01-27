library(ggplot2)
library(readr)
library(reshape2)
library(tibble)
library(svglite)

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
path <- "exchangeFluxes.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
colNames <- names(data)

# read experimental data
pathExperiment <- "../../ExperimentalData/chemostatExperiment.txt"
experiment <- read_delim(pathExperiment, "\t", col_names = TRUE, comment = "#")
experiment <- experiment[c(1, 3:7, 9)]

# define colours 
colours <- c("#297A4A", "#29526D", "#708FA3", "#AA4E39", "#AA7839", "#252525",
            "#7AB793", "#636363", "#D47F6A", "#D4A66A")

# take only the one used here
coloursHere <- colours[c(2, 1, 5, 4, 10, 9)]

# plot (exchange fluxes)
meltedData <- melt(data[c(1:7)], id = "growth rate")
meltedExperiment <- melt(experiment, id = "D (1/h)")
fig = ggplot() +
  geom_point(data = meltedExperiment, aes(x = `D (1/h)`, y = value, 
                                          color = variable), size = 4) +
  geom_line(data = meltedData, aes(x = `growth rate`, y = value, 
                                   color = variable), size = 1.5) +
  scale_color_manual(breaks = colNames[2:7], labels = colNames[2:7], 
                     values = c(coloursHere, coloursHere)) +
  xlab(expression("growth rate [h"^"-1"~"]")) +
  ylab(expression("fluxes [mmol (gDW h)"^"-1"~"]")) +
  ylim(0, 20) +
  guides(colour = guide_legend(nrow = 3)) +
  my_theme
ggsave(file = "exchangeFluxes.svg", plot = fig)

# read damage data
path <- "damage.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
colNames <- names(data)

# take only the one used here
coloursHere <- colours[c(1:8)]

# plot (damage)
meltedData <- melt(data[c(1:9)], id = "growth rate")
fig <- ggplot() +
  geom_line(data = meltedData, aes(x = `growth rate`, y = value, 
                                   color = variable), size = 1.5) +
  scale_color_manual(breaks = colNames[2:9], labels = colNames[2:9], 
                     values = c(coloursHere)) + 
  xlab(expression("growth rate [h"^"-1"~"]")) +
  ylab(expression("fluxes [mmol (gDW h)"^"-1"~", h"^"-1"~"]")) +
  guides(colour = guide_legend(nrow = 5)) +
  my_theme
ggsave(file = "damage.svg", plot = fig)

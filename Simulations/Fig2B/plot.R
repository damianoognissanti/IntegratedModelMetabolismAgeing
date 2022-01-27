library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)

# set our theme
my_theme <- theme_light(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(size = 20)) + 
  theme(rect = element_rect(fill = "transparent")) +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))

# read data
path <- "life.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
colNames <- names(data)

# multiply mass by general dry mass of cell
data$M <- data$M * 2.9e-11

# plot (growth) 
fig <- ggplot(data = data, aes(x = time, y = M)) +
  my_theme + 
  geom_line(colour = "#636363", size = 2) +
  xlab("time") +
  ylab("cell mass [gDW]") +
  ylim(0, 5e-11) +
  theme(legend.position = "none")
ggsave(file = "mass.svg", plot = fig)

# plot (boolean inputs) 
meltedData <- melt(data[, c("time", "glucose(bool)", "H2O2(bool)", "Trx1/2(bool)")], 
                   id = "time")
fig <- ggplot(meltedData, aes(x = time, y = variable, fill = value)) + 
  my_theme + 
  geom_tile() + 
  scale_fill_gradient(low = "grey", high = "#4C996B", 
                      limits = c(0, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 3/16)
ggsave(file = "booleanInputs.svg", plot = fig)

# plot (protein content) 
fig <- ggplot(data, aes(x = time)) +
  my_theme + 
  geom_line(aes(y = P), color = "#4C996B", size = 2) +
  geom_line(aes(y = D), color = "#252525", size = 2) +
  xlab("time") +
  ylab(expression("protein content [g (gDW)"^"-1"~"]"))
ggsave(file = "proteinContent.svg", plot = fig)

# plot (exchange fluxes normalized) 
fig <- ggplot(data, aes(x = time)) +
  my_theme + 
  geom_hline(yintercept = -1, size = 1.5, linetype = "dashed", 
             color = "#636363") +
  geom_line(aes(y = - O2 / glucose), color = "#552F00", size = 2) +
  geom_line(aes(y = CO2 / glucose), color = "#D4A66A", size = 2) +
  geom_line(aes(y = ethanol / glucose), color = "#AA4E39", 
            size = 2) +
  geom_line(aes(y = pyruvate / glucose), color = "#D47F6A", 
            size = 2) +
  geom_line(aes(y = acetate / glucose), color = "#AA7839", 
            size = 2) +
  xlab("time") +
  ylab(expression("normalised fluxes"))
ggsave(file = "exchangeFluxes_norm.svg", plot = fig)

# plot (growth/damage formation) 
fig <- ggplot(data, aes(x = time)) +
  my_theme + 
  geom_line(aes(y = growth), color = "#708FA3", size = 2) +
  geom_line(aes(y = enzymeUsage), color = "#4C996B", size = 2) +
  geom_line(aes(y = `H2O2 production` / glucose), color = "#AA7839", size = 2) +
  geom_line(aes(y = Dformation / glucose), color = "#252525", size = 2) +
  xlab("time") +
  ylab(expression("fluxes")) +
  guides(colour = guide_legend(nrow = 2))
ggsave(file = "outputFluxes.svg", plot = fig)
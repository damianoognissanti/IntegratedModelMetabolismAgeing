library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)
library(scales)

# set theme
my_theme <- theme_light(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(size = 20)) + 
  theme(rect = element_rect(fill = "transparent")) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))

# read data
path <- "f0r0Grid.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
data <- subset(data, (rls < 50) & (status != "alive"))

# plot (rls) 
fig <- ggplot(data, aes(x = `damage formation`, y = `damage repair`, 
                        fill = rls)) +
  geom_tile() +
  facet_grid(`regulation factor` ~ .) +
  xlab("non-metabolic damage formation") +
  ylab("damage repair") +
  scale_fill_gradient(low = "lightgrey", high = "#032236",
                      limits = c(12,32)) +
  my_theme
ggsave(file = "rls.svg", plot = fig)


# plot (generation time) 
fig <- ggplot(data, aes(x = `damage formation`,
                        y = `damage repair`,
                        fill = `average generation time`)) +
  geom_tile() +
  facet_grid(`regulation factor` ~ .) +
  xlab("non-metabolic damage formation") +
  ylab("damage repair") +
  scale_fill_gradient(low = "lightgrey", high = "#003D19") +
  my_theme
ggsave(file = "genTime.svg", plot = fig)

# plot (time phase 1) 
fig <- ggplot(data, aes(x = `rls`,
                        y = `time phase 1`)) +
  geom_point(colour = "#636363", size = 4) +
  facet_grid(`regulation factor` ~ .) +
  xlab("rls") +
  ylab("time in phase 1") +
  my_theme
ggsave(file = "timePhase1.svg", plot = fig)

# plot (damage at end of phase 1) 
fig <- ggplot(data, aes(x = `damage formation`,
                        y = `damage repair`,
                        fill = `damage phase 1`)) +
  geom_tile() +
  facet_grid(`regulation factor` ~ .) +
  xlab("non-metabolic damage formation") +
  ylab("damage repair") +
  scale_fill_gradient(low = "lightgrey", high = "black") +
  my_theme
ggsave(file = "damagePhase1.svg", plot = fig)
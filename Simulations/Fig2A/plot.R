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
data <- subset(data, `damage formation` > 0.0)
 
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
library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)
library(forcats)

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
path <- "enzymeDeletions.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
colNames <- names(data)

# plot (divisions in phases)
meltedData <- melt(data[c(4,5,10,13,16)], id = c("deletion span"))
meltedData <- subset(meltedData, !is.na(`deletion span`))
fig <- ggplot(meltedData, aes(x = variable, y = value,  
                              fill = `deletion span`,
                              alpha = variable)) +
  my_theme +
  geom_hline(yintercept = c(1, 6, 16, 23), color = "grey", size = 1) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = `deletion span`), size = 0.8) +
  facet_grid(`deletion span` ~.) +
  xlab("") +
  ylab("number of divisions") +
  ggtitle("deletion span") +
  scale_fill_manual(values = c("#29526D", "#4C996B", "#D4A66A", "#D47F6A"),
                    breaks = c("complete", "phase 1", "phase 2", "phase 3")) +
  scale_color_manual(values = c("#29526D", "#4C996B", "#D4A66A", "#D47F6A"),
                     breaks = c("complete", "phase 1", "phase 2", "phase 3")) +
  scale_alpha_discrete(range = c(0.8, 0.1)) +
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
  scale_x_discrete(limits = rev, labels = c("in III", "in II", "in I", "total")) +
  theme(legend.position = "none") +
  coord_flip()
ggsave(file = "rls.svg", plot = fig)

# plot (enzyme usages)
subset1 <- subset(data, `enzyme change` > 1)
fig <- ggplot(data, aes(x = `rls change`, y = `enzyme change`, 
                        color = `deletion span`)) +
  my_theme +
  geom_hline(yintercept = 1, size = 1.5, colour = "grey") +
  geom_vline(xintercept = 1, size = 1.5, colour = "grey") +
  geom_point(size = 4, alpha = 0.4) +
  # geom_text(data = subset1, aes(label = `enzyme(s)`), hjust = 0, nudge_x = 0.05) +
  scale_color_manual(values = c("#29526D", "#4C996B", "#D4A66A", "#D47F6A"), 
                     breaks = c("complete", "phase 1", "phase 2", "phase 3")) +
  xlab("relative change in replicative lifespan") +
  ylab("relative change in total enzyme usage")
ggsave(file = "rlsVsUsage.svg", plot = fig) 

# plot (enzymes that change rls)
subset1 <- subset(data, rls - 23 != 0) 
fig <- ggplot(subset1, aes(x = `enzyme(s)`,
                           y = `deletion span`, fill = `rls` - 23)) +
  my_theme +
  geom_tile() +
  ylab("deletion span") +
  xlab("") +
  facet_grid(. ~ pathway, scale = "free_x", space = "free_x") +
  scale_fill_gradientn(limits = c(-23, 5),
                       colors = c("grey", "#29526D", "white", "#003D19"), 
                       values = c(0.0, 18/28, 23/28, 1.0), name = "∆  rls") + 
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5)) +
  theme(legend.title = element_text()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(panel.grid.major = element_blank())
ggsave(file = "data.svg", plot = fig)
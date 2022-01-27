library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)
library(forcats)
library(stringr)

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
path <- "proteinKnockouts.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
colNames <- names(data)

# plot (rls and generation time change) 
coeff <- 20/2
order <- levels(as.factor(data$knockout))
medianGenTime = median(data$`average generation time`[data$knockout == "wildtype"] * coeff)
fig <- ggplot(data, aes(x = knockout)) +
  my_theme +
  geom_hline(yintercept = data$rls[data$knockout == "wildtype"], 
             color = "grey", size = 1) +
  geom_hline(yintercept = medianGenTime, color = "grey", size = 1) +
  geom_boxplot(aes(y = rls), outlier.shape = NA, fill = "#708FA3", alpha = 0.5) +
  geom_jitter(aes(y = rls), width = 0.1, size = 0.8, color = "#708FA3") +
  geom_boxplot(aes(y = `average generation time` * coeff), outlier.shape = NA, 
               fill = "#4C996B", alpha = 0.5) +
  geom_jitter(aes(y = `average generation time` * coeff), width = 0.1, 
              size = 0.8, color = "#4C996B") +
  xlab("knockout") +
  scale_x_discrete(limits = order[c(11, 1:10)]) +
  scale_y_continuous(name = "rls", 
                     sec.axis = sec_axis(~./coeff, 
                                         name = "average generation time")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 12)) +
  theme(legend.position = "none")
ggsave(file = "deletedProteins.svg", plot = fig)

# plot (rls in phases) 
meltedData = melt(data[c(1, 2, 4, 5, 6)], id = "knockout")
fig <- ggplot(meltedData, aes(x = knockout)) +
  my_theme +
  geom_boxplot(aes(y = value, color = variable, fill = variable), 
               outlier.shape = NA, alpha = 0.5, position = "identity") +
  geom_jitter(aes(y = value, color = variable, fill = variable), 
              width = 0.1, size = 0.8) +
  scale_fill_manual(values = c("#29526D", "#4C996B", "#D4A66A", "#D47F6A"),
                    labels = c("total", "phase 1", "phase 2", "phase 3"),
                    breaks = c("rls", "rls phase 1", "rls phase 2", "rls phase 3")) +
  scale_color_manual(values = c("#29526D", "#4C996B", "#D4A66A", "#D47F6A"),
                     labels = c("total", "phase 1", "phase 2", "phase 3"),
                     breaks = c("rls", "rls phase 1", "rls phase 2", "rls phase 3")) +
  xlab("knockout") +
  ylab("number of divisions") +
  scale_x_discrete(limits = order[c(11, 1:10)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 12))
ggsave(file = "phases.svg", plot = fig)

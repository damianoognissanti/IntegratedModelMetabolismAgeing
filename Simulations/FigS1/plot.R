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

# go through all result folders
paths <- list.files(path = getwd(), pattern = "glc") 
for (p in paths) {
  
  # read in the data
  data <- read_delim(paste(p, "/activity.txt", sep = ""), 
                     "\t", col_names = FALSE, comment = "#")
  n <- ncol(data)
  colNames <- names(data)
  
  # get pathways 
  pathways <- unique(data[2])
  pathways <- pathways[!is.na(pathways)]
  
  # prepare x labels
  iterationLabels <- colNames[3:n]
  for (j in 3:n) {
    iterationLabels[j-2] <- substring(colNames[j], 2)
  }
  
  # find the point when the switch happens
  # set switch to -2.5 to account that the first few columns are not
  # the actual iteration data
  switch <- -2.5
  for (j in 4:n) {
    comparison <- data[j] == data[j-1]
    if (all(comparison[!is.na(comparison)])) {
        switch <- switch + j
        break
    }
  }

  # heatmap
  for(i in 1:length(pathways)) {
    
    # take subset of data
    tmpData <- subset(data, data[2] == pathways[i])
    nProteins <- nrow(tmpData)
    
    # reformat to have three columns with coordinates of data points
    meltedData <- melt(tmpData, id.vars = colNames[1],
                       measure.vars = colNames[3:ncol(data)])
    meltedData <- subset(meltedData, meltedData[3] != "nothing")
    
    # convert String to Bool to Int in order to plot
    meltedData$value <- as.integer(as.logical(meltedData$value))
    
    # plot and save
    fig <- ggplot(meltedData, aes(x = variable, 
                                  y = factor(X1, levels = rev(levels(factor(X1)))))) + 
      geom_tile(aes(fill = value)) + 
      scale_fill_gradient(low = "grey", high = "#4C996B", limits = c(0, 1)) +
      geom_vline(xintercept = switch, colour = "#252525", size = 1.5) +
      my_theme +
      theme(legend.position="none") +
      xlab("iteration") +
      scale_x_discrete(labels = labels, expand = c(0.2, 0)) +
      ylab(NULL) +
      theme(aspect.ratio = nProteins/16)
    ggsave(file = paste(p, "/", pathways[i], ".svg", sep = ""), plot = fig)
  }
}



library(tidyverse)

samples = readr::read_csv("samples.csv") %>%
    rename(sample_no = 1)


p = ggplot(samples) +
    geom_line(aes(sample_no, lowrank), color = "blue") +
    geom_line(aes(sample_no, mra), color = "red") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size=20)) +
    ylab(expression(sigma["w"]^2)) +
    xlab("sample #")

ggsave(filename="samples.pdf", width=400, height=100, units="mm")


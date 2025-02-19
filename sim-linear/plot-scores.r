library(tidyverse)

NITER = 10
DATA_FILE = sprintf("sum_of_scores.csv")

crps = read_csv(DATA_FILE) %>%
    rename(time = 1) %>%
    mutate(mra = mra / exact, lr = lr / exact, exact = exact / exact)

p = ggplot( crps ) +
    geom_line(aes(time, mra), color = "red") +
    geom_line(aes(time, lr), color = "blue") +
    geom_line(aes(time, exact), color = "black") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size=20)) +
    ylab("CRPS ratio")

ggsave(filename="CRPS.pdf", width=400, height=100, units="mm")

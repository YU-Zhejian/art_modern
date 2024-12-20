library(dplyr)
library(ggplot2)

df <- readr::read_tsv("correlation.tsv")

p <- ggplot(df) +
    geom_line(aes(x=pos, y=meanqual, color=profile_2)) +
    scale_y_continuous(limits=c(0, 45)) +
    facet_grid(profile_1~.) +
    theme_bw()
ggsave("fig/meanqual.png", p)
#
# geom_abline(intercept=0.05, slope=0) +

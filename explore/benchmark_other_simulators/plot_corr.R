library(dplyr)
library(ggplot2)

df <- readr::read_tsv("correlation.tsv")

p <- ggplot(df) +
    geom_line(aes(x=pos, y=meanqual, color=profile_2)) +
    scale_y_continuous(limits=c(0, 45)) +
    facet_grid(profile_1~rlen, scales="free_x") +
    theme_minimal()
ggsave("fig/meanqual.pdf", p, width=8, height=4)
p <- ggplot(df) +
    geom_line(aes(x=pos, y=correlation, color=profile_2)) +
    facet_grid(profile_1~rlen, scales="free_x") +
    theme_minimal()
ggsave("fig/correlation.pdf", p, width=8, height=4)
#
# geom_abline(intercept=0.05, slope=0) +

library(dplyr)
library(ggplot2)

replacement_list <- list(
  "pirs" = "pIRS",
  "dwgsim" = "DWGSIM",
  "wgsim" = "wgsim",
  "art" = "Original ART",
  "art_modern" = "art_modern",
  "real_1" = "Read 1",
  "real_2" = "Read 2",
  "real" = "Input Data"
)

df <- readr::read_tsv("correlation.tsv") %>%
dplyr::mutate(
  read_id = sapply(profile_1, function(x)
    replacement_list[[x]]),
  SOFTWARE = sapply(profile_2, function(x)
    replacement_list[[x]])
) %>%
  dplyr::filter(SOFTWARE != "wgsim" & SOFTWARE != "DWGSIM")

p <- ggplot(df) +
    geom_line(aes(x=pos, y=meanqual, color=factor(
      SOFTWARE,
      levels = c("wgsim", "DWGSIM", "pIRS", "art_modern", "Original ART", "Input Data")
    ))) +
    scale_y_continuous(limits=c(10, 45)) +
    scale_color_discrete("Software") +
    facet_grid(read_id~rlen, scales="free_x") +
    theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("fig/meanqual.pdf", p, width=7, height=4)
p <- ggplot(df) +
    geom_line(aes(x=pos, y=correlation, color=factor(
      SOFTWARE,
      levels = c("wgsim", "DWGSIM", "pIRS", "art_modern", "Original ART")
    ))) +
    facet_grid(profile_1~rlen, scales="free_x") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("fig/correlation.pdf", p, width=8, height=4)
#
# geom_abline(intercept=0.05, slope=0) +

library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

df <- readr::read_tsv("time_complexity.tsv")

df_defaults <- df %>%
  dplyr::filter(name %in% c("BamReadOutput_l=4_t=1", "HeadlessBamReadOutput_l=4_t=1", "SamReadOutput_l=4_t=1", "FastqReadOutput", "PwaReadOutput", "FastaReadOutput", "DumbReadOutput", "EmptyImplicitLFIOReadOutput", "EmptyLFIOReadOutput"))

g <- ggplot(df_defaults) +
  geom_boxplot(aes(y=time_complexity, x=as.factor(threads))) +
  scale_y_continuous(limits=c(0, NA)) +
  scale_x_discrete("#parallel driver threads") +
  facet_wrap(~name) +
                                        theme_bw()
ggsave("boxplot_formats.pdf", g, width = 10, height = 10)

quit()

df_bam_levels <- df %>%
  dplyr::filter(name %in% grep("^BamReadOutput_l=._t=1", unique(name), value = TRUE))

g <- ggplot(df_bam_levels) +
  geom_boxplot(aes(y=time_complexity, x=as.factor(threads))) +
  scale_y_continuous(limits=c(0, NA)) +
  scale_x_discrete("#parallel driver threads") +
  facet_wrap(~name) +
                     theme_bw()
ggsave("boxplot_bam_levels.pdf", g, width = 10, height = 10)

df_bam_threads <- df %>%
  dplyr::filter(name %in% grep("^BamReadOutput_l=4_t=.", unique(name), value = TRUE))

g <- ggplot(df_bam_threads) +
  geom_boxplot(aes(y=time_complexity, x=as.factor(threads))) +
  scale_y_continuous(limits=c(0, NA)) +
  scale_x_discrete("#parallel driver threads") +
  facet_wrap(~name) +
  theme_bw()
ggsave("boxplot_bam_threads.pdf", g, width = 10, height = 10)


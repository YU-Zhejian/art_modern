library(ggplot2)
library(dplyr)
df <- readr::read_tsv("time.tsv")

p <- ggplot(df) +
    geom_boxplot(aes(y=TEST_CASE, x=WALL_CLOCK)) +
    theme_bw()
ggsave("fig/wall_clock.png", p)

p <- ggplot(df) +
    geom_boxplot(aes(y=TEST_CASE, x=SYSTEM+USER)) +
    theme_bw()
ggsave("fig/cpu_time.png", p)

p <- ggplot(df) +
    geom_boxplot(aes(y=TEST_CASE, x=RSS)) +
    theme_bw()
ggsave("fig/residential_mem.png", p)

p <- df %>% dplyr::select(TEST_CASE, MAJ_PG_F, MIN_PG_F, VOL_CTX_S, IV_CTX_S) %>%
    tidyr::pivot_long(names_from=TEST_CASE) %>%
    ggplot() +
    geom_boxplot(aes(y=TEST_CASE, x=RSS)) +
    theme_bw()
ggsave("fig/page_faults.png", p)


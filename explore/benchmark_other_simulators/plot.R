library(ggplot2)
library(dplyr)
df <- readr::read_tsv("time.tsv") %>%
    dplyr::mutate(CPU_TIME=SYSTEM+USER)

replacement_list <- list(
    "CPU_TIME" = "CPU Time",
    "WALL_CLOCK" = "Wall Clock Time",
    "RSS" = "Residential Memory",
    "MAJ_PG_F" = "Major Page Faults",
    "MIN_PG_F" = "Minor Page Faults",
    "VOL_CTX_S" = "Voluntary Context Switches",
    "IV_CTX_S" = "Involuntary Context Switches"
)

p <- df %>%
    dplyr::select(CPU_TIME, WALL_CLOCK, RSS, TEST_CASE) %>%
    dplyr::mutate(DATA=ifelse(stringr::str_count(TEST_CASE, "genome") != 0, "GENOME", "TRANSCRIPTOME"))%>%
    tidyr::pivot_longer(
        cols=c("CPU_TIME", "WALL_CLOCK", "RSS"),
        names_to="ASPECTS",
        values_to="VALUES"
    ) %>%
    dplyr::mutate(
        ASPECTS = sapply(ASPECTS,  function(x) replacement_list[[x]])
    )%>%
    ggplot() +
    geom_boxplot(aes(y=TEST_CASE, x=VALUES)) +
    facet_grid(DATA~ASPECTS, scales="free") +
    theme_bw()
ggsave("fig/time_memory.png", p, width=12, height=4)

p <- df %>%
    dplyr::select(MAJ_PG_F, MIN_PG_F, VOL_CTX_S, IV_CTX_S, TEST_CASE) %>%
    tidyr::pivot_longer(
      cols=c("MAJ_PG_F", "MIN_PG_F", "VOL_CTX_S", "IV_CTX_S"),
      names_to="PFCS_TYPE",
      values_to="PFCS"
    ) %>%
    dplyr::mutate(
        PFCS_TYPE = sapply(PFCS_TYPE,  function(x) replacement_list[[x]])
    )%>%
    ggplot() +
    geom_boxplot(aes(y=TEST_CASE, x=PFCS)) +
    facet_wrap(~PFCS_TYPE, scales="free") +
    theme_bw()
ggsave("fig/page_faults.png", p)

sessionInfo()

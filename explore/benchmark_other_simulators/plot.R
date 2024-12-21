library(ggplot2)
library(dplyr)
df <- readr::read_tsv("time.tsv") %>%
    dplyr::mutate(CPU_TIME = SYSTEM + USER) %>%
    dplyr::mutate(
        DATA = ifelse(stringr::str_count(TEST_CASE, "genome") != 0, "GENOME", "TRANSCRIPTOME"),
        SOFTWARE = stringr::str_extract(TEST_CASE, "^[^-]+"),
        RLEN = ifelse(stringr::str_count(TEST_CASE, "300") != 0, "300", "100")
    )

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
    dplyr::select(CPU_TIME, WALL_CLOCK, RSS, DATA, SOFTWARE, RLEN) %>%
    tidyr::pivot_longer(
        cols = c("CPU_TIME", "WALL_CLOCK", "RSS"),
        names_to = "ASPECTS",
        values_to = "VALUES"
    ) %>%
    dplyr::mutate(
        ASPECTS = sapply(ASPECTS, function(x) replacement_list[[x]])
    ) %>%
    ggplot() +
    geom_boxplot(aes(y = SOFTWARE, x = VALUES, fill=RLEN, color=RLEN)) +
    scale_x_continuous(trans="log10", labels=scales::label_number(scale_cut = scales::cut_si(""))) +
    facet_grid(DATA ~ ASPECTS, scales = "free") +
    theme_bw()
ggsave("fig/time_memory.pdf", p, width = 12, height = 4)

p <- df %>%
    dplyr::select(MAJ_PG_F, MIN_PG_F, VOL_CTX_S, IV_CTX_S, DATA, SOFTWARE, RLEN) %>%
    tidyr::pivot_longer(
        cols = c("MAJ_PG_F", "MIN_PG_F", "VOL_CTX_S", "IV_CTX_S"),
        names_to = "PFCS_TYPE",
        values_to = "PFCS"
    ) %>%
    dplyr::mutate(
        PFCS_TYPE = sapply(PFCS_TYPE, function(x) replacement_list[[x]])
    ) %>%
    ggplot() +
    geom_boxplot(aes(y = SOFTWARE, x = PFCS, fill=RLEN, color=RLEN)) +
    scale_x_continuous(trans="log10", labels=scales::label_number(scale_cut = scales::cut_si(""))) +
    facet_grid(DATA ~ PFCS_TYPE, scales = "free") +
    theme_bw()
ggsave("fig/page_faults.pdf", p, width = 12, height = 4)

sessionInfo()

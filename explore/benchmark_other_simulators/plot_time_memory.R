library(ggplot2)
library(dplyr)
df <- readr::read_tsv("time.tsv") %>%
  dplyr::mutate(CPU_TIME = SYSTEM + USER) %>%
  dplyr::mutate(
    DATA = ifelse(
      stringr::str_count(TEST_CASE, "genome") != 0,
      "GENOME",
      "TRANSCRIPTOME"
    ),
    SOFTWARE = stringr::str_extract(TEST_CASE, "^[^-]+"),
    RLEN = ifelse(stringr::str_count(TEST_CASE, "300") != 0, "300", "100")
  )

df_2 <- df %>%
  dplyr::filter(SOFTWARE %in% c("art_original", "art_modern")) %>%
  dplyr::select(WALL_CLOCK, SOFTWARE, RLEN, DATA) %>%
  tidyr::pivot_wider(names_from=SOFTWARE, values_from=WALL_CLOCK, values_fn =mean)


replacement_list <- list(
  "CPU_TIME" = "CPU Time (s)",
  "WALL_CLOCK" = "Clock Time (s)",
  "RSS" = "Residential Memory (MB)",
  "MAJ_PG_F" = "Major Page Faults",
  "MIN_PG_F" = "Minor Page Faults",
  "VOL_CTX_S" = "Voluntary Context Switches",
  "IV_CTX_S" = "Involuntary Context Switches",
  "pirs" = "pIRS",
  "dwgsim" = "DWGSIM",
  "wgsim" = "wgsim",
  "art_original" = "Original ART",
  "art_modern" = "art_modern (HEAD)",
  "art_modern_prev_ver" = "art_modern (master)",
  "art_modern_gcc" = "art_modern (GCC)",
  "art_modern_jemalloc" = "art_modern (HEAD/jemalloc)",
  "art_modern_asio" = "art_modern (HEAD/ASIO)",
  "art_modern_mimalloc" = "art_modern (HEAD/mi-malloc)"
)
levels <- c("wgsim", "DWGSIM", "pIRS", "art_modern (HEAD)", "art_modern (master)", "art_modern (GCC)", "art_modern (HEAD/jemalloc)","art_modern (HEAD/mi-malloc)","art_modern (HEAD/ASIO)",  "Original ART")
p <- df %>%
  dplyr::select(CPU_TIME, WALL_CLOCK, RSS, DATA, SOFTWARE, RLEN) %>%
  dplyr::mutate(RSS = RSS / 1024) %>%
  tidyr::pivot_longer(
    cols = c("CPU_TIME", "WALL_CLOCK", "RSS"),
    names_to = "ASPECTS",
    values_to = "VALUES"
  ) %>%
  dplyr::filter(!(SOFTWARE=="pirs" & DATA == "TRANSCRIPTOME" & ASPECTS=="CPU_TIME")) %>%
  dplyr::mutate(
    ASPECTS = sapply(ASPECTS, function(x)
      replacement_list[[x]]),
    SOFTWARE = sapply(SOFTWARE, function(x)
      replacement_list[[x]])
  ) %>%
  ggplot() +
  geom_boxplot(aes(
    y = factor(
      SOFTWARE,
      levels = levels
    ),
    x = VALUES,
    fill = RLEN,
    color = RLEN
  )) +
  scale_x_continuous(
    labels = scales::label_number(scale_cut = scales::cut_si("")),
    breaks = scales::breaks_pretty(n = 5),
    limits = c(0, NA)
  ) +
  scale_y_discrete(
    "Software"
  ) +
  scale_fill_discrete("Read length") +
  scale_color_discrete("Read length") +
  facet_grid(DATA ~ ASPECTS, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("fig/time_memory.pdf", p, width = 7, height = 4)

p <- df %>%
  dplyr::select(MAJ_PG_F, MIN_PG_F, VOL_CTX_S, IV_CTX_S, DATA, SOFTWARE, RLEN) %>%
  tidyr::pivot_longer(
    cols = c("MAJ_PG_F", "MIN_PG_F", "VOL_CTX_S", "IV_CTX_S"),
    names_to = "PFCS_TYPE",
    values_to = "PFCS"
  ) %>%
  dplyr::mutate(
    PFCS_TYPE = sapply(PFCS_TYPE, function(x)
      replacement_list[[x]]),
    SOFTWARE = sapply(SOFTWARE, function(x)
      replacement_list[[x]])
  ) %>%
  ggplot() +
  geom_boxplot(aes(
    y = factor(
      SOFTWARE,
      levels = levels
    ),
    x = PFCS,
    fill = RLEN,
    color = RLEN
  )) +
  scale_x_continuous(trans = "log10", labels = scales::label_number()) +
  facet_grid(DATA ~ PFCS_TYPE, scales = "free") +
  theme_bw()
ggsave("fig/page_faults.pdf", p, width = 12, height = 4)

sessionInfo()

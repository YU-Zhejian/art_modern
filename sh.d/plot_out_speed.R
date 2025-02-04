library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

df <- readr::read_tsv("cmake-build-relwithdebinfo-llvm_intel/time_complexity.tsv")

df_defaults <- df %>%
  dplyr::filter(name %in% c("BamReadOutput_l=4_t=1", "HeadlessBamReadOutput_l=4_t=1", "SamReadOutput_l=4_t=1", "FastqReadOutput", "PwaReadOutput", "FastaReadOutput", "DumbReadOutput", "EmptyLockedIOReadOutput", "EmptyLFIOReadOutput"))

ggplot(df_defaults) +
  geom_boxplot(aes(y=time_complexity, x=as.factor(threads))) +
  facet_wrap(~name, scales = "free_y")

df_bam_levels <- df %>%
  dplyr::filter(name %in% grep("BamReadOutput_l=._t=1", unique(name), value = TRUE))

ggplot(df_bam_levels) +
  geom_boxplot(aes(y=time_complexity, x=as.factor(threads))) +
  facet_wrap(~name)

df_bam_threads <- df %>%
  dplyr::filter(name %in% grep("BamReadOutput_l=4_t=.", unique(name), value = TRUE))

ggplot(df_bam_threads) +
  geom_boxplot(aes(y=time_complexity, x=as.factor(threads))) +
  facet_wrap(~name)


library(dplyr)
library(ggplot2)

all_depth <- data.frame()

for (software in c("art_modern", "art", "pirs", "wgsim", "dwgsim")){
  for(rlen in c("100", "300")){
    depth_path <- sprintf("data_out/depth/pe%s_%s.depth", rlen, software)
    depth_table <- readr::read_tsv(depth_path, col_names = c("chr", "pos", "depth"))%>%
      dplyr::mutate(software = software, rlen = rlen)
    all_depth <- rbind(all_depth, depth_table)
  }
}

p <- ggplot(all_depth) +
  geom_line(aes(x = pos, y = depth, color = software)) +
  scale_x_continuous(limits = c(0, 1000)) +
  theme_bw() +
  facet_grid(rlen ~ .)

ggsave("fig/depth_first1k.pdf", p, width = 8, height = 4)

p <- ggplot(all_depth) +
  geom_line(aes(x = pos, y = depth, color = software)) +
  scale_x_continuous(limits = c(max(all_depth$pos) - 1000, max(all_depth$pos))) +
  theme_bw() +
  facet_grid(rlen ~ .)

ggsave("fig/depth_last1k.pdf", p, width = 8, height = 4)

p <- ggplot(all_depth) +
  geom_violin(aes(y = software, x = depth, fill=rlen)) +
  geom_vline(xintercept = 100) +
  theme_bw()
ggsave("fig/depth_box.pdf", p, width = 5, height = 4)

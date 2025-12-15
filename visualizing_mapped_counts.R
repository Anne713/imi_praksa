setwd("/home/bioinfo/Spiroplasma/host_mapping/img")

library(data.table)
library(tidyverse)
library(patchwork)

df_full <- fread("../counts_per_window.sum.bed", col.names = c("chrom","start","end","count"))
df_full <- df_full[count > 10]

df <- as_tibble(df_full %>%
  filter(grepl("^NC_0000[0-2][0-9]|^NC_012920", chrom)))

df <- df %>%
  mutate(
    pos = ((start + end) / 2) / 1e6,
    chrom = factor(chrom, levels = unique(chrom))
  )

# compute cumulative offsets for concatenating chromosomes
chr_lengths <- df %>%
  group_by(chrom) %>%
  summarise(chr_len = max(end)/1e6, .groups = "drop") %>%
  mutate(offset = cumsum(lag(chr_len, default = 0)))

df <- df %>%
  left_join(chr_lengths %>% select(chrom, offset), by = "chrom") %>%
  mutate(pos_global = pos + offset)

# scatterplot position vs read counts
p <- ggplot(df, aes(x=pos_global, y=count)) +
  geom_point(size=0.3, alpha = 0.5) +
  theme_minimal() +
  theme(legend.position="none") +
  labs(
    title="Distribution of host reads across genome",
    x="Genomic position (Mb, concatenated chromosomes)",
    y="Read counts per 10kb window (logarithmic)"
  ) +
  scale_y_log10() 

ggsave("host_reads_distribution.png", p, bg = "white", width=12, height=6, dpi=300)

y_min <- min(df$count, na.rm = TRUE)
y_max <- max(df$count, na.rm = TRUE)

# per chromosome counts, linear scale 
plot_list <- list()

for (chr in unique(df$chrom)[1:24]) {
  p <- df %>%
    filter(chrom == chr) %>%
    ggplot(aes(x = pos, y = count)) + 
    geom_line() +
    theme_minimal() +
    labs(
      title = paste0("Read coverage across ", chr),
      x = "Position (Mb)",
      y = "Read counts (per 10kb window)"
    ) + 
    ylim(y_min, y_max)  
  
  plot_list[[chr]] <- p
  
  filename <- paste0(gsub("\\.", "_", chr), "_read_coverage.png")
  ggsave(filename, p, bg = "white", width = 12, height = 8, dpi = 300)
}

combined_plot <- wrap_plots(plot_list, ncol = 6, nrow = 4) &
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)
  )

ggsave("all_chrom_read_coverage.png", combined_plot, bg = "white", width=24, height=12, dpi=300)


# per chromosome counts, log scale 
plot_list <- list()

for (chr in unique(df$chrom)[1:24]) {
  p <- df %>%
    filter(chrom == chr) %>%
    ggplot(aes(x = pos, y = count)) + 
    geom_point() +
    theme_minimal() +
    labs(
      title = paste0("Read coverage across ", chr),
      x = "Position (Mb)",
      y = "Read counts (per 10kb window)"
    ) +
    scale_y_log10(limits = c(y_min, y_max))  
  
  plot_list[[chr]] <- p
  
  filename <- paste0(gsub("\\.", "_", chr), "_read_coverage_log.png")
  ggsave(filename, p, bg = "white", width = 12, height = 8, dpi = 300)
}

combined_plot <- wrap_plots(plot_list, ncol = 6, nrow = 4) &
  update_geom_defaults("point", list(size = 0.3, alpha = 0.5)) &  
  theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)
  )

ggsave("all_chrom_read_coverage_log.png", combined_plot, bg = "white", width=24, height=12, dpi=300)



setwd("C:/Users/okzme/OneDrive/Desktop/Alja/Spiroplasma")

library(tidyverse)
library(readxl)
library(writexl)
library(reshape2)
library(viridis)

original_data <- read_excel("rezultati.xlsx")
df <- original_data %>%
  select(-pro_id) %>%
  mutate(ngs_id = str_sub(ngs_id, -4, -1)) %>%
  arrange(ngs_id)

###############################################################################
# preprocessing

df <- df %>%
  # normalize to reads per million
  mutate(across(
    c(c_pav_s:k_silva),  
    ~ .x * 1000000 / c_pav_total,  
    .names = "{.col}_rpm"
  )) %>%
  pivot_longer(
    cols = ends_with("_rpm"),
    names_to = "pipeline",
    values_to = "rpm"
  ) %>%
  mutate(
    # collapse host + type into one column
    host_type = paste(host, type, sep = "_")
  ) %>%
  select(ngs_id, run, host_type, pipeline, rpm)

old_names <- df$pipeline %>% as.factor() %>% levels()
new_names <- c(
  "centrifuge_pavian_species", 
  "centrifuge_genus",
  "centrifuge_species",
  "kraken2_silva", 
  "nextflow_hostfree_kraken2_ncbi_genus", 
  "nextflow_kraken2_ncbi_genus",
  "nextflow_kraken2_ncbi_species",
  "nextflow_kraken2_silva_genus", 
  "nextflow_minimap2_ncbi_species", 
  "nextflow_minimap2_silva_genus"
)
name_map <- setNames(new_names, old_names)

df <- df %>%
  mutate(pipeline = recode(pipeline, !!!name_map))

write_tsv(df, "results/df_normalized.tsv")

###############################################################################
# comparisons between pipelines

# pipelines
pipeline_summary <- df %>%
  group_by(pipeline) %>%
  summarize(
    mean_rpm = mean(rpm, na.rm = TRUE),
    median_rpm = median(rpm, na.rm = TRUE),
    sd_rpm = sd(rpm, na.rm = TRUE),
    max_rpm = max(rpm, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(pipeline_summary, "results/pipeline_summary.tsv")

# host (and sample type) 
host_summary <- df %>%
  group_by(host_type) %>%
  summarize(
    mean_rpm = mean(rpm, na.rm = TRUE),
    median_rpm = median(rpm, na.rm = TRUE),
    sd_rpm = sd(rpm, na.rm = TRUE),
    max_rpm = max(rpm, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(host_summary, "results/host_summary.tsv")

# host and pipeline
host_pipeline_summary <- df %>%
  group_by(host_type, pipeline) %>%
  summarize(
    mean_rpm = mean(rpm, na.rm = TRUE),
    median_rpm = median(rpm, na.rm = TRUE),
    sd_rpm = sd(rpm, na.rm = TRUE),
    max_rpm = max(rpm, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(host_pipeline_summary, "results/host_pipeline_summary.tsv")

# boxplot per host_type & pipeline 
p_box <- ggplot(df, aes(x = pipeline, y = rpm)) +
  geom_boxplot() +
  facet_wrap(~host_type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(title = "Spiroplasma RPM per pipeline by sample type", y = "RPM (log scale)", x = "Pipeline") +
  scale_y_log10()

ggsave("results/boxplot_pipelines.png", plot = p_box, bg = 'white', width = 10, height = 6, dpi = 300)

# correlations between pipelines 
df_wide <- df %>%
  select(ngs_id, pipeline, rpm) %>%
  pivot_wider(names_from = pipeline, values_from = rpm)

cor_matrix <- df_wide %>%
  select(-ngs_id) %>%
  cor(use = "pairwise.complete.obs", method = "spearman")

# visualize correlations as a heatmap
cor_long <- melt(cor_matrix, varnames = c("Pipeline1", "Pipeline2"), value.name = "Spearman_rho")

p_heat <- ggplot(cor_long, aes(x = Pipeline1, y = Pipeline2, fill = Spearman_rho)) +
  geom_tile() +
  scale_fill_viridis_c(option = "cividis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Spearman correlation between pipelines", fill = "rho")

ggsave("results/pipeline_correlations.png", plot = p_heat, bg = 'white', width = 8, height = 5, dpi = 300)

###############################################################################
# comparison between runs 2 and 3

# wilcoxon rank-sum test 
wilcox_results <- df %>%
  group_by(pipeline) %>%
  summarise(
    n_run2 = sum(run == 2, na.rm = TRUE),
    n_run3 = sum(run == 3, na.rm = TRUE),
    med_run1 = median(rpm[run == 2], na.rm = TRUE),
    med_run2 = median(rpm[run == 3], na.rm = TRUE),
    p_value = if (!is.na(n_run2) & !is.na(n_run3)) {
      wilcox.test(
        x = rpm[run == 2],
        y = rpm[run == 3]
      )$p.value
    } else {
      NA
    },
    .groups = "drop"
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) # multiple testing correction

write_tsv(wilcox_results, "results/wilcoxon_run_comparison.tsv")

# bins for low but not zero reads 
df_bins <- df %>%
  filter(!is.na(rpm)) %>%
  mutate(read_bin = case_when(
    rpm == 0 ~ "zero",
    rpm > 0 & rpm < 10 ~ "low",
    rpm >= 10 ~ "pos"
  )) %>%
  mutate(read_bin = factor(read_bin, levels = c("zero", "low", "pos")))

# contingency table of run vs bins
bin_table <- df_bins %>%
  count(run, read_bin) %>%
  pivot_wider(names_from = read_bin, values_from = n, values_fill = 0) %>%
  rowwise() %>%
  mutate(total = sum(c_across(zero:pos))) %>%
  mutate(across(zero:pos, ~ .x / total)) 

write_tsv(bin_table, "results/run_vs_bin_table.tsv")

# run 2 vs run 3
df_low <- df %>%
  filter(run %in% c(2,3)) %>%
  mutate(is_low = if_else(rpm > 0 & rpm < 10, 1, 0))

prop_res <- prop.test(table(df_low$is_low, df_low$run))
capture.output(prop_res, file = "results/prop_test_low_per_run.txt")

# visualize proportion of low reads per run
ggplot(df_bins, aes(x = factor(run), fill = read_bin)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Run", y = "Proportion", title = "Read count categories by run") +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "cividis")

ggsave("results/run_vs_bin_barplot.png", bg = 'white', width = 6, height = 4)

###############################################################################
# genus level hits

df_subset <- df %>%
  filter(pipeline %in% c("centrifuge_species", "centrifuge_genus")) %>%
  mutate(host = str_extract(host_type, "^[^_]+"),
         level = str_extract(pipeline, "(?<=_).*")) %>% 
  select(-host_type, -pipeline)

summary_genus <- df_subset %>%
  group_by(host, level) %>%
  summarise(
    mean_rpm = mean(rpm, na.rm = TRUE),
    median_rpm = median(rpm, na.rm = TRUE),
    sd_rpm = sd(rpm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(host, level)

write_tsv(summary_genus, "results/genus_summary.tsv")

# boxplot
p_genus_box <- ggplot(df_subset, aes(x = level, y = rpm, fill = host)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_y_log10() +
  scale_fill_manual(values = c("human" = "#4c5073", "tick" = "#edd51c")) +
  labs(
    title = "Centrifuge maximum taxonomic level matched by host",
    x = "Taxonomic level",
    y = "RPM (log scale)",
    fill = "Host"
  ) +
  theme_minimal() 

ggsave("results/genus_boxplot.png", plot = p_genus_box, bg = 'white', width = 8, height = 5)

###############################################################################
# plasma vs blood

paired_ids <- original_data %>%
  filter(host == "human") %>%
  mutate(main_type = str_extract(type, "^[^_]+")) %>%
  group_by(pro_id) %>%
  summarize(types_present = list(unique(main_type)), .groups = "drop") %>%
  filter(map_lgl(types_present, ~ all(c("blood", "plasma") %in% .x))) %>%
  pull(pro_id)

df_with_id <- df %>%
  left_join(
    original_data %>%
      select(ngs_id, pro_id, type) %>%
      mutate(ngs_id = str_sub(ngs_id, -4, -1),
             type = str_extract(type, "^[^_]+")),
    by = "ngs_id"
  )

df_paired <- df_with_id %>%
  filter(pro_id %in% paired_ids, 
         pipeline != "kraken2_silva") %>%
  group_by(pro_id, type, pipeline) %>%
  summarize(rpm = mean(rpm, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = rpm)

# fold changes 
df_fc <- df_paired %>%
  mutate(
    fc = (plasma + 1) / (blood + 1),
  ) %>%
  group_by(pipeline) %>%
  summarize(
    mean_fc = mean(fc, na.rm = TRUE),
    median_fc = median(fc, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(df_fc, "results/fold_changes_plasma_vs_blood.tsv")

# scatter plot
bloody_p <- ggplot(df_paired, aes(x = blood, y = plasma)) +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Blood RPM",
    y = "Plasma RPM",
    title = "Blood vs. plasma Spiroplasma reads"
  ) +
  theme_minimal()

ggsave("results/scatter_blood_vs_plasma.png", plot = bloody_p, bg = 'white', width = 9, height = 5)

# paired line plot
df_long <- df_paired %>%
  pivot_longer(c(blood, plasma), names_to = "type", values_to = "rpm")

paired_blood <- ggplot(df_long, aes(x = type, y = rpm, group = pro_id)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  facet_wrap(~pipeline, scales = "free_y") +
  scale_y_log10() +
  labs(
    title = "Blood vs. plasma Spiroplasma reads",
    x = "Sample type",
    y = "RPM"
  ) +
  theme_minimal()

ggsave("results/paired_lines_blood_vs_plasma.png", plot = paired_blood, bg = 'white', width = 8, height = 6)

###############################################################################
# hostfree human vs tick

df_host_compare <- df %>%
  filter(pipeline %in% c(
    "nextflow_hostfree_kraken2_ncbi_genus", 
    "nextflow_kraken2_ncbi_genus"
  )) %>%
  mutate(host = str_extract(host_type, "^[^_]+"),
    host = factor(host)) %>% 
  select(-host_type, -run) 

# fold change
df_fc_host <- df_host_compare %>%
  pivot_wider(
    names_from = pipeline, 
    values_from = rpm
  ) %>%
  mutate(
    fc = (`nextflow_kraken2_ncbi_genus` + 1) / 
      (`nextflow_hostfree_kraken2_ncbi_genus` + 1),
  )

# summary per host
summary_host <- df_fc_host %>%
  group_by(host) %>%
  summarize(
    mean__fc = mean(fc, na.rm = TRUE),
    median__fc = median(fc, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(summary_host, "results/fold_changes_hostfree.tsv")

# boxplot of fold changes by host
hosts <- ggplot(df_fc_host, aes(x = host, y = fc, fill = host)) +
  geom_boxplot() +
  labs(
    title = "Comparison of Kraken2 pipelines by host",
    x = "Host",
    y = "Fold change (standard / hostfree)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("human" = "#4c5073", "tick" = "#edd51c")) +
  scale_y_log10() 

ggsave("results/fold_change_hostfree_plot.png", plot = hosts, bg = 'white', width = 8, height = 6)


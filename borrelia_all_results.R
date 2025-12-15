setwd("/home/bioinfo/Spiroplasma")

library(tidyverse)
library(readxl)
library(writexl)

target_species = "Borrelia miyamotoi"
target_species_full = "Borrelia_miyamotoi"
target_genus = "Borrelia"

data <- read_excel("zbrani_vzorci/ids_all.xlsx")[, 1:3]

# Centrifuge 
###############################################################################
# tsv tabele s stevilom readov na takson
data$c_tab_s <- NA
data$c_tab_g <- NA

for(i in seq_len(nrow(data))) {
  sample_id <- data$ngs_id[i]
  file_path <- paste0("klasifikacija_centrifuge_nt/tabele/", sample_id, "_report.tsv")
  
  if(file.exists(file_path)) {
    tsv_data <- read_tsv(file_path, col_types = cols())
    
    species_reads <- tsv_data$numReads[tsv_data$name == target_species]
    genus_reads   <- tsv_data$numReads[tsv_data$name == target_genus]
    # default to 0 if not found
    data$c_tab_s[i] <- if(length(species_reads) > 0) species_reads else 0
    data$c_tab_g[i] <- if(length(genus_reads) > 0) genus_reads else 0
    
  } else {
    warning(paste("File not found: ", file_path))
  }
}

# summary od Paviana iz kreportov
pav_species <- read_tsv("klasifikacija_centrifuge_nt/b_miyamotoi_all.tsv", col_types = cols())
pav_species <- pav_species %>%
  mutate(
    ngs_id = gsub("^X|\\.cladeReads$", "", name), 
    ngs_id = gsub("\\.", "-", ngs_id)
  ) %>%
  select(ngs_id, c_pav_s = `Borrelia miyamotoi`) %>%
  mutate(c_pav_s = replace_na(c_pav_s, 0))

data <- left_join(data, pav_species, by = "ngs_id")

pav_total <- read_excel("klasifikacija_centrifuge_nt/summary_pavian_all.xlsx")[,1:2] %>%
  select(ngs_id = Name, c_pav_total = `Number of raw reads`)

data <- left_join(data, pav_total, by = "ngs_id")

# Nextflow 
###############################################################################
# results mape za vsak vzorec s tsv filom s st readov na takson (zapisan cel lineage)
extract_nextflow_count <- function(sample_id, dir, tax_lvl) {
  file <- file.path(dir, paste0(sample_id, "_results/abundance_table_", tax_lvl, ".tsv"))
  target <- if (tax_lvl == "species") target_species else target_genus
  
  if (!file.exists(file)) return(NA)
  
  tsv <- read_tsv(file, col_types = cols())
  row <- tsv$total[grepl(target, tsv$tax)]
  
  if (length(row) > 0) row else 0
}

data <- data %>%
  rowwise() %>%
  mutate(
    n_k2_silva_g = extract_nextflow_count(id, "xwf-16s_results/results_20250822_wf_kraken2_silva_138_1/", "genus"),
    n_k2_ncbi_s = extract_nextflow_count(id, "xwf-16s_results/results_20250822_wf_kraken2_ncbi_16s_18s/", "species"),
    n_mm_silva_g = extract_nextflow_count(id, "xwf-16s_results/results_20250822_wf_minimap_SILVA_138_1/", "genus"),
    n_mm_ncbi_s = extract_nextflow_count(id, "xwf-16s_results/results_20250822_wf_minimap_ncbi_16s_18s/", "species"), 
    n_k2_ncbi_g = extract_nextflow_count(id, "wf-16s/results/", "genus"), 
    n_hostfree_k2_ncbi_g = extract_nextflow_count(id, "wf-16s/results_host_g/", "genus") 
  ) %>%
  ungroup()

# Kraken
###############################################################################
# run 2: log file s stats headerjem in tabelo numreads, covbases ... na ref genom
read_log_table <- function(file) {
  header_line <- read_lines(file) %>%
    str_which("^#rname") 
  
  if (length(header_line) == 0) {
    stop(paste("No table header (#rname) found in ", file))
  }
  
  read_tsv(file, skip = header_line - 1, col_types = cols()) 
}

data$k_silva <- NA

for(i in seq_len(nrow(data))) {
  sample_id <- data$id[i]
  file_path <- paste0("xrezultati_kraken2_silva_mappiranje/05_16s_2025_08_14/02_mappiranje_spiroplazme/", sample_id, "_stats.log")
  
  if(file.exists(file_path)) {
    tsv_data <- read_log_table(file_path)
    
    data$k_silva[i] <- sum(tsv_data$numreads[grepl(target_species_full, tsv_data$`#rname`)], na.rm = TRUE)
    
  } else {
    warning(paste("File not found: ", file_path))
  }
}

###############################################################################

write_xlsx(data, "zbrani_vzorci/borrelia_results_comparison.xlsx")

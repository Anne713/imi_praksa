# skripta za extraction best hita in best target hita iz rezultatov mapiranja 

library(tidyverse)

min_cov = 90
target_org = "Spiroplasma"

# extract table iz vsakega log fila, omitting header in vrstice z 0 reads 
read_log_table <- function(file) {
  header_line <- read_lines(file) %>%
    str_which("^#rname") 
  
  if (length(header_line) == 0) {
    stop(paste("No table header (#rname) found in ", file))
  }
  
  read_tsv(
    file,
    skip = header_line - 1,
    col_types = cols()
  ) %>% 
    filter(numreads > 0) %>%
    mutate(sample = str_remove(basename(file), "_stats\\.log$")) %>%
    relocate(sample, .before = 1) %>%
    arrange(desc(numreads), desc(covbases), desc(meandepth))
}

files <- list.files(pattern = "\\.log$")

# criteria za best hit
find_best <- function(tbl) {
  if (any(tbl$coverage > min_cov)) {
    tbl %>%
      filter(coverage > min_cov) %>%
      arrange(desc(numreads)) %>%
      slice(1)
  } else {
    tbl %>%
      arrange(desc(coverage)) %>%
      slice(1)
  }
}

# najdi najboljsega in najboljsi target, ce sta enaka vrne samo enkrat
pick_best_hits <- function(tbl) {
  best_overall <- find_best(tbl)
  
  if (any(str_detect(tbl$`#rname`, target_org))) {
    best_target <- find_best(tbl %>% filter(str_detect(`#rname`, target_org)))
  } else {
    best_target <- NULL
  }
  
  if (best_overall$`#rname` == best_target$`#rname`) {
    best_target <- NULL
  }
  
  best_overall <- best_overall %>% mutate(hittype = "overall")
  if (!is.null(best_target)) {
    best_target <- best_target %>% mutate(hittype = "target")
    return(bind_rows(best_overall, best_target))
  } else {
    return(best_overall)
  }
}

# choose best, zapisi v tabelo best hitov 
best_hits <- map_dfr(files, function(f) {
  tbl <- read_log_table(f)
  if (nrow(tbl) == 0) return(NULL) 
  pick_best_hits(tbl)
})

write_tsv(best_hits, "best_hits_per_sample.tsv")
message("Done! Best hits saved to best_hits_per_sample.tsv") 

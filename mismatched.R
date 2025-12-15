setwd("/home/bioinfo/Spiroplasma")

library(tidyverse)
library(readxl)
library(writexl)

data <- read_excel("zbrani_vzorci/Spiroplasma_ixodetis_bolniki_0827_manj.xlsx")

ggplot(data, aes(x = `23S`, y = log10(`16S`))) +
  geom_point() +
  labs(
    x = "ct 23S PCR S.ixodetis",
    y = "16S S.ixodetis reads (log scale)",
    title = "Scatterplot of 23S S.ixodetis reads vs. 16S S.ixodetis reads"
  )

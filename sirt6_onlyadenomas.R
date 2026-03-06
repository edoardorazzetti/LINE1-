library(tidyverse)
library(readxl)
library(readxl)
library(dplyr)

library(readxl)
library(dplyr)
library(tidyr)
library(biomaRt)
library(ggpubr)

library(tidyverse)
library(biomaRt)
library(ggpubr)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

data <- read_csv("fulldataset.csv")

# -----------------------------------------------------------------------------
# 2. Extract sample metadata from annotation rows
# -----------------------------------------------------------------------------

tissue_info <- data.frame(
  Sample = colnames(data)[-1],
  Tissue = unlist(data[2, -1], use.names = FALSE),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 3. Filter samples: Adenoma and Adenoma adjacent only
# -----------------------------------------------------------------------------

campioni_adenoma <- tissue_info$Sample[
  tissue_info$Tissue %in% c("Adenoma", "Adenoma adjacent")
]

data_adenoma    <- dplyr::select(data, ID, all_of(campioni_adenoma))
metadata_adenoma <- tissue_info %>%
  filter(Tissue %in% c("Adenoma", "Adenoma adjacent"))

cat("Campioni selezionati:", length(campioni_adenoma), "\n")

# -----------------------------------------------------------------------------
# 4. Retrieve SIRT6 transcript IDs from Ensembl via biomaRt
# -----------------------------------------------------------------------------

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

sirt6_transcripts <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters    = "external_gene_name",
  values     = "SIRT6",
  mart       = mart,
  useCache   = FALSE
)

cat("Trascritti SIRT6 da Ensembl:", nrow(sirt6_transcripts), "\n")

# -----------------------------------------------------------------------------
# 5. Filter SIRT6 transcripts from expression matrix
#    (strip version suffix before matching, e.g. ENST00000305232.10 → ENST00000305232)
# -----------------------------------------------------------------------------

sirt6_expr <- data_adenoma %>%
  filter(gsub("\\..*", "", ID) %in% sirt6_transcripts$ensembl_transcript_id) %>%
  mutate(across(-ID, as.numeric))

cat("Trascritti SIRT6 trovati nel dataset:", nrow(sirt6_expr), "\n")

# -----------------------------------------------------------------------------
# 6. Compute per-sample median expression across SIRT6 transcripts
# -----------------------------------------------------------------------------

sirt6_mediana_paziente <- sirt6_expr %>%
  summarise(across(-ID, \(x) median(x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Mediana_SIRT6") %>%
  left_join(metadata_adenoma, by = "Sample")

print(sirt6_mediana_paziente)

# Median per tissue group
sirt6_mediana_paziente %>%
  group_by(Tissue) %>%
  summarise(Mediana = median(Mediana_SIRT6)) %>%
  print()

# -----------------------------------------------------------------------------
# 7. Boxplot with Wilcoxon test
# -----------------------------------------------------------------------------

ggplot(sirt6_mediana_paziente, aes(x = Tissue, y = Mediana_SIRT6, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  theme_classic(base_size = 13) +
  labs(
    title = "Espressione di SIRT6 — Adenoma vs Adenoma Adjacent",
    x     = "Gruppo tissutale",
    y     = "Mediana espressione SIRT6"
  ) +
  theme(legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Adenoma", "Adenoma adjacent")),
    method      = "wilcox.test",
    label       = "p.signif"
  )

# -----------------------------------------------------------------------------
# 8. Classify patients as High / Low SIRT6 expression (relative to global median)
# -----------------------------------------------------------------------------

mediana_globale <- median(sirt6_mediana_paziente$Mediana_SIRT6)

sirt6_classificato <- sirt6_mediana_paziente %>%
  mutate(Gruppo_mediana = ifelse(Mediana_SIRT6 >= mediana_globale, "Alto", "Basso"))

casistica <- sirt6_classificato %>%
  group_by(Tissue, Gruppo_mediana) %>%
  summarise(N_pazienti = n(), .groups = "drop") %>%
  arrange(Tissue, Gruppo_mediana)

print(casistica)



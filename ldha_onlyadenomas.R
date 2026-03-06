# =============================================================================
# LDHA Expression Analysis — Adenoma vs Adenoma Adjacent
# =============================================================================
# Analyzes LDHA transcript-level expression in adenoma and adjacent tissue
# using RNA-seq data. Includes median expression per patient, group comparison,
# and high/low expression classification.
#
# Input:  fulldataset.csv  — transcript-level expression matrix
#                            (rows = transcripts, cols = samples)
#                            rows 1–2 contain "Type" and "Tissue" annotations
# Output: boxplot + casistica table
# =============================================================================

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

data_adenoma     <- dplyr::select(data, ID, all_of(campioni_adenoma))
metadata_adenoma <- tissue_info %>%
  filter(Tissue %in% c("Adenoma", "Adenoma adjacent"))

cat("Campioni selezionati:", length(campioni_adenoma), "\n")

# -----------------------------------------------------------------------------
# 4. Retrieve LDHA transcript IDs from Ensembl via biomaRt
# -----------------------------------------------------------------------------

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ldha_transcripts <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters    = "external_gene_name",
  values     = "LDHA",
  mart       = mart,
  useCache   = FALSE
)

cat("Trascritti LDHA da Ensembl:", nrow(ldha_transcripts), "\n")

# -----------------------------------------------------------------------------
# 5. Filter LDHA transcripts from expression matrix
#    (strip version suffix before matching, e.g. ENST00000305232.10 → ENST00000305232)
# -----------------------------------------------------------------------------

ldha_expr <- data_adenoma %>%
  filter(gsub("\\..*", "", ID) %in% ldha_transcripts$ensembl_transcript_id) %>%
  mutate(across(-ID, as.numeric))

cat("Trascritti LDHA trovati nel dataset:", nrow(ldha_expr), "\n")

# -----------------------------------------------------------------------------
# 6. Compute per-sample median expression across LDHA transcripts
# -----------------------------------------------------------------------------

ldha_mediana_paziente <- ldha_expr %>%
  summarise(across(-ID, \(x) median(x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Mediana_LDHA") %>%
  left_join(metadata_adenoma, by = "Sample")

print(ldha_mediana_paziente)

# Median per tissue group
ldha_mediana_paziente %>%
  group_by(Tissue) %>%
  summarise(Mediana = median(Mediana_LDHA)) %>%
  print()

# -----------------------------------------------------------------------------
# 7. Boxplot with Wilcoxon test
# -----------------------------------------------------------------------------

ggplot(ldha_mediana_paziente, aes(x = Tissue, y = Mediana_LDHA, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  theme_classic(base_size = 13) +
  labs(
    title = "Espressione di LDHA — Adenoma vs Adenoma Adjacent",
    x     = "Gruppo tissutale",
    y     = "Mediana espressione LDHA"
  ) +
  theme(legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Adenoma", "Adenoma adjacent")),
    method      = "wilcox.test",
    label       = "p.signif"
  )

# -----------------------------------------------------------------------------
# 8. Classify patients as High / Low LDHA expression (relative to global median)
# -----------------------------------------------------------------------------

mediana_globale <- median(ldha_mediana_paziente$Mediana_LDHA)

ldha_classificato <- ldha_mediana_paziente %>%
  mutate(Gruppo_mediana = ifelse(Mediana_LDHA >= mediana_globale, "Alto", "Basso"))

casistica <- ldha_classificato %>%
  group_by(Tissue, Gruppo_mediana) %>%
  summarise(N_pazienti = n(), .groups = "drop") %>%
  arrange(Tissue, Gruppo_mediana)

print(casistica)

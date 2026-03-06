library(readxl)
library(dplyr)
metadata <- read_excel("metadata.xlsx")
all(metadata$ID %in% colnames(livelli_di_espressione))
type_row <- c("Type", metadata$Type[match(colnames(livelli_di_espressione)[-1], metadata$ID)])
tissue_row <- c("Tissue", metadata$Tissue[match(colnames(livelli_di_espressione)[-1], metadata$ID)])
livelli_annotato <- rbind(type_row, tissue_row, livelli_di_espressione)

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

sirt6_transcripts <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters = "external_gene_name",
  values = "SIRT6",
  mart = mart,
  useCache = FALSE
)

print(sirt6_transcripts)
# Filtra i trascritti di SIRT6 nel tuo dataset
sirt6_expr <- livelli_di_espressione %>%
  filter(grepl(paste(sirt6_transcripts$ensembl_transcript_id, collapse="|"), ID))
cat("Trascritti SIRT6 trovati:", nrow(sirt6_expr), "\n")
print(sirt6_expr$ID)
library(tidyr)
sirt6_mediana <- sirt6_expr %>%
  summarise(across(-ID, median)) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Mediana_SIRT6") %>%
  left_join(metadata %>% select(ID, Type, Tissue), by = c("Sample" = "ID"))
print(sirt6_mediana)
sirt6_mediana_tissue <- sirt6_expr %>%
  summarise(across(-ID, median)) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Mediana_SIRT6") %>%
  left_join(metadata %>% select(ID, Tissue), by = c("Sample" = "ID")) %>%
  group_by(Tissue) %>%
  summarise(Mediana = median(Mediana_SIRT6))

print(sirt6_mediana_tissue)




library(ggpubr)

sirt6_mediana_paziente <- sirt6_expr %>%
  summarise(across(-ID, median)) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Mediana_SIRT6") %>%
  left_join(metadata %>% select(ID, Tissue), by = c("Sample" = "ID"))




ggplot(sirt6_mediana_paziente, aes(x = Tissue, y = Mediana_SIRT6, fill = Tissue)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_classic() +
  labs(title = "Espressione di SIRT6 per gruppo tissutale",
       x = "Gruppo",
       y = "Mediana espressione SIRT6") +
  theme(legend.position = "none") +
  stat_compare_means(
    comparisons = list(
      c("CRC", "CRC adjacent"),
      c("Adenoma", "Adenoma adjacent"),
      c("CRC", "Adenoma")
    ),
    method = "wilcox.test",
    label = "p.signif"  # mostra *, **, *** oppure metti label = "p.format" per il valore esatto
  )

mediana_globale <- median(sirt6_mediana_paziente$Mediana_SIRT6)
sirt6_classificato <- sirt6_mediana_paziente %>%
  mutate(Gruppo_mediana = ifelse(Mediana_SIRT6 >= mediana_globale, "Alto", "Basso"))

casistica <- sirt6_classificato %>%
  group_by(Tissue, Gruppo_mediana) %>%
  summarise(N_pazienti = n(), .groups = "drop") %>%
  arrange(Tissue, Gruppo_mediana)

print(casistica)

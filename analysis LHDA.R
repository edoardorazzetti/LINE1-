# Trova i trascritti di LDHA
ldha_transcripts <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters = "external_gene_name",
  values = "LDHA",
  mart = mart,
  useCache = FALSE
)

# Filtra nel dataset (con grepl per gestire i suffissi)
ldha_expr <- livelli_di_espressione %>%
  filter(grepl(paste(ldha_transcripts$ensembl_transcript_id, collapse="|"), ID))

cat("Trascritti LDHA trovati:", nrow(ldha_expr), "\n")

# Calcola la mediana per paziente
ldha_mediana_paziente <- ldha_expr %>%
  summarise(across(-ID, median)) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Mediana_LDHA") %>%
  left_join(metadata %>% select(ID, Tissue), by = c("Sample" = "ID"))

# Boxplot con significatività
ggplot(ldha_mediana_paziente, aes(x = Tissue, y = Mediana_LDHA, fill = Tissue)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_classic() +
  labs(title = "Espressione di LDHA per gruppo tissutale",
       x = "Gruppo",
       y = "Mediana espressione LDHA") +
  theme(legend.position = "none") +
  stat_compare_means(
    comparisons = list(
      c("CRC", "CRC adjacent"),
      c("Adenoma", "Adenoma adjacent"),
      c("CRC", "Adenoma")
    ),
    method = "wilcox.test",
    label = "p.signif"
  )

# Casistica sopra/sotto mediana
mediana_globale_ldha <- median(ldha_mediana_paziente$Mediana_LDHA)

ldha_classificato <- ldha_mediana_paziente %>%
  mutate(Gruppo_mediana = ifelse(Mediana_LDHA >= mediana_globale_ldha, "Alto", "Basso"))

casistica_ldha <- ldha_classificato %>%
  group_by(Tissue, Gruppo_mediana) %>%
  summarise(N_pazienti = n(), .groups = "drop") %>%
  arrange(Tissue, Gruppo_mediana)

print(casistica_ldha)
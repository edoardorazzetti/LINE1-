setwd("C:/Users/edino/OneDrive/Desktop/tesi/00_Input")
getwd()
genexpression<-read.delim("TPM_genes_CRCadj.tsv")
df <- df %>%
  left_join(genexpression, by = "tss_tss_gene_id")
library(dplyr)

# 1. Prepar0 un sotto-dataframe solo con le colonne utili e rinomino per il join
genexpression_clean <- genexpression %>%
  select(gene_name, gene_type) %>%
  rename(tss_tss_gene_id = gene_name)

# 2. Faccio il join col tuo dataframe df
df <- df %>%
  left_join(genexpression_clean, by = "tss_tss_gene_id")

# 3. Creo una colonna  is_coding (TRUE se protein_coding, FALSE altrimenti)
df <- df %>%
  mutate(is_coding = gene_type == "protein_coding")
##faccio con Biomart
# Installa se serve: BiocManager::install("biomaRt")
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Recupera tabella di mapping
n
n
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mart
) %>%
  rename(tss_tss_gene_id = ensembl_gene_id,
         gene_name = hgnc_symbol)
# 1. Unisci df con i nomi dei geni
df_with_names <- df %>%
  left_join(gene_annotations, by = "tss_tss_gene_id")

# 2. Poi fai il join con genexpression
df_final <- df_with_names %>%
  left_join(genexpression %>% select(gene_name, gene_type), by = "gene_name") %>%
  mutate(is_coding = gene_type == "protein_coding")


#Confrontare i livelli di espressione dei geni con LINE1 (presenti in df) con quelli senza LINE1
genexpression$has_LINE1 <- genexpression$gene_id %in% df$tss_tss_gene_id
expression_levels<-read.delim("Expression_Level_TPM_CRC_Adenoma.txt")
#aggiungo lenghts
toadd <- lenghts %>% filter(Transcript.stable.ID %in% expression$tss_tss_transcript_id)
toadd <- toadd[,c(3,5)]
#File con livelli di espressione
livelli_di_espressione<-read.delim("Expression_Level_TPM_CRC_Adenoma.txt")
#copio codice dalla console
getwd()
"C:/Users/edino/OneDrive/Desktop/tesi/00_Input"
genexpression<-read.delim("TPM_genes_CRCadj.tsv")
df <- df %>%
  +   left_join(genexpression, by = "tss_tss_gene_id")

View(genexpression)
expression<-read.delim("overlaps_same_opposite.txt")
lenghts<-read.delim("mart_export (1).txt")
View(lenghts)
lenghts<-read.delim("mart_export (1).txt")
merged_data <- merge(expression, lenghts, by = "isoform_id")
#aggiungo lenghts
lenghts %>% filter(Transcript.stable.ID %in% expression$tss_tss_transcript_id)

#aggiungo lenghts
toadd <- lenghts %>% filter(Transcript.stable.ID %in% expression$tss_tss_transcript_id)
toadd<-toadd %>% rename("tss_tss_transcript_id"=Transcript.stable.ID)
colnames(expression)
expression <- merge(expression, toadd, by = "tss_tss_transcript_id", all.x = TRUE)
livelli_di_espressione<-read.delim("Expression_Level_TPM_CRC_Adenoma.txt")
View(livelli_di_espressione)
library(tidyverse)
livelli_di_espressione<-data.frame(livelli_di_espressione)
livelli_di_espressione <- livelli_di_espressione %>%
  select(-VT437S, -VT448K, -VT462S)
livelli_di_espressione<-t(livelli_di_espressione)
livelli_di_espressione<-column_to_rownames(livelli_di_espressione,"ID")
livelli_di_espressione$class<-tissues$Tissue
tissues <- read_excel("metadata.xlsx")
livelli_di_espressione <- livelli_di_espressione %>%
  dplyr::select(-VT437S, -VT448K, -VT462S)
livelli_di_espressione<-t(livelli_di_espressione)
livelli_di_espressione <- livelli_di_espressione %>%
  select(ID, everything())
livelli_di_espressione <- livelli_di_espressione %>%
  dplyr::select(-c("-VT437S", "-VT448K", "-VT462S"))
righe_da_eliminare_ids <- c("-VT437S", "-VT448K", "-VT462S")
livelli_di_espressione <- livelli_di_espressione %>%
  rownames_to_column(var = "temp_id_col") %>%
  filter(!temp_id_col %in% righe_da_eliminare_ids) %>%
  column_to_rownames(var = "temp_id_col")  
rownames(livelli_di_espressione)=NULL
install.packages("readxl")
library(readxl)
metadata<-read_excel("metadata.xlsx")
livelli_di_espressione$sample<-metadata$Tissue



##mediana for each class
livelli_di_espressione <-
livelli_di_espressione%>%gather(isoforms,level,-sample)%>%group_by(sample,isoforms)%>%summarize(median=median(level))


isoformsL1<-read.delim("overlaps_same_opposite.txt")
#rimuovo versione isoforma

names(livelli_di_espressione)=gsub("\\..", "", names(livelli_di_espressione))
#trovo le isoforme che voglio tenere
id_trascritti_da_mantenere <- unique(isoformsL1$tss_tss_transcript_id)
livelli_di_espressione <- livelli_di_espressione %>%
  filter(tss_tss_transcript_id %in% id_trascritti_da_mantenere)
names(livelli_di_espressione)
names(isoformsL1)
ids_da_cercare <- unique(isoformsL1$tss_tss_transcript_id)
nomi_colonne_livelli <- names(livelli_di_espressione)
colonne_da_mantenere <- intersect(nomi_colonne_livelli, ids_da_cercare)
livelli_di_espressione <- livelli_di_espressione %>%
  select(all_of(colonne_da_mantenere)) 
#calcolo mediana
mediane_per_categoria_isoforma <- livelli_di_espressione %>%
  pivot_longer(
    cols = -sample, # Impila tutte le colonne eccetto 'sample'
    names_to = "isoforms", # Le nuove colonne impilate verranno chiamate 'isoforms'
    values_to = "level"    # I valori corrispondenti verranno chiamati 'level'
  ) %>%
  group_by(sample, isoforms) %>% # Raggruppa per la categoria del campione e per ogni isoforma
  summarize(median_level = median(level, na.rm = TRUE), .groups = "drop") # Calcola la mediana, ignorando i valori NA

mediane_wider <- mediane_per_categoria_isoforma %>%
  pivot_wider(
    names_from = sample,  # Prendi i valori della colonna 'sample' e usali come nomi delle nuove colonne
    values_from = median_level # Prendi i valori della colonna 'median_level' per popolare le nuove colonne
  )

#Plotto
library(ggplot2)
ggplot(mediane_per_categoria_isoforma, aes(x = sample, y = median_level, fill = sample)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.6, width = 0.2) + # Aggiunge i punti per visualizzare i singoli dati
  labs(
    title = "Distribuzione dei Livelli di Espressione Medi di Isoforme con L1",
    x = "Categoria del Campione",
    y = "Mediana del Livello di Espressione (log10+1)"
  ) +
  theme_bw() + # Un tema pulito per il grafico
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#cambio nomi
lunghezze<-lenghts
rm(lenghts)
rm(mediane_per_categoria_isoforma)
espressione_generale<-livelli_di_espressione
rm(livelli_di_espressione)
espressioneL1<-mediane_wider
rm(mediane_wider)

#unisco mediane
isoformsL1_con_mediane <- isoformsL1 %>%
  left_join(espressioneL1, by = c("tss_tss_transcript_id" = "isoforms"))

#Faccio boxplot con facet
ggplot(isoformsL1_con_mediane, aes(x = Sample_Category, y = Median_Expression_log10, fill = Sample_Category)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.6, width = 0.2) +
  labs(
    title = "Distribuzione dei Livelli di Espressione Medi (log10) per Categoria e Overlap Class",
    x = "Categoria del Campione",
    y = "Mediana del Livello di Espressione (log10 + 1)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Overlap_class)


#questo lo fa gemini
data_for_boxplot <- isoformsL1_con_mediane %>%
  pivot_longer(
    cols = all_of(category_cols), # Impila solo le colonne delle categorie specificate
    names_to = "Sample_Category", # La nuova colonna per i nomi delle categorie (e.g., "Adenoma")
    values_to = "Median_Expression_log10" # La nuova colonna per i valori mediani log10
  )

# 2. Crea il boxplot con facet_wrap
ggplot(data_for_boxplot, aes(x = Sample_Category, y = Median_Expression_log10, fill = Sample_Category)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.6, width = 0.2) + # Aggiunge i punti per visualizzare i singoli dati
  labs(
    title = "Distribuzione dei Livelli di Espressione Medi (log10) per Categoria e Overlap Class",
    x = "Categoria del Campione",
    y = "Mediana del Livello di Espressione (log10 + 1)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Overlap_class) # Questa è la parte per dividere i grafici per Overlap_class


#Cambio il formato come richiesto dalla library
category_cols <- c("Adenoma", "Adenoma adjacent", "CRC", "CRC adjacent")

data_for_boxplot <- isoformsL1_con_mediane %>%
  pivot_longer(
    cols = all_of(category_cols), # Impila solo queste colonne
    names_to = "Sample_Category", # La nuova colonna con i nomi delle categorie
    values_to = "Median_Expression_log10" # La nuova colonna con i valori mediani log10
  )

# --- 3. Crea il boxplot con facet_wrap ---
# Ora 'data_for_boxplot' dovrebbe contenere 'Sample_Category', 'Median_Expression_log10' e 'Overlap_class'
ggplot(data_for_boxplot, aes(x = Sample_Category, y = Median_Expression_log10, fill = Sample_Category)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.6, width = 0.2) +
  labs(
    title = "Distribuzione dei Livelli di Espressione Medi (log10) per Categoria e Overlap Class",
    x = "Categoria del Campione",
    y = "Mediana del Livello di Espressione (log10 + 1)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Overlap_class) # Dividi per la colonna Overlap_class


# Assicurati che dplyr, tidyr e ggplot2 siano caricati
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 0. Assicurati che i tuoi dataset iniziali siano caricati nell'ambiente ---
# (Come 'isoformsL1' e 'espressioneL1' come da tue indicazioni)

# --- 1. Assicurati che 'isoformsL1_con_mediane_raw' esista e sia corretto ---
# Questo è il dataset risultante dall'unione di isoformsL1 ed espressioneL1.
# Se non l'hai già creato, o se vuoi ricalcolarlo per sicurezza:
isoformsL1_con_mediane_raw <- isoformsL1 %>%
  left_join(espressioneL1, by = c("tss_tss_transcript_id" = "isoforms"))

# --- 2. Applica la normalizzazione log10 alle colonne delle mediane in 'isoformsL1_con_mediane_raw' ---
# Identifica le colonne delle categorie di campione che contengono i livelli mediani.
category_cols <- c("Adenoma", "Adenoma adjacent", "CRC", "CRC adjacent")

# Applica log10(x + 1) solo a queste colonne
isoformsL1_con_mediane_log10 <- isoformsL1_con_mediane_raw %>%
  mutate(across(all_of(category_cols), ~ log10(.x + 1), .names = "{.col}_log10"))

category_cols_log10 <- paste0(category_cols, "_log10")

# --- 3. Trasforma 'isoformsL1_con_mediane_log10' in formato "LONG" per il boxplot ---
# Questo creerà 'data_for_boxplot'
data_for_boxplot <- isoformsL1_con_mediane_log10 %>%
  pivot_longer(
    cols = all_of(category_cols_log10), # Impila le colonne normalizzate log10
    names_to = "Sample_Category",        # La nuova colonna per i nomi delle categorie
    values_to = "Median_Expression_log10" # La nuova colonna per i valori mediani log10
  ) %>%
  # --- NUOVO PASSAGGIO: Filtra i valori di mediana pari a 0 (o log10(1)) ---
  filter(Median_Expression_log10 > 0) # Esclude le righe dove la mediana log10 è esattamente 0
# Che corrisponde a un valore mediano originale di 0.

# --- 4. Crea il boxplot con facet_wrap ---
ggplot(data_for_boxplot, aes(x = Sample_Category, y = Median_Expression_log10, fill = Sample_Category)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.5, width = 0.2) +
  labs(
    title = "Distribuzione dei Livelli di Espressione Medi (log10) per Categoria e Overlap Class (Escluse le Mediane Zero)",
    x = "Categoria del Campione",
    y = "Mediana del Livello di Espressione (log10 + 1)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Overlap_class) 

#Estraggo nome geni
nome_geni <- unique(isoformsL1_con_mediane$tss_tss_gene_id)
print(nome_geni)

coherency<-read.delim("overlaps_same_opposite.txt")

#Unisco il nome del gene al numero del trascritto
library(dplyr)
head(coherency)
head(genexpression)
coherency <- left_join(coherency, 
                       genexpression, 
                       by = c("tss_tss_gene_id" = "gene_id"))

#Salvo come tabella excel
library(openxlsx)
write.xlsx(coherency, "nomegeni_coerenzaoverlap_codificaono.xlsx")

#Ottengo nome geni che hanno coerenza overlap
same_sense_genes <- coherency %>%
  filter(Overlap_class == "Same_sense") %>%
  dplyr::select(gene_name)
#Rimuovo gli NA
same_sense_genes_no_na <- coherency %>%
  filter(Overlap_class == "Same_sense") %>%
  dplyr::select(gene_name) %>%
  filter(!is.na(gene_name))

gene_names_vector <- same_sense_genes_no_na$gene_name
genes_for_copy_paste <- paste(gene_names_vector, collapse = ",")
print(genes_for_copy_paste)

#prendo quelli con overlap non coerente
opposite_sense_genes_vector <- coherency %>%
  filter(Overlap_class == "Opposite_sense") %>%
  pull(gene_name) %>%
  na.omit()
cat(opposite_sense_genes_vector, sep = "\n")


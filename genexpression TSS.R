library(tidyverse)
library(pheatmap)
df <- read_delim("overlaps_same_strand.csv", delim = ";")
head(df)
heatmap_data <- df %>%
  group_by(tss_seqnames, tss_strand) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = tss_strand, values_from = count, values_fill = 0)
heatmap_data <- heatmap_data %>%
  arrange(factor(tss_seqnames, levels = mixedsort(unique(tss_seqnames))))
library(readr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
# Installa il pacchetto (solo la prima volta)
install.packages("gtools")

# Carica la libreria
library(gtools)
heatmap_data <- heatmap_data %>%
  arrange(factor(tss_seqnames, levels = mixedsort(unique(tss_seqnames))))
Heatmap(mat,
        name = "Overlap count",
        col = colorRamp2(c(0, max(mat)), c("white", "firebrick")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "TSS strand",
        row_title = "Chromosome",
        heatmap_legend_param = list(title_position = "topcenter"),
        rect_gp = gpar(col = "black"))
##distanza tra TSS e LINE
library(readr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Leggi il file
df <- read_delim("overlaps_same_strand.csv", delim = ";")

# Calcola i centri
df <- df %>%
  mutate(
    tss_center = tss_TSS_original_coord,
    line1_center = (line1_start + line1_end) / 2,
    distance = line1_center - tss_center
  )

# Binning delle distanze (es. per heatmap)
df$distance_bin <- cut(df$distance, breaks = seq(-10000, 10000, by = 500))

# Tabella: count di overlap per distanza (binned) e strand
dist_table <- df %>%
  group_by(distance_bin, tss_strand) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = tss_strand, values_from = count, values_fill = 0)

# Matrice
dist_mat <- as.matrix(column_to_rownames(dist_table, var = "distance_bin"))

# Heatmap
Heatmap(dist_mat,
        name = "N overlaps",
        col = colorRamp2(c(0, max(dist_mat)), c("white", "navyblue")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        column_title = "TSS strand",
        row_title = "Distance TSS–LINE1 (binned)",
        heatmap_legend_param = list(title_position = "topcenter"),
        rect_gp = gpar(col = "black"))
# Calcola i centri
df <- df %>%
  mutate(
    tss_center = tss_TSS_original_coord,
    line1_center = (line1_start + line1_end) / 2,
    distance = line1_center - tss_center
  )

# Binning delle distanze (es. per heatmap)
df$distance_bin <- cut(df$distance, breaks = seq(-10000, 10000, by = 500))

# Tabella: count di overlap per distanza (binned) e strand
dist_table <- df %>%
  group_by(distance_bin, tss_strand) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = tss_strand, values_from = count, values_fill = 0)

# Matrice
dist_mat <- as.matrix(column_to_rownames(dist_table, var = "distance_bin"))

# Heatmap
Heatmap(dist_mat,
        name = "N overlaps",
        col = colorRamp2(c(0, max(dist_mat)), c("white", "navyblue")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        column_title = "TSS strand",
        row_title = "Distance TSS–LINE1 (binned)",
        heatmap_legend_param = list(title_position = "topcenter"),
        rect_gp = gpar(col = "black"))
##
line1_table <- df %>%
  group_by(tss_seqnames, line1_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = line1_name, values_from = count, values_fill = 0)

line1_mat <- as.matrix(column_to_rownames(line1_table, var = "tss_seqnames"))

Heatmap(line1_mat,
        name = "LINE1 counts",
        col = colorRamp2(c(0, max(line1_mat)), c("white", "darkgreen")),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "LINE1 subfamily",
        row_title = "Chromosome",
        show_column_names = FALSE,
        heatmap_legend_param = list(title_position = "topcenter"),
        rect_gp = gpar(col = "grey50"))
##
gene_line1 <- df %>%
  group_by(tss_tss_gene_id, line1_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = line1_name, values_from = count, values_fill = 0)

gene_line1_mat <- as.matrix(column_to_rownames(gene_line1, var = "tss_tss_gene_id"))

Heatmap(gene_line1_mat,
        name = "Overlap count",
        col = colorRamp2(c(0, max(gene_line1_mat)), c("white", "purple4")),
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        heatmap_legend_param = list(title_position = "topcenter"))
##frequenza LINE1 per cromosoma ma log10
library(readr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# Crea tabella count LINE1 per cromosoma
line1_table <- df %>% 
  group_by(tss_seqnames, line1_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = line1_name, values_from = count, values_fill = 0)

# Converti in matrice
line1_mat <- as.matrix(column_to_rownames(line1_table, var = "tss_seqnames"))

# Applic0 log10(count + 1)
line1_mat_log <- log10(line1_mat + 1)

# Heatmap log-scalata
Heatmap(line1_mat_log,
        name = "log10(count + 1)",
        col = colorRamp2(c(0, max(line1_mat_log)), c("white", "darkgreen")),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "LINE1 subfamily",
        row_title = "Chromosome",
        show_column_names = FALSE,
        heatmap_legend_param = list(title_position = "topcenter"),
        rect_gp = gpar(col = "grey50"))

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
library(biomaRt)
BiocManager::install("biomaRt")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mart
) %>%
  rename(tss_tss_gene_id = ensembl_gene_id,
         gene_name = hgnc_symbol)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mart
) %>%
  rename(tss_tss_gene_id = ensembl_gene_id,
         gene_name = hgnc_symbol)
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
table(genexpression$has_LINE1)
boxplot(median_adj ~ has_LINE1, data = genexpression,
        main = "Espressione genica con/senza LINE1",
        xlab = "Presenza di LINE1 nel TSS", ylab = "Espressione (median_adj)",
        col = c("lightblue", "lightgreen"))
genexpression_clean$log10_median_adj <- log10(genexpression_clean$median_adj + 1)
boxplot(log10_median_adj ~ has_LINE1, data = genexpression_clean,
        main = "Espressione log10(median_adj): LINE1 vs no LINE1",
        xlab = "Presenza LINE1 nel TSS", ylab = "log10(Espressione + 1)",
        col = c("lightgray", "lightblue"))

wilcox.test(log10_median_adj ~ has_LINE1, data = genexpression_clean)
###
genexpression_clean$log10_crc <- log10(genexpression_clean$median_crc + 1)
genexpression_clean$log10_adj <- log10(genexpression_clean$median_adj + 1)
boxplot(log10_crc ~ has_LINE1, data = genexpression_clean,
        main = "CRC (log10): LINE1 vs no LINE1",
        xlab = "LINE1 nel TSS", ylab = "log10(median_crc + 1)",
        col = c("gray", "skyblue"))

wilcox.test(log10_crc ~ has_LINE1, data = genexpression_clean)
boxplot(median_crc ~ has_LINE1, data = genexpression_clean,
        main = "Espressione CRC: LINE1 vs no LINE1",
        xlab = "LINE1 nel TSS", ylab = "Espressione CRC (mediana)",
        col = c("lightgray", "lightblue"))

# Assicurati di avere i pacchetti necessari installati e caricati
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(dplyr)
library(tidyr)
library(ggplot2)

# --- Carica il tuo dataframe 'genexpression_clean' ---
# Sostituisci "percorso/al/tuo/file.csv" con il percorso effettivo del tuo file.
# Se 'genexpression_clean' è già caricato nell'ambiente, puoi saltare questa riga.
# genexpression_clean <- read.csv("percorso/al/tuo/file.csv", header = TRUE, stringsAsFactors = FALSE)

# --- Inizio del codice di analisi reale ---

# Verifica il tipo della colonna has_LINE1 e convertila a logica se necessario
# Dallo screenshot, sembra che 'FALSE' sia già un booleano, ma è meglio essere sicuri.
# Se `as.logical()` dovesse fallire ancora, usa la conversione da stringa:
# genexpression_clean$has_LINE1 <- as.logical(genexpression_clean$has_LINE1)
# O, più robusto se sono stringhe "TRUE"/"FALSE":
genexpression_clean$has_LINE1 <- (genexpression_clean$has_LINE1 == TRUE)
# Se i tuoi dati sono stati caricati e "FALSE" è già un booleano, la riga sopra lo manterrà così com'è.
# Se sono stringhe "FALSE", li convertirà correttamente.

# 1. Prepara i dati per il plotting: trasforma da wide a long per le mediane
# Prenderemo le colonne 'median_adj' e 'median_crc' e le metteremo in un'unica colonna 'MedianExpression'
# con una nuova colonna 'SampleType' che indichi se è 'adj' o 'crc'.
data_for_plotting <- genexpression_clean %>%
  pivot_longer(
    cols = c(median_adj, median_crc),
    names_to = "SampleTypeRaw",
    values_to = "MedianExpression"
  ) %>%
  mutate(
    # Creiamo la colonna SampleType pulita (CRC vs ADJ)
    SampleType = case_when(
      SampleTypeRaw == "median_adj" ~ "ADJ",
      SampleTypeRaw == "median_crc" ~ "CRC",
      TRUE ~ NA_character_ # Per gestire eventuali altri valori, anche se non attesi qui
    ),
    # Applica la trasformazione log10(x+1). Dallo screenshot, sembra tu abbia 'log10_median_adj',
    # ma è più sicuro ricalcolarlo qui da 'MedianExpression' per coerenza tra ADJ e CRC.
    Log10MedianExpression = log10(MedianExpression + 1)
  ) %>%
  select(gene_id, gene_name, has_LINE1, SampleType, Log10MedianExpression)

# Converti 'has_LINE1' in un fattore per etichette migliori nei plot
data_for_plotting$has_LINE1_label <- factor(data_for_plotting$has_LINE1,
                                            levels = c(FALSE, TRUE),
                                            labels = c("No LINE-1 in TSS", "LINE-1 in TSS"))

# Ordina i livelli di SampleType per un ordine consistente nel plot
data_for_plotting$SampleType <- factor(data_for_plotting$SampleType, levels = c("ADJ", "CRC"))

# Verifica il dataframe preparato
head(data_for_plotting)
summary(data_for_plotting)


# --- Creazione dei Boxplot ---

# Boxplot 1: Geni SENZA LINE-1 nel TSS
plot_no_line1 <- data_for_plotting %>%
  filter(has_LINE1 == FALSE) %>% # Filtra per has_LINE1 è FALSE
  ggplot(aes(x = SampleType, y = Log10MedianExpression, fill = SampleType)) +
  geom_boxplot() +
  labs(
    title = "Without LINE-1 in TSS",
    x = "Sample Type",
    y = "Log10(Median Expression + 1)"
  ) +
  scale_fill_manual(values = c("ADJ" = "lightblue", "CRC" = "lightcoral")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  ylim(c(0,5.5))

print(plot_no_line1)

# Boxplot 2: Geni CON LINE-1 nel TSS
plot_with_line1 <- data_for_plotting %>%
  filter(has_LINE1 == TRUE) %>% # Filtra per has_LINE1 è TRUE
  ggplot(aes(x = SampleType, y = Log10MedianExpression, fill = SampleType)) +
  geom_boxplot() +
  labs(
    title = "With LINE-1 in TSS",
    x = "Sample Type",
    y = "Log10(Median Expression + 1)"
  ) +
  scale_fill_manual(values = c("ADJ" = "lightblue", "CRC" = "lightcoral")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  ylim(c(0,5.5))

print(plot_with_line1)

# Per visualizzarli fianco a fianco (richiede il pacchetto 'patchwork' o 'cowplot')
 #install.packages("patchwork") # Se non l'hai già
 library(patchwork)
 plot_no_line1 + plot_with_line1 +
   plot_layout(guides = 'collect') + # Raccoglie le legende se multiple
   plot_annotation(title = 'Comparison of Gene Expression Medians by LINE-1 Presence') &
   theme(plot.title = element_text(hjust = 0.5))
 
getwd()
overlaps<-read.delim("overlaps_same_opposite.txt") 
expression_levels<-read.delim("Expression_Level_TPM_CRC_Adenoma.txt")
##
library(dplyr)
ggplot(overlaps, aes(x = Overlap_class, y = expression_level)) + # SOSTITUISCI 'expression_level' se il nome è diverso
  geom_boxplot(notch = TRUE, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "blue", color = "blue") + # Aggiunge la media come un punto
  labs(
    title = "Mediana dei Livelli di Espressione per Coerenza Overlap",
    x = "Coerenza Overlap (Classe)",
    y = "Livello di Espressione" # SOSTITUISCI con l'etichetta appropriata per i tuoi dati
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
##
# Carica le librerie necessarie
library(dplyr)
library(tidyr) # Per la funzione pivot_longer()
library(ggplot2) # Per il boxplot

# --- Esegui l'unione dei due dataset ---
# 'overlaps' è il tuo dataset principale
# 'expression_levels' è il dataset con i livelli di espressione

overlaps_merged <- left_join(
  overlaps, # Usiamo il nome corretto del tuo dataset
  expression_levels,
  by = c("tss_tss_transcript_id" = "ID") # Unisci basandoti sugli ID delle isoforme
)

# Ora 'overlaps_merged' conterrà tutte le colonne di 'overlaps'
# più le colonne di 'expression_levels' (inclusi i tuoi livelli di espressione)
# per le righe che hanno un ID isoforma corrispondente.

# Puoi controllare le prime righe del nuovo dataset unito per vedere le nuove colonne
head(overlaps_merged)

# E anche i nomi delle colonne per identificare quelle dei livelli di espressione
colnames(overlaps_merged)

# --- Trasforma il dataset unito in formato "long" per i livelli di espressione ---
# Identifichiamo le colonne di espressione dal dataset originale `expression_levels`
# (tutte tranne la prima, che è 'ID')
expression_cols_to_pivot <- names(expression_levels)[-1]

# Ora, applica pivot_longer al dataset *unito* 'overlaps_merged'
overlaps_long <- overlaps_merged %>%
  pivot_longer(
    cols = all_of(expression_cols_to_pivot), # Seleziona le colonne di espressione da pivotare
    names_to = "Sample_ID",               # Nuova colonna per il nome del campione/condizione
    values_to = "Expression_Value"        # Nuova colonna per il valore di espressione
  )

# Controlla il nuovo dataset in formato long
head(overlaps_long)
colnames(overlaps_long)
#Boxplot
# --- Crea il Boxplot ---
ggplot(overlaps_long, aes(x = Overlap_class, y = Expression_Value)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "blue", color = "blue") + # Aggiunge la media come un punto
  labs(
    title = "Mediana dei Livelli di Espressione per Coerenza Overlap",
    x = "Coerenza Overlap (Classe)",
    y = "Livello di Espressione"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
merged_data <- merge(expression, lenghts, by = "isoform_id")


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
 getwd()
 "C:/Users/edino/OneDrive/Desktop/tesi/00_Input"
 expression<-read.delim("overlaps_same_opposite.txt")
 View(expression)
 View(genexpression)
 View(expression)
 View(isoforms)
 View(expression)
 lenghts<-read.delim("mart_export (1).txt")
 View(lenghts)
 lenghts<-read.delim("mart_export (1).txt")
 gc()
used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells  4859294 259.6    7932538 423.7  5728704 306.0
Vcells 36134465 275.7   67879505 517.9 67878581 517.9
 merged_data <- merge(expression, lenghts, by = "isoform_id")
Errore in fix.by(by.x, x) : 
  'by' deve specificare unicamente una colonna valida

 #aggiungo lenghts
 lenghts %>% filter(Transcript.stable.ID %in% expression$tss_tss_transcript_id)

 #aggiungo lenghts
 toadd <- lenghts %>% filter(Transcript.stable.ID %in% expression$tss_tss_transcript_id)
 View(toadd)
 gc()
used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells  4990681 266.6    7932538 423.7  7932538 423.7
Vcells 34700546 264.8   67879505 517.9 67878581 517.9
 toadd <- toadd[,c(3,5)]
 View(toadd)
 toadd<-toadd %>% rename("tss_tss_transcript_id"=Transcript.stable.ID)
 colnames(expression)
[1] "tss_seqnames"           "tss_start"              "tss_end"               
[4] "tss_width"              "tss_strand"             "tss_TSS_original_coord"
[7] "tss_tss_gene_id"        "tss_tss_transcript_id"  "line1_seqnames"        
[10] "line1_start"            "line1_end"              "line1_width"           
[13] "line1_strand"           "line1_name"             "line1_line_family"     
[16] "line1_line_type"        "Overlap_class"         
 colnames(toadd)
[1] "tss_tss_transcript_id"                     
[2] "Transcript.length..including.UTRs.and.CDS."
 expression <- merge(expression, toadd, by = "tss_tss_transcript_id", all.x = TRUE)
 gc()
used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells  5238858 279.8    7932538 423.7  7932538 423.7
Vcells 35336861 269.6   67879505 517.9 67878581 517.9
 gc()
used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells  5238727 279.8    7932538 423.7  7932538 423.7
Vcells 35336663 269.6   67879505 517.9 67878581 517.9
 gc()
used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells  5238727 279.8    7932538 423.7  7932538 423.7
Vcells 35336663 269.6   67879505 517.9 67878581 517.9
> gc()
used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells  5238766 279.8    7932538 423.7  7932538 423.7
Vcells 35341276 269.7   67879505 517.9 67878581 517.9
> gc()
used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells  5238784 279.8    7932538 423.7  7932538 423.7
Vcells 35341306 269.7   67879505 517.9 67878581 517.9
 #File con livelli di espressione
   livelli_di_espressione<-read.delim("Expression_Level_TPM_CRC_Adenoma.txt")
 View(livelli_di_espressione)
> gc()
used  (Mb) gc trigger   (Mb)  max used   (Mb)
Ncells  5383202 287.5    7932538  423.7   7932538  423.7
Vcells 90947673 693.9  141184380 1077.2 141182290 1077.2
> rm(data)
> rm(df)
> rm(dist_mat)
> rm(dist_table)
> rm(gene_line1)
> rm(gene_line1_mat)
> rm(toadd)
> library(tidyverse)
> patient_cols <- names(df)[grepl("^VT", names(df))]
> df_long <- df %>%
  +     select(ID, all_of(patient_cols)) %>% # Seleziona l'ID e tutte le colonne dei pazienti
  +     pivot_longer(
    +         cols = -ID, # Tutte le colonne tranne 'ID'
    +         names_to = "Patient_ID", # Crea una nuova colonna 'Patient_ID' con i nomi originali delle colonne
    +         values_to = "Expression_Level" # Crea una nuova colonna 'Expression_Level' con i valori
    +     ) %>%
  +     # Creare una nuova colonna 'Type' (Sano o Tumore) basandosi sull'ID del paziente
  +     mutate(
    +         Type = ifelse(grepl("K$", Patient_ID), "Tumore", "Sano") # Se finisce con 'K', è 'Tumore', altrimenti 'Sano'
    +     )
Errore in UseMethod("select") : 
  su un oggetto di classe "select" è stato usato un metodo non applicabile per "function"

> df_long <- df %>%
  +     dplyr::select(ID, all_of(patient_cols)) %>% # Usa esplicitamente select da dplyr
  +     pivot_longer(
    +         cols = -ID, # Tutte le colonne tranne 'ID'
    +         names_to = "Patient_ID", # Crea una nuova colonna 'Patient_ID' con i nomi originali delle colonne
    +         values_to = "Expression_Level" # Crea una nuova colonna 'Expression_Level' con i valori
    +     ) %>%
  +     # Creare una nuova colonna 'Type' (Sano o Tumore) basandosi sull'ID del paziente
  +     mutate(
    +         Type = ifelse(grepl("K$", Patient_ID), "Tumore", "Sano") # Se finisce con 'K', è 'Tumore', altrimenti 'Sano'
    +     )
Errore in UseMethod("select") : 
  su un oggetto di classe "select" è stato usato un metodo non applicabile per "function"

> df_long <- df %>%
  +     dplyr::select(ID, all_of(patient_cols)) %>% # Usa esplicitamente select da dplyr
  +     pivot_longer(
    +         cols = -ID, # Tutte le colonne tranne 'ID'
    +         names_to = "Patient_ID", # Crea una nuova colonna 'Patient_ID' con i nomi originali delle colonne
    +         values_to = "Expression_Level" # Crea una nuova colonna 'Expression_Level' con i valori
    +     ) %>%
  +     # Creare una nuova colonna 'Type' (Sano o Tumore) basandosi sull'ID del paziente
  +     mutate(
    +         Type = ifelse(grepl("K$", Patient_ID), "Tumore", "Sano") # Se finisce con 'K', è 'Tumore', altrimenti 'Sano'
    +     )
Errore in UseMethod("select") : 
  su un oggetto di classe "select" è stato usato un metodo non applicabile per "function"

> df_long <- df %>%
  +     dplyr::select(ID, all_of(patient_cols)) %>% # Usa esplicitamente select da dplyr
  +     pivot_longer(
    +         cols = -ID, # Tutte le colonne tranne 'ID'
    +         names_to = "Patient_ID", # Crea una nuova colonna 'Patient_ID' con i nomi originali delle colonne
    +         values_to = "Expression_Level" # Crea una nuova colonna 'Expression_Level' con i valori
    +     ) %>%
  +     # Creare una nuova colonna 'Type' (Sano o Tumore) basandosi sull'ID del paziente
  +     mutate(
    +         Type = ifelse(grepl("K$", Patient_ID), "Tumore", "Sano") # Se finisce con 'K', è 'Tumore', altrimenti 'Sano'
    +     )
Errore in UseMethod("select") : 
  su un oggetto di classe "select" è stato usato un metodo non applicabile per "function"

> class(select)
[1] "function"
> 
  > 
  > exists("select", mode = "function")
[1] TRUE
> 
  > 
  > library(tidyverse)
> ggplot(df_long, aes(x = Type, y = Expression_Level, fill = Type)) +
  +     geom_boxplot(outlier.shape = NA) +
  +     geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  +     scale_fill_manual(values = c("Tumore" = "firebrick", "Sano" = "forestgreen")) +
  +     labs(
    +         title = "Livelli di Espressione delle Isoforme Geniche",
    +         subtitle = "Confronto tra Pazienti Sani e con Tumore",
    +         x = "Tipo di Paziente",
    +         y = "Livello di Espressione"
    +     ) +
  +     theme_minimal() +
  +     theme(
    +         plot.title = element_text(hjust = 0.5),
    +         plot.subtitle = element_text(hjust = 0.5),
    +         legend.position = "none"
    +     )
Errore in ggplot(df_long, aes(x = Type, y = Expression_Level, fill = Type)) : 
  oggetto 'df_long' non trovato

> # --- 2. Preparare i dati per il boxplot ---
  > # Identificare le colonne dei pazienti
  > # Usiamo 'livelli_di_espressione' invece di 'df'
  > patient_cols <- names(livelli_di_espressione)[grepl("^VT", names(livelli_di_espressione))]
> 
  > # Trasformare il dataset da formato "wide" a formato "long"
  > # Usiamo 'livelli_di_espressione' invece di 'df'
  > df_long <- livelli_di_espressione %>%
  +     dplyr::select(ID, all_of(patient_cols)) %>% # Usa esplicitamente select da dplyr
  +     pivot_longer(
    +         cols = -ID, # Tutte le colonne tranne 'ID'
    +         names_to = "Patient_ID", # Crea una nuova colonna 'Patient_ID' con i nomi originali delle colonne
    +         values_to = "Expression_Level" # Crea una nuova colonna 'Expression_Level' con i valori
    +     ) %>%
  +     # Creare una nuova colonna 'Type' (Sano o Tumore) basandosi sull'ID del paziente
  +     mutate(
    +         Type = ifelse(grepl("K$", Patient_ID), "Tumore", "Sano") # Se finisce con 'K', è 'Tumore', altrimenti 'Sano'
    +     )
install.packages("readxl")
library("readxl")
tissues <- read_excel("metadata.xlsx")
livelli_di_espressione <- livelli_di_espressione %>%
     select(-VT437S, -VT448K, -VT462S)
livelli_di_espressione<-t(livelli_di_espressione)


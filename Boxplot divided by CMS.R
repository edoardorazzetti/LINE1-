install.packages("readxl")
library("readxl")
library("tidyverse")
getwd()
setwd("C:/Users/edino/OneDrive/Desktop/tesi/00_Input")
#metadata generale
metadata<-read_excel("metadata.xlsx")
#carico livelli di espressione relativi a met/L1
metadata_met<-read.delim("expression.txt")
#unisco i due dataset
metadata <- metadata_met %>%
  inner_join(metadata %>% select(ID, CMS), by = "ID")
#metto in formato long
metadata <- metadata %>%
  pivot_longer(
    cols = !c(ID, CMS),
    names_to = "Gene", # Il nome della colonna per i nomi dei geni
    values_to = "Livello_Espressione" # Il nome della colonna per i valori di espressione
  )  
#boxplot
ggplot(dati_long, aes(x = Gene, 
                      y = Livello_Espressione, 
                      fill = CMS)) +
  
  # Aggiungi il box plot
  ggplot(metadata, aes(x = Gene, y = Expression_Levels, fill = CMS)) +
  geom_boxplot() +
  theme_minimal() +
  
  # RUOTA le etichette dell'asse X (dove si trovano i geni)
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # 45 gradi è diagonale
  )
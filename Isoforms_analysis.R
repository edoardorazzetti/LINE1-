#Import library
library(tidyverse)

setwd("C:/Users/edino/OneDrive/Desktop/tesi/00_Input")


metadata <- read.delim("metadata.txt", row.names = 1)
countdata<-read.delim("Expression_Level_TPM_CRC_Adenoma.txt", row.names = 1)
met <- c("ENST00000397752.8","ENST00000318493.11","ENST00000422097.1","ENST00000456159.1","ENST00000454623.1","ENST00000436117.2","ENST00000495962.1")
#Select MET isoforms
countdata<-countdata[met,]
#Check only what's present
countdata<-countdata[,row.names(metadata)]
metadata<-metadata[names(countdata),]

# data setting
data<-data.frame(metadata[3],t(countdata),check.names = F)
data=data %>% gather(isoforme, levelexpression, -Tissue) 
%>% group_by(Tissue, isoforme)
%>% summarise(med=median(levelexpression))
data$Tissue <- factor(data$Tissue, levels = c("Adenoma adjacent","CRC adjacent","Adenoma","CRC"))

# plot
ggplot(data, aes(x=Tissue, y=log10(levelexpression),col=Tissue))+
  geom_boxplot( outliers=F)+
  #geom_jitter(alpha=0.4,width = 0.1)+
  #geom_line(aes(group=isoforme),col="grey70")+
  scale_color_manual(values = c("orange1","red2","orange3","red4"))+
  #geom_jitter(size=2, alpha=0.1,width = 0.1)+
  geom_pwc(method = "wilcox.test",label="p.signif",hide.ns = T, vjust =0.5)+
  #scale_color_manual(values = c("orange","red4"))+
  theme_bw()+
  facet_wrap(~isoforme,nrow = 1,scales = "free_x")+labs(y="log10(median levels)")+theme(axis.text.x = element_text(angle = 45,hjust = 1))

# plot
ggplot(data, aes(x = Tissue, y = log10(levelexpression), col = Tissue)) +
  geom_boxplot(outlier.shape = NA) + # Removed outliers for a cleaner look (optional)
  # geom_jitter(alpha = 0.4, width = 0.1) + # This line was removed
  #geom_line(aes(group = isoforme), col = "grey70") +
  scale_color_manual(values = c("orange1", "red2", "orange3", "red4")) +
  #geom_jitter(size = 2, alpha = 0.1, width = 0.1) + # This line was also removed (commented out)
 geom_pwc(method = "wilcox.test", label = "p.signif")+ 
  # scale_color_manual(values = c("orange", "red4")) +
  theme_bw() +
  facet_wrap(~isoforme, nrow = 1, scales = "free_x") +
  labs(y = "log10(median levels)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#Divido i boxplot per ogni singola isoforma
library(ggplot2)
library(dplyr)
unique_tissues <- unique(data$Tissue)
unique_isoforms <- unique(data$Isoforme)
colnames(data)  # Lists all column names
head(data)      # Shows the first few rows of the data frame
str(data)       # Shows the structure of the data frame, including column names and types
unique_isoforms <- unique(data$Isoform)
unique_isoforms <- unique(data$isoform)
unique_isoforms <- unique(data$isoforme)
# Assuming the correct isoform column name is 'transcript_id'

unique_tissues <- unique(data$Tissue)
unique_isoforms <- unique(data$transcript_id) # Corrected column name

for (tissue in unique_tissues) {
  for (isoform in unique_isoforms) {
    subset_data <- data %>%
      filter(Tissue == tissue, transcript_id == isoform) # Corrected column name
    
    # ... rest of your plotting code ...
  }
}
for (tissue in unique_tissues) {
  for (isoform in unique_isoforms) {
    subset_data <- data %>%
      filter(Tissue == tissue, Isoform == isoform)
    
    if (nrow(subset_data) > 0) {
      p <- ggplot(subset_data, aes(y = Expression)) +
        geom_boxplot() +
        labs(
          title = paste("Tissue:", tissue, "\nIsoform:", isoform),
          y = "log10(median levels)"
        ) +
        theme_bw()
      
      print(p)
      # You can save each plot to a file if needed:
      # ggsave(filename = paste("boxplot_", gsub(" ", "_", tissue), "_", isoform, ".png", sep = ""), plot = p)
    } else {
      cat(paste("No data for Tissue:", tissue, "and Isoform:", isoform, "\n"))
    }
  }
}
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming your data frame is called 'data' and has columns like 'Gene', 'Isoform',
# and then columns for each tissue: 'CRC adjacent', 'CRC', 'Adenoma adjacent', 'Adenoma'

# 1. Reshape the data from wide to long format
data_long <- data %>%
  pivot_longer(
    cols = c("CRC adjacent", "CRC", "Adenoma adjacent", "Adenoma"),
    names_to = "Tissue",
    values_to = "Expression"
  )

# 2. Create a separate plot for each tissue and each isoform
unique_tissues <- unique(data_long$Tissue)
unique_isoforms <- unique(data_long$Isoform)

for (tissue in unique_tissues) {
  for (isoform in unique_isoforms) {
    subset_data <- data_long %>%
      filter(Tissue == tissue, Isoform == isoform)
    
    if (nrow(subset_data) > 0) {
      p <- ggplot(subset_data, aes(y = Expression)) +
        geom_boxplot() +
        labs(
          title = paste("Tissue:", tissue, "\nIsoform:", isoform),
          y = "Expression Level" # Adjust y-axis label if needed
        ) +
        theme_bw()
      
      print(p)
      # Optional: Save each plot
      # ggsave(filename = paste("boxplot_", gsub(" ", "_", tissue), "_", isoform, ".png", sep = ""), plot = p)
    } else {
      cat(paste("No data for Tissue:", tissue, "and Isoform:", isoform, "\n"))
    }
  }
}
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming your data frame is called 'data' and has columns like 'Gene', 'Isoform',
# and then columns for each tissue: 'CRC adjacent', 'CRC', 'Adenoma adjacent', 'Adenoma'

# 1. Reshape the data from wide to long format
data_long <- data %>%
  pivot_longer(
    cols = c("CRC adjacent", "CRC", "Adenoma adjacent", "Adenoma"),
    names_to = "Tissue",
    values_to = "Expression"
  )

# 2. Create a separate plot for each tissue and each isoform
unique_tissues <- unique(data_long$Tissue)
unique_isoforms <- unique(data_long$Isoform)

for (tissue in unique_tissues) {
  for (isoform in unique_isoforms) {
    subset_data <- data_long %>%
      filter(Tissue == tissue, Isoform == isoform)
    
    if (nrow(subset_data) > 0) {
      p <- ggplot(subset_data, aes(y = Expression)) +
        geom_boxplot() +
        labs(
          title = paste("Tissue:", tissue, "\nIsoform:", isoform),
          y = "Expression Level" # Adjust y-axis label if needed
        ) +
        theme_bw()
      
      print(p)
      # Optional: Save each plot
      # ggsave(filename = paste("boxplot_", gsub(" ", "_", tissue), "_", isoform, ".png", sep = ""), plot = p)
    } else {
      cat(paste("No data for Tissue:", tissue, "and Isoform:", isoform, "\n"))
    }
  }
}
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming your data frame is called 'data' and has columns like 'Gene', 'Isoform',
# and then columns for each tissue: 'CRC adjacent', 'CRC', 'Adenoma adjacent', 'Adenoma'

# 1. Reshape the data from wide to long format
data_long <- data %>%
  pivot_longer(
    cols = c("CRC adjacent", "CRC", "Adenoma adjacent", "Adenoma"),
    names_to = "Tissue",
    values_to = "Expression"
  )

# 2. Create a separate plot for each tissue and each isoform
unique_tissues <- unique(data_long$Tissue)
unique_isoforms <- unique(data_long$Isoform)

for (tissue in unique_tissues) {
  for (isoform in unique_isoforms) {
    subset_data <- data_long %>%
      filter(Tissue == tissue, Isoform == isoform)
    
    if (nrow(subset_data) > 0) {
      p <- ggplot(subset_data, aes(y = Expression)) +
        geom_boxplot() +
        labs(
          title = paste("Tissue:", tissue, "\nIsoform:", isoform),
          y = "Expression Level" # Adjust y-axis label if needed
        ) +
        theme_bw()
      
      print(p)
      # Optional: Save each plot
      # ggsave(filename = paste("boxplot_", gsub(" ", "_", tissue), "_", isoform, ".png", sep = ""), plot = p)
    } else {
      cat(paste("No data for Tissue:", tissue, "and Isoform:", isoform, "\n"))
    }
  }
}
unique_tissues <- unique(data_long$Tissue)
head(data_long)
str(data_long)
colnames(data_long)
data_long <- data %>%
  pivot_longer(
    cols = c("CRC adjacent", "CRC", "Adenoma adjacent", "Adenoma"),
    names_to = "Tissue",
    values_to = "Expression"
  )
data$Tissue <- factor(data$Tissue, levels = c("Adenoma adjacent", "CRC adjacent", "Adenoma", "CRC"))
ggplot(data, aes(x = Tissue, y = log10(LevelExpression), col = Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = Isoforme), col = "grey70") +
  scale_color_manual(values = c("orange1", "red2", "orange3", "red4")) +
  theme_bw() +
  facet_wrap(~Isoforme, nrow = 1, scales = "free_x") +
  labs(y = "log10(median levels)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
install(ggpubr)

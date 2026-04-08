library("tidyverse")
library("readxl")
library(openxlsx)
metadata<-read_excel("C:/Users/edino/OneDrive/Desktop/tesi/00_Input/metadata.xlsx")
filtered_metadata<-metadata %>% filter(Tissue%in% c("Adenoma","Adenoma adjacent"))
write.xlsx(filtered_metadata, "metadata_for_paper.xlsx", sheetName = "Sheet1")
install.packages("officer")
install.packages("magrittr")
library(officer)
library(magrittr)
miss<-read_excel("file miss.xlsx")  
library(readxl)
library(magrittr)
miss <- read_excel("file miss.xlsx")
miss <- subset(miss, !is.na(`FPO NGS VAF`))
miss$FA  <- as.numeric(miss$`FA [%]`)
miss$NGS <- as.numeric(miss$`FPO NGS VAF`)
modello <- lm(NGS ~ FA, data = miss)
summary(modello)
correlazione <- cor(miss$FA, miss$NGS, use = "complete.obs")
fit <- lm(NGS ~ FA, data = miss)
summary(fit)
library(ggplot2)
ggplot(miss, aes(x = FA, y = NGS)) +
  geom_point(aes(color = Target), size = 3) +  
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Confronto dPCR vs NGS",
       subtitle = paste("R-squared =", round(summary(fit)$r.squared, 3)),
       x = "dPCR - FA [%]",
       y = "NGS - VAF [%]") +
  theme_minimal()

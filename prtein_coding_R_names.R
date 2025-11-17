file<-read_excel("nomegeni_coerenzaoverlap_codificaono.xlsx")
library(readxl)
install.packages("readxl")
file[file$Overlap_class == "same_sense", "gene_name"]
file[file$Overlap_class == "Same_sense" & !is.na(file$gene_name), "gene_name"]
cat(file[file$Overlap_class == "Same_sense" & !is.na(file$gene_name), "gene_name"], sep = "\n")
cat(file[file$Overlap_class == "Same_sense" & !is.na(file$gene_name), "gene_name", drop = TRUE], sep = "\n")
cat(file$gene_name[file$Overlap_class == "Opposite_sense" & !is.na(file$gene_name)], sep = "\n")
cat(file$gene_name[file$gene_type == "protein_coding" & !is.na(file$gene_name)], sep = "\n")

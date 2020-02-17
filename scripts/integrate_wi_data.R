#!/usr/bin/env Rscript
library(data.table)
fastqs <- data.table::fread("../data/WI_FASTQs.tsv", col.names = c("strain", "id", "library", "fq1", "fq2", "seq_folder"))
df <- data.table::fread("https://docs.google.com/spreadsheets/d/1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI/pub?output=tsv",
                        select = c("strain", "reference_strain", "isotype"))

# Filter strains that have not yet been classified as having a reference strain
df <- df[!is.na(reference_strain), ]

df <- fastqs[df, on = "strain"][,c("strain", "reference_strain", "isotype", "id", "library", "fq1", "fq2", "seq_folder")]

# Filter strains without sequence data
df <- df[!is.na(fq1),]

df[, reference_strain := reference_strain == 1]

data.table::fwrite(df, "../WI_sample_sheet.tsv", sep ="\t")
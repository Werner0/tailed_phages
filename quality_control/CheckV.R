#Script to check the quality of assembled virus genomes

library(data.table)
library(dplyr)
library(strex)

#Parse source files
dt50 <- fread("checkv_family.txt")
dt50$source <- "family"
dt70 <- fread("checkv_genus.txt")
dt70$source <- "genus"
dt95 <- fread("checkv_species.txt")
dt95$source <- "species"
dtraw <- fread("checkv_strain.txt")
dtraw$source <- "strain"
dt <- rbind(dt50, dt70, dt95, dtraw)
dt <- dt[,c("contig_id", "provirus", "miuvig_quality", "checkv_quality", "completeness", "source")]
dt$contig_id <- str_before_first(dt$contig_id, "_")

#Binary quality classification
dt[dt$miuvig_quality=="High-quality", quality:="high"]
dt[dt$checkv_quality=="High-quality", quality:="high"]
dt[dt$checkv_quality=="Complete", quality:="high"]
dt[dt$checkv_quality=="Medium-quality", quality:="low"]
dt[dt$checkv_quality=="Not-determined", quality:="low"]
dt[dt$checkv_quality=="Low-quality", quality:="low"]
dt[dt$miuvig_quality=="Genome-fragment", quality:="low"]

#Aggregate results
x <- count(dt, contig_id, provirus, source, quality)

#Export results
#fwrite(x, file = "results.csv")
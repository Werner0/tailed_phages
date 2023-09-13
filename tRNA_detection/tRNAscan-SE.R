#Script to parse tRNAscan-SE output
#Result set 1: Count by topology, source and type of tRNA
#Result set 2: Proportional abundance of detected tRNA (Coda)
#Result set 3: Ranked proportional abundance of detected tRNA

library(data.table)
library(tidyr)
library(dplyr)
library(compositions)

#Parse function
parse_function <- function (input) {
  input <- input[,c(1,5,6,9,11:13)]
  colnames(input) <- c("genome", "trna", "anticodon", "prim_score", "sec_score", "pseudo_flag", "source")
  input$prim_score <- as.double(input$prim_score)
  input$sec_score <- as.double(input$sec_score)
  output <- separate(input, genome, into = c("topology", "seq_number"), sep = "_")
}

#Parse input
dt_strain <- fread("strain_trna.txt", skip = 3)
dt_species <- fread("species_trna.txt", skip = 3)
dt_genus <- fread("genus_trna.txt", skip = 3)
dt_family <- fread("family_trna.txt", skip = 3)
dt_strain$source <- "strain"
dt_species$source <- "species"
dt_genus$source <- "genus"
dt_family$source <- "family"
dt_strain <- parse_function(dt_strain)
dt_species <- parse_function(dt_species)
dt_genus <- parse_function(dt_genus)
dt_family <- parse_function(dt_family)

#Combine results from strain, species, genus and family sets
x <- as.data.table(rbind(dt_strain, dt_species, dt_genus, dt_family))
count_x_1 <- count(x, topology, source, pseudo_flag)

#Export first set of results
#fwrite(count_x_1, file = "results_1.csv")

#Determine tRNA proportional abundances
count_x_2 <- count(x, topology, trna, source)
aggregate_x <- as.data.table(pivot_wider(count_x_2, values_from = "n", names_from = c("topology","source")))
coda_x <- t(as.data.table(acomp(t(aggregate_x[,-1]))))
colnames(coda_x) <- colnames(aggregate_x[,c(2:9)])
coda_x <- as.data.table(round(coda_x*100,2))
coda_x <- cbind(aggregate_x$trna,coda_x)
colnames(coda_x)[1] <- "trna"
coda_x <- coda_x[,c(1,5,4,3,2,9,8,7,6)]

#Export second set of results
#fwrite(coda_x, file = "results_2.csv")

#Rank tRNA proportional abundances
rank_x <- as.data.table(apply(-coda_x[,2:9], MARGIN=2, FUN=rank, ties="first"))
rank_x <- cbind(aggregate_x$trna,rank_x)
colnames(rank_x)[1] <- "trna"
rank_x <- rank_x[order(rank_x$Circular_strain),] #Replace rank_x$Circular_strain with applicable column

#Export third set of results
#fwrite(rank_x, file = "results_3.csv")

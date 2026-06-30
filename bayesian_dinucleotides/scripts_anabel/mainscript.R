#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

setwd("/users/genomics/anabel/genomes/")
library(seqinr)
library(dplyr)
library(tidyverse)
library(stringi)

source('../scripts/AllFunctions_Cluster.R')
organism = args[1]
seq = args[2]
window_size = as.numeric(args[3])
output_folder = args[4]


print(window_size)


if (!file.exists(output_folder)){
  dir.create(file.path(output_folder))
}

sequence <- read.fasta(paste0(organism, "/", seq))
# path cluster: /users/genomics/anabel/genomes/

dimer <- c("AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG", "CA", "CT", "CC", "CG", "GA", "GT", "GC", "GG")
hexamer <- paste0("CA", dimer, "TG")

n_windows <- ceiling(length(sequence[[1]]) / (window_size -4))

all_windows <- list()

start <- 0

for (n in 1:n_windows) {
  end <- start + window_size
  
  if(end > getLength(sequence) ){
    end <- getLength(sequence)
  }
  
  all_windows[[n]] <- ebox_ratio_window(seq = sequence[[1]], start = start, end = end, specie = organism,  chr = 1)
  
  start <- end - 4
}

bind_rows(all_windows) -> total2

write.table(total2, file = paste0(output_folder,"/",organism,"_",seq,".txt"), sep="\t",row.names = FALSE, quote=FALSE)



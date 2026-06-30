ebox <- function(seq) {
  # name_seq <- names(seq)
  seq <- seq %>% paste0(collapse = "")
  seq<- seq %>% toupper() # ho converteixo amb majúscules per evitar problemes
  
  eboxes <- str_extract_all(seq, regex_ebox) %>% unlist()
  
  if (length(eboxes) > 0) { # aquesta condició serveix per a que no doni error quan no hi ha E-boxes en una seqüència
    
    positions <- str_locate_all(seq, regex_ebox) %>% as.data.frame()
    
    data.frame(chr = args[2],
               # sequence = name_seq,
               start = positions$start,
               end = positions$end,
               ebox = sapply(eboxes, function(x){substring(x, first = 3, last = 8)}),
               internal = sapply(eboxes, function(x){substring(x, first = 5, last = 6)}),
               flank5 = sapply(eboxes, function(x){substring(x, first = 1, last = 2)}),
               flank3 = sapply(eboxes, function(x){substring(x, first = 9, last = 10)}))
  }
}

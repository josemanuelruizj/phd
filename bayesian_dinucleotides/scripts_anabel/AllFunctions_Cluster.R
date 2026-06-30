markov_hexamer <- function(hexamer, window) {
  nucleotide <- seqinr::count(window, 1, freq = T, alphabet = s2c("ACTG")) #absolute freq
  dinucleotide <- seqinr::count(window, 2, freq = T, alphabet = s2c("ACTG"))
  
  c1 <- matrix(dinucleotide, 4, 4, byrow = TRUE, dimnames = list(c("A", "C", "G", "T"), c("A", "C", "G", "T")))
  mt <- c1[,]/(c1[,1]+c1[,2]+c1[,3]+c1[,4])
  probstotals <- c()
  
  for (i in 1:6) {
    if(i==1){
      pos1 = substring(hexamer, 1, 1)
      prob1 = nucleotide[pos1]
      probstotals <- c(probstotals, prob1)
    }
    
    else{
      posi = substring(hexamer, i, i)
      posIprev = substring(hexamer,i-1, i-1)
      probstotals = c(probstotals, mt[posIprev, posi])
    }
  }
  
  prod(probstotals)
} 

ebox <- function(seq, chr, actual_window = 1, window_size) {
  # name_seq <- names(seq)
  seq <- seq %>% paste0(collapse = "")
  seq<- seq %>% toupper()

  eboxes <- stri_match_all_regex(seq,"(?=(CA[A,T,G,C]{2}TG))")[[1]][,2]
  
  if (!is.na(eboxes[1])) { 
    positions <- stri_locate_all_regex(seq, "(?=(CA[A,T,G,C]{2}TG))")[[1]] %>% as.data.frame()
    positions$start + actual_window * window_size

    data.frame(
      id = paste0(chr,"-", positions$start, "-", positions$start + 5),
      chr = chr,
      # sequence = name_seq,
      start = positions$start,
      end = positions$start + 5,
      ebox = eboxes) # sapply(eboxes, function(x){substring(x, first = 3, last = 8)})
    # internal = sapply(eboxes, function(x){substring(x, first = 5, last = 6)}),
    # flank5 = sapply(eboxes, function(x){substring(x, first = 1, last = 2)}), # tal como está puesto no descomentar, da información falsa, los dos primeros nucleotidos de la ebox
    # flank3 = sapply(eboxes, function(x){substring(x, first = 9, last = 10)}))
  }
}

count_hex <- function(window) { 
  windowchar <- window %>% paste0(collapse = "")
  windowchar<- windowchar %>% toupper()
  possible_hex <- stri_match_all_regex(windowchar, "(?=([ACTG]{6}))")[[1]][,2]
  x <- table(possible_hex)
  
  sum(x)
  # 
  # data.frame(n_different_hexamer = length(names(x)),
  #            n_hexamer = sum(x))
}

comparative <- function(window, masqued = F) {
  
  if (masqued == T) {
    window[grepl("[a-z]", window)] <- "N"
  } else {
    window <- toupper(window)
  }
  ebox_table <- ebox(seq = window, chr = "chr1")
  
  observed <- table(ebox_table$ebox)
  
  expected <- c()
  
  # we calculate the expected proportion of each ebox
  for (ebox in names(observed)) {
    expected[ebox] <- markov_hexamer(hexamer = ebox, window = window)
  }
  
  comparativa <- data.frame(observed = observed,
                            expected = expected)
  
  comparativa$observed.Var1 <- NULL
  colnames(comparativa) <- gsub(".Freq", "", colnames(comparativa))
  
  
  # Convert from observed to relative frequency
  n_hexamers <- count_hex(window)
  
  comparativa$observed <- comparativa$observed / n_hexamers
  
  comparativa$ratio <- comparativa$observed / comparativa$expected
  
  return(comparativa)
  
}


ebox_ratio_window <- function(seq, start, end, specie, chr) {
  
  table <- data.frame(id = paste0(specie, "-", chr,"-", start, "-", end),
                      specie = specie,
                      chr = chr,
                      # sequence = name_seq,
                      start = start,
                      end = end,
                      ebox = hexamer)
  
  window = seq[start:end] %>% toupper()
  windowrev <- rev(comp(window, forceToLower = FALSE))
  fullwindow <- c(window, rep("N",6), windowrev)
  window_colaps =  fullwindow %>% paste(collapse = "")
  
  
  
  # observed hexamer
  
  obs_hexamers <- stri_match_all_regex(window_colaps, "(?=(CA[A,T,G,C]{2}TG))")[[1]][,2] %>%
    table() %>% as.data.frame()
  
  if(nrow(obs_hexamers)>0){
    colnames(obs_hexamers) <- c("ebox", "observed")
    
    table <- table %>% 
      left_join(obs_hexamers, by = "ebox")
   
  } else {
    table$observed <- 0
  }
  
  # expected hexamer
  expected <- sapply(hexamer, markov_hexamer, window = fullwindow) * count_hex(fullwindow)
  expected <- data.frame(expected)
  expected$ebox <- rownames(expected)
  
  
  table <- table %>% 
    left_join(expected, by = "ebox")
  
  table$ratio <-  table$observed /  table$expected
  
  # observed dimer
  table$dimer <- substring(table$ebox, 3, 4)
  # browser()
  observed_dimer <- seqinr::count(fullwindow, 2, freq = T, alphabet = s2c("ACTG"))
  observed_dimer <- data.frame(observed_dimer)
  colnames(observed_dimer) <- gsub("Var1", "dimer", colnames(observed_dimer))
  
  table <- table %>% left_join(observed_dimer, by = "dimer")
  colnames(table) <- gsub("Freq", "observed_dimer", colnames(table))

  # expected dimer
  frecuencias_nucleotidos <- seqinr::count(fullwindow, 1, freq = T, alphabet = s2c("ACTG"))
  dimeros_split <- strsplit(table$dimer, "")
  
  expected_dimer <- data.frame(dimer = table$dimer,
                               expected_dimer = sapply(dimeros_split, function(x) {
                                 prod(frecuencias_nucleotidos[x])
                               }))
  
  table <- table %>% left_join(expected_dimer, by = "dimer")
  
  table$ratio_dimer <- table$observed_dimer / table$expected_dimer
  
  table$expected_dimer <- table$expected_dimer %>% round(2)
  table$ratio_dimer <- table$ratio_dimer %>% round(2)
  table$ratio <- table$ratio %>% round(2)
  
  
  table[is.na(table)] <- 0
  return(table)
}

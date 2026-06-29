library(rotl)
library(ape)

species <- c(
  'Otolemur garnettii',
  'Microcebus murinus',
  'Tarsius syrichta',
  'Saimiri boliviensis',
  'Callithrix jacchus',
  'Rhinopithecus roxellana',
  'Nasalis larvatus',
  'Papio hamadryas',
  'Papio anubis',
  'Macaca mulatta',
  'Macaca fascicularis',
  'Chlorocebus sabaeus',
  'Nomascus leucogenys',
  'Pongo pygmaeus abelii',
  'Gorilla gorilla' ,
  'Pan paniscus',
  'Pan troglodytes',
  'Homo sapiens'

)

resolved <- tnrs_match_names(species)
tree <- tol_induced_subtree(ott_ids = resolved$ott_id)

# ---- limpiar tip labels (quita ott#### y texto raro) ----
tree$tip.label <- gsub("_ott[0-9]+$", "", tree$tip.label)
tree$tip.label <- gsub(" \\(.*\\)", "", tree$tip.label)
tree$tip.label <- gsub("_", " ", tree$tip.label)

# ---- humano arriba: forzar que sea la primera tip ----
human <- "Homo sapiens"
desired_order <- c(human, setdiff(tree$tip.label, human))
tree2 <- rotateConstr(tree, desired_order)

plot(tree2, cex = 1.1)

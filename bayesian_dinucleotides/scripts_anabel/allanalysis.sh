declare -a StringArray=( 'amellifera' 'celegans' 'chicken' 'chimpanzee' 'dog' 'drosophmelanogaster' 'epsteinbarr' 'gorilla' 'human' 'macac' 'mouse' 'opossum' 'orangutan' 'pig' 'scerevisiae' 'xtropicalis' 'zebrafish')

module load R/4.1.2

for species in ${StringArray[@]}; do
  echo $species
  chrs=$(ls ../genomes/"$species"/ | grep '.fa$\|.fasta$')

  for chr in ${chrs[@]}; do
    echo $chr
    fname="/users/genomics/anabel/genomes/Results/"$species"_"$chr".txt"
    
    echo $fname
    if test -f "$fname"; then
       echo "$fname exists"
    else 
       echo "launching"
       sbatch -J "$species""$chr" \
         --partition=normal,long,short,bigmem \
         --mem=50GB \
         --time=6:00:00 \
         --nodes=1 \
         --export=ALL \
         --ntasks=1 \
         --output="/users/genomics/anabel/genomes/Results/""$species""$chr"".log.txt" \
         --wrap="Rscript  /projects_eg/projects/anabel/mainscript.R $species $chr 50000 Results"
    fi
  done
done

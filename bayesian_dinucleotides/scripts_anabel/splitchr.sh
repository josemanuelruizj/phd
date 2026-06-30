#! /bin/bash
#SBATCH --partition=normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1 --nodes=1
#SBATCH --output=/users/genomics/anabel/logs/split-%j.out
#SBATCH --error=/users/genomics/anabel/logs/split-%j.e

genome=$1
folder=$2

cd /users/genomics/anabel/genomes/${folder}

zcat ${genome} | awk '{
if (substr($0, 1, 1)==">") {filename=(substr($0,2)".fasta")}
print $0 >> filename
close(filename)
}'

echo "I am happy to inform you that your genome ${genomes} has been splitted"

#!/bin/bash
#
# 1. batch blastn to local copy of nt database 
#
for file in *.fa
do
{
echo "blasting file $file"
blastn -db /Applications/1_genomics/ncbi-blast-2.2.30+/db/nt -query $file -out ./Out$file.txt -num_threads 3
}
done
wait
#
# 2. Batch process and LCA in Megan
#
#
for file in *.fa
do
{
echo "import blastFile='/Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/blast_v3/Out$file.txt' fastaFile='/Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/blast_v3/$file' meganFile='/Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/blast_v3/megan2_$file.rma' minScore=200 topPercent=5 minSupport=1 blastFormat=BlastN">comands.txt
echo "export what=paths file= '/Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/blast_v3/Taxonomy_$file.txt'">>comands.txt
echo "quit">>comands.txt
echo "">>comands.txt
echo "Processing file $file in Megan"
/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub -g -E < comands.txt
}
done
#
# 3. Merging of all samples was done manually in megan afterwards
#
#
# 4. Blast against unite
#
##blastn -db /Users/ferninfm/Desktop/ANTONIA/0_databases/Unite_INDSC/UNITE_public_02.03.2015 -query /Users/ferninfm/Desktop/ANTONIA/1_Collapsed/2_Meganv2/Splidata_A032.fas -out /Users/ferninfm/Desktop/ANTONIA/1_Collapsed//3_UNITE_Blast/UNITE_Splidata_A032.fas.txt -num_threads 3
for file in ls /Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/1_Input_unique_filtered_fasta/*.fa
do
{
echo "blasting file $file to UNITE"
blastn -db /Users/ferninfm/Desktop/ANTONIA/0_databases/Unite_INDSC/UNITE_public_02.03.2015 -query $file -out /Users/ferninfm/Desktop/ANTONIA/2_Not_collapsed/4B_Blast_UNITE/UNITE_$(basename $file).txt -num_threads 3 -outfmt 10
}
done

#
#
# Bash to run ITSx
#
#
./ITSx -i /Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/4C_ITSx/All_renamed.fa -o /Users/Fernando/Desktop/ANTONIA/2_Not_collapsed/4C_ITSx/All_its.fa -t f --cpu 4 --fasta T --summary T --graphical T --preserve T --save_regions ITS1
#
#
# Bash to run Swarm
#
#
swarm -f -t 4 -w OTU_represen -i internal -l its1.otus.log -o its1.otus.fas -s stats.txt -u uclust.txt its1.fas
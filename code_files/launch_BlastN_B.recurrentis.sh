#!/usr/bin/env bash
#
#SBATCH -o BLASTN_hits.out -e BLASTN_hits.err
#SBATCH --qos=normal
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
#SBATCH --qos=normal
#SBATCH --job-name BLASTN_hits
#SBATCH -p common
#SBATCH --array=1-16%16

module load apptainer

wdir=$(echo "${HOME}/01_B.recurrentis")
nthread=$(echo "30")

cd $wdir
echo $wdir
echo $sample


sample=$(cat samples.ids | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 )

echo "4 Run BlastN analyses of hits" ; \
 echo "4.1 Launch BlastN in NCBI format" ; \
  blastn -task blastn-short \
   -db /pasteur/zeus/projets/p02/Hotpaleo/common_data/db/ncbi/nt/nt \
   -max_target_seqs 25 \
   -evalue 1E-05 \
   -num_threads ${nthread} \
   -outfmt 11  \
   -query 02_blast_hits/${sample}.fasta \
   -out 02_blast_hits/${sample}_hits_vs_nt.blastn ; \
 echo "4.2 BlastN in text format" ; \
  blast_formatter \
   -archive 02_blast_hits/${sample}_hits_vs_nt.blastn \
   -outfmt 0 \
   -out 02_blast_hits/${sample}_hits_vs_nt_txt.blastn \
   -num_descriptions 25 \
   -num_alignments 3 ; \
 echo "4.3 BlastN in table format" ; \
  blast_formatter \
   -archive 02_blast_hits/${sample}_hits_vs_nt.blastn \
   -max_target_seqs 1 \
   -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen length pident evalue score staxids stitle' \
   -out 02_blast_hits/${sample}_hits_vs_nt_table.blastn.tsv ; \
echo "Job for $sample is done"

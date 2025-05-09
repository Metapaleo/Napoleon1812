module load apptainer
module load tabix

wdir=$(echo "${HOME}/01_B.recurrentis")

gatk=$(echo "singularity exec --bind /pasteur --home $HOME:/home/$USER ${HOME}/VMs/singularity/broadinstitute-gatk-4.1.9.0.img gatk")

ref=$(echo "${HOME}/01_B.recurrentis/00_genome/ncbi_dataset/data/GCF_000019705.1/GCF_000019705.1_ASM1970v1_genomic.fna")

datemsa=$(echo "Apr 17")
daym=$(echo "$datemsa" | awk '{print $2}'); [[ ${#daym} -eq 1 ]] && datemsa=$(echo "$datemsa" | sed 's/ /  /')
msadir=ae/90a3c86daa95ca4fa238940ac9c1bf
bed=work/d1/018b393f3fcc7082b0f3b4cf952c58/genomic.gff_genes.bed
datetree=$(echo "Apr 17")
dayt=$(echo "$datetree" | awk '{print $2}'); [[ ${#dayt} -eq 1 ]] && datetree=$(echo "$datetree" | sed 's/ /  /')
treedir=aa/0b6bf461f2542e2781e0e33c2527cf
treemodel=$(ls -lh ${wdir}/work/${treedir}/ml_tree.raxml.bestTree | rev | cut -f 1 -d ' ' | rev | sed 's/ml_tree.raxml.bestTree/ml_tree.raxml.bestModel/')

cat 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt-ex_all_borrelia.txt | cut -f 1 > 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt_Borrelia_hits.ids
cat 02_blast_hits/Las_Gobas12.fasta | grep '>' | tr -d '>' | grep -Fwf 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt_Borrelia_hits.ids > 02_blast_hits/Las_Gobas12_Borrelia_hits.ids
cat 02_blast_hits/C11907.fasta | grep '>' | tr -d '>' | grep -Fwf 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt_Borrelia_hits.ids > 02_blast_hits/C11907_Borrelia_hits.ids
cat 02_blast_hits/YYY093A_hits_vs_nt_txt-ex_all_Borrelia.txt | cut -f 1 > 02_blast_hits/YYY093A_Borrelia_hits.ids
cat 02_blast_hits/YYY092B_mapped_vs_nt_txt-ex_all_borrelia.txt | cut -f 1 > 02_blast_hits/YYY092B_Borrelia_hits.ids

for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do samtools view -h -N ${sample}_Borrelia_hits.ids ../results/bams/${sample}_merged.rmdup.bam | samtools sort -o ${sample}_Borrelia_hits.bam ; done

for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do \
file=$(echo "${sample}_Borrelia_hits") ; \
samtools index ${file}.bam ; \
$gatk --java-options "-Xmx40G -Djava.io.tmpdir=$TMPDIR" HaplotypeCaller -R ${ref} -I ${file}.bam --min-base-quality-score 20 --sample-ploidy 1 --emit-ref-confidence BP_RESOLUTION --output-mode EMIT_VARIANTS_ONLY -O ${file}.haplotyper.vcf ; \
bgzip ${file}.haplotyper.vcf ; \
tabix -p vcf ${file}.haplotyper.vcf.gz ; \
$gatk --java-options "-Xmx40G -Djava.io.tmpdir=$TMPDIR" GenotypeGVCFs -R ${ref} --variant ${file}.haplotyper.vcf.gz --sample-ploidy 1 --include-non-variant-sites 	 -O ${file}.haplotyper.genotypegvcfs.vcf ; \
awk -F'\t' '$0 ~ /^#/ || (length($4) == 1 && length($5) == 1)' ${file}.haplotyper.genotypegvcfs.vcf > ${file}.haplotyper.genotypegvcfs_clean.vcf ; \
mv ${file}.haplotyper.genotypegvcfs_clean.vcf ${file}.haplotyper.genotypegvcfs.vcf ; \
bgzip ${file}.haplotyper.genotypegvcfs.vcf ; \
tabix -p vcf ${file}.haplotyper.genotypegvcfs.vcf.gz ; \
done

for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do
file=$(echo "${sample}_Borrelia_hits") ; \
bcftools index ${file}.haplotyper.genotypegvcfs.vcf.gz -f ; \
bcftools consensus -M N -a N -e "INFO/DP='.'" -f ${ref} ${file}.haplotyper.genotypegvcfs.vcf.gz > ${file}_genome_from_hits.fasta ; \
done

# Generate EPA directories
mkdir -p ${wdir}/01_epa_Borrelia_genus/01_MEGAN_reads
epadir=$(echo "${wdir}/01_epa_Borrelia_genus/01_MEGAN_reads")
mkdir ${epadir}/00_files
mkdir ${epadir}/01_genes
mkdir ${epadir}/02_tree
mkdir ${epadir}/03_results

# Copy files to the EPA directory eliminating annotation
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do cat ${wdir}/02_blast_hits/${sample}_Borrelia_hits_genome_from_hits.fasta | cut -f 1 -d ' ' > ${epadir}/00_files/${sample}_genome_from_hits.fasta ; done

# Copy the necessary files to generate the msa

cp ${wdir}/work/${msadir}/msa.fasta \
${wdir}/work/${msadir}/include.txt \
${wdir}/work/${msadir}/exclude.txt \
${wdir}/work/${msadir}/genes.txt \
${epadir}/01_genes/

cp ${bed} ${epadir}/00_files/

cp ${HOME}/scripts/command_used_to_msa_to_place.sh ${epadir}/01_genes/

# Generate the genes_filtered files
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do 
bedtools getfasta -fi ${epadir}/00_files/${sample}_genome_from_hits.fasta -bed ${bed} -fo stdout -s -fullHeader > ${epadir}/01_genes/${sample}_genes_filtered.fasta ; \
done

# Get the modified version of the script to generate another msa file to be placed:
cd ${epadir}/01_genes/ ; perl command_used_to_msa_to_place.sh ; cd ${wdir}

# Check that the msa and msa_to_place.fasta have identical lenghts
plgs=$(${HOME}/scripts/count_seq_lengths.sh ${epadir}/01_genes/msa_to_place.fasta | head -n 2 | tail -n 1)
msags=$(${HOME}/scripts/count_seq_lengths.sh ${epadir}/01_genes/msa.fasta | head -n 2 | tail -n 1)
echo $plgs $msags | awk '{if ($1 == $2) {print "Genome sizes are compatible"} else {print "Incompatible Genome sizes"}}'

# Copy the latest bootstrapped tree file and bestModel from the work directory (change grep 'Mar 24' by the right date)
cp ${wdir}/results/tree/ml_tree.raxml_with_support.nw ${epadir}/02_tree/
cp ${treemodel} ${epadir}/02_tree/

# Execute epa-ng
echo "source /opt/gensoft/adm/etc/profile.d/modules.sh ; module load epa-ng/0.3.8 ; cd ${epadir}/03_results; epa-ng --ref-msa ${epadir}/01_genes/msa.fasta --tree ${epadir}/02_tree/ml_tree.raxml_with_support.nw --query ${epadir}/01_genes/msa_to_place.fasta --model ${epadir}/02_tree/ml_tree.raxml.bestModel -T 20" > ${epadir}/03_results/launch_epa.sh

srun --cpus-per-task=20 --mem=100G -q ultrafast -p dedicated sh ${epadir}/03_results/launch_epa.sh

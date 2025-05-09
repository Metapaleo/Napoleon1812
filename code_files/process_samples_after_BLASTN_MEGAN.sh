# Set working directory
# wdir=$(echo "${HOME}/01_B.recurrentis")
wdir=$(echo "${HOME}/02_S.enterica/000_new_amphy_run_correcting092")
cd ${wdir}/02_blast_hits/

# Set directories and variables that are needed

datemsa=$(ls -l --time-style=+"%b %e" ${wdir}/results/msa/msa.fasta | awk '{print $6, $7}')
daym=$(echo "$datemsa" | awk '{print $2}'); [[ ${#daym} -eq 1 ]] && datemsa=$(echo "$datemsa" | sed 's/ /  /')
msadir=$(ls -l --time-style=+"%b %e" ${wdir}/work/*/*/msa.fasta | grep "${datemsa}" | grep -v '\->' | rev | cut -f 2,3 -d '/' | rev)
bed=$(ls -lh ${wdir}/work/*/*/*.gff_genes.bed | grep -v '\->' | rev | cut -f 1 -d ' ' | rev)
datetree=$(ls -l --time-style=+"%b %e" ${wdir}/results/tree/ml_tree.raxml_with_support.nw | awk '{print $6, $7}')
dayt=$(echo "$datetree" | awk '{print $2}'); [[ ${#dayt} -eq 1 ]] && datetree=$(echo "$datetree" | sed 's/ /  /')
treedir=$(ls -l --time-style=+"%b %e" ${wdir}/work/*/*/ml_tree.raxml_with_support.nw | grep "${datetree}" | grep -v '\->' | rev | cut -f 2,3 -d '/' | rev)
treemodel=$(ls -lh ${wdir}/work/${treedir}/ml_tree.raxml.bestTree | rev | cut -f 1 -d ' ' | rev | sed 's/ml_tree.raxml.bestTree/ml_tree.raxml.bestModel/')

gatk=$(echo "singularity exec --bind /pasteur --home $HOME:/home/$USER ${HOME}/VMs/singularity/broadinstitute-gatk-4.1.9.0.img gatk")

#ref=$(echo "${HOME}/01_B.recurrentis/00_genome/ncbi_dataset/data/GCF_000019705.1/GCF_000019705.1_ASM1970v1_genomic.fna")
ref=$(echo "${HOME}/02_S.enterica/00_genome/ncbi_dataset/data/GCF_000018385.1/GCF_000018385.1_ASM1838v1_genomic.fna")

# BLAST the reads that map on the targeted genome
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do
samtools view results/bams/${sample}_merged.rmdup.rescaled.bam | cut -f 1,10 | awk '{print ">"$1"\n"$2}' > 02_blast_hits/${sample}.fasta
done

# Parse BLAST output to get imported to MEGAN. Import to MEGAN and extract the IDs of the reads summarized at the targeted species. 
cat 02_blast_hits/Gobas-C11907/tmp_1_hits_vs_nt_txt.blastn | head -n 15 > 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt.blastn
for i in {1..15026}; do cat 02_blast_hits/Gobas-C11907/tmp_${i}_hits_vs_nt_txt.blastn | head -n -11 | tail -n +16 ; done >> 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt.blastn
cat 02_blast_hits/Gobas-C11907/tmp_1_hits_vs_nt_txt.blastn | tail -n 11 >> 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt.blastn

cat 02_blast_hits/YYY092B/tmp_1_hits_vs_nt_txt.blastn | head -n 15 > 02_blast_hits/YYY092B_mapped_vs_nt_txt.blastn
for i in {1..322}; do cat 02_blast_hits/YYY092B/tmp_${i}_hits_vs_nt_txt.blastn | head -n -11 | tail -n +16 ; done >> 02_blast_hits/YYY092B_mapped_vs_nt_txt.blastn
cat 02_blast_hits/YYY092B/tmp_1_hits_vs_nt_txt.blastn | tail -n 11 >> 02_blast_hits/YYY092B_mapped_vs_nt_txt.blastn


# Get ids of BLASTN+MEGAN hits
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do grep -Fwf all_together-ex_S.enterica.ids ${sample}.fasta | tr -d '>' > ${sample}_hits.ids ; done
# cat 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt-ex.txt | cut -f 1 > 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt_hits.ids
# cat 02_blast_hits/Las_Gobas12.fasta | grep '>' | tr -d '>' | grep -Fwf 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt_hits.ids > 02_blast_hits/Las_Gobas12_hits.ids
# cat 02_blast_hits/C11907.fasta | grep '>' | tr -d '>' | grep -Fwf 02_blast_hits/Gobas-C11907_mapped_vs_nt_txt_hits.ids > 02_blast_hits/C11907_hits.ids

# Get BAMs from the hit ids
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do samtools view -h -N ${sample}_hits.ids ../results/bams/${sample}_merged.rmdup.bam | samtools sort -o ${sample}_hits.bam ; done
# for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do samtools view -h -N ${sample}_hits.ids ../results/bams/${sample}_merged.rmdup.bam | samtools sort -o ${sample}_hits.bam ; done

# Get BAMs from the hit ids with ED ≤ 1) to see if YYY092B gets corrected
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do samtools view -H ../results/bams/${sample}_merged.rmdup.bam > header.tmp ; samtools view -h -N ${sample}_hits.ids ../results/bams/${sample}_merged.rmdup.bam | grep -E 'NM:i:0|NM:i:1' | cat header.tmp - | samtools view -b - | samtools sort -o ${sample}_hits_ED1.bam ; done
rm header.tmp

# load necessary modules
module load apptainer
module load tabix

# generate variable to run the software

gunzip -c ${ref}.gz > ${ref}
samtools faidx ${ref}
$gatk CreateSequenceDictionary -R ${ref}

# Run GATK on the selected sample
#for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do \
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do \
file=$(echo "${sample}_hits") ; \
#file=$(echo "${sample}_hits_ED1") ; \
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


# Generate fasta from VCF (any position mapped by at least 1 read is genotyped)
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do \
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do
file=$(echo "${sample}_hits") ; \
#file=$(echo "${sample}_hits_ED1") ; \
bcftools index ${file}.haplotyper.genotypegvcfs.vcf.gz -f ; \
#bcftools consensus -M N -a N -e "INFO/DP='.'" -f ${ref} ${file}.haplotyper.genotypegvcfs.vcf.gz > ${sample}_genome_from_hits_ED1.fasta ; \
bcftools consensus -M N -a N -e "INFO/DP='.'" -f ${ref} ${file}.haplotyper.genotypegvcfs.vcf.gz > ${sample}_genome_from_hits.fasta ; \
done

# Generate EPA directories
mkdir ${wdir}/01_epa/01_MEGAN_reads
epadir=$(echo "${wdir}/01_epa/01_MEGAN_reads")
#mkdir -p ${wdir}/01_epa_ED1/01_MEGAN_reads
#epadir=$(echo "${wdir}/01_epa_ED1/01_MEGAN_reads")
mkdir ${epadir}/00_files
mkdir ${epadir}/01_genes
mkdir ${epadir}/02_tree
mkdir ${epadir}/03_results

# Copy files to the EPA directory eliminating annotation
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do cat ${wdir}/02_blast_hits/${sample}_genome_from_hits.fasta | cut -f 1 -d ' ' > ${epadir}/00_files/${sample}_genome_from_hits.fasta ; done
#for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do cat ${wdir}/02_blast_hits/${sample}_genome_from_hits_ED1.fasta | cut -f 1 -d ' ' > ${epadir}/00_files/${sample}_genome_from_hits_ED1.fasta ; done
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do cat ${wdir}/02_blast_hits/${sample}_genome_from_hits.fasta | cut -f 1 -d ' ' > ${epadir}/00_files/${sample}_genome_from_hits.fasta ; done

# Copy the necessary files to generate the msa

cp ${wdir}/work/${msadir}/msa.fasta \
${wdir}/work/${msadir}/include.txt \
${wdir}/work/${msadir}/exclude.txt \
${wdir}/work/${msadir}/genes.txt \
${epadir}/01_genes/

cp ${bed} ${epadir}/00_files/

cp ${HOME}/scripts/command_used_to_msa_to_place.sh ${epadir}/01_genes/

# Generate the genes_filtered files
#for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do \
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do 
bedtools getfasta -fi ${epadir}/00_files/${sample}_genome_from_hits.fasta -bed ${bed} -fo stdout -s -fullHeader > ${epadir}/01_genes/${sample}_genes_filtered.fasta ; \
#bedtools getfasta -fi ${epadir}/00_files/${sample}_genome_from_hits_ED1.fasta -bed ${bed} -fo stdout -s -fullHeader > ${epadir}/01_genes/${sample}_ED1_genes_filtered.fasta ; \
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


####### Generate genotyping stats and genotype tables for heatmaps #################

# Count variable sites per strain within the genes used in the MSA (i.e., informative sites used in phylogenetic analyses)
## Get a list of positions used in the MSA:
cat ${wdir}/work/${msadir}/genes.txt | tr ':' '\t' | tr '-' '\t' | cut -f 1 -d '(' |  awk '{for(i = $2; i < $3; i++) {print $1 "\t" i}}' > positions_used_in_msa.txt

# Get samples.ids list For Borrelia recurrentis:
grep '>' ../results/msa/msa.fasta | tr -d '>' | grep -v 'reference' > samples.ids
cat ${wdir}/01_epa*/01_MEGAN_reads/01_genes/msa_to_place.fasta | grep '>' | tr -d '>' | sort -u >> samples.ids


# Get variant positions per strain that overlap with the positions in the MSA:
for sample in `cat samples.ids`
do
zcat ${wdir}/results/vcfs/${sample}*.vcf.gz | grep -v '#' | awk '$5 != "."' | cut -f 1,2,4,5 | grep -Fwf positions_used_in_msa.txt | awk -v sample=${sample} '{print sample"\t"$1"_"$2"_"$3"-"$4}' > ${sample}_var_sites.tsv
done

# Get variant positions after applying filters (BLASTN+MEGAN hits, ED ≤ 1 reads, etc.)
## For S. enterica: (Se-BMh: Salmonella enterica BLASTN+MEGAN hits)
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do zcat ../02_blast_hits/${sample}_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$5 != "."' | cut -f 1,2,4,5 | grep -Fwf positions_used_in_msa.txt | awk -v sample=${sample}_Se-BMh '{print sample"\t"$1"_"$2"_"$3"-"$4}' > ${sample}_Se-BMh_var_sites.tsv; done
## For S. enterica: (Se-BMh: Salmonella enterica BLASTN+MEGAN hits, only reads with ED1 ≤ 1)
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do zcat ../02_blast_hits/${sample}_hits_ED1.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$5 != "."' | cut -f 1,2,4,5 | grep -Fwf positions_used_in_msa.txt | awk -v sample=${sample}_Se-BMh_ED1 '{print sample"\t"$1"_"$2"_"$3"-"$4}' > ${sample}_Se-BMh_ED1_var_sites.tsv; done

## For B. recurrentis: (Br-BMh: Borrelia recurrentis BLASTN+MEGAN hits)
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do zcat ../02_blast_hits/${sample}_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$5 != "."' | cut -f 1,2,4,5 | grep -Fwf positions_used_in_msa.txt | awk -v sample=${sample}_Br-BMh '{print sample"\t"$1"_"$2"_"$3"-"$4}' > ${sample}_Br-BMh_var_sites.tsv; done
## For B. recurrentis: (Br-BMh: Borrelia recurrentis BLASTN+MEGAN hits, including also reads classified as Borrelia genus)
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do zcat ../02_blast_hits/${sample}_Borrelia_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$5 != "."' | cut -f 1,2,4,5 | grep -Fwf positions_used_in_msa.txt | awk -v sample=${sample}_Br-BMh_Bg '{print sample"\t"$1"_"$2"_"$3"-"$4}' > ${sample}_Br-BMh_inclBorrelia_genus_var_sites.tsv; done

# Get list of variant sites found in all samples
cat *var_sites.tsv | cut -f 2 | sort -u | cut -f 1-3 -d '_' | sort -u | rev | sed 's/_/@/' | rev | tr '@' '\t' > variant_sites.tsv

# Generate a genotype table for all samples
for sample in `cat samples.ids`
do
zcat ${wdir}/results/vcfs/${sample}_merged*.vcf.gz \
| grep -v '^#' \
| awk -F'\t' -v sample="$sample" '
    NR==FNR {
        sites[$1"_"$2]=1
        next
    }
    # Only single-base REF & ALT and in our site list
    length($4)==1 && length($5)==1 && ($1"_"$2) in sites {
        # extract the GT from the 10th column
        split($10, a, ":")
        gt = a[1]

        code = (gt == "."   ? 0 :      # no genotype => no read => code 0
               (gt == "0") ? 1 :      # genotype 0 => matches REF => code 1
               2)                     # anything else => has ALT => code 2

        print sample, $1"_"$2, code
    }
' variant_sites.tsv - \
> ${sample}_genotypes_in_var_sites.tsv
done

# Generate genotype tables for low-coverage samples with applied filters
## For S. enterica:
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do zcat ${wdir}/02_blast_hits/${sample}_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '^#' | awk -v sample=${sample}_Se-BMh -F'\t' 'NR==FNR{sites[$1"_"$2]; next} length($4)==1 && length($5)==1 && ($1"_"$2 in sites) {split($10,a,":"); gt=a[1]; code=(gt=="."?0:(gt=="0"?1:2)); print sample, $1"_"$2, code}' variant_sites.tsv - > ${sample}_Se-BMh_genotypes_in_var_sites.tsv; done

for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do zcat ${wdir}/02_blast_hits/${sample}_hits_ED1.haplotyper.genotypegvcfs.vcf.gz | grep -v '^#' | awk -v sample=${sample}_Se-BMh_ED1 -F'\t' 'NR==FNR{sites[$1"_"$2]; next} length($4)==1 && length($5)==1 && ($1"_"$2 in sites) {split($10,a,":"); gt=a[1]; code=(gt=="."?0:(gt=="0"?1:2)); print sample, $1"_"$2, code}' variant_sites.tsv - > ${sample}_Se-BMh_ED1_genotypes_in_var_sites.tsv; done


## For B.recurrentis:
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do zcat ../02_blast_hits/${sample}_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk -v sample=${sample}_Br-BMh -F'\t' 'NR==FNR{sites[$1"_"$2]; next} length($4)==1 && length($5)==1 && ($1"_"$2 in sites) {split($10,a,":"); gt=a[1]; code=(gt=="."?0:(gt=="0"?1:2)); print sample, $1"_"$2, code}' variant_sites.tsv - > ${sample}_Br-BMh_genotypes_in_var_sites.tsv; done

for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do zcat ../02_blast_hits/${sample}_Borrelia_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk -v sample=${sample}_Br-BMh_Bg -F'\t' 'NR==FNR{sites[$1"_"$2]; next} length($4)==1 && length($5)==1 && ($1"_"$2 in sites) {split($10,a,":"); gt=a[1]; code=(gt=="."?0:(gt=="0"?1:2)); print sample, $1"_"$2, code}' variant_sites.tsv - > ${sample}_Br-BMh_Bg_genotypes_in_var_sites.tsv; done


# Get the total list of samples ordered phylogenetically
module load gotree
gotree labels -i ml_tree.raxml.bestTree_rooted.nwk > strains_tree_order.txt
cat strains_tree_order.txt | grep -v 'reference' > samples_table.ids
cat *_genotypes_in_var_sites.tsv | cut -f 1 -d ' ' | sort -u | grep -vFwf strains_tree_order.txt >> samples_table.ids

# Generate genotype tables for each ancient strain, each table containing only the variable positions covered by at least one read in the concerned ancient strain. In the output table the columns (strains) are ordered by their phylogenetic positioning and rows (variant genomic positions) are clustered hierarchically
## For S. enterica:
for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do cat ${sample}_genotypes_in_var_sites.tsv | awk '$3 != 0' | cut -f 2 -d ' ' > ${sample}_mapped_var_pos.ids ; (echo "Sample CHROM_POS GT" ; cat *_genotypes_in_var_sites.tsv | grep -Fwf ${sample}_mapped_var_pos.ids) | tr ' ' '\t' | ~/scripts/traspose_as_otu.pl samples_table.ids - 1 2 | awk -v order="$(cat samples_table.ids)" 'BEGIN{FS=OFS="\t"; n=split(order,a,"\n")} NR==1{for(i=1;i<=NF;i++) c[$i]=i; printf "%s", $1; for(i=1;i<=n;i++) if(a[i] in c) printf "%s%s", OFS, a[i]; print ""; next} {printf "%s", $1; for(i=1;i<=n;i++) if(a[i] in c) printf "%s%s", OFS, $(c[a[i]]); print ""}' | ./cluster_variants.py > table_gt_var_sites_mapped_${sample}.tsv
done

## For B. recurrentis:
for sample in {C11907,Las_Gobas12,YYY092B,YYY093A}; do cat ${sample}_genotypes_in_var_sites.tsv | awk '$3 != 0' | cut -f 2 -d ' ' > ${sample}_mapped_var_pos.ids ; (echo "Sample CHROM_POS GT" ; cat *_genotypes_in_var_sites.tsv | grep -Fwf ${sample}_mapped_var_pos.ids) | tr ' ' '\t' | ./traspose_as_otu.pl samples_table.ids - 1 2 | awk -v order="$(cat samples_table.ids)" 'BEGIN{FS=OFS="\t"; n=split(order,a,"\n")} NR==1{for(i=1;i<=NF;i++) c[$i]=i; printf "%s", $1; for(i=1;i<=n;i++) if(a[i] in c) printf "%s%s", OFS, a[i]; print ""; next} {printf "%s", $1; for(i=1;i<=n;i++) if(a[i] in c) printf "%s%s", OFS, $(c[a[i]]); print ""}' | ./cluster_variants.py > table_gt_var_sites_mapped_${sample}.tsv
done


cat YYY0*_genotypes_in_var_sites.tsv | awk '$3 != "0"' | cut -f 2 -d ' ' | rev | sed 's/_/@/' | rev | tr '@' '\t' | sort -u > mapped_var_positions_YYY_samples.tsv

cat YYY0*_genotypes_in_var_sites.tsv | awk '$3 != "0"' | cut -f 2 -d ' ' | sort -u > mapped_var_positions_YYY_samples.ids







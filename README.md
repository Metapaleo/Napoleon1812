# Paratyphoid Fever and Relapsing Fever in 1812 Napoleon's Devastated Army
Scripts, code and results from R. Barbieri, et al.2025 study

# Salmonella enterica & Borrelia recurrentis Analysis Pipeline

This repository contains all shell and R commands used to:

1. Extract reads mapping to targeted species by BLASTN & MEGAN
2. Call variants with GATK
3. Generate consensus genomes
4. Build multiple sequence alignments (MSA) and phylogenetic trees
5. Place new genomes onto existing trees using EPA-ng
6. Produce genotype heatmaps for low-coverage samples

---

## Directory Structure

```
├── 00_genome/               # Reference genomes
├── 01_epa/                  # EPA-ng placement inputs & outputs
├── 02_blast_hits/           # FASTA & hit IDs from BLASTN
├── 02_blast_hits_Bg/        # BLAST+MEGAN hits at Borrelia genus level
├── 02_epa_Borrelia_genus/   # EPA-ng placement for Borrelia genus
├── results/                 # MSA, trees, VCFs, etc.
│   ├── msa/                 # Multiple sequence alignments
│   ├── tree/                # RAxML-NG & EPA-ng trees
│   └── vcfs/                # GATK & bcftools variants
├── work/                    # Nextflow intermediate files
└── genotype_tables/         # Genotype tables & heatmaps
```

---

## Shell Commands

### 1. Set Working Directory & Variables

```bash
wdir="$HOME/02_S.enterica/000_new_amphy_run_correcting092"
cd $wdir/02_blast_hits

# Identify MSA and tree working directories based on file modification dates
datemsa=$(date -r $wdir/results/msa/msa.fasta +"%b %e")
msadir=$(find $wdir/work -type f -name msa.fasta -newermt "$datemsa" \
         -printf '%h\n' | sed "s|$wdir/work/||")

datetree=$(date -r $wdir/results/tree/ml_tree.raxml_with_support.nw +"%b %e")
treedir=$(find $wdir/work -type f -name ml_tree.raxml_with_support.nw \
          -newermt "$datetree" -printf '%h\n' | sed "s|$wdir/work/||")
treemodel=$(basename $wdir/work/$treedir/ml_tree.raxml.bestTree \
           | sed 's/bestTree/bestModel/')

gatk="singularity exec --bind /pasteur --home $HOME:/home/$USER \
     $HOME/VMs/singularity/broadinstitute-gatk-4.1.9.0.img gatk"

ref="$HOME/02_S.enterica/00_genome/ncbi_dataset/.../genomic.fna"
```

### 2. Extract Reads for BLASTN

```bash
for sample in C11907 Las_Gobas12 YYY092B YYY093A; do
  samtools view results/bams/${sample}_merged.rmdup.rescaled.bam \
  | cut -f1,10 \
  | awk '{print ">"$1"\n"$2}' \
  > ${sample}.fasta
done
```

### 3. Parse BLAST Output & Extract Hit IDs

```bash
# Example for sample YYY092B
cat YYY092B/tmp_1_hits_vs_nt_txt.blastn | head -n15 > YYY092B_mapped.blastn
for i in {1..322}; do
  sed -n '16,$p;1,11p' YYY092B/tmp_${i}_hits_vs_nt_txt.blastn \
    >> YYY092B_mapped.blastn
done

# Extract read IDs from MEGAN summaries
grep -Fwf all_together-ex_S.enterica.ids YYY092B.fasta \
  | tr -d '>' > YYY092B_hits.ids
```

### 4. Subset BAM Files & Call Variants

```bash
# Subset by hit IDs
for sample in YYY087A YYY092B YYY095A YYY097B; do
  samtools view -h -N ${sample}_hits.ids \
    results/bams/${sample}_merged.rmdup.bam \
  | samtools sort -o ${sample}_hits.bam
done

# Keep only reads with edit distance ≤ 1
for sample in YYY087A YYY092B YYY095A YYY097B; do
  samtools view -H results/bams/${sample}_merged.rmdup.bam > header.tmp
  samtools view -h -N ${sample}_hits.ids \
    results/bams/${sample}_merged.rmdup.bam \
  | grep -E 'NM:i:0|NM:i:1' \
  | cat header.tmp - \
  | samtools view -b - \
  | samtools sort -o ${sample}_hits_ED1.bam
done
rm header.tmp
```

### 5. Initialize GATK Reference

```bash
module load apptainer tabix

gunzip -c ${ref}.gz > ${ref}
samtools faidx ${ref}
$gatk CreateSequenceDictionary -R ${ref}
```

### 6. GATK Variant Calling & Genotyping

```bash
for sample in C11907 Las_Gobas12 YYY092B YYY093A; do
  file=${sample}_hits
  samtools index ${file}.bam

  $gatk --java-options "-Xmx40G" HaplotypeCaller \
    -R ${ref} -I ${file}.bam \
    --min-base-quality-score 20 --sample-ploidy 1 \
    --emit-ref-confidence BP_RESOLUTION \
    -O ${file}.haplotyper.vcf

  bgzip ${file}.haplotyper.vcf
  tabix -p vcf ${file}.haplotyper.vcf.gz

  $gatk --java-options "-Xmx40G" GenotypeGVCFs \
    -R ${ref} --variant ${file}.haplotyper.vcf.gz \
    --sample-ploidy 1 --include-non-variant-sites \
    -O ${file}.haplotyper.genotypegvcfs.vcf

  # Keep only SNPs
  awk -F'\t' '$0 ~ /^#/ || (length($4)==1 && length($5)==1)' \
    ${file}.haplotyper.genotypegvcfs.vcf \
    > tmp.vcf
  mv tmp.vcf ${file}.haplotyper.genotypegvcfs.vcf
  bgzip ${file}.haplotyper.genotypegvcfs.vcf
  tabix -p vcf ${file}.haplotyper.genotypegvcfs.vcf.gz
done
```

### 7. Generate Consensus FASTA

```bash
for sample in C11907 Las_Gobas12 YYY092B YYY093A; do
  file=${sample}_hits
  bcftools index ${file}.haplotyper.genotypegvcfs.vcf.gz -f
  bcftools consensus -M N -a N -e "INFO/DP='.'" \
    -f ${ref} ${file}.haplotyper.genotypegvcfs.vcf.gz \
    > ${file}_genome_from_hits.fasta
done
```

### 8. EPA-ng Placement

```bash
epadir=$wdir/01_epa/01_MEGAN_reads
mkdir -p ${epadir}/{00_files,01_genes,02_tree,03_results}

# Copy consensus FASTA without annotations
for sample in C11907 Las_Gobas12 YYY092B YYY093A; do
  cut -f1 -d' ' ${file}_genome_from_hits.fasta \
    > ${epadir}/00_files/${sample}_genome_from_hits.fasta
done

# Copy MSA inputs
cp $wdir/work/$msadir/{msa.fasta,include.txt,exclude.txt,genes.txt} \
   ${epadir}/01_genes/
cp $wdir/$bed ${epadir}/00_files/

# Extract genes
for sample in C11907 Las_Gobas12 YYY092B YYY093A; do
  bedtools getfasta -fi ${epadir}/00_files/${sample}_genome_from_hits.fasta \
    -bed $bed -fo ${epadir}/01_genes/${sample}_genes_filtered.fasta \
    -s -fullHeader
done

# Generate MSA-to-place
cd ${epadir}/01_genes
perl $HOME/scripts/command_used_to_msa_to_place.sh
cd $wdir

# Copy tree & model
cp results/tree/ml_tree.raxml_with_support.nw ${epadir}/02_tree/
cp work/$treedir/$treemodel ${epadir}/02_tree/

# Launch EPA-ng
cat << EOF > ${epadir}/03_results/launch_epa.sh
module load epa-ng
epa-ng \
  --ref-msa ${epadir}/01_genes/msa.fasta \
  --tree ${epadir}/02_tree/ml_tree.raxml_with_support.nw \
  --query ${epadir}/01_genes/msa_to_place.fasta \
  --model ${epadir}/02_tree/ml_tree.raxml.bestModel \
  -T 20
EOF
srun --cpus-per-task=20 --mem=100G sh ${epadir}/03_results/launch_epa.sh
```

---

## R Commands: Genotype Heatmaps

```r
setwd("~/02_S.enterica/000_new_amphy_run_correcting092")
library(gplots)
library(RColorBrewer)

samples <- c("YYY087A","YYY092B","YYY095A","YYY097B")
for (i in samples) {
  tfile <- paste0("genotype_tables/table_gt_var_sites_mapped_", i, ".tsv")
  mat   <- as.matrix(read.table(tfile, sep="\t", header=TRUE, row.names=1))
  
  # Recode 1→3, 2→4
  cols <- grep(i, colnames(mat))
  mat[,cols] <- ifelse(mat[,cols]==1,3,
                 ifelse(mat[,cols]==2,4,mat[,cols]))
  
  # Choose palette for values {0,1,2,3,4}
  pal <- colorRampPalette(c("white","lightsteelblue1","tomato3",
                            "lightblue","darkred"))(100)
  
  pdf(paste0("genotype_tables/heatmap_var_sites_mapped_",i,".pdf"),
      height=10,width=10)
  heatmap.2(t(mat), dendrogram="none", Colv=FALSE, Rowv=FALSE,
            col=pal, trace="none", key=FALSE, scale="none",
            cexRow=1, cexCol=1, margins=c(6,6))
  dev.off()
}
```

---

## Borrelia Genus-Level Variation

Steps differ only in:

1. Generating hit IDs from `*mapped_vs_nt_txt-ex_all_borrelia.txt`
2. Subsetting BAM with `${sample}_Borrelia_hits.ids`
3. Naming outputs `${sample}_Borrelia_hits.*`
4. EPA directory: `01_epa_Borrelia_genus/`
5. Consensus FASTA from `${file}_Borrelia_hits`

All other commands remain the same.

---

### Contact

For questions, please open an issue in this repository or send an email to "adna [AT] pasteur.fr"

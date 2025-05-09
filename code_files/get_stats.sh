echo -e "Sample\tStat Name\tValue" > stats_table_raw.tsv

for sample in `cat samples.ids`; do
echo "De-duplicated mapped reads (counts)" 
ncounts=$(samtools view results/bams/${sample}_merged.rmdup.bam | wc -l)
echo -e "${sample}\tDe-duplicated mapped reads\t${ncounts}" >> stats_table_raw.tsv

echo "Coverage mean (X)"
echo "Average Coverage on mapped positions in BAM and VCF"
covmapbam=$(samtools depth results/bams/${sample}_merged.rmdup.bam | awk '{if ($1 == "NC_012125.1") {sum+=$3; count++}} END {print sum/count}')
echo -e "${sample}\tAverage Coverage on mapped positions in BAM\t${covmapbam}" >> stats_table_raw.tsv

covmapvcf=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk '{s+=$1 ; n++}END{print s/n}')
echo -e "${sample}\tAverage Coverage on mapped positions in VCF\t${covmapvcf}" >> stats_table_raw.tsv

echo "Average Coverage over the whole genome in BAM and VCF"
covWGbam=$(samtools coverage results/bams/${sample}_merged.rmdup.bam | grep 'NC_012125.1' | cut -f 7)
echo -e "${sample}\tAverage Coverage over the whole genome in BAM\t${covWGbam}" >> stats_table_raw.tsv

gl=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep '^#' | grep 'contig=' | cut -f 4 -d '=' | tr -d '>' | head -n 1)
covWGvcf=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk -v gl=${gl} '{s+=$1}END{print s/gl}')
echo -e "${sample}\tAverage Coverage over the whole genome in VCF\t${covWGvcf}" >> stats_table_raw.tsv

echo "Breath of Coverage in BAM and VCF"
BoCbam=$(samtools coverage results/bams/${sample}_merged.rmdup.bam | grep 'NC_012125.1' | cut -f 6)
BoCvcf=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | grep -c 'DP=' | awk -v gl=${gl} '{print $1/gl*100}')
echo -e "${sample}\tBreath of Coverage in BAM\t${BoCbam}" >> stats_table_raw.tsv
echo -e "${sample}\tBreath of Coverage in VCF\t${BoCvcf}" >> stats_table_raw.tsv

echo "SNPs (counts) in VCF"
snp=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | awk -F'\t' '$5 != "."'  | wc -l)
echo -e "${sample}\tSNPs counts\t${snp}" >> stats_table_raw.tsv

echo "SNPs average depth (X) in VCF"
snpDP=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | awk -F'\t' '$5 != "."'  | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk '{s+=$1 ; n++}END{print s/n}')
echo -e "${sample}\tSNPs Average Depth\t${snpDP}" >> stats_table_raw.tsv

echo "SNPs >=4X "
snpDP8x=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | awk -F'\t' '$5 != "."'  | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk '$1 >= 4' | wc -l)
echo -e "${sample}\tSNPs with Depth ≥8\t${snpDP8x}" >> stats_table_raw.tsv

echo "SNPs <4X "
snpDPu7x=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | awk -F'\t' '$5 != "."'  | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk '$1 < 4' | wc -l)
echo -e "${sample}\tSNPs with Depth <8\t${snpDPu7x}" >> stats_table_raw.tsv

echo "Het sites >8X"
het=$(zcat results/vcfs/${sample}_*.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk -F'\t' '$5 != "."' | ./identify_heterozygous.sh | grep -cFw 'LikelyHeterozygous')
echo -e "${sample}\tHeterozygous sites in positions with DP > 8x\t${het}" >> stats_table_raw.tsv

echo "Get stats in MSA"
./faCount results/msa/msa.fasta | awk -v sample=$sample '{if ($1 == sample) {print sample"\tMSA length\t"$2"\n"sample"\tNumber of N-positions in MSA\t"$7"\n"sample"\tNon-N positions in MSA\t"$2-$7}}' >> stats_table_raw.tsv

done


for sample in {YYY087A,YYY092B,YYY095A,YYY097B}; do
echo "Unique mapped reads classified as S. enterica by BLASTN+MEGAN (Se-BMh) (counts)" 
ncounts=$(samtools view 02_blast_hits/${sample}_hits.bam | wc -l)
echo -e "${sample}_Se-BMh\tDe-duplicated mapped reads\t${ncounts}" >> stats_table_raw.tsv

echo "Unique mapped reads classified as S. enterica by BLASTN+MEGAN (counts) with ED ≤ 1" 
ncountsed1=$(samtools view 02_blast_hits/${sample}_hits_ED1.bam | wc -l)
echo -e "${sample}_Se-BMh_ED1\tDe-duplicated mapped reads\t${ncountsed1}" >> stats_table_raw.tsv

echo "Coverage mean (X)"
echo "Average Coverage on mapped positions in BAM and VCF from Se-BMh"
covmapbam=$(samtools depth 02_blast_hits/${sample}_hits.bam | awk '{if ($1 == "NC_012125.1") {sum+=$3; count++}} END {print sum/count}')
echo -e "${sample}_Se-BMh\tAverage Coverage on mapped positions in BAM\t${covmapbam}" >> stats_table_raw.tsv

covmapvcf=$(zcat 02_blast_hits/${sample}_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk '{s+=$1 ; n++}END{print s/n}')
echo -e "${sample}_Se-BMh\tAverage Coverage on mapped positions in VCF\t${covmapvcf}" >> stats_table_raw.tsv

echo "Average Coverage on mapped positions in BAM and VCF from Se-BMh with ED ≤ 1"
covmapbamed1=$(samtools depth 02_blast_hits/${sample}_hits_ED1.bam | awk '{if ($1 == "NC_012125.1") {sum+=$3; count++}} END {print sum/count}')
echo -e "${sample}_Se-BMh_ED1\tAverage Coverage on mapped positions in BAM\t${covmapbamed1}" >> stats_table_raw.tsv

covmapvcfed1=$(zcat 02_blast_hits/${sample}_hits_ED1.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk '{s+=$1 ; n++}END{print s/n}')
echo -e "${sample}_Se-BMh_ED1\tAverage Coverage on mapped positions in VCF\t${covmapvcfed1}" >> stats_table_raw.tsv

echo "Breath of Coverage in BAM and VCF reads classified as S. enterica by BLASTN+MEGAN (Se-BMh) "
BoCbam=$(samtools coverage 02_blast_hits/${sample}_hits.bam | tail -n 1 | cut -f 6)
BoCvcf=$(zcat 02_blast_hits/${sample}_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | grep -c 'DP=' | awk -v gl=${gl} '{print $1/gl*100}')
echo -e "${sample}_Se-BMh\tBreath of Coverage in BAM\t${BoCbam}" >> stats_table_raw.tsv
echo -e "${sample}_Se-BMh\tBreath of Coverage in VCF\t${BoCvcf}" >> stats_table_raw.tsv

echo "Breath of Coverage in BAM and VCF reads classified as S. enterica by BLASTN+MEGAN (Se-BMh) with ED ≤ 1"
BoCbamed1=$(samtools coverage 02_blast_hits/${sample}_hits_ED1.bam | tail -n 1 | cut -f 6)
BoCvcfed1=$(zcat 02_blast_hits/${sample}_hits_ED1.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | grep -c 'DP=' | awk -v gl=${gl} '{print $1/gl*100}')
echo -e "${sample}_Se-BMh_ED1\tBreath of Coverage in BAM\t${BoCbam}" >> stats_table_raw.tsv
echo -e "${sample}_Se-BMh_ED1\tBreath of Coverage in VCF\t${BoCvcf}" >> stats_table_raw.tsv

echo "SNPs (counts) in VCF"
snp=$(zcat 02_blast_hits/${sample}_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | awk -F'\t' '$5 != "."'  | wc -l)
snped1=$(zcat 02_blast_hits/${sample}_hits_ED1.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | awk -F'\t' '$5 != "."'  | wc -l)
echo -e "${sample}_Se-BMh\tSNPs counts\t${snp}" >> stats_table_raw.tsv
echo -e "${sample}_Se-BMh_ED1\tSNPs counts\t${snp}" >> stats_table_raw.tsv

echo "SNPs average depth (X) in VCF"
snpDP=$(zcat 02_blast_hits/${sample}_hits.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | awk -F'\t' '$5 != "."'  | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk '{s+=$1 ; n++}END{print s/n}')
snpDPed1=$(zcat 02_blast_hits/${sample}_hits_ED1.haplotyper.genotypegvcfs.vcf.gz | grep -v '#' | awk '$1 == "NC_012125.1"' | awk -F'\t' '$5 != "."'  | grep 'DP=' | tr -d '@' | sed 's/DP=/@DP=/' | cut -f 2 -d '@' | cut -f 1 -d ';' | cut -f 1 | cut -f 2 -d '=' | awk '{s+=$1 ; n++}END{print s/n}')
echo -e "${sample}_Se-BMh\tSNPs Average Depth\t${snpDP}" >> stats_table_raw.tsv
echo -e "${sample}_Se-BMh_ED1\tSNPs Average Depth\t${snpDP}" >> stats_table_raw.tsv

echo "Get stats in MSA"
./faCount 01_epa_first_submission/01_MEGAN_reads/01_genes/msa_to_place.fasta | awk -v sample=${sample} '{if ($1 == sample) {print sample"_Se-BMh\tMSA length\t"$2"\n"sample"_Se-BMh\tNumber of N-positions in MSA\t"$7"\n"sample"_Se-BMh\tNon-N positions in MSA\t"$2-$7}}' >> stats_table_raw.tsv

./faCount 01_epa_ED1/01_MEGAN_reads/01_genes/msa_to_place.fasta | awk -v sample=${sample}_ED1 '{if ($1 == sample) {print sample"_Se-BMh_ED1\tMSA length\t"$2"\n"sample"_Se-BMh_ED1\tNumber of N-positions in MSA\t"$7"\n"sample"_Se-BMh_ED1\tNon-N positions in MSA\t"$2-$7}}' >> stats_table_raw.tsv

done

cat stats_table_raw.tsv | cut -f 1 | sort -u | grep -v '^Sample' > samples_stats.ids

~/scripts/traspose_as_otu.pl samples_stats.ids stats_table_raw.tsv 1 2 > stats_table.tsv

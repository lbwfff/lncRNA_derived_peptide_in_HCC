#!/bin/bash
#SBATCH -p compute
#SBATCH -c 60
#SBATCH -t 24:00:00

##reads mapping

source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate rna

#rm rrna
index_file="/scratch/lb4489/bioindex/homo_rRNA"

for ID in SRR25115687 SRR25115688 SRR25115689;do #SRR25115686
    trim_galore --quality 20 -j 120 -length 20 ./raw/"$ID".fastq
done

for ID in SRR27466630 SRR27466631 SRR29424511 SRR29424512 SRR29424513 SRR29424516 SRR29424517 SRR29424519;do
	trim_galore --paired --quality 20 -j 120 --length 20 ./raw/"$ID"_1.fastq  ./raw/"$ID"_2.fastq 
done

trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135324.fastq ./raw/SRR8135325.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135326.fastq ./raw/SRR8135327.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135328.fastq ./raw/SRR8135329.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135330.fastq ./raw/SRR8135331.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135332.fastq ./raw/SRR8135333.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135334.fastq ./raw/SRR8135335.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135336.fastq ./raw/SRR8135337.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135338.fastq ./raw/SRR8135339.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135340.fastq ./raw/SRR8135341.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135342.fastq ./raw/SRR8135343.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135344.fastq ./raw/SRR8135345.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3
trim_galore --paired --quality 20 -j 60 --length 20 ./raw/SRR8135346.fastq ./raw/SRR8135347.fastq -a AGATCGGAAGAGCACACGTCT --three_prime_clip_R1 3

for ID in SRR8135325 SRR8135327 SRR8135329 SRR8135331 SRR8135333 SRR8135335 SRR8135337 SRR8135339 SRR8135341 SRR8135343 SRR8135345 SRR8135347;do
    bowtie -x "$index_file" --un ./derRNA/"$ID".derRNA.fq -q ./"$ID"_val_2.fq -p 60 --norc -S "$ID".rRNA.mapped.sam
done

for ID in SRR25115686 SRR25115687 SRR25115688 SRR25115689;do
    bowtie -x "$index_file" --un ./derRNA/"$ID".derRNA.fq -q ./"$ID"_trimmed.fq -p 60 --norc -S "$ID".rRNA.mapped.sam
done

for ID in SRR8135325 SRR8135327 SRR8135329 SRR8135331 SRR8135333 SRR8135335 SRR8135337 SRR8135339 SRR8135341 SRR8135343 SRR8135345 SRR8135347 SRR25115686 SRR25115687 SRR25115688 SRR25115689; do
	STAR --runThreadN 60 --outSAMattributes All --genomeDir /scratch/lb4489/project/liver_ribo/index/GRCh38 --readFilesIn ./derRNA/"$ID".derRNA.fq --outFilterType BySJout --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --outFileNamePrefix "$ID" --outReadsUnmapped None
done

##Ribo-seq QC

conda activate ribocode

for id in Hepato_1 Hepato_2 Hepato_3 Hepato_4 Hepato_5 HepG2_Ribo HepG2_TIARKO_Ribo LC001_normal_RPF LC001_tumor_RPF LC033_normal_RPF LC033_tumor_RPF LC034_normal_RPF LC034_tumor_RPF LC501_normal_RPF LC501_tumor_RPF LC502_normal_RPF LC502_tumor_RPF LC505_normal_RPF LC505_tumor_RPF LC506_normal_RPF LC506_tumor_RPF LC507_normal_RPF LC507_tumor_RPF LC508_normal_RPF LC508_tumor_RPF LC509_normal_RPF LC509_tumor_RPF SRR8135325 SRR8135327 SRR8135329 SRR8135331 SRR8135333 SRR8135335 SRR8135337 SRR8135339 SRR8135341 SRR8135343 SRR8135345 SRR8135347 SRR25115686 SRR25115687 SRR25115688 SRR25115689; do
	metaplots -a /scratch/lb4489/project/liver_ribo/index/ribocode/ribocode -m 27 -M 31 -r /scratch/lb4489/project/liver_ribo/mapping/rna/"$id"Aligned.toTranscriptome.out.bam -o "$id"
done

GTF_FILE="/scratch/lb4489/project/liver_ribo/index/GENCODE_NONCODE_merged.gtf"

for id in Hepato_1 Hepato_2 Hepato_3 Hepato_4 Hepato_5 HepG2_Ribo HepG2_TIARKO_Ribo LC001_normal_RPF LC001_tumor_RPF LC033_normal_RPF LC033_tumor_RPF LC034_normal_RPF LC034_tumor_RPF LC501_normal_RPF LC501_tumor_RPF LC502_normal_RPF LC502_tumor_RPF LC505_normal_RPF LC505_tumor_RPF LC506_normal_RPF LC506_tumor_RPF LC507_normal_RPF LC507_tumor_RPF LC508_normal_RPF LC508_tumor_RPF LC509_normal_RPF LC509_tumor_RPF SRR8135325 SRR8135327 SRR8135329 SRR8135331 SRR8135333 SRR8135335 SRR8135337 SRR8135339 SRR8135341 SRR8135343 SRR8135345 SRR8135347 SRR25115686 SRR25115687 SRR25115688 SRR25115689; do
    ribotish quality -b /scratch/lb4489/project/liver_ribo/mapping/ribotish/"$id"Aligned.sortedByCoord.out.bam -g "$GTF_FILE" -p 10 -l 26,31

done

conda activate ribocode

RiboCode -a /scratch/lb4489/project/liver_ribo/index/ribocode/ribocode -c ribocode_config.txt -l no -g -b -o ribocode_liver_sep -A TTG,CTG,GTG

INPUT_DIR="/scratch/lb4489/project/liver_ribo/mapping/ribotish/Hepato_1Aligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/Hepato_2Aligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/Hepato_3Aligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/Hepato_4Aligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/Hepato_5Aligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/HepG2_TIARKO_RiboAligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/LC001_tumor_RPFAligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/LC033_tumor_RPFAligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/LC034_normal_RPFAligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/LC501_normal_RPFAligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/LC501_tumor_RPFAligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/LC502_tumor_RPFAligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/LC505_tumor_RPFAligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/SRR25115689Aligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/SRR25115688Aligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/SRR25115687Aligned.sortedByCoord.out.bam,/scratch/lb4489/project/liver_ribo/mapping/ribotish/SRR25115686Aligned.sortedByCoord.out.bam"

ribotish predict -b "$INPUT_DIR" --seq --altcodons CTG,GTG,TTG  -g /scratch/lb4489/project/liver_ribo/index/GENCODE_NONCODE_merged.gtf -f /scratch/lb4489/project/liver_ribo/index/GRCh38.p14.genome.fa --framebest -p 20  -o ribotish_liver_sep.txt --minaalen 5 

ribotish predict -b "$INPUT_DIR" --seq --altcodons CTG,GTG,TTG  -g /scratch/lb4489/project/liver_ribo/index/GENCODE_NONCODE_merged.gtf -f /scratch/lb4489/project/liver_ribo/index/GRCh38.p14.genome.fa --longest -p 20  -o ribotish_liver_longest.txt --minaalen 5

conda activate ribotricer

input_file="config_remain.txt"

while IFS=$'\t' read -r sample_name file_name col3 read_lengths psite_offsets
do

sample_base=$(echo "$sample_name" | sed 's/Aligned.toTranscriptome.out//')
bam_file="/scratch/lb4489/project/liver_ribo/mapping/ribotish/${sample_base}Aligned.sortedByCoord.out.bam"

echo "Processing sample: $sample_base with BAM file: $bam_file"

ribotricer detect-orfs \
            --bam "$bam_file" \
            --ribotricer_index /scratch/lb4489/project/liver_ribo/index/ribotricer/ribotricer_candidate_orfs.tsv \
            --prefix "$sample_base" \
            --read_lengths "$read_lengths" \
	    --psite_offsets "$psite_offsets" --min_valid_codons_ratio  0.4
            
done < "$input_file"

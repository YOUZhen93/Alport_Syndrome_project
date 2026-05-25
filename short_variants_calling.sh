# Short variants calling pipeline 
# author Jiahui Zhang
# notes: those codes were used in Alport Syndrome cohort analysis
#!/bin/bash

# 
set -euo pipefail

# ------------------------------
# global variables for tools that were used
# ------------------------------
# 
TRIMGALORE="trim_galore"          
BWA="bwa"
SAMTOOLS="samtools"
PICARD="picard"                   
GATK="gatk"
ANNOVAR="table_annovar.pl"
BCFTOOLS="bcftools"
# global variables for ref genome (GRCh38.p14)
REF_FASTA="/path/to/hg38.genome.fa"

# dbsnp file with gnomAD, ExAC allele frequencies info (common snp with MAF > 0.001)
KNOWN_SITES_DBSNP="/path/to/dbsnp.vcf"
FILTER_DBSNP="/path/to/dbsnp_maf0001.vcf"

# ANNOVAR database
ANNOVAR_DB="/path/to/humandb"

# WES target panel regions
TARGET_REGIONS="wes_targets.bed"

# AnnotSV fetched database
annodb="/path/to/AnnotSV/"
export ANNOTSV=/path/to/AnnotSV

# in/out path
IN_DIR="/path/to/fastq"
OUT_DIR="/path/to/results"


# sample name
SAMPLE=$1
# ------------------------------

# timer starts
echo "========== Starting pipeline for $SAMPLE at $(date) =========="

# tmp path
TMP_DIR="${OUT_DIR}/tmp_${SAMPLE}"
mkdir -p "$TMP_DIR"

# -------------------------------------------------------------------
# TrimGalore QC
# -------------------------------------------------------------------
echo "Trimming adapters and low-quality bases with TrimGalore..."
$TRIMGALORE \
    --paired \
    --quality 20 \
    --length 50 \
    --output_dir "$OUT_DIR" \
    "${IN_DIR}/${SAMPLE}_1.fq.gz" \
    "${IN_DIR}/${SAMPLE}_2.fq.gz"

TRIMMED_R1="${OUT_DIR}/${SAMPLE}_1_val_1.fq.gz"
TRIMMED_R2="${OUT_DIR}/${SAMPLE}_2_val_2.fq.gz"

# -------------------------------------------------------------------
# BWA-MEM align to GRCh38.p14
# -------------------------------------------------------------------
echo "Aligning with BWA-MEM..."
$BWA mem \
    -t 8 \
    -M \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
    "$REF_FASTA" \
    "$TRIMMED_R1" \
    "$TRIMMED_R2" \
    > "$OUT_DIR/${SAMPLE}.sam"

# -------------------------------------------------------------------
# SAMtools sort 
# -------------------------------------------------------------------
echo "Sorting SAM file with SAMtools..."
$SAMTOOLS sort \
    -@ 8 \
    -o "$OUT_DIR/${SAMPLE}.sorted.bam" \
    "$OUT_DIR/${SAMPLE}.sam"
rm "$OUT_DIR/${SAMPLE}.sam"   

# -------------------------------------------------------------------
# Picard MarkDuplicates
# -------------------------------------------------------------------
echo "Marking duplicates with Picard..."
$PICARD MarkDuplicates \
    INPUT="$OUT_DIR/${SAMPLE}.sorted.bam" \
    OUTPUT="$OUT_DIR/${SAMPLE}.dedup.bam" \
    METRICS_FILE="$OUT_DIR/${SAMPLE}.dup_metrics.txt" \
    TMP_DIR="$TMP_DIR" \
    REMOVE_DUPLICATES=false \
    CREATE_INDEX=true


$SAMTOOLS index "$OUT_DIR/${SAMPLE}.dedup.bam"

# -------------------------------------------------------------------
# GATK HaplotypeCaller 
# -------------------------------------------------------------------
echo "Calling variants with GATK HaplotypeCaller..."
$GATK HaplotypeCaller \
    -R "$REF_FASTA" \
    -I "$OUT_DIR/${SAMPLE}.dedup.bam" \
    -O "$OUT_DIR/${SAMPLE}.raw_variants.vcf.gz" \
    --dbsnp "$KNOWN_SITES_DBSNP" \
    --standard-min-confidence-threshold-for-calling 30.0


# -------------------------------------------------------------------
# keep rare mutation
# -------------------------------------------------------------------
echo "Removing variants with high MAF..."
$BCFTOOLS isec -C "$OUT_DIR/${SAMPLE}.raw_variants.vcf.gz" "$FILTER_DBSNP" -O z -o $OUT_DIR/${SAMPLE}.rare_variants.vcf.gz


# -------------------------------------------------------------------
# ANNOVAR annotation: pre-fetched databases: avsnp150, dbnsfp47a, and dbscsnv11
# -------------------------------------------------------------------
echo "Annotating variants with ANNOVAR..."
perl "$ANNOVAR" \
    "$OUT_DIR/${SAMPLE}.rare_variants.vcf.gz" \
    "$ANNOVAR_DB" \
    --buildver hg38 \
    --out "$OUT_DIR/${SAMPLE}_annovar" \
    --remove \
    --protocol refGene,cytoBand,clinvar_20220320,gnomad211_exome,avsnp150,dbnsfp47a,dbscsnv11 \
    --operation g,r,f,f,f,f,f \
    --nastring . \
    --vcfinput \
    --otherinfo


# -------------------------------------------------------------------
# spliceai annotation
# -------------------------------------------------------------------
echo "Annotating splicing variants with spliceai..."
spliceai -I "$OUT_DIR/${SAMPLE}.rare_variants.vcf.gz" -O "${ID}.rare.spliceai.vcf" -R "${REF_FASTA}" -A grch38


# -------------------------------------------------------------------
# copy number variations by cnvkit
# -------------------------------------------------------------------
echo "Copy number calling with cnvkit..."
cnvkit.py coverage "$OUT_DIR/${SAMPLE}.dedup.bam" "$TARGET_REGIONS" -o "$OUT_DIR/${SAMPLE}.coverage.cnn"

# -------------------------------------------------------------------
# copy number variations annotation by AnnotSV
# -------------------------------------------------------------------
echo "Copy number annotations..."
AnnotSV -bedFile "$OUT_DIR/${SAMPLE}.coverage.cnn" -annotationsDir ${annodb} -outputDir "$OUT_DIR" -genomeBuild GRCh38 

# -------------------------------------------------------------------
# clean up
# -------------------------------------------------------------------
rm -rf "$TMP_DIR"

echo "========== All steps completed for $SAMPLE at $(date) =========="

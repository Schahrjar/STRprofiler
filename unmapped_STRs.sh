#!/bin/bash



# This script is to find novel STR loci absent from the hg38 reference genome assembly by
# screening unmapped sequence read from sequence data (BAM files).
# 
# Shahryar Alavi
# University of Isfahan, Isfahan, Iran
# UCL Institute of Neurology, London, UK
# January 2026



# Set working directory
WrkDir=/path/to/working/directory



# ===== ===== =====
# Indevidual-level
# unmapped reads processing
# ===== ===== =====
# Retreive and preprocess unmapped reads of exomes
for DirBAM in $(cat $WrkDir/dir-BAMs.txt); do
    # Retreive sample ID:
    sampleID=$(basename $DirBAM | cut -d "." -f 1)
    # Make required directories
    mkdir $WrkDir/UBAM/${sampleID}
    mkdir $WrkDir/UBAM/${sampleID}/pileups

    # Capture unmapped reads
    # Retreive unmapped reads (both segments of a pair are not mapped) from a BAM file:
    samtools view -f 12 -F 256 -b -o $WrkDir/UBAM/${sampleID}/pileups/${sampleID}_unmapped-pairs.bam -@ 8 $DirBAM
    # Unmapped segments where the mate is mapped:
    samtools view -f 4 -F 264 -b -o $WrkDir/UBAM/${sampleID}/pileups/${sampleID}_unmapped-mappedpair.bam -@ 8 $DirBAM
    # Mapped segments of the unmapped ones:
    samtools view -f 8 -F 260 -b -o $WrkDir/UBAM/${sampleID}/pileups/${sampleID}_mapped-unmappedpair.bam -@ 8 $DirBAM

    # Aggregate unmapped reads
    # Merge the BAM files:
    samtools merge -f $WrkDir/UBAM/${sampleID}/pileups/${sampleID}_unmapped-merged.bam -@ 2 $WrkDir/UBAM/${sampleID}/pileups/${sampleID}_unmapped-pairs.bam $WrkDir/UBAM/${sampleID}/pileups/${sampleID}_unmapped-mappedpair.bam $WrkDir/UBAM/${sampleID}/pileups/${sampleID}_mapped-unmappedpair.bam
    # Sort the merged BAM file:
    samtools sort -n $WrkDir/UBAM/${sampleID}/pileups/${sampleID}_unmapped-merged.bam -o $WrkDir/UBAM/${sampleID}/${sampleID}_unmapped-reads.bam -@ 4

    # Create FASTQ
    # Convert unmapped BAM to paired FASTQ files:
    samtools fastq -1 $WrkDir/UBAM/${sampleID}/${sampleID}_unmapped-reads_R1.fastq.gz -2 $WrkDir/UBAM/${sampleID}/${sampleID}_unmapped-reads_R2.fastq.gz -@ 4 $WrkDir/UBAM/${sampleID}/${sampleID}_unmapped-reads.bam

done
# ===== ===== =====



# ===== ===== =====
# Cohort-level
# unmapped reads processing
# ===== ===== =====
# Aggregating all unmapped reads FASTQ files
mkdir $WrkDir/cohort
# Forward strand sequence reads:
zcat $WrkDir/UBAM/*/*_unmapped-reads_trimmed_R1.fastq.gz > $WrkDir/cohort/unmapped-reads_cohort_R1.fastq
gzip $WrkDir/cohort/unmapped-reads_cohort_R1.fastq
# Reverse strand sequence reads:
zcat $WrkDir/UBAM/*/*_unmapped-reads_trimmed_R2.fastq.gz > $WrkDir/cohort/unmapped-reads_cohort_R2.fastq
gzip $WrkDir/cohort/unmapped-reads_cohort_R2.fastq

# a) Aligning cohort unmapped raw reads (paired FASTQ) to the T2T-CHM13 human genome assembly
mkdir $WrkDir/cohort/BWA
bwa mem -M -t 8 -R "@RG\tID:ABC\tSM:cohort-unmapped\tLB:ABCD\tPL:Illumina" $WrkDir/chm13v2.0.fa $WrkDir/cohort/unmapped-reads_cohort_R?.fastq.gz | samtools view -Sb - > $WrkDir/cohort/BWA/cohort_unmapped-reads_T2T-map.bam
# Sort BAM file
samtools sort $WrkDir/cohort/BWA/cohort_unmapped-reads_T2T-map.bam -o $WrkDir/cohort/BWA/cohort_unmapped-reads_T2T-map_sorted.bam -@ 8
# Index BAM file
samtools index $WrkDir/cohort/BWA/cohort_unmapped-reads_T2T-map_sorted.bam
# Remove excess BAM files
rm $WrkDir/cohort/BWA/cohort_unmapped-reads_T2T-map.bam

# b) de novo assembly of the cohort unmapped reads FASTQ data
megahit \
        -1 $WrkDir/cohort/unmapped-reads_cohort_R1.fastq.gz \
        -2 $WrkDir/cohort/unmapped-reads_cohort_R2.fastq.gz \
        --out-dir $WrkDir/cohort/MEGAHIT --out-prefix cohort_unmapped-reads
# Using CD-HIT for clustering high identity contigs
mkdir $WrkDir/cohort/CD-HIT
cd-hit-est -i $WrkDir/cohort/MEGAHIT/cohort_unmapped-reads.contigs.fa -o $WrkDir/cohort/CD-HIT/cohort_clustered-contigs.fa -c 0.95 -n 10
# Aligning clustered unmapped contigs against human genome assembly T2T-CHM13
mkdir $WrkDir/cohort/minimap2
minimap2 -x asm5 $WrkDir/chm13v2.0.fa $WrkDir/cohort/CD-HIT/cohort_clustered-contigs.fa > $WrkDir/cohort/minimap2/cohort_clustered-contigs_T2T-map.paf


# c) Screening for STRs in cohort-level contigs mapped to the hs1 T2T-CHM13 genome.
mkdir $WrkDir/cohort/TRF
cd $WrkDir/cohort/TRF
# Run TRF on the cohort-level clustered contigs
trf409.linux64 $WrkDir/cohort/CD-HIT/cohort_clustered-contigs.fa 2 7 7 80 10 50 500 -d
# Screen TRF ouput for STRs with this criteria:
    # $4 >= 5               # repeat length greater than 5
    # $6 >= 80              # motifs match % (for repeat sequence perfection)
    # $3 >= 2 && $3 <= 6    # motif length between 2 and 6 bp
    # Save the output in this format:
    # contig  start  end  motif_length  repeats  percent_match  motif
awk '
        /^Sequence:/ { contig=$2 }
        $1 ~ /^[0-9]/ &&
        $3 >= 2 && $3 <= 6 &&
        $4 >= 5 &&
        $6 >= 80 {
        print contig "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $14
        }
' $WrkDir/cohort/TRF/cohort_clustered-contigs.fa*.dat > $WrkDir/cohort/TRF/cohort_clustered-contigs_STR.tsv
# ===== ===== =====

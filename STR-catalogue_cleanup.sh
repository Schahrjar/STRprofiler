#!/bin/bash



# This script provides major steps for creating a STRs catalogue for reliable genotyping.
# It begins with merging BED files of multiple variant catalogues, followed by removing
# loci with low sequence read depth or low quality mapped sequence reads. Finally it
# scans the STRs for imperfect repeats where the reference genome has multiple interuptions.
# Removing loci which are not well-covered from sequence data and those loci with imperfect repeats
# reduces the rate of false genotyping and incorrect interpretations in later steps.
# 
# Shahryar Alavi
# University of Isfahan, Isfahan, Iran
# UCL Institute of Neurology, London, UK
# February 2025



# Define paths to tools and reference data
WrkDir=/path/to/working/directory
BEDtools=/path/to/bedtools2/bin
RefGenome=/path/to/hg38.fasta

# Merge BED files of STR catalogues
cat $WrkDir/known-clinical_EH.bed $WrkDir/Illumina-polymorphic.bed > $WrkDir/STRs_exonic_catalogue.bed

# Some STR loci are not well covered in exome datasets. Use samtools to calculate depth and mapping quality of each STR region on a randomly selected list of samples.
shuf -n 200 $WrkDir/path-to-BAMs_list.txt > $WrkDir/selected_BAMs.txt

# Define the headers of the output file:
echo -e "Chromosome\tStart\tEnd\tmotif\tGene\tMeanDepth\tMeanMAPQ" > $WrkDir/STRs-exonic_coverage.tsv

# calculate depth and mapping quality
while read BEDregion; do
    BEDcoordinate=$(echo "$BEDregion" | awk -F '\t' '{print $1 ":" $2 "-" $3}')
    echo -e "Calculating coverage of region:\n$BEDcoordinate"
    CoverageRegion=$(samtools coverage -r $BEDcoordinate --bam-list $WrkDir/selected_BAMs.txt | tail +2 | awk -F '\t' '{print $7 "\t" $9}')
    echo -e "$BEDregion\t$CoverageRegion" >> $WrkDir/STRs-exonic_coverage.tsv
done < $WrkDir/STRs_exonic_catalogue.bed
# Some STR loci have low MAPQ (mean MAPQ < 50). Exclude loci with low MAPQ and those with mean depth < 5 (equals to depth < 1,000 for 200 samples).
# Save the retained STRs to a file named STRs-exonic_coverage-filtered.bed (the header line should be removed).

# There might be some imperfect repeats. Revise the list of STRs by finding imperfect repeats:
$BEDtools/fastaFromBed -fi $RefGenome -bed $WrkDir/STRs-exonic_coverage_filtered.bed -bedOut > $WrkDir/STRs-exonic_coverage_filtered_fastaSeq.bed
# This return the reference sequence at each STR loci, which can be inspected for STRs harbouring interuptions at their repeat tracks.
# Save the rateined STRs in STRs-exonic_coverage_filter_fastaSeq_revised.bed file.

# Finally, clean the BED file of the STR catalogue.
awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' $WrkDir/STRs-exonic_coverage_filter_fastaSeq_revised.bed > $WrkDir/STRs-exonic_catalogue_cleaned.bed
# Using STRipy (https://stripy.org/expansionhunter-catalog-creator), this BED file can be converted to JSON format for use in tools such as ExpansionHunter.

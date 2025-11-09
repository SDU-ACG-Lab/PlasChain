#!/bin/bash

common_simulation_path="./simulation/sim1"
common_plasmid_path="./plasmid"
database_simulation_path="/home/fengshi/src/database/Mock1"
Ref_path="$common_plasmid_path/Refseq"
Home_path=/home/fengshi/src/bacteria


plasmid_path="$common_plasmid_path"
simulation_path="$common_simulation_path"

outfile="$Ref_path/all_Refer.fasta"

> "$outfile" 

for file in "$plasmid_path"/*.fasta;do
        cat "$file" >> "$outfile"
done

bwa index "$outfile"
bwa mem "$outfile" "$simulation_path"/sim1_R1.fastq "$simulation_path"/sim1_R2.fastq > "$Ref_path/sim1.sam"
samtools view -b "$Ref_path/sim1.sam" > "$Ref_path/sim1.bam"
samtools sort "$Ref_path/sim1.bam" -o "$Ref_path/sim1_sorted.bam"

bedtools genomecov -ibam "$Ref_path/sim1_sorted.bam" -bg > "$Ref_path/coverage.bedgraph"
awk '{covered_length[$1]+=$3-$2} END {for (ref in covered_length) print ref, covered_length[ref]}' "$Ref_path/coverage.bedgraph" > "$Ref_path/covered_length.txt"
awk 'NR==FNR{a[$1]=$2;next}{b=$1; if (b in a) print b"\t"($2/a[b])*100"\t"a[b]}' $Home_path/plasmid_size.txt "$Ref_path/covered_length.txt" > "$Ref_path/cov_percentage.tsv"
awk '$2>95{print $0}'  $Ref_path/cov_percentage.tsv > $Ref_path/gold_Ref.txt

rm $Ref_path/sim1.sam $Ref_path/sim1.bam $Ref_path/sim1_sorted.bam $Ref_path/coverage.bedgraph $Ref_path/covered_length.txt
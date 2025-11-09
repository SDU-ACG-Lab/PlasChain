#!/bin/bash

SIM_PATH="./simulation"
cat "$SIM_PATH"/sim1_chromosome_R1.fastq "$SIM_PATH"/sim1_plasmid_R1.fastq > "$SIM_PATH"/sim1_R1.fastq
cat "$SIM_PATH"/sim1_chromosome_R2.fastq "$SIM_PATH"/sim1_plasmid_R2.fastq > "$SIM_PATH"/sim1_R2.fastq

rm $SIM_PATH"/sim1_chromosome_R1.fastq" $SIM_PATH"/sim1_plasmid_R1.fastq" $SIM_PATH"/sim1_chromosome_R2.fastq" $SIM_PATH"/sim1_plasmid_R2.fastq"
#!/bin/bash
SIMULATION_INFO=./simulation_info.txt
inter_path="./inter"
CHROMOSOME_READ_NUM=$1
echo "chromosome_reads	$CHROMOSOME_READ_NUM" 
echo "chromosome_reads	$CHROMOSOME_READ_NUM" >> $SIMULATION_INFO


#根据指定数量模拟染色体的read
iss generate --cpus 32 --genomes chromosome/*.fasta --n_reads $CHROMOSOME_READ_NUM --model hiseq --output simulation/sim1_chromosome

#利用python脚本计算待模拟的质粒的read数和丰度
python utilty.py -r $CHROMOSOME_READ_NUM -ab simulation/sim1_chromosome_abundance.txt -as association.txt -cz chromosome_size.txt -pz plasmid_size.txt -cc "$inter_path/"ccn.txt -ua "$inter_path"/ua.txt -oa "$inter_path"/abundance.txt -or "$inter_path"/read_number_file

PLASMID_READ_NUM=$(awk '{print $2}' ./inter/read_number_file)

echo “即将模拟$PLASMID_READ_NUM reads\n”

#开始模拟质粒
iss generate --cpus 32 --genomes plasmid/*.cycle --abundance_file "$inter_path"/abundance.txt --n_reads $PLASMID_READ_NUM --model hiseq --output simulation/sim1_plasmid

#将染色体模拟结果和质粒模拟结果拼接起来
./merge.sh


echo "plasmid reads: $PLASMID_READ_NUM " >> $SIMULATION_INFO

chr_size=$(grep "chromosome total size" "$SIMULATION_INFO" | awk '{print $4}')
pla_size=$(grep "plasmid total size" "$SIMULATION_INFO" | awk '{print $4}')
chr_reads=$(grep "chromosome_reads" "$SIMULATION_INFO" | awk '{print $2}')
pla_reads=$(grep "plasmid reads" "$SIMULATION_INFO" | awk '{print $3}')
read_len=$(grep "read length" "$SIMULATION_INFO" | awk '{print $3}')

# 计算总长度和总数据量
total_size=$((chr_size + pla_size))
total_reads=$((chr_reads + pla_reads))
total_bases=$((total_reads * read_len))

# 计算测序深度
depth=$(echo "scale=2; $total_bases / $total_size" | bc)

echo "Sequencing depth: ${depth}×"
echo "Sequencing depth: ${depth}×" >> $SIMULATION_INFO

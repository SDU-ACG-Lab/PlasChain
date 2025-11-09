#!/bin/bash

READ_SIZE=126
SIMULATION_INFO=./simulation_info.txt
#清空原有文件
> chromosome_size.txt
> plasmid_size.txt

#计算染色体的长度
for file in ./chromosome/*.fasta;do
	#获取第一行的样品名
	chromosome_name=$(head -n 1 "$file"| awk -F '[>, ]' '{print $2}')
	
	#计算去除第一行后的字符长度和
	size=$(tail -n +2 "$file"| awk '{total += length($0)}END{print total }')
	#输出
	echo "$chromosome_name	$size" >> chromosome_size.txt
done

#计算质粒的长度
for file in ./plasmid/*.fasta;do
	#获取第一行的样品名
	plasmid_name=$(head -n 1 "$file"| awk -F '[>, ]' '{print $2}')
	
	#计算去除第一行后的字符长度和
	size=$(tail -n +2 "$file"| awk '{total += length($0)}END{print total }')
	#输出
	echo "$plasmid_name	$size" >> plasmid_size.txt
done

total_chr_size=$(awk '{sum+=$2}END{print sum}' chromosome_size.txt)

total_plasmid_size=$(awk '{sum+=$2}END{print sum}' plasmid_size.txt)
echo "read length	$READ_SIZE" >> $SIMULATION_INFO
echo "chromosome total size	$total_chr_size" >> $SIMULATION_INFO
echo "plasmid total size	$total_plasmid_size" >> $SIMULATION_INFO

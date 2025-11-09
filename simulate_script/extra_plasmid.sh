#!/bin/bash

SIMULATION_INFO=./simulation_info.txt

#删除上一次模拟的数据
rm ./chromosome/*.fasta
rm ./plasmid/*.fasta
rm ./plasmid/*.cycle

# 从源文件中剥离染色体和质粒，并构建基本联系
> association.txt
> extra_log.txt
> selected_fna.txt

for file in ./bacteria/GC*;do
	echo "$file" >> selected_fna.txt
	echo "Dealing with $file" >> extra_log.txt
	./extra_plasmid.awk "$file"
done

# 将染色体和质粒分别放在不同文件夹
mv *.plasmid.fasta ./plasmid
mv *.chromosome.fasta ./chromosome

# 随机选择短质粒加入质粒文件，并初始化它们的联系
./select_short_plasmid.sh $1

# 环化质粒并计算其长度(原始长度)

for file in ./plasmid/*.fasta;do
	awk 'NR>1' "$file" > "${file}.tmp" &&\
	cat "$file" "${file}.tmp" > "${file}.cycle"&&\
	rm "${file}.tmp"
done

./calculate_size.sh

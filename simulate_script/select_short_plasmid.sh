#!/bin/bash

# 文件夹路径
folder_path="./src/short_plasmid"

# 要选择的文件数
n=$1

# 获取文件夹中的所有文件列表
files=("$folder_path"/*)

# 获取文件总数
total_files=${#files[@]}

# 确保 n 不超过文件总数
if ((n > total_files)); then
    echo "Error: n is greater than the total number of files in the folder."
    exit 1
fi

# 随机选择 n 个文件
selected_files=($(shuf -e "${files[@]}" | head -n "$n"))

echo "random select ${n} short plasmids. "
echo "random select ${n} short plasmids. " >> ./simulation_info.txt

# 将质粒移动到指定文件夹，再补充association文件
for file in "${selected_files[@]}";do
	echo $file >> selected_fna.txt
	cp "$file" ./plasmid
	file_name=$(basename "$file")
	echo "${file_name%.plasmid.fasta}"$'\t'"null" >> ./association.txt      
done


#!/bin/bash

SIMULATION_INFO=./simulation_info.txt

# 定义源文件夹和目标文件夹
src_dir="./src/all_bacteria"
complement_src_dir="./src/refseq"

dst_dir="./bacteria"

# 检查参数数量
if [ $# -ne 2 ]; then
    echo "Usage2: $0 <num of bacteria(no matter it has plasmid or not)> <num of plasmids randomly selected>"
    exit 1
fi


# 删除目标目录中的所有文件
rm -f ./bacteria/*.fna

# 指定要复制的文件数量  
num_files=$1

files=("$src_dir"/*.fna)
exisit_files=${#files[@]}

if [ "$num_files" -gt $exisit_files ]; then
    echo "Randomly selected ${num_files} bacteria" > "${SIMULATION_INFO}"
    echo "Randomly selected ${exisit_files} bacteria from scapp reference"
    find "$src_dir" -type f -name "*.fna" | shuf | head -n $exisit_files | xargs -I {} cp {} "$dst_dir"

    # 计算剩余需要复制的文件数
    let remaining=num_files-exisit_files

    echo "Randomly selected $remaining bacteria from refseq"

    find "$complement_src_dir" -type f -name "*.fna" | shuf | head -n $remaining | xargs -I {} cp {} "$dst_dir"

else
    echo "Randomly selected ${num_files} bacteria" > "${SIMULATION_INFO}"
    echo "Randomly selected ${num_files} bacteria"
    find "$src_dir" -type f -name "*.fna" | shuf | head -n "$num_files" | xargs -I {} cp {} "$dst_dir"
fi


./extra_plasmid.sh $2

#!/usr/bin/env python 

import numpy as np
import math,argparse
import random

READ_SIZE = 126



def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'find gold standard plasmid through bedgragh.'
        )
    parser.add_argument('-r','--read',
     help='input reads of chromosome will be simulated',
     required=True, type=int
     )

    parser.add_argument('-ab','--abundance',
     help='input abundance file containing the abundance of chromosome simulated',
     required=True, type=str
     )
    parser.add_argument('-as', '--association',
      help=' input association file reveling the association between plasmid and chromosome ',
      required=True, type=str
     )
    parser.add_argument('-cz','--chromosome_size',
     help='chromosome size file',
     required=True, type=str
     )
    parser.add_argument('-pz','--plasmid_size',
     help='plasmid size file',
     required=True, type=str
     )
    parser.add_argument('-cc','--chromosome_copy_number',
     help='plasmid size file',
     required=True, type=str
     )
    parser.add_argument('-ua','--updated_association',
     help='plasmid size file',
     required=True, type=str
     )
    parser.add_argument('-oa','--output_abundance',
     help='output abundance file containing the abundance of plasmid will be simulated ',
     required=True, type=str
     )
    parser.add_argument('-or','--output_read',
     help='output read file containing the read of plasmid will be simulated ',
     required=True, type=str
     )
    return parser.parse_args()







#获得细菌DNA的复制数
def genome_copy_number(abundance_file, bacteria_size_file, output_dir,total_size):
    bacteria_copy_number = {}

    # 将abundance file按chromosome name 排序
    with open(abundance_file, 'r') as file:
        abundance_lines = [line.strip().split('\t') for line in file]
    abundance_lines.sort(key=lambda x: x[0])

    # 将chromosome size file按chromosome name 排序
    with open(bacteria_size_file, 'r') as file:
        bacteria_size_lines = [line.strip().split('\t') for line in file]
    bacteria_size_lines.sort(key=lambda x: x[0])

    # copy_number =  abundance * n_read * read_size / chromosome_size  
    for i in range(len(abundance_lines)):
        bacteria_copy_number[abundance_lines[i][0]] = [total_size * float(abundance_lines[i][1]) / int(bacteria_size_lines[i][1]), int(bacteria_size_lines[i][1])]
    with open(output_dir,'w') as f:
        for bacteria, copy_number in bacteria_copy_number.items():
            f.write(f"{bacteria} {copy_number[0]:.0f} {copy_number[1]}\n")

#获得质粒模拟的read数
def read_of_plasmid(updated_according_file,bacteria_copy_number_file, plasmid_size_file, abundance_file, read_file): 

    #构建质粒和细菌DNA的对应关系
    according_bacterias = {}
    with open(updated_according_file, 'r') as f:
        for line in f:
            plasmid, bacteria = line.strip().split()[:2]
            according_bacterias[plasmid] = bacteria
    
    #构建细菌DNA和copy number的关系
    bacteria_copy_number = {}

    with open(bacteria_copy_number_file, 'r') as f0:
        for line in f0:
            bacteria, copy_number = line.strip().split()[:2]
            bacteria_copy_number[bacteria] = int(copy_number)
    
    plasmid_size = {}
    with open(plasmid_size_file, 'r') as f2:
        for line in f2:
            plasmid, size = line.strip().split()[:2]
            plasmid_size[plasmid] = int(size)

    plasmid_read = {}

    plasmid_abundance = {}
    total_plasmid_size = 0


    for plasmid, bacteria in according_bacterias.items():
        plasmid_abundance[plasmid] = {'abundacne':0, 'copy':0}
        copy_number = bacteria_copy_number[bacteria] * calculate_p(plasmid_size[plasmid])
        plasmid_relatve_abundance = plasmid_size[plasmid] * copy_number
        plasmid_abundance[plasmid]['abundance'] = plasmid_relatve_abundance
        total_plasmid_size += plasmid_relatve_abundance

    simulate_read = total_plasmid_size / READ_SIZE
    
    for plasmid in according_bacterias.keys():
        plasmid_abundance[plasmid]['abundance']/=total_plasmid_size
        plasmid_abundance[plasmid]['copy'] = int(total_plasmid_size * plasmid_abundance[plasmid]['abundance']  / READ_SIZE)
    
    with open(abundance_file, 'w') as f3:
         for plasmid in plasmid_abundance.keys():
             f3.write(f"{plasmid}\t\t{plasmid_abundance[plasmid]['abundance']:.12f}\t\t{plasmid_abundance[plasmid]['copy']}\n")
    with open(read_file, 'w') as f4:
        f4.write(f"read_number  {simulate_read:.0f}")
    return simulate_read
    
#获得质粒和细菌DNA之间的对应关系
def according(bacteria_file,according_file,updated_according_file):
    # 读取bacteria文件中的所有DNA名称
    bacterias = {}
    with open(bacteria_file, 'r') as f:
        for line in f:
            bacteria, len = line.strip('\t').split()[:2]
            bacterias[bacteria] = len

    # 初始化一个空列表来存储处理后的according文件行
    updated_according = []

    # 读取according文件并处理每一行
    bacteria_list = list(bacterias.keys())
    with open(according_file, 'r') as acc_f:
        for line in acc_f:
            plasmid, bacteria = line.strip().split()
            if bacteria.lower() == 'null':
                # 当细菌DNA为null时，从bacteria_list中随机选择一个替换
                new_bacteria = random.choice(bacteria_list)
                updated_according.append(f"{plasmid}\t\t{new_bacteria}\n")
            else:
                # 若不为null，则保持不变
                updated_according.append(line)

    # 将处理后的数据写入新的文件
    with open(updated_according_file, 'w') as out_f:
        out_f.writelines(updated_according)

#计算拷贝数
def calculate_p(len):
    if(1000<= len <10000):
        p = math.log(len,10)/30
    elif(len<100000):
        p = math.log(len,10)/20
    elif(len<1000000):
        p = math.log(len,10)/10
    else:
        p = 1
    copy_number = np.random.geometric(p)
    return copy_number


def main():
    args = parse_user_input()
    total_size = args.read * READ_SIZE
    _abundance_file = args.abundance
    according_file = args.association
    bacteria_size_file = args.chromosome_size
    plasmid_size_file = args.plasmid_size
    bacteria_copy_number_file = args.chromosome_copy_number
    updated_according_file = args.updated_association
    abundance_file = args.output_abundance
    read_file = args.output_read
    genome_copy_number(_abundance_file, bacteria_size_file, bacteria_copy_number_file,total_size)
    according(bacteria_size_file,according_file,updated_according_file)
    read =  read_of_plasmid(updated_according_file,bacteria_copy_number_file, plasmid_size_file, abundance_file, read_file)
    print(f"simulated {read:.0f} reads\n")

if __name__=='__main__':
    main()




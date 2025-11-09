#!/usr/bin/awk -f

BEGIN {
    plasmid_name = ""
    chromosome_name = ""
    plasmid = 0
    association_file = "association.txt"
    plasmid_num=0
    chromosome_num=0
    logfile="extra_log.txt"
}

#匹配到标注的首行
/^>/ {
    if (tolower($0) ~ /plasmid/) {
        plasmid_name = $1
	plasmid_num = plasmid_num + 1
	print plasmid_name "detected;" >> logfile
	#gsub全局替换,sub仅替换第一个匹配到的
        # 去掉质粒名字前面的">"
	gsub(">", "", plasmid_name)
	#sub (regular expression, substitution string):
	#sub (regular expression, substitution string, target string)
       	# 获取第一个空格前面的部分作为质粒名字
        sub(" .*", "", plasmid_name) 
        plasmid_file = plasmid_name".plasmid.fasta"
	plasmid = 1
	print plasmid_name"\t"chromosome_name >> association_file
	print $0 > plasmid_file
        next
    }else{
    	chromosome_num = chromosome_num + 1
    	chromosome_name = $1
	print chromosome_name "detected;" >> logfile
	gsub(">", "", chromosome_name)  # 去掉质粒名字前面的">"
        sub(" .*", "", chromosome_name) # 获取第一个空格前面的部分作为质粒名字
       	chromosome_file = chromosome_name ".chromosome.fasta"
	plasmid = 0
        print $0 > chromosome_file
        next

    }


}
#将对应的碱基行加入对应的序列文件中
{
    if (plasmid_name != "" && plasmid == 1) {
        print $0 >> plasmid_file
    }
    if (chromosome_name !="" && plasmid ==0){
	    print $0 >> chromosome_file
	}
}
END{
    print plasmid_num " plasmids and "chromosome_num " chromosomes detected." >> logfile
}


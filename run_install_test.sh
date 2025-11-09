#!/bin/bash
if [ -z "$1" ]; then
    echo "Usage: $0 <platon db path>"
    exit 1
fi

echo "It will take a few minutes"

python scapp.py -g ./test/test.fastg -b ./test/test.bam -path ./test/test.paths -o ./test/out -db $1 > /dev/null 2>&1

if [[ "$(grep '7501' ./test/out/test.confident_cycs.fasta | awk -F '_' '{print $6}')" == "11.62325" && \
      "$(grep '6099' ./test/out/test.confident_cycs.fasta | awk -F '_' '{print $6}')" == "12.46750" && \
      "$(grep '1302' ./test/out/test.confident_cycs.fasta | awk -F '_' '{print $6}')" == "8.42703" ]]; then
    echo "Test passed"
    rm -rf ./test/out
else
    echo "Test failed"
fi

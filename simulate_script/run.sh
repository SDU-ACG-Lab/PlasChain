#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <num of bacteria> <num of supplementary short plasmid> <num of reads>"
    exit 1
fi

SIMULATION_INFO=./simulation_info.txt

[ -d "bacteria" ] || mkdir bacteria
[ -d "plasmid" ] || mkdir plasmid
[ -d "chromosome" ] || mkdir chromosome
[ -d "simulation" ] || mkdir simulation
[ -d "inter" ] || mkdir inter


start=$(date +%s)


./random_select.sh $1 $2

./simulate.sh $3

end=$(date +%s)

echo "Simulation completed in $((end - start)) seconds"
echo "Simulation info saved to $SIMULATION_INFO"
echo "Simulated reads are saved to ./simulation"

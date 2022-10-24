#!/bin/bash

SEEDS=(1541540285 1927685433 1801974995 4529662466 2113415321 1344297384 7713334490 1637529655 1087798589 1888785167)

INPUT=$1 #input base gnome, such as: /global/cfs/projectdirs/m1641/inputs/large_genomes/human_chr8/human_chr8.fa
OUTPUT=$2 #output of the reads sampled, such as: /global/cfs/projectdirs/m1641/inputs/large_genomes/human_chr8/reads_100x
BINARY=$3 #art executable, such as: /global/homes/p/pghosh/dbg/art_bin_MountRainier/art_illumina

for i in {1..10}; do

    PREFIX=chr_illumina_10x_part_$i
    seednum=$((i-1))
    $BINARY -p -m 1000 -s 10 -i $INPUT -l 100 -f 10 -nf 1 -na -rs ${SEEDS[$seednum]} -ss HS25 -o $OUTPUT/$PREFIX -d EAS139:136:FC706VJ:2:5:1000: > $OUTPUT/stats$i.log &

done

echo Waiting
wait

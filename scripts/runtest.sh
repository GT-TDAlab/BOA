#!/bin/bash

# Update me
indir=/project/jing/boa/StringBOA_LFS/art_short_paired_reads # short reads dir: genome_1.fa, genome_2.fa for paired input
output=/project/jing/boa/debug # output dir for storing partitioned reads
programdir=./ # where the executatbles are at
exe1=pbucketing-code-release # one of the two executables to compare with
exe2=pbucketing-varible-length # one of the two executables to compare with

runinstans () {
    program=$1
    mpiprocess=$2
    parts=$3
    genome=$4
    outdir=$5

    echo "mpiexec -n $mpiprocess $programdir/$program -f $indir/${genome}_1.fa -b 100000000 -r 100 -p $parts -o $outdir/${genome}+$parts+$mpiprocess+Zoltan+333937177 -t Zoltan -d 333937177"
    mpiexec -n $mpiprocess $programdir/$program -f $indir/${genome}_1.fa -b 100000000 -r 100 -p $parts -o $outdir/${genome}+$parts+$mpiprocess+Zoltan+333937177 -t Zoltan -d 333937177
    echo "mpiexec -n $mpiprocess $programdir/$program -f $indir/${genome}_1.fa -m $indir/${genome}_2.fa -b 100000000 -r 100 -p $parts -o $outdir/${genome}+$parts+$mpiprocess+Zoltan+333937177+paired -t Zoltan -d 333937177"
    mpiexec -n $mpiprocess $programdir/$program -f $indir/${genome}_1.fa -m $indir/${genome}_2.fa -b 100000000 -r 100 -p $parts -o $outdir/${genome}+$parts+$mpiprocess+Zoltan+333937177+paired -t Zoltan -d 333937177
    echo "mpiexec -n $mpiprocess $programdir/$program -f $indir/${genome}_1.fa -b 100000000 -r 100 -p $parts -o $outdir/${genome}+$parts+$mpiprocess+ParMetis+333937177 -t ParMetis -d 333937177"
    mpiexec -n $mpiprocess $programdir/$program -f $indir/${genome}_1.fa -b 100000000 -r 100 -p $parts -o $outdir/${genome}+$parts+$mpiprocess+ParMetis+333937177 -t ParMetis -d 333937177
    echo "mpiexec -n $mpiprocess $programdir/$program -f $indir/${genome}_1.fa -m $indir/${genome}_2.fa -b 100000000 -r 100 -p $parts -o $outdir/${genome}+$parts+$mpiprocess+ParMetis+333937177+paired -t ParMetis -d 333937177"
    mpiexec -n $mpiprocess $programdir/$program -f $indir/${genome}_1.fa -m $indir/${genome}_2.fa -b 100000000 -r 100 -p $parts -o $outdir/${genome}+$parts+$mpiprocess+ParMetis+333937177+paired -t ParMetis -d 333937177
}

runexp () {
    outdir=/project/jing/boa/debug/$1
    mkdir -p $outdir
    runinstans $1 32 20 fruitfly $outdir
    runinstans $1 3 20 test $outdir

    pushd $outdir
        echo "" > $outdir.txt
        for folder in ./*; do
            for file in $folder in $folder/*; do
                md5sum $file >> $outdir.txt
            done
        done
    popd
}

runexp $exe1
runexp $exe2

# if the md5sums of this two file are the same then the partitioning results of them are the same
md5sum $output/$exe1.txt
md5sum $output/$exe2.txt

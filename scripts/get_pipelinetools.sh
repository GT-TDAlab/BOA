#!/bin/bash

pushd $(dirname "$0")
scriptdir="$(pwd)"
popd

withArt=1
withMinia=1
withQuast=1

download_dir=${scriptdir}/../thirdparty

#parse input
if [ "$#" -eq 3 ]; then
    withArt=$1
    withMinia=$2
    withQuast=$3
elif [ "$#" -eq 4 ]; then
    withArt=$1
    withMinia=$2
    withQuast=$3
    download_dir=$4
elif [ "$#" -ne 0 ]; then 
    echo "Usage: ./scripts/get_piplinetools.sh 'whether install art (0/1)' 'whether install Minia (0/1)' 'whether install Quast (0/1)' "
    echo "    'whether install art', 0 (No) or 1 (Yes), is default as 1"
    echo "    'whether install Minia', 0 (No) or 1 (Yes), is default as 1"
    echo "    'whether install Quast', 0 (No) or 1 (Yes), is default as 1"
    echo "    'where thirdparties tools should be download at' is default as './thirdparty'"
    exit 1
fi

source $scriptdir/download_lib.sh

#download art (no need to install)
if [[ "$withArt" -eq "1" ]]; then
    tool=art
    tool_v=2016.06.05
    tool_tarfile=artbinmountrainier${tool_v}linux64.tgz
    tool_link=https://www.niehs.nih.gov/research/resources/assets/docs/${tool_tarfile} #linux 64-bit
    tool_dir=$download_dir/art_bin_MountRainier
    download_lib $tool $tool_v $tool_tarfile $tool_link $tool_dir
fi

#download minia
if [[ "$withArt" -eq "1" ]]; then
    tool=minia
    tool_v=v0.0.102
    tool_tarfile=minia-${too_v}-bin-Linux.tar.gz
    tool_link=https://github.com/GATB/minia/releases/download/v0.0.102/${too_tarfile} #linux
    tool_dir=$download_dir/minia-v0.0.102-bin-Linux
    download_lib $tool $tool_v $tool_tarfile $tool_link $tool_dir
fi

if [[ "$withQuast" -eq "1" ]]; then
    clone_lib quast https://github.com/ablab/quast.git $download_dir/quast
fi

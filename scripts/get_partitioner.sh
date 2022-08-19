#!/bin/bash

pushd $(dirname "$0")
scriptdir="$(pwd)"
popd

withParMetis=1
download_dir=${scriptdir}/../thirdparty
config_file=${scriptdir}/../config.txt

#parse input
if [ "$#" -eq 1 ]; then
    withParMetis=$1
elif [ "$#" -eq 3 ]; then
    withParMetis=$1
    download_dir=$2
    config_file=$3
elif [ "$#" -ne 0 ]; then 
    echo "Usage: ./scripts/get_thirdparties.sh 'whether install Zoltan with ParMetis(0/1)' 'where the thirdparties should be download at' 'where config.txt is'"
    echo "    'whether Zoltan install with ParMetis', 0 (No) or 1 (Yes), is default as 1"
    echo "    'where thirdparties (Zoltan, ParMetis) should be download at' is default as './thirdparty'"
    echo "    'where config.txt is at' is default as './config.txt'"
    exit 1
fi

source $scriptdir/download_lib.sh

install_parmetis () {
    dir=$1
    pm_build_dir=$2
    pm_inc_dir=$3
    pm_lib_dir=$4
    if [ ! -d $pm_build_dir ]; then
        mkdir -p $pm_build_dir
    fi
    if [ -z "$(ls -A $pm_build_dir)" ]; then
        pushd $dir

        #update idx type width to 64 bits
        headerfile=./metis/include/metis.h
        awk '/#define IDXTYPEWIDTH/{gsub(/32/, "64")};{print}' ${headerfile} > ${headerfile}.temp
        mv ${headerfile}.temp ${headerfile}

        mkdir -p $pm_build_dir
        echo "--Installing parmetis to $pm_build_dir"
        make config prefix=${pm_build_dir} > installlog
        make install >> installlog

        #since zoltan needs both metis and parmetis, so we are copying them together
        cp $headerfile $pm_inc_dir
        cp $pm_build_dir/Linux-x86_64/libmetis/libmetis.a  $pm_lib_dir

        rm installlog
        popd
    else
        echo "--ParMetis's build folder is not empty, so skip building and installing"
        echo "    $pm_build_dir"
    fi
}

pm=parmetis
pm_v=4.0.3
pm_tarfile=parmetis-${pm_v}.tar.gz
pm_link=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/$pm_tarfile
pm_dir=$download_dir/${pm}-$pm_v
pm_build_dir=$pm_dir/build
pm_inc_dir=$pm_build_dir/include
pm_lib_dir=$pm_build_dir/lib
if [[ "$withParMetis" -eq "1" ]]; then
    #download ParMetis
    download_lib $pm $pm_v $pm_tarfile $pm_link $pm_dir
    #install ParMetis
    install_parmetis $pm_dir $pm_build_dir $pm_inc_dir $pm_lib_dir
fi

#download zoltan
z=Zoltan
z_v=3.83
z_tarfile=v${z_v}.tar.gz
z_link=https://github.com/sandialabs/Zoltan/archive/$z_tarfile
z_dir=$download_dir/Zoltan-$z_v
download_lib $z $z_v $z_tarfile $z_link $z_dir

#install zoltan
z_build_dir=$z_dir/build
if [ ! -d $z_build_dir ]; then
    mkdir -p $z_build_dir
    pushd $z_build_dir
    echo "--Installing zoltan to $z_build_dir"
    if [[ "$withParMetis" -eq "1" ]]; then
        ../configure --prefix=${z_build_dir} --with-id-type=ulong --with-parmetis --with-parmetis-incdir=$pm_inc_dir --with-parmetis-libdir=$pm_lib_dir > installlog
    else
        ../configure --prefix=${z_build_dir} --with-id-type=ulong > installlog
    fi
    make everything >> installlog
    make install >> installlog
    rm installlog
    popd
else
    echo "--Zoltan's built folder exists already, so skip building and installing:"
    echo "    $z_build_dir"
fi

#update config.txt
pm_install_var=PARMETIS_INSTALL
z_install_var=ZOLTAN_INSTALL
z_home_var=ZOLTAN_HOME
echo "--Update the config file at $config_file"
if [ -f $config_file ]; then
    if [[ "$withParMetis" -eq "1" ]]; then
        awk -v var=${pm_install_var} '$0!~var' $config_file > ${config_file}.temp
        mv  ${config_file}.temp ${config_file}
    fi
    awk -v var=${z_install_var} '$0!~var' $config_file > ${config_file}.temp
    awk -v var=${z_home_var} '$0!~var' ${config_file}.temp > ${config_file}
    rm ${config_file}.temp
fi

echo "$pm_install_var=$pm_build_dir" >> $config_file
echo "$z_install_var=$z_build_dir" >> $config_file
echo "$z_home_var=$z_dir" >> $config_file

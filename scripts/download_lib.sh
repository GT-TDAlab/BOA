#!/bin/bash

#download and untar function
download_lib () {
    libname=$1
    v=$2
    tarfile=$3
    link=$4
    dir=$5

    if [ ! -d $dir ]; then
        echo "--downloading $libname v$v to $download_dir"
        mkdir -p $download_dir
        pushd $download_dir
        wget $link > log
        echo "--Unzipping the package"
        tar -xvf ./$tarfile >> log
        rm log
        rm ./$tarfile
        popd
        if [ ! -d $dir ]; then
            echo "$dir doesn't exist after downloading"
            exit 1
        fi
    else
        echo "--$libname is downloaded already, skip downloading"
    fi
}

clone_lib () {
    libname=$1
    link=$2
    dir=$3
    if [ ! -d $dir ]; then
        echo "--downloading $libname latest version from git source"
        mkdir -p $download_dir
        pushd $download_dir
        git clone $link > log
        popd
        if [ ! -d $dir ]; then
            echo "$dir doesn't exist after downloading"
            exit 1
        fi
    else
        echo "--$libname is downloaded already, skip downloading"
    fi
}

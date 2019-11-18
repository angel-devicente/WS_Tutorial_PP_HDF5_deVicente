#!/bin/bash

tar -zxf libpng-1.6.37.tar.gz
cd libpng-1.6.37

module purge
module load gnu/7.2.0
module load szip/gnu

./configure --prefix=$PWD/../libpng
make
make install
